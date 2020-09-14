#include <cmath>   
#include <iostream>

using namespace std;

#include "state.hh"
#include "timers.hh"

State::State(const Details &details, const Data &data, const Model &model, const ObservationModel &obsmodel) : comp(model.comp), trans(model.trans), param(model.param), details(details), data(data), model(model), obsmodel(obsmodel)
{
	disc_spline.resize(model.spline.size());

	transev.resize(details.ndivision); 
	indev.resize(data.popsize);

	Qmap.resize(details.ndivision);
	for(auto sett = 0u; sett < details.ndivision; sett++){
		Qmap[sett].resize(data.narage); for(auto v = 0u; v < data.narage; v++) Qmap[sett][v] = 0;
	}
	
	popw.resize(data.nardp);                                      
	lambda.resize(data.nsettardp);
}


/// This simulates from the compartmental model and generates an event list for individual i
/// The individual states in state c at time t 
void State::simulate_compartmental_transitions(unsigned int i, unsigned int c, double t)
{
	vector <Event> &evlist = indev[i];

	auto timep = 0u;                                       // Works out which time period t is in
	while(timep < model.ntime_period && t > model.time_period[timep].tend) timep++;

	Event ev;
	ev.ind = i; ev.timep = timep; 
		
	auto a = data.democatpos[data.ind[i].dp][0];           // Sets the age of the individual
		
	unsigned int tra; 
	if(c == 0){                                            // If an infection then adds first 
		evlist.clear();	
		tra = 0;
		ev.trans = tra; ev.t = t; 
		evlist.push_back(ev);
		c = trans[tra].to;
	}
	
	while(comp[c].trans.size() > 0){                       // Generates the event sequence as a function of time
		auto tra = select_branch(c,a);                       // Selects which branch to go down
		
		auto dt = sample_duration(tra);                      // Samples the time for the transition
		
		t += dt;
		
		while(timep < model.ntime_period-1 && model.time_period[timep].tend < t){    // Adds in changes in time period
			ev.trans = comp[c].transtimep; ev.t = model.time_period[timep].tend;
			evlist.push_back(ev);
			timep++;
			ev.timep = timep; 
		}
		
		ev.trans = tra; ev.t = t;
		evlist.push_back(ev);

		c = trans[tra].to; 
	}
}

/// Selects which branch to follow
unsigned int State::select_branch(unsigned int c, unsigned int a)
{
	unsigned int tra, kmax = comp[c].trans.size();
	if(kmax == 1) return comp[c].trans[0];
	else{                                                // Uses the branching probability to select branch
		double z = ran(); auto k = 0u; while(k < kmax && z > comptransprob[c].probsum[a][k]) k++;
		if(k == kmax) emsgEC("Model",2);
		tra = comp[c].trans[k];
	}
	return tra;
}
		
/// Samples the time taken for a particular transition to happen
double State::sample_duration(unsigned int tra)
{
	double dt;
	
	switch(trans[tra].type){                             // Samples duration spent in compartment
	case EXP_DIST:
		dt = exp_sample_time(paramval[trans[tra].param_mean]);
		break;
	
	case GAMMA_DIST:
		{
			auto mean = paramval[trans[tra].param_mean]; auto sd = paramval[trans[tra].param_cv]*mean;
			dt = gamma_sample(mean*mean/(sd*sd),mean/(sd*sd));
		}
		break;
		
	case LOGNORM_DIST:
		{
			auto mean_ns = paramval[trans[tra].param_mean], cv_ns = paramval[trans[tra].param_cv];
			auto sd = sqrt(log((1+cv_ns*cv_ns))), mean = log(mean_ns) - sd*sd/2;
			dt = lognormal_sample(mean,sd); 
		}
		break;
		
	default: emsgEC("Model",3); break;
	}

	if(dt < TINY) dt = TINY;
	
	return dt;
}

/// Calculates the latent process likelihood and puts the results in Lev
void State::set_process_likelihood()
{		
	for(auto c = 0u; c < data.narea; c++){
		for(auto dp = 0u; dp < data.ndemocatpos; dp++){
			auto w = c*data.ndemocatpos + dp;
			popw[w] = data.area[c].ind[dp].size();
		}
	}		
			
	Lev = 0.0;
	auto t = 0.0; 
	auto n = 0u;
	for(auto sett = 0u; sett < details.ndivision; sett++){
		auto phi = disc_spline[model.phispline_ref][sett]; 
		auto beta = disc_spline[model.betaspline_ref][sett];
	
		auto tmax = details.division_time[sett+1];
		
		for(auto c = 0u; c < data.narea; c++){
			auto fac = beta*areafactor[c];
			for(auto dp = 0u; dp < data.ndemocatpos; dp++){
				auto w = c*data.ndemocatpos + dp;
				auto v = c*data.nage + data.democatpos[dp][0];
				lambda[w] = susceptibility[dp]*(fac*Qmap[sett][v] + phi);		
				if(lambda[w] < 0) emsgEC("State",46);
				
				Lev -= lambda[w]*popw[w]*(tmax-t);
			}
		}
		
		while(n < infev.size()){
			auto i = infev[n].ind;
			auto ev = indev[i][infev[n].e];
			auto tt = ev.t;
			if(tt >= tmax) break;
	
			t = tt;
			
			auto c = data.ind[i].area;
			auto w = c*data.ndemocatpos + data.ind[i].dp;
			Lev += log(lambda[w]);
			popw[w]--;
			n++;
			
			Lev += lambda[w]*(tmax-t);
		}
		t = tmax;
	} 
	if(std::isnan(Lev)) emsgEC("State",47);
}


/// Checks the latent process likelihood and prior are correct 
void State::check_Lev_and_Pr()
{
	auto Levst = Lev;
	set_process_likelihood();
	double dd = Levst - Lev; if(dd*dd > TINY) emsgEC("State",49);
	
	dd = Pr - model.prior(paramval); if(sqrt(dd*dd) > TINY) emsgEC("State",51);
}
	
	
/// Clears all the infections from the state
void State::clear()
{
	for(const auto& i : infev) indev[i.ind].clear();
	
	infev.clear();
	transev.clear(); transev.resize(details.ndivision);
}


/// Copies the state from another state
void State::copy_state(const State &from)
{
	paramval = from.paramval;
	L = from.L;
	EF = from.EF;
	Pr = from.Pr;
	
	for(const auto& iev : infev) indev[iev.ind].clear();
	for(const auto& iev : from.infev) indev[iev.ind] = from.indev[iev.ind];
		
	infev = from.infev;
	transev = from.transev;
	Qmap = from.Qmap;
		
	disc_spline = from.disc_spline;
	susceptibility = from.susceptibility;
	areafactor = from.areafactor;
	comptransprob = from.comptransprob;
}


/// Checks quantities in the state are correct
void State::check() const
{
	for(auto j = 0u; j < infev.size(); j++){                         // Checks bounds
		auto i = infev[j].ind, e = infev[j].e;
		if(i >= indev.size()) emsgEC("State",16);
		if(e >= indev[i].size()) emsgEC("State",17);
	}
	
	for(auto j = 0u; j < infev.size()-1; j++){                       // Checks order of infection events
		if(get_infection_time(j) > get_infection_time(j+1)) emsgEC("State",18);
	}
	
	for(const auto& iev : indev){                                    // Looks at individual event sequences 
		auto emax = iev.size();
		if(emax > 0){
			auto c = 0u; 
			double t = 0;
			for(const auto& ev : iev){                         
				auto tt = ev.t; if(tt <= t) emsgEC("State",19);            // CHecks time ordering
				auto tra = ev.trans;
				if(trans[tra].from != c) emsgEC("State",20);               // Checks for consistency
				c = trans[tra].to; t = tt;
			}
			if(comp[c].trans.size() != 0) emsgEC("State",21);
			
			for(auto timep = 0u; timep < model.ntime_period; timep++){   // Checks time periods are correct
				t = model.time_period[timep].tend;
				if(t > iev[0].t && t < iev[emax-1].t){
					unsigned int e;
					for(e = 0u; e < emax; e++) if(iev[e].t == t) break;
					if(timep < model.ntime_period-1){
						if(e == emax) emsgEC("State",22);
					}
					else{
						if(e != emax) emsgEC("State",23);
					}
				}
			}
		}
	}
	
	vector < vector <int> > done;
	done.resize(indev.size());
	auto num = 0u;
	for(auto i = 0u; i < indev.size(); i++){                               // Check that infev correspond to indev
		if(indev[i].size() > 0){
			num++;
			done[i].resize(indev[i].size());
			for(auto e = 0u; e < indev[i].size(); e++) done[i][e] = 0;
		}
	}
	if(num != infev.size()) emsgEC("State",40);

	for(auto sett = 0u; sett < details.ndivision; sett++){                 // Chechs transev consistent with indev
		for(auto &trev : transev[sett]){
			auto i = trev.ind, e = trev.e;
			if(e >= indev[i].size()) emsgEC("State",41);
			
			auto se = (unsigned int)(details.ndivision*indev[i][e].t/details.period); 
			if(se != sett) emsgEC("State",42);
			if(done[i][e] != 0) emsgEC("State",43);
			done[i][e] = 1;
		}
	}
	
	for(auto i = 0u; i < indev.size(); i++){
		for(auto e = 0u; e < indev[i].size(); e++){
			if(indev[i][e].t < details.period){
				if(done[i][e] != 1) emsgEC("State",44);
			}
			else{
				if(done[i][e] != 0) emsgEC("State",45);
			}
		}
	}
	
	double dd;                                                        // Checks observation likelihood and prior correct
	dd = L - obsmodel.observation_likelihood(transev,indev); if(sqrt(dd*dd) > TINY) emsgEC("State",1);
	dd = Pr - model.prior(paramval); if(sqrt(dd*dd) > TINY) emsgEC("State",2);
}


/// Sets up the state with a specified set of parameters
/// If the values for the branching probabilities not possible then returns fail.
Status State::set_param(const vector <double> &paramv)
{
	paramval = paramv;
	
	if(model.create_comptransprob(comptransprob,paramval) == 1) return FAIL;
	for(auto sp = 0u; sp < model.spline.size(); sp++) disc_spline[sp] = model.create_disc_spline(sp,paramval);
	susceptibility = model.create_susceptibility(paramval);    
	areafactor = model.create_areafactor(paramval);    
	
	return SUCCESS;
}


/// Sets the values of phi and beta (this is used to speed up computation)
void State::set_beta_and_phi(unsigned int sett)
{	
	phi = disc_spline[model.phispline_ref][sett]; 
	beta = disc_spline[model.betaspline_ref][sett]; 
}


/// Sets the observation likelihood and prior
void State::set_L_and_Pr()
{
	L = obsmodel.observation_likelihood(transev,indev);
	Pr = model.prior(paramval);
}


/// Sets the error function
void State::set_EF()
{
	EF = obsmodel.observation_likelihood(transev,indev);
}


/// Gets the time of an infection event
double State::get_infection_time(unsigned int n) const
{
	if(n == infev.size()) return LARGE;
	return indev[infev[n].ind][infev[n].e].t;
}


/// Adds an individual event sequence defined in infev[i] to transev
void State::add_indev(unsigned int i)
{
	unsigned int e, emax, se;
	EventRef evref;
	
	auto iev = indev[i];
	emax = iev.size(); if(emax == 0) return;
	
	auto fac = double(details.ndivision)/details.period;
	
	evref.ind = i; evref.e = 0;
	infev.push_back(evref);
	for(e = 0; e < emax; e++){
		evref.e = e;
		se = (unsigned int)(fac*iev[e].t); 
		if(se < details.ndivision) transev[se].push_back(evref);
	}
}


/// For a given time sett, sets Qmap by using the Qmap from another state and the difference given by dQmap
void State::set_Qmap_using_dQ(unsigned int sett, const State &state, const vector <double> &dQmap)
{
	auto& stQmap = state.Qmap[sett];
	auto& neQmap = Qmap[sett];
	for(auto v = 0u; v < data.narage; v++){
		double val = stQmap[v] + dQmap[v];
		if(val < -TINY){ cout << val << "val\n"; emsgEC("State",1);}
		if(val < 0) val = 0;	
		neQmap[v] = val;
	}
}


/// Sets the infection map Qmap as a function of time based in transition events transev 
void State::set_Qmap(unsigned int check)
{
	vector <double> Qma(data.narage);
	
 	for(auto& Qm : Qma) Qm = 0;

	auto nage = data.nage;
	for(auto sett = 0u; sett < details.ndivision; sett++){
		for(auto v = 0u; v < data.narage; v++){
			auto val = Qma[v];
			if(check == 1){
				if(val < -TINY) emsgEC("State",2);
				if(val < Qmap[sett][v]-TINY || val > Qmap[sett][v]+TINY) emsgEC("State",3);
			}
			if(val < 0){ val = 0; Qma[v] = 0;}	
			
			Qmap[sett][v] = val;
		}
		
		for(const auto& tre : transev[sett]){
			auto i = tre.ind;
			Event fev = indev[i][tre.e];

			auto v = data.ind[i].area*data.nage+data.democatpos[data.ind[i].dp][0];
			auto dq = trans[fev.trans].DQ[fev.timep];
			if(dq != UNSET){
				for(auto loop = 0u; loop < 2; loop++){
					auto q = model.DQ[dq].q[loop];
					if(q != UNSET){
						auto fac = model.DQ[dq].fac[loop];
						
						auto& genQten = data.genQ.Qten[data.Q[q].Qtenref];
						auto kmax = genQten.ntof[v];
						auto& cref = genQten.tof[v];
						auto& valref = genQten.valf[v];
						if(nage == 1){
							for(auto k = 0u; k < kmax; k++){
								Qma[cref[k]*nage] += fac*valref[k][0];
							}
						}
						else{
							for(auto k = 0u; k < kmax; k++){
								auto vv = cref[k]*nage;	
								for(auto a = 0u; a < nage; a++){
									Qma[vv] += fac*valref[k][a];
									vv++;
								}
							}
						}
					}
				}
			}
		}
	}
}
	

/// Used for time ordering event references	
static bool compEventRefTime(EventRefTime lhs, EventRefTime rhs)
{
	return lhs.t < rhs.t;
};


/// Time orders infection events infev
void State::sort_infev()
{
	vector <EventRefTime> xt;
	for(const auto& iev : infev){
		EventRefTime evreft;
		evreft.ind = iev.ind; evreft.e = iev.e;
		if(indev[iev.ind].size() == 0) emsgEC("State",54);
		
		evreft.t = indev[iev.ind][iev.e].t;
		xt.push_back(evreft);	
	}
	sort(xt.begin(),xt.end(),compEventRefTime);

	for(auto i = 0u; i < infev.size(); i++){
		infev[i].ind = xt[i].ind;	infev[i].e = xt[i].e;
	}
}


/// Initialises the state based on a particle (used for abc methods)
void State::initialise_from_particle(const Particle &part)
{
	paramval = part.paramval;
	
	clear();                                               // Removes the existing initial sequence 

	vector <int> indlist;
	for(const auto& ev : part.ev){                         // Decompresses the events stored in the particles
		int i = ev.ind;
		if(indev[i].size() == 0) indlist.push_back(i);
		indev[i].push_back(ev);
	}	
	
	for(auto i : indlist) add_indev(i);
	
	sort_infev();
	
	set_Qmap(0);

	EF = part.EF;
	Pr = model.prior(paramval);
	
	if(EF != obsmodel.observation_likelihood(transev,indev)) emsg("Observation does not agree");
}
		
		
/// Generates a particle from the state (used for abc methos)
void State::generate_particle(Particle &part) const
{
	part.EF = EF;
	part.paramval = paramval;
	
	vector <Event> store;
	for(const auto& inde : indev){                           // Compresses the events to take up as little memory as possible 
		for(const auto& ev : inde) store.push_back(ev);
	}
	part.ev = store;
}


/// STANDARD PARAMETER PROPOSALS 

/// The functions below all refer to changes in parameter with fixed event sequence
void State::standard_parameter_prop(Jump &jump)
{	
	stand_param_betaphi_prop(jump);
	stand_param_area_prop(jump);
	stand_param_compparam_prop(jump);

	if(checkon == true) check_Lev_and_Pr();
}


/// Makes proposal to beta and phi
void State::stand_param_betaphi_prop(Jump &jump)
{	
	timers.timebetaphiinit -= clock();
	PrecalcBetaPhi precalc;
	
	likelihood_beta_phi_initialise(precalc);
	
	vector <unsigned int> parampos;
	for(const auto& spl : model.spline[model.betaspline_ref]) parampos.push_back(spl.param);
	for(const auto& spl : model.spline[model.phispline_ref]) parampos.push_back(spl.param);
		
	timers.timebetaphiinit += clock();
		
	unsigned int loopmax = 12/parampos.size(); if(loopmax == 0) loopmax = 1;
	
	timers.timebetaphi -= clock();
	for(auto loop = 0u; loop < loopmax; loop++){
		for(auto th : parampos){
			if(param[th].min != param[th].max){
				auto param_store = paramval[th];	
				paramval[th] += normal_sample(0,jump.stand[th]);               // Makes a change to a parameter

				vector < vector <double> > disc_spline_prop;

				double al=0, Lev_prop=Lev, Pr_prop=Pr;
				if(model.inbounds(paramval[th],th) == true){
					for(auto sp = 0u; sp < model.spline.size(); sp++) disc_spline_prop.push_back(model.create_disc_spline(sp,paramval));

					Lev_prop = likelihood_beta_phi(disc_spline_prop,precalc);
					
					if(smooth_spline == 1) Pr_prop = model.prior(paramval);
					al = exp(Pr_prop - Pr + Lev_prop-Lev);
				}
			
				if(ran() < al){
					Pr = Pr_prop;
					Lev = Lev_prop;
					disc_spline = disc_spline_prop;
					jump.stand_accept(th);
				}
				else{
					paramval[th] = param_store;
					jump.stand_reject(th);					
				}
			}
		}
	}
	timers.timebetaphi += clock();
}


/// Initialise quantities for fast likelihood calculation
void State::likelihood_beta_phi_initialise(PrecalcBetaPhi &precalc)
{
	vector <unsigned int> map;
	map.resize(data.nardp); for(auto& ma : map) ma = 0;
	
	popw.resize(data.nardp);                                        
	for(auto c = 0u; c < data.narea; c++){
		for(auto dp = 0u; dp < data.ndemocatpos; dp++){
			auto w = c*data.ndemocatpos + dp;
			popw[w] = data.area[c].ind[dp].size();
		}
	}		
	
	precalc.betafac.resize(details.ndivision); precalc.phifac.resize(details.ndivision);
	
	auto t = 0.0; 
	auto n = 0u;
	
	for(auto sett = 0u; sett < details.ndivision; sett++){
		auto tmax = details.division_time[sett+1];
		vector <BetaPhiFactors> lcontlist;
		
		auto betasum = 0.0, phisum = 0.0;
		for(auto c = 0u; c < data.narea; c++){
			auto fac = areafactor[c];
			for(auto dp = 0u; dp < data.ndemocatpos; dp++){
				auto w = c*data.ndemocatpos + dp;
				auto v = c*data.nage + data.democatpos[dp][0];	
				betasum -= fac*susceptibility[dp]*Qmap[sett][v]*popw[w]*(tmax-t);
				phisum -= susceptibility[dp]*popw[w]*(tmax-t);
			}
		}
		
		while(n < infev.size()){
			auto i = infev[n].ind;
			Event ev = indev[i][infev[n].e];
			auto tt = ev.t;
			if(tt >= tmax) break;
	
			t = tt;
			
			auto c = data.ind[i].area;
			auto fac = areafactor[c];
			
			auto dp = data.ind[i].dp;
			auto w = c*data.ndemocatpos + dp;
			auto v = c*data.nage + data.democatpos[dp][0]; 
			
			if(map[w] == 0){
				BetaPhiFactors lcont;
				lcont.w = w; lcont.betafac = fac*susceptibility[dp]*Qmap[sett][v]; lcont.phifac = susceptibility[dp]; lcont.num = UNSET;
				lcontlist.push_back(lcont);
			}
			map[w]++;
		
			betasum += fac*susceptibility[dp]*Qmap[sett][v]*(tmax-t);
			phisum += susceptibility[dp]*(tmax-t);
			popw[w]--;
			n++;
		}
		
		precalc.betafac[sett] = betasum; precalc.phifac[sett] = phisum;
					
		for(auto& lcont : lcontlist){
			auto w = lcont.w;
			lcont.num = map[w];
			map[w] = 0;
		}
					
		precalc.lc.push_back(lcontlist);
		
		t = tmax;
	}
}


/// Calculates the latent likelihood using pre-calculated quantities
double State::likelihood_beta_phi(const vector < vector <double> > &disc_spline, const PrecalcBetaPhi &precalc) const 
{
	auto Lev = 0.0;
	for(auto sett = 0u; sett < details.ndivision; sett++){
		auto phi = disc_spline[model.phispline_ref][sett]; 
		auto beta = disc_spline[model.betaspline_ref][sett];

		Lev += precalc.betafac[sett]*beta + precalc.phifac[sett]*phi;
		
		for(const auto& l : precalc.lc[sett]) Lev += l.num*log(l.betafac*beta + l.phifac*phi);
		if(std::isnan(Lev)) emsgEC("State",52);
	}
	
	return Lev;
}

/// Makes proposal to change factors affecting transmission rates in areas
void State::stand_param_area_prop(Jump &jump)
{	
	PrecalcArea precalc;
	
	likelihood_area_initialise(precalc);

	timers.timecovar -= clock();
	
	vector <unsigned int> thlist;
	for(auto th : model.covariate_param) thlist.push_back(th); 
	
	if(model.region_effect > 0){
		for(auto th : model.regioneffect_param) thlist.push_back(th); 
	}
	
	if(model.region_effect == 1) thlist.push_back(model.sigma_param);
		
	unsigned int loopmax = 12/thlist.size(); if(loopmax == 0) loopmax = 1;
	
	for(auto loop = 0u; loop < loopmax; loop++){
		for(auto th : thlist){
			if(param[th].min != param[th].max){
				auto param_store = paramval[th];	
				paramval[th] += normal_sample(0,jump.stand[th]);               // Makes a change to a parameter

				double al=0, Lev_prop=Lev, Pr_prop=Pr;
				vector <double> areafactor_prop;
				if(model.inbounds(paramval[th],th) == true){
					areafactor_prop = model.create_areafactor(paramval); 
				
					Lev_prop = likelihood_area(areafactor_prop,precalc);
					Pr_prop = model.prior(paramval);
					
					al = exp(Pr_prop-Pr + Lev_prop-Lev);
				}

				if(ran() < al){
					Lev = Lev_prop;
					Pr = Pr_prop;
					areafactor = areafactor_prop;
					
					jump.stand_accept(th);
				}
				else{
					paramval[th] = param_store;
					jump.stand_reject(th);
				}
			}
		}
	}
	timers.timecovar += clock();
}


/// Initialise quantities for fast likelihood calculation
void State::likelihood_area_initialise(PrecalcArea &precalc)
{
	timers.timecovarinit -= clock();

	vector <double> lambdaareafactor, lambdaphifac;
	lambdaareafactor.resize(data.nardp);
	lambdaphifac.resize(data.nardp);

	precalc.mult.resize(data.narea);
	precalc.add.resize(data.narea);
	precalc.areasum.resize(data.narea);
	
	for(auto c = 0u; c < data.narea; c++){
		precalc.areasum[c] = 0;
		for(auto dp = 0u; dp < data.ndemocatpos; dp++){
			auto w = c*data.ndemocatpos + dp;
			popw[w] = data.area[c].ind[dp].size();
		}
	}		
	
	auto L0 = 0.0, t = 0.0;
	auto n = 0u;
	for(auto sett = 0u; sett < details.ndivision; sett++){
		auto phi = disc_spline[model.phispline_ref][sett]; 
		auto beta = disc_spline[model.betaspline_ref][sett];
	
		auto tmax = details.division_time[sett+1];
		
		for(auto c = 0u; c < data.narea; c++){
			for(auto dp = 0u; dp < data.ndemocatpos; dp++){
				auto w = c*data.ndemocatpos + dp;
				auto v = c*data.nage + data.democatpos[dp][0];
				lambdaareafactor[w] = susceptibility[dp]*beta*Qmap[sett][v];
				lambdaphifac[w] = susceptibility[dp]*phi;
				
				precalc.areasum[c] -= lambdaareafactor[w]*popw[w]*(tmax-t);
				L0 -= lambdaphifac[w] *popw[w]*(tmax-t);
			}
		}
		
		while(n < infev.size()){
			auto i = infev[n].ind;
			Event ev = indev[i][infev[n].e];
			auto tt = ev.t;
			if(tt >= tmax) break;
	
			t = tt;
			
			auto c = data.ind[i].area;
			auto w = c*data.ndemocatpos + data.ind[i].dp;
			
			precalc.mult[c].push_back(lambdaareafactor[w]);
			precalc.add[c].push_back(lambdaphifac[w]);
		
			popw[w]--;
			n++;
			
 			precalc.areasum[c] += lambdaareafactor[w]*(tmax-t);
			L0 += lambdaphifac[w]*(tmax-t);
		}
		
		t = tmax;
	} 
	precalc.L0_store = L0;
	
	timers.timecovarinit += clock();
}


/// Calculates the latent likelihood using pre-calculated quantities
double State::likelihood_area(const vector <double> &areafactor, const PrecalcArea &precalc) const
{
	double Lev = precalc.L0_store;
	for(auto c = 0u; c < data.narea; c++){
		auto fac = areafactor[c];
		Lev += precalc.areasum[c]*fac;
		auto kmax = precalc.mult[c].size();
		for(auto k = 0u; k < kmax; k++) Lev += log(precalc.mult[c][k]*fac + precalc.add[c][k]);
	}
	if(std::isnan(Lev)) emsgEC("State",53);
	
	return Lev;
}
					
					
/// Makes proposals to compartmental parameters (branching probabilities and transition distributions)
void State::stand_param_compparam_prop(Jump &jump)
{	
	timers.timecompparam -= clock();

	vector <PrecalcCompParam> precalc(trans.size());
	
	for(auto tra = 0u; tra < trans.size(); tra++){
		if(trans[tra].istimep == false){
			precalc[tra].num.resize(data.nage); for(auto& num : precalc[tra].num) num = 0;
		}
		precalc[tra].numvisittot = 0; precalc[tra].dtsum = 0; precalc[tra].dtlist.clear();
	}

	for(const auto& iev : infev){                  // Extracts values based on the event sequence
		auto i = iev.ind;
		auto dp = data.ind[i].dp;
		auto a = data.democatpos[dp][0];
				
		auto t = 0.0;
		for(const auto& ev : indev[i]){
			auto tra = ev.trans;
			if(trans[tra].istimep == false){
				precalc[tra].num[a]++;
	
				auto dt = ev.t-t;
				if(dt == 0) emsgEC("Model",10);
				t += dt;
				
				precalc[tra].numvisittot++;
				switch(trans[tra].type){
				case EXP_DIST: precalc[tra].dtsum += dt; break;
				case GAMMA_DIST: case LOGNORM_DIST: precalc[tra].dtlist.push_back(dt); break;
				case INFECTION_DIST: break;
				default: emsgEC("model",99); break;
				}
			}
		}
	}

	auto Li_dt = likelihood_dt(precalc,paramval);
	auto Li_prob = likelihood_prob(precalc,comptransprob);

	for(auto th = 0u; th < param.size(); th++){
		if((param[th].type == DISTVAL_PARAM || param[th].type == BRANCHPROB_PARAM) && param[th].min != param[th].max){
			vector <double> param_store = paramval;	
		
			paramval[th] += normal_sample(0,jump.stand[th]);               // Makes a change to a parameter
			
			auto dL = 0.0, Lp_prob = Li_prob, Pr_prop = Pr;
			auto flag=0u;

			vector <CompTransProb> comptransprobp;
			if(model.inbounds(paramval[th],th) == true){
				if(param[th].type == BRANCHPROB_PARAM){
					if(model.create_comptransprob(comptransprobp,paramval) == 1) Lp_prob = -LARGE;
					else Lp_prob = likelihood_prob(precalc,comptransprobp);
				}
				else Lp_prob = Li_prob;
				
				dL = dlikelihood_dt(precalc,param_store,paramval);
				Pr_prop = model.prior(paramval);
				
				auto al = exp(dL+ Lp_prob - Li_prob + Pr_prop - Pr);
				if(ran() < al) flag = 1; else flag = 0;
			}
			
			if(flag == 1){
				Li_dt += dL;
				Li_prob = Lp_prob;
				Pr = Pr_prop;
				if(param[th].type == BRANCHPROB_PARAM) comptransprob = comptransprobp;
				
				jump.stand_accept(th);
			}
			else{
				paramval = param_store;
				jump.stand_reject(th);
			}	
		}
	}
	
	if(checkon == true){
		double dd;
		dd = likelihood_dt(precalc,paramval)-Li_dt; if(dd*dd > TINY) emsgEC("Model",13);
		dd = likelihood_prob(precalc,comptransprob)-Li_prob; if(dd*dd > TINY) emsgEC("Model",14);
	}
	
	timers.timecompparam += clock();
}


/// Calculates the likelihood relating to branching probabilities
double State::likelihood_prob(vector <PrecalcCompParam> &precalc, vector <CompTransProb> &comptransprob) const
{
	auto L = 0.0;
	for(auto c = 0u; c < comp.size(); c++){	
		auto kmax = comp[c].trans.size();
		if(kmax > 1){
			for(auto k = 0u; k < kmax; k++){
				for(auto a = 0u; a < data.nage; a++){
					auto num = precalc[comp[c].trans[k]].num[a];
					if(num > 0) L += num*log(comptransprob[c].prob[a][k]);
				}
			}
		}
	}
	return L;
}
			
			
/// Calculates the likelihood for transition distributions
double State::likelihood_dt(vector <PrecalcCompParam> &precalc, vector <double> &paramv) const
{
	auto L = 0.0;
	for(auto tra = 0u; tra < trans.size(); tra++){
		switch(trans[tra].type){
		case EXP_DIST:
			{
				auto r = 1.0/paramv[trans[tra].param_mean];
				L += precalc[tra].numvisittot*log(r) - r*precalc[tra].dtsum;
			}
			break;
		
		case GAMMA_DIST:
			{
				auto mean = paramv[trans[tra].param_mean], sd = paramv[trans[tra].param_cv]*mean;
				for(auto dt : precalc[tra].dtlist) L += gamma_probability(dt,mean*mean/(sd*sd),mean/(sd*sd));
			}
			break;
			
		case LOGNORM_DIST:
			{
				auto mean_ns = paramv[trans[tra].param_mean], cv_ns = paramv[trans[tra].param_cv];
				auto sd = sqrt(log(1+cv_ns*cv_ns)), mean = log(mean_ns) - sd*sd/2;
				for(auto dt : precalc[tra].dtlist) L += lognormal_probability(dt,mean,sd*sd);
			}
			break;
		}
	}
	
	return L;
}


/// Calculates the change in likelihood for a given change in parameters
double State::dlikelihood_dt(vector <PrecalcCompParam> &precalc, vector <double> &paramvi, vector <double> &paramvf) const
{
	auto L = 0.0;
	for(auto tra = 0u; tra < trans.size(); tra++){
		switch(trans[tra].type){
		case EXP_DIST:
			{
				auto ri = 1.0/paramvi[trans[tra].param_mean], rf = 1.0/paramvf[trans[tra].param_mean];
				if(ri != rf){
					L += precalc[tra].numvisittot*log(rf) - rf*precalc[tra].dtsum;
					L -= precalc[tra].numvisittot*log(ri) - ri*precalc[tra].dtsum;
				}
			}
			break;
		
		case GAMMA_DIST:
			{
				auto meani = paramvi[trans[tra].param_mean], sdi = paramvi[trans[tra].param_cv]*meani;
				auto meanf = paramvf[trans[tra].param_mean], sdf = paramvf[trans[tra].param_cv]*meanf;
		
				if(meani != meanf || sdi != sdf){
					for(auto dt : precalc[tra].dtlist){
						L += gamma_probability(dt,meanf*meanf/(sdf*sdf),meanf/(sdf*sdf));
						L -= gamma_probability(dt,meani*meani/(sdi*sdi),meani/(sdi*sdi));
					}
				}				
			}
			break;
			
		case LOGNORM_DIST:
			{
				auto mean_nsi = paramvi[trans[tra].param_mean], cv_nsi = paramvi[trans[tra].param_cv];
				auto mean_nsf = paramvf[trans[tra].param_mean], cv_nsf = paramvf[trans[tra].param_cv];
				if(mean_nsi != mean_nsf || cv_nsi != cv_nsf){
					auto sdi = sqrt(log(1+cv_nsi*cv_nsi)), meani = log(mean_nsi) - sdi*sdi/2;
					auto sdf = sqrt(log(1+cv_nsf*cv_nsf)), meanf = log(mean_nsf) - sdf*sdf/2;
		
					for(auto dt : precalc[tra].dtlist){
						L += lognormal_probability(dt,meanf,sdf*sdf);
						L -= lognormal_probability(dt,meani,sdi*sdi);
					}
				}
			}
			break;
		}
	}
	
	return L;
}
