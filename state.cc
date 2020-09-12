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
	
	popw.resize(data.nardp);                                        // Used for event based changes	
	lambda.resize(data.nsettardp);
}

/// This simulates from the model and generates an event list
void State::simulate_compartmental_transitions(unsigned int i, unsigned int c, double t)
{
	vector <Event> &evlist = indev[i];

	auto timep = 0u; while(timep < model.ntime_period && t > model.time_period[timep].tend) timep++;

	Event ev;
	ev.ind = i; ev.timep = timep; 
		
	auto a = data.democatpos[data.ind[i].dp][0];
		
	unsigned int tra;
	if(c == 0){
		evlist.clear();	
		tra = 0;
		ev.trans = tra; ev.ind = i; ev.t = t; 
		evlist.push_back(ev);
	
		c = trans[tra].to;
	}
	
	double dt, z;
	do{
		auto kmax = comp[c].trans.size();
		if(kmax == 0) break;
		
		if(kmax == 1) tra = comp[c].trans[0];
		else{
			z = ran(); auto k = 0u; while(k < kmax && z > comptransprob[c].probsum[a][k]) k++;
			if(k == kmax) emsgEC("Model",2);
			tra = comp[c].trans[k];
		}
		
		switch(trans[tra].type){
		case EXP_DIST:
			dt = -log(ran())*paramval[trans[tra].param_mean];
			break;
		
		case GAMMA_DIST:
			{
				auto mean = paramval[trans[tra].param_mean]; auto sd = paramval[trans[tra].param_cv]*mean;
				dt = gammasamp(mean*mean/(sd*sd),mean/(sd*sd));
			}
			break;
			
		case LOGNORM_DIST:
			{
				auto mean_ns = paramval[trans[tra].param_mean], cv_ns = paramval[trans[tra].param_cv];
				auto sd = sqrt(log((1+cv_ns*cv_ns))), mean = log(mean_ns) - sd*sd/2;
				dt = exp(normal(mean,sd));
			}
			break;
			
		default: emsgEC("Model",3); break;
		}

		if(dt < TINY) dt = TINY;
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
	}while(1 == 1);
}

/// Calculates the latent process likelihood 
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
				if(lambda[w] < 0) emsgEC("Chain",46);
				
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
			if(std::isnan(L)) emsgEC("Chain",47);
			popw[w]--;
			n++;
			
			Lev += lambda[w]*(tmax-t);
		}
		t = tmax;
	} 
}

/// Checks the likelihood and prior ire correct 
void State::check_Lev_and_Pr()
{
	auto Levst = Lev;
	set_process_likelihood();
	double dd = Levst - Lev; if(dd*dd > TINY) emsgEC("State",49);
	
	dd = Pr - model.prior(paramval); if(sqrt(dd*dd) > TINY) emsgEC("State",51);
}
	
/// Clears the state
void State::clear()
{
	for(const auto& i : infev) indev[i.ind].clear();
	
	infev.clear();
	transev.clear(); transev.resize(details.ndivision);
}

/// Copies from another state
void State::copy_state(const State &from)
{
	paramval = from.paramval;
	L = from.L;
	EF = from.EF;
	Pr = from.Pr;
	
	for(const auto& iev : infev) indev[iev.ind].clear();
	for(const auto& iev : from.infev) indev[iev.ind] = from.indev[iev.ind];
	//indev = from.indev;
		
	infev = from.infev;
	transev = from.transev;
	Qmap = from.Qmap;
		
	disc_spline = from.disc_spline;
	susceptibility = from.susceptibility;
	areafactor = from.areafactor;
	comptransprob = from.comptransprob;
}

/// Checks quanties in the state are correct
void State::check() const
{
	for(auto j = 0u; j < infev.size(); j++){ // Checks order
		auto i = infev[j].ind, e = infev[j].e;
		if(i >= indev.size()) emsgEC("Chain",16);
		if(e >= indev[i].size()) emsgEC("Chain",17);
		if(j < infev.size()-1){
			if(indev[i][e].t > indev[infev[j+1].ind][infev[j+1].e].t) emsgEC("State",18);
		}
	}
	
	for(const auto& iev : indev){
		auto emax = iev.size();
		if(emax > 0){
			auto c = 0u; 
			double tt = 0;
			for(const auto& ev : iev){
				auto ttt = ev.t; if(ttt <= tt) emsgEC("Chain",19);
				auto tra = ev.trans;
				if(trans[tra].from != c) emsgEC("Chain",20);
				c = trans[tra].to; tt = ttt;
			}
			if(comp[c].trans.size() != 0) emsgEC("Chain",21);
			
			for(auto timep = 0u; timep < model.ntime_period; timep++){
				tt = model.time_period[timep].tend;
				if(tt > iev[0].t && tt < iev[emax-1].t){
					unsigned int e;
					for(e = 0u; e < emax; e++) if(iev[e].t == tt) break;
					if(timep <  model.ntime_period-1){
						if(e == emax) emsgEC("Chain",22);
					}
					else{
						if(e != emax) emsgEC("Chain",23);
					}
				}
			}
		}
	}
	
	for(auto j = 0u; j < infev.size(); j++){
		auto i = infev[j].ind, e = infev[j].e;
		if(i >= indev.size()) emsgEC("Chain",37);
		if(e >= indev[i].size()) emsgEC("Chain",38);
		if(j < infev.size()-1){
			if(indev[i][e].t > indev[infev[j+1].ind][infev[j+1].e].t) emsgEC("Chain",39);
		}
	}

	vector < vector <int> > done;
	done.resize(indev.size());
	auto num = 0u;
	for(auto i = 0u; i < indev.size(); i++){
		if(indev[i].size() > 0){
			num++;
			done[i].resize(indev[i].size());
			for(auto e = 0u; e < indev[i].size(); e++) done[i][e] = 0;
		}
	}
	if(num != infev.size()) emsgEC("Chain",40);

	for(auto sett = 0u; sett < details.ndivision; sett++){
		for(auto j = 0u; j < transev[sett].size(); j++){
			auto i = transev[sett][j].ind, e = transev[sett][j].e;
			if(e >= indev[i].size()) emsgEC("Chain",41);
			
			auto se = (unsigned int)(details.ndivision*indev[i][e].t/details.period); 
			if(se != sett) emsgEC("Chain",42);
			if(done[i][e] != 0) emsgEC("Chain",43);
			done[i][e] = 1;
		}
	}
	
	for(auto i = 0u; i < indev.size(); i++){
		for(auto e = 0u; e < indev[i].size(); e++){
			if(indev[i][e].t < details.period){
				if(done[i][e] != 1) emsgEC("Chain",44);
			}
			else{
				if(done[i][e] != 0) emsgEC("Chain",45);
			}
		}
	}
	
	double dd;
	dd = L - obsmodel.observation_likelihood(transev,indev); if(sqrt(dd*dd) > TINY) emsgEC("State",1);
	dd = Pr - model.prior(paramval); if(sqrt(dd*dd) > TINY) emsgEC("State",2);
}

/// Sets up the state with a specifies set of parameters
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

/// Sets the initial observation likelihood and prior
void State::set_L_and_Pr()
{
	L = obsmodel.observation_likelihood(transev,indev);
	Pr = model.prior(paramval);
}

/// Gets the time of an infection event
double State::get_infection_time(unsigned int n) const
{
	if(n == infev.size()) return LARGE;
	return indev[infev[n].ind][infev[n].e].t;
}

/// Adds an individual event sequence
void State::add_indev(unsigned int i)
{
	unsigned int e, emax, se;
	EventRef evref;
	
	emax = indev[i].size();
	if(emax == 0) return;
	
	evref.ind = i; evref.e = 0;
	infev.push_back(evref);
	for(e = 0; e < emax; e++){
		evref.e = e;
		se = (unsigned int)(details.ndivision*indev[i][e].t/details.period); 
		if(se < details.ndivision) transev[se].push_back(evref);
	}
}

/// For a given time sett, sets Qmap by using the Qmao from another state and the difference given by dQmap
void State::set_Qmap_using_dQ(unsigned int sett, const State &state, const vector <double> &dQmap)
{
	for(auto v = 0u; v < data.narage; v++){
		double val = state.Qmap[sett][v] + dQmap[v];
		if(val < -TINY){ cout << val << "val\n"; emsgEC("Chain",1);}
		if(val < 0) val = 0;	
		Qmap[sett][v] = val;
	}
}

/// Sets Qmap
void State::set_Qmap(unsigned int check)
{
	vector <double> Qma(data.narage);
	
 	for(auto& Qm : Qma) Qm = 0;

	auto nage = data.nage;
	for(auto sett = 0u; sett < details.ndivision; sett++){
		for(auto v = 0u; v < data.narage; v++){
			auto val = Qma[v];
			if(check == 1){
				if(val < -TINY) emsgEC("Chain",2);
				if(val < Qmap[sett][v]-TINY || val > Qmap[sett][v]+TINY) emsgEC("Chain",3);
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
						
						auto qt = data.Q[q].Qtenref;
						auto kmax = data.genQ.Qten[qt].ntof[v];
						auto& cref = data.genQ.Qten[qt].tof[v];
						auto& valref = data.genQ.Qten[qt].valf[v];
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

struct EventRefT {                
	unsigned int ind;                   
	unsigned int e;	              
	double t;	 
};

static bool compEventRefT(EventRefT lhs, EventRefT rhs)
{
	return lhs.t < rhs.t;
};

/// Time orders x
void State::sort_infev()
{
	vector <EventRefT> xt;
	for(const auto& iev : infev){
		EventRefT evreft;
		evreft.ind = iev.ind; evreft.e = iev.e;
		if(indev[iev.ind].size() == 0) emsgEC("Chain",54);
		
		evreft.t = indev[iev.ind][iev.e].t;
		xt.push_back(evreft);	
	}
	sort(xt.begin(),xt.end(),compEventRefT);

	for(auto i = 0u; i < infev.size(); i++){
		infev[i].ind = xt[i].ind;	infev[i].e = xt[i].e;
	}
}

/// Initialises the chain based on a particle (used for abcmbp)
void State::initialise_from_particle(const Particle &part)
{
	paramval = part.paramval;
	
	for(const auto& iev : infev) indev[iev.ind].clear();   // Removes the existing initial sequence 
	infev.clear();
	
	transev.clear(); transev.resize(details.ndivision); 
	
	vector <int> indlist;
	for(const auto& ev : part.ev){
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
		
/// STANDARD PARAMETER PROPOSALS

/// This incorporates standard proposals which adds and removes events as well as changes parameters
void State::standard_parameter_prop(Jump &jump)
{	
	stand_param_betaphi_prop(jump);
	stand_param_area_prop(jump);
	stand_param_compparam_prop(jump);

	if(checkon == 1) check_Lev_and_Pr();
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
				paramval[th] += normal(0,jump.stand[th]);               // Makes a change to a parameter

				vector < vector <double> > disc_spline_prop;

				double al, Lev_prop=Lev, Pr_prop=Pr;
				if(paramval[th] < param[th].min || paramval[th] > param[th].max) al = 0;
				else{
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
		if(std::isnan(Lev)) emsgEC("Chain",52);
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
				paramval[th] += normal(0,jump.stand[th]);               // Makes a change to a parameter

				double al, Lev_prop=Lev, Pr_prop=Pr;
				vector <double> areafactor_prop;
				if(paramval[th] < param[th].min || paramval[th] > param[th].max) al = 0;
				else{
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
	if(std::isnan(Lev)) emsgEC("Chain",53);
	
	return Lev;
}
					
/// Makes proposal to compartmental parameters
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
		
			paramval[th] += normal(0,jump.stand[th]);               // Makes a change to a parameter
			
			auto Lp_prob = Li_prob;
			auto dL = 0.0;
			auto Pr_prop = Pr;
			auto flag=0u;
	
			vector <CompTransProb> comptransprobp;
			if(paramval[th] < param[th].min || paramval[th] > param[th].max) flag = 0;
			else{	
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
	
	if(checkon == 1){
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
			
/// Calculates the likelihood for the timings of the transitions
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
				for(auto dt : precalc[tra].dtlist) L += gammaprob(dt,mean*mean/(sd*sd),mean/(sd*sd));
			}
			break;
			
		case LOGNORM_DIST:
			{
				auto mean_ns = paramv[trans[tra].param_mean], cv_ns = paramv[trans[tra].param_cv];
				auto sd = sqrt(log(1+cv_ns*cv_ns)), mean = log(mean_ns) - sd*sd/2;
				for(auto dt : precalc[tra].dtlist) L += lognormprob(dt,mean,sd*sd);
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
						L += gammaprob(dt,meanf*meanf/(sdf*sdf),meanf/(sdf*sdf));
						L -= gammaprob(dt,meani*meani/(sdi*sdi),meani/(sdi*sdi));
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
						L += lognormprob(dt,meanf,sdf*sdf);
						L -= lognormprob(dt,meani,sdi*sdi);
					}
				}
			}
			break;
		}
	}
	
	return L;
}
