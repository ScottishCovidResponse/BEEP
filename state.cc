#include <cmath>   
#include <iostream>

using namespace std;

#include "state.hh"
#include "timers.hh"

State::State(const Details &details, const DATA &data, const MODEL &model, const Obsmodel &obsmodel) : comp(model.comp), trans(model.trans), param(model.param), details(details), data(data), model(model), obsmodel(obsmodel)
{
	disc_spline.resize(model.spline.size());

	trev.resize(details.nsettime); 
	indev.resize(data.popsize);

	Qmap.resize(details.nsettime);
	for(auto sett = 0u; sett < details.nsettime; sett++){
		Qmap[sett].resize(data.narage); for(auto v = 0u; v < data.narage; v++) Qmap[sett][v] = 0;
	}
	
	popw.resize(data.nardp);                                        // Used for event based changes	
	lam.resize(data.nsettardp);
}

/// This simulates from the model and generates an event list
void State::simulate_compartmental_transitions(unsigned int i, unsigned int c, double t)
{
	vector <FEV> &evlist = indev[i];

	auto timep = 0u; while(timep < model.ntimeperiod && t > model.timeperiod[timep].tend) timep++;

	FEV ev;
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
			z = ran(); auto k = 0u; while(k < kmax && z > comptrans[c].probsum[a][k]) k++;
			if(k == kmax) emsgEC("Model",2);
			tra = comp[c].trans[k];
		}
		
		switch(trans[tra].type){
		case exp_dist:
			dt = -log(ran())*paramval[trans[tra].param_mean];
			break;
		
		case gamma_dist:
			{
				auto mean = paramval[trans[tra].param_mean]; auto sd = paramval[trans[tra].param_cv]*mean;
				dt = gammasamp(mean*mean/(sd*sd),mean/(sd*sd));
			}
			break;
			
		case lognorm_dist:
			{
				auto mean_ns = paramval[trans[tra].param_mean], cv_ns = paramval[trans[tra].param_cv];
				auto sd = sqrt(log((1+cv_ns*cv_ns))), mean = log(mean_ns) - sd*sd/2;
				dt = exp(normal(mean,sd));
			}
			break;
			
		default: emsgEC("MODEL",3); break;
		}

		if(dt < tiny) dt = tiny;
		t += dt;
		
		while(timep < model.ntimeperiod-1 && model.timeperiod[timep].tend < t){    // Adds in changes in time period
			ev.trans = comp[c].transtimep; ev.t = model.timeperiod[timep].tend;
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
void State::set_likelihood()
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
	for(auto sett = 0u; sett < details.nsettime; sett++){
		auto phi = disc_spline[model.phispline_ref][sett]; 
		auto beta = disc_spline[model.betaspline_ref][sett];
	
		auto tmax = details.settime[sett+1];
		
		for(auto c = 0u; c < data.narea; c++){
			auto fac = beta*areafac[c];
			for(auto dp = 0u; dp < data.ndemocatpos; dp++){
				auto w = c*data.ndemocatpos + dp;
				auto v = c*data.nage + data.democatpos[dp][0];
				lam[w] = sus[dp]*(fac*Qmap[sett][v] + phi);		
				if(lam[w] < 0) emsgEC("Chain",46);
				
				Lev -= lam[w]*popw[w]*(tmax-t);
			}
		}
		
		while(n < x.size()){
			auto i = x[n].ind;
			auto ev = indev[i][x[n].e];
			auto tt = ev.t;
			if(tt >= tmax) break;
	
			t = tt;
			
			auto c = data.ind[i].area;
			auto w = c*data.ndemocatpos + data.ind[i].dp;
			Lev += log(lam[w]);
			if(std::isnan(L)) emsgEC("Chain",47);
			popw[w]--;
			n++;
			
			Lev += lam[w]*(tmax-t);
		}
		t = tmax;
	} 
}

/// Checks the likelihood and prior ire correct 
void State::check_LevPr() 
{
	auto Levst = Lev;
	set_likelihood();
	double dd = Levst - Lev; if(dd*dd > tiny) emsgEC("State",49);
	
	dd = Pr - model.prior(paramval); if(sqrt(dd*dd) > tiny) emsgEC("State",51);
}
	
/// Clears the state
void State::clear()
{
	for(const auto& i : x) indev[i.ind].clear();
	
	x.clear();
	trev.clear(); trev.resize(details.nsettime);
}

/// Copies from another state
void State::copy(const State &from)
{
	paramval = from.paramval;
	L = from.L;
	EF = from.EF;
	Pr = from.Pr;
	
	for(const auto& xx : x) indev[xx.ind].clear();
	for(const auto& xx : from.x) indev[xx.ind] = from.indev[xx.ind];
	//indev = from.indev;
		
	x = from.x;
	trev = from.trev;
	Qmap = from.Qmap;
		
	disc_spline = from.disc_spline;
	sus = from.sus;
	areafac = from.areafac;
	comptrans = from.comptrans;
}

/// Checks quanties in the state are correct
void State::check() const
{
	for(auto j = 0u; j < x.size(); j++){ // Checks order
		auto i = x[j].ind, e = x[j].e;
		if(i >= indev.size()) emsgEC("Chain",16);
		if(e >= indev[i].size()) emsgEC("Chain",17);
		if(j < x.size()-1){
			if(indev[i][e].t > indev[x[j+1].ind][x[j+1].e].t) emsgEC("State",18);
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
			
			for(auto timep = 0u; timep < model.ntimeperiod; timep++){
				tt = model.timeperiod[timep].tend;
				if(tt > iev[0].t && tt < iev[emax-1].t){
					unsigned int e;
					for(e = 0u; e < emax; e++) if(iev[e].t == tt) break;
					if(timep <  model.ntimeperiod-1){
						if(e == emax) emsgEC("Chain",22);
					}
					else{
						if(e != emax) emsgEC("Chain",23);
					}
				}
			}
		}
	}
	
	for(auto j = 0u; j < x.size(); j++){
		auto i = x[j].ind, e = x[j].e;
		if(i >= indev.size()) emsgEC("Chain",37);
		if(e >= indev[i].size()) emsgEC("Chain",38);
		if(j < x.size()-1){
			if(indev[i][e].t > indev[x[j+1].ind][x[j+1].e].t) emsgEC("Chain",39);
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
	if(num != x.size()) emsgEC("Chain",40);

	for(auto sett = 0u; sett < details.nsettime; sett++){
		for(auto j = 0u; j < trev[sett].size(); j++){
			auto i = trev[sett][j].ind, e = trev[sett][j].e;
			if(e >= indev[i].size()) emsgEC("Chain",41);
			
			auto se = (unsigned int)(details.nsettime*indev[i][e].t/details.period); 
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
	dd = L - obsmodel.Lobs(trev,indev); if(sqrt(dd*dd) > tiny) emsgEC("State",1);
	dd = Pr - model.prior(paramval); if(sqrt(dd*dd) > tiny) emsgEC("State",2);
}

/// Sets up the state with a specifies set of parameters
Status State::set_param(const vector <double> &paramv)
{
	paramval = paramv;
	if(model.create_comptrans(comptrans,paramval) == 1) return fail;
	for(auto sp = 0u; sp < model.spline.size(); sp++) disc_spline[sp] = model.create_disc_spline(sp,paramval);
	sus = model.create_sus(paramval);    
	areafac = model.create_areafac(paramval);    
	return success;
}

/// Sets the values of phi and beta (this is used to speed up computation)
void State::set_betaphi(unsigned int sett)
{	
	phi = disc_spline[model.phispline_ref][sett]; 
	beta = disc_spline[model.betaspline_ref][sett]; 
}

/// Sets the initial observation likelihood and prior
void State::set_LPr()
{
	L = obsmodel.Lobs(trev,indev);
	Pr = model.prior(paramval);
}

/// Gets the time of an infection event
double State::get_infection_time(unsigned int n) const
{
	if(n == x.size()) return large;
	return indev[x[n].ind][x[n].e].t;
}

/// Adds an individual event sequence
void State::add_indev(unsigned int i)
{
	unsigned int e, emax, se;
	EVREF evref;
	
	emax = indev[i].size();
	if(emax == 0) return;
	
	evref.ind = i; evref.e = 0;
	x.push_back(evref);
	for(e = 0; e < emax; e++){
		evref.e = e;
		se = (unsigned int)(details.nsettime*indev[i][e].t/details.period); 
		if(se < details.nsettime) trev[se].push_back(evref);
	}
}

/// For a given time sett, sets Qmap by using the Qmao from another state and the difference given by dQmap
void State::set_Qmap_using_dQ(unsigned int sett, const State &state, const vector <double> &dQmap)
{
	for(auto v = 0u; v < data.narage; v++){
		double val = state.Qmap[sett][v] + dQmap[v];
		if(val < -tiny){ cout << val << "val\n"; emsgEC("Chain",1);}
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
	for(auto sett = 0u; sett < details.nsettime; sett++){
		for(auto v = 0u; v < data.narage; v++){
			auto val = Qma[v];
			if(check == 1){
				if(val < -tiny) emsgEC("Chain",2);
				if(val < Qmap[sett][v]-tiny || val > Qmap[sett][v]+tiny) emsgEC("Chain",3);
			}
			if(val < 0){ val = 0; Qma[v] = 0;}	
			
			Qmap[sett][v] = val;
		}
		
		for(const auto& tre : trev[sett]){
			auto i = tre.ind;
			FEV fev = indev[i][tre.e];

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

struct EVREFT {                
	unsigned int ind;                   
	unsigned int e;	              
	double t;	 
};

static bool compEVREFT(EVREFT lhs, EVREFT rhs)
{
	return lhs.t < rhs.t;
};

/// Time orders x
void State::sort_x()
{
	vector <EVREFT> xt;
	for(const auto& xx : x){
		EVREFT evreft;
		evreft.ind = xx.ind; evreft.e = xx.e;
		if(indev[xx.ind].size() == 0) emsgEC("Chain",54);
		
		evreft.t = indev[xx.ind][xx.e].t;
		xt.push_back(evreft);	
	}
	sort(xt.begin(),xt.end(),compEVREFT);

	for(auto i = 0u; i < x.size(); i++){
		x[i].ind = xt[i].ind;	x[i].e = xt[i].e;
	}
}

/// Initialises the chain based on a particle (used for abcmbp)
void State::initialise_from_particle(const Particle &part)
{
	paramval = part.paramval;
	
	for(const auto& xx : x) indev[xx.ind].clear();   // Removes the existing initial sequence 
	x.clear();
	
	trev.clear(); trev.resize(details.nsettime); 
	
	vector <int> indlist;
	for(const auto& ev : part.ev){
		int i = ev.ind;
		if(indev[i].size() == 0) indlist.push_back(i);
		indev[i].push_back(ev);
	}	
	
	for(auto i : indlist) add_indev(i);
	
	sort_x();
	
	set_Qmap(0);

	EF = part.EF;
	Pr = model.prior(paramval);
	
	if(EF != obsmodel.Lobs(trev,indev)) emsg("Observation does not agree");
}
		
/// STANDARD PARAMETER PROPOSALS

/// This incorporates standard proposals which adds and removes events as well as changes parameters
void State::standard_parameter_prop(Jump &jump)
{
	timers.timestandard -= clock();
	
	timers.timembptemp -= clock();
	set_likelihood();
	timers.timembptemp += clock();
	
	timers.timeparam -= clock();
	stand_param_betaphi_prop(jump);
	stand_param_area_prop(jump);
	stand_param_compparam_prop(jump);
	
	timers.timeparam += clock();
		
	if(checkon == 1) check_LevPr();

	timers.timestandard += clock();
}

struct LCONT {                
	unsigned int w;       
	unsigned int num;       	
	double betafac;	              
	double phifac;	 
};

/// Makes proposal to beta and phi
void State::stand_param_betaphi_prop(Jump &jump)
{	
	timers.timebetaphiinit -= clock();
	
	vector <unsigned int> map;
	map.resize(data.nardp); for(auto& ma : map) ma = 0;
	
	for(auto c = 0u; c < data.narea; c++){
		for(auto dp = 0u; dp < data.ndemocatpos; dp++){
			auto w = c*data.ndemocatpos + dp;
			popw[w] = data.area[c].ind[dp].size();
		}
	}		
	
	vector <double> betafac, phifac;
	betafac.resize(details.nsettime); phifac.resize(details.nsettime);
	
	auto t = 0.0; 
	auto n = 0u;
	
	vector<	vector <LCONT> > lc;
	for(auto sett = 0u; sett < details.nsettime; sett++){
		auto tmax = details.settime[sett+1];
		vector <LCONT> lcontlist;
		
		auto betasum = 0.0, phisum = 0.0;
		for(auto c = 0u; c < data.narea; c++){
			auto fac = areafac[c];
			for(auto dp = 0u; dp < data.ndemocatpos; dp++){
				auto w = c*data.ndemocatpos + dp;
				auto v = c*data.nage + data.democatpos[dp][0];	
				betasum -= fac*sus[dp]*Qmap[sett][v]*popw[w]*(tmax-t);
				phisum -= sus[dp]*popw[w]*(tmax-t);
			}
		}
		
		while(n < x.size()){
			auto i = x[n].ind;
			FEV ev = indev[i][x[n].e];
			auto tt = ev.t;
			if(tt >= tmax) break;
	
			t = tt;
			
			auto c = data.ind[i].area;
			auto fac = areafac[c];
			
			auto dp = data.ind[i].dp;
			auto w = c*data.ndemocatpos + dp;
			auto v = c*data.nage + data.democatpos[dp][0]; 
			
			if(map[w] == 0){
				LCONT lcont;
				lcont.w = w; lcont.betafac = fac*sus[dp]*Qmap[sett][v]; lcont.phifac = sus[dp]; lcont.num = UNSET;
				lcontlist.push_back(lcont);
			}
			map[w]++;
		
			betasum += fac*sus[dp]*Qmap[sett][v]*(tmax-t);
			phisum += sus[dp]*(tmax-t);
			popw[w]--;
			n++;
		}
		
		betafac[sett] = betasum; phifac[sett] = phisum;
					
		for(auto& lcont : lcontlist){
			auto w = lcont.w;
			lcont.num = map[w];
			map[w] = 0;
		}
					
		lc.push_back(lcontlist);
		
		t = tmax;
	} 
	
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

					Lev_prop = 0;
					for(auto sett = 0u; sett < details.nsettime; sett++){
						auto phi = disc_spline_prop[model.phispline_ref][sett]; 
						auto beta = disc_spline_prop[model.betaspline_ref][sett];
	
						
						Lev_prop += betafac[sett]*beta + phifac[sett]*phi;
						
						for(const auto& l : lc[sett]) Lev_prop += l.num*log(l.betafac*beta + l.phifac*phi);
						if(std::isnan(Lev_prop)) emsgEC("Chain",52);
					}
					
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

/// Makes proposal to change factors affecting transmission rates in areas
void State::stand_param_area_prop(Jump &jump)
{	
	timers.timecovarinit -= clock();

	vector <double> lamareafac, lamphifac;
	lamareafac.resize(data.nardp);
	lamphifac.resize(data.nardp);
	
	vector < vector <double> >	mult, add;
	mult.resize(data.narea);
	add.resize(data.narea);
	
	vector <double> areasum;
	areasum.resize(data.narea);
	
	for(auto c = 0u; c < data.narea; c++){
		areasum[c] = 0;
		for(auto dp = 0u; dp < data.ndemocatpos; dp++){
			auto w = c*data.ndemocatpos + dp;
			popw[w] = data.area[c].ind[dp].size();
		}
	}		
	
	auto L0 = 0.0, t = 0.0;
	auto n = 0u;
	for(auto sett = 0u; sett < details.nsettime; sett++){
		auto phi = disc_spline[model.phispline_ref][sett]; 
		auto beta = disc_spline[model.betaspline_ref][sett];
	
		auto tmax = details.settime[sett+1];
		
		for(auto c = 0u; c < data.narea; c++){
			for(auto dp = 0u; dp < data.ndemocatpos; dp++){
				auto w = c*data.ndemocatpos + dp;
				auto v = c*data.nage + data.democatpos[dp][0];
				lamareafac[w] = sus[dp]*beta*Qmap[sett][v];
				lamphifac[w] = sus[dp]*phi;
				
				areasum[c] -= lamareafac[w]*popw[w]*(tmax-t);
				L0 -= lamphifac[w] *popw[w]*(tmax-t);
			}
		}
		
		while(n < x.size()){
			auto i = x[n].ind;
			FEV ev = indev[i][x[n].e];
			auto tt = ev.t;
			if(tt >= tmax) break;
	
			t = tt;
			
			auto c = data.ind[i].area;
			auto w = c*data.ndemocatpos + data.ind[i].dp;
			
			mult[c].push_back(lamareafac[w]);
			add[c].push_back(lamphifac[w]);
		
			popw[w]--;
			n++;
			
 			areasum[c] += lamareafac[w]*(tmax-t);
			L0 += lamphifac[w]*(tmax-t);
		}
		
		t = tmax;
	} 
	
	timers.timecovarinit += clock();
		
	timers.timecovar -= clock();
	
	vector <unsigned int> thlist;
	for(auto th : model.covar_param) thlist.push_back(th); 
	
	if(model.regioneffect > 0){
		for(auto th : model.regioneff_param) thlist.push_back(th); 
	}
	
	if(model.regioneffect == 1) thlist.push_back(model.sigma_param);
		
	unsigned int loopmax = 12/thlist.size(); if(loopmax == 0) loopmax = 1;
	
	for(auto loop = 0u; loop < loopmax; loop++){
		for(auto th : thlist){
			if(param[th].min != param[th].max){
				auto param_store = paramval[th];	
				paramval[th] += normal(0,jump.stand[th]);               // Makes a change to a parameter

				double al, Lev_prop=Lev, Pr_prop=Pr;
				vector <double> areafac_prop;
				if(paramval[th] < param[th].min || paramval[th] > param[th].max) al = 0;
				else{
					areafac_prop = model.create_areafac(paramval); 
				
					Lev_prop = L0;
					for(auto c = 0u; c < data.narea; c++){
						auto fac = areafac_prop[c];
						Lev_prop += areasum[c]*fac;
						auto kmax = mult[c].size();
						for(auto k = 0u; k < kmax; k++) Lev_prop += log(mult[c][k]*fac + add[c][k]);
					}
					if(std::isnan(Lev_prop)) emsgEC("Chain",53);
				
					Pr_prop = model.prior(paramval);
					al = exp(Pr_prop-Pr + Lev_prop-Lev);
				}

				//ntrstand[th]++;
				if(ran() < al){
					Lev = Lev_prop;
					Pr = Pr_prop;
					areafac = areafac_prop;
					
					jump.stand_accept(th);
					//nacstand[th]++;
					//if(samp < burnin){ if(samp < 50) paramjumpstand[th] *= 1.05; else paramjumpstand[th] *= 1.01;}
				}
				else{
					paramval[th] = param_store;
					jump.stand_reject(th);
					//if(samp < burnin){ if(samp < 50) paramjumpstand[th] *= 0.975; else paramjumpstand[th] *= 0.995;}
				}
			}
		}
	}
	timers.timecovar += clock();
}

/// Makes proposal to compartmental paramters
void State::stand_param_compparam_prop(Jump &jump)
{	
	timers.timecompparam -= clock();

	vector <TransInfo> transinfo(trans.size());
	
	for(auto tra = 0u; tra < trans.size(); tra++){
		if(trans[tra].istimep == 0){
			transinfo[tra].num.resize(data.nage); for(auto& num : transinfo[tra].num) num = 0;
		}
		transinfo[tra].numvisittot = 0; transinfo[tra].dtsum = 0; transinfo[tra].dtlist.clear();
	}

	for(const auto& xx : x){                  // Extracts values based on the event sequence
		auto i = xx.ind;
		auto dp = data.ind[i].dp;
		auto a = data.democatpos[dp][0];
				
		auto t = 0.0;
		for(const auto& ev : indev[i]){
			auto tra = ev.trans;
			if(trans[tra].istimep == 0){
				transinfo[tra].num[a]++;
	
				auto dt = ev.t-t;
				if(dt == 0) emsgEC("Model",10);
				t += dt;
				
				transinfo[tra].numvisittot++;
				switch(trans[tra].type){
				case exp_dist: transinfo[tra].dtsum += dt; break;
				case gamma_dist: case lognorm_dist: transinfo[tra].dtlist.push_back(dt); break;
				case infection_dist: break;
				default: emsgEC("model",99); break;
				}
			}
		}
	}

	auto Li_dt = likelihood_dt(transinfo,paramval);
	auto Li_prob = likelihood_prob(transinfo,comptrans);

	for(auto th = 0u; th < param.size(); th++){
		if((param[th].type == distval_paramtype || param[th].type == branchprob_paramtype) && param[th].min != param[th].max){
			vector <double> param_store = paramval;	
		
			paramval[th] += normal(0,jump.stand[th]);               // Makes a change to a parameter
			
			auto Lp_prob = Li_prob;
			auto dL = 0.0;
			auto Pr_prop = Pr;
			auto flag=0u;
	
			vector <CompTrans> comptransp;
			if(paramval[th] < param[th].min || paramval[th] > param[th].max) flag = 0;
			else{	
				if(param[th].type == branchprob_paramtype){
					if(model.create_comptrans(comptransp,paramval) == 1) Lp_prob = -large;
					else Lp_prob = likelihood_prob(transinfo,comptransp);
				}
				else Lp_prob = Li_prob;
				
				dL = dlikelihood_dt(transinfo,param_store,paramval);
				Pr_prop = model.prior(paramval);
				
				auto al = exp(dL+ Lp_prob - Li_prob + Pr_prop - Pr);
				if(ran() < al) flag = 1; else flag = 0;
			}
			
			if(flag == 1){
				Li_dt += dL;
				Li_prob = Lp_prob;
				Pr = Pr_prop;
				if(param[th].type == branchprob_paramtype) comptrans = comptransp;
				
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
		dd = likelihood_dt(transinfo,paramval)-Li_dt; if(dd*dd > tiny) emsgEC("Model",13);
		dd = likelihood_prob(transinfo,comptrans)-Li_prob; if(dd*dd > tiny) emsgEC("Model",14);
	}
	
	timers.timecompparam += clock();
}

/// Calculates the likelihood relating to branching probabilities
double State::likelihood_prob(vector <TransInfo> &transinfo, vector <CompTrans> &comptrans) const
{
	auto L = 0.0;
	for(auto c = 0u; c < comp.size(); c++){	
		auto kmax = comp[c].trans.size();
		if(kmax > 1){
			for(auto k = 0u; k < kmax; k++){
				for(auto a = 0u; a < data.nage; a++){
					auto num = transinfo[comp[c].trans[k]].num[a];
					if(num > 0) L += num*log(comptrans[c].prob[a][k]);
				}
			}
		}
	}
	return L;
}
			
/// Calculates the likelihood for the timings of the transitions
double State::likelihood_dt(vector <TransInfo> &transinfo, vector <double> &paramv) const
{
	auto L = 0.0;
	for(auto tra = 0u; tra < trans.size(); tra++){
		switch(trans[tra].type){
		case exp_dist:
			{
				auto r = 1.0/paramv[trans[tra].param_mean];
				L += transinfo[tra].numvisittot*log(r) - r*transinfo[tra].dtsum;
			}
			break;
		
		case gamma_dist:
			{
				auto mean = paramv[trans[tra].param_mean], sd = paramv[trans[tra].param_cv]*mean;
				for(auto dt : transinfo[tra].dtlist) L += gammaprob(dt,mean*mean/(sd*sd),mean/(sd*sd));
			}
			break;
			
		case lognorm_dist:
			{
				auto mean_ns = paramv[trans[tra].param_mean], cv_ns = paramv[trans[tra].param_cv];
				auto sd = sqrt(log(1+cv_ns*cv_ns)), mean = log(mean_ns) - sd*sd/2;
				for(auto dt : transinfo[tra].dtlist) L += lognormprob(dt,mean,sd*sd);
			}
			break;
		}
	}
	
	return L;
}

/// Calculates the change in likelihood for a given change in parameters
double State::dlikelihood_dt(vector <TransInfo> &transinfo, vector <double> &paramvi, vector <double> &paramvf) const
{
	auto L = 0.0;
	for(auto tra = 0u; tra < trans.size(); tra++){
		switch(trans[tra].type){
		case exp_dist:
			{
				auto ri = 1.0/paramvi[trans[tra].param_mean], rf = 1.0/paramvf[trans[tra].param_mean];
				if(ri != rf){
					L += transinfo[tra].numvisittot*log(rf) - rf*transinfo[tra].dtsum;
					L -= transinfo[tra].numvisittot*log(ri) - ri*transinfo[tra].dtsum;
				}
			}
			break;
		
		case gamma_dist:
			{
				auto meani = paramvi[trans[tra].param_mean], sdi = paramvi[trans[tra].param_cv]*meani;
				auto meanf = paramvf[trans[tra].param_mean], sdf = paramvf[trans[tra].param_cv]*meanf;
		
				if(meani != meanf || sdi != sdf){
					for(auto dt : transinfo[tra].dtlist){
						L += gammaprob(dt,meanf*meanf/(sdf*sdf),meanf/(sdf*sdf));
						L -= gammaprob(dt,meani*meani/(sdi*sdi),meani/(sdi*sdi));
					}
				}				
			}
			break;
			
		case lognorm_dist:
			{
				auto mean_nsi = paramvi[trans[tra].param_mean], cv_nsi = paramvi[trans[tra].param_cv];
				auto mean_nsf = paramvf[trans[tra].param_mean], cv_nsf = paramvf[trans[tra].param_cv];
				if(mean_nsi != mean_nsf || cv_nsi != cv_nsf){
					auto sdi = sqrt(log(1+cv_nsi*cv_nsi)), meani = log(mean_nsi) - sdi*sdi/2;
					auto sdf = sqrt(log(1+cv_nsf*cv_nsf)), meanf = log(mean_nsf) - sdf*sdf/2;
		
					for(auto dt : transinfo[tra].dtlist){
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
