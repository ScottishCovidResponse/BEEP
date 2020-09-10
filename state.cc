#include <cmath>   
#include <iostream>

using namespace std;

#include "state.hh"
#include "timers.hh"

State::State(const Details &details, const DATA &data, const MODEL &model, const Obsmodel &obsmodel) : comp(model.comp), trans(model.trans), details(details), data(data), model(model), obsmodel(obsmodel)
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

/// Calculates the latent process likelihood 
void State::setlikelihood()
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
void State::checkLevPr() 
{
	auto Levst = Lev;
	setlikelihood();
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
Status State::setparam(const vector <double> &paramv)
{
	paramval = paramv;
	if(model.create_comptrans(comptrans,paramval) == 1) return fail;
	for(auto sp = 0u; sp < model.spline.size(); sp++) disc_spline[sp] = model.create_disc_spline(sp,paramval);
	sus = model.create_sus(paramval);    
	areafac = model.create_areafac(paramval);    
	return success;
}

/// Sets the values of phi and beta (this is used to speed up computation)
void State::setbetaphi(unsigned int sett)
{	
	phi = disc_spline[model.phispline_ref][sett]; 
	beta = disc_spline[model.betaspline_ref][sett]; 
}

/// Sets the initial observation likelihood and prior
void State::setLPr()
{
	L = obsmodel.Lobs(trev,indev);
	Pr = model.prior(paramval);
}

/// Gets an infection event
FEV State::getinfev(unsigned int n) const
{
	return indev[x[n].ind][x[n].e];
}

/// Adds an individual event sequence
void State::addindev(unsigned int i)
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
void State::setQmapUsingdQ(unsigned int sett, const State &state, const vector <double> &dQmap)
{
	for(auto v = 0u; v < data.narage; v++){
		double val = state.Qmap[sett][v] + dQmap[v];
		if(val < -tiny){ cout << val << "val\n"; emsgEC("Chain",1);}
		if(val < 0) val = 0;	
		Qmap[sett][v] = val;
	}
}

/// Sets Qmap
void State::setQmap(unsigned int check)
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
void State::sortx()
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
	
	for(auto i : indlist) addindev(i);
	
	sortx();
	
	setQmap(0);

	EF = part.EF;
	Pr = model.prior(paramval);
	
	if(EF != obsmodel.Lobs(trev,indev)) emsg("Observation does not agree");
}
		
/// STANDARD PROPOSALS


struct LCONT {                
	unsigned int w;       
	unsigned int num;       	
	double betafac;	              
	double phifac;	 
};

/// Makes proposal to beta and phi
void State::betaphi_prop(unsigned int samp, unsigned int burnin, vector <float> &paramjumpstand, vector <unsigned int> &ntrstand, 	vector <unsigned int> &nacstand)
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
			if(model.param[th].min != model.param[th].max){
				auto valst = paramval[th];	
				paramval[th] += normal(0,paramjumpstand[th]);               // Makes a change to a parameter

				vector < vector <double> > disc_spline_prop;

				double al, Lev_prop=Lev, Pr_prop=Pr;
				if(paramval[th] < model.param[th].min || paramval[th] > model.param[th].max) al = 0;
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
			
				ntrstand[th]++;
				if(ran() < al){
					Pr = Pr_prop;
					Lev = Lev_prop;
					disc_spline = disc_spline_prop;
					nacstand[th]++;
					if(samp < burnin){ if(samp < 50) paramjumpstand[th] *= 1.05; else paramjumpstand[th] *= 1.01;}
				}
				else{
					paramval[th] = valst;
					if(samp < burnin){ if(samp < 50) paramjumpstand[th] *= 0.975; else  paramjumpstand[th] *= 0.995;}
				}
			}
		}
	}
	timers.timebetaphi += clock();
}


/// Makes proposal to change factors affecting transmission rates in areas
void State::area_prop(unsigned int samp, unsigned int burnin, vector <float> &paramjumpstand, vector <unsigned int> &ntrstand, 	vector <unsigned int> &nacstand)
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
	
	auto num = model.covar_param.size(); if(model.regioneffect == 1) num += data.nregion;
	unsigned int loopmax = 12/num; if(loopmax == 0) loopmax = 1;
	
	for(auto loop = 0u; loop < loopmax; loop++){
		for(auto th : model.covar_param){ 
			if(model.param[th].min != model.param[th].max) area_prop2(samp,burnin,th,L0,areasum,mult,add,paramjumpstand,ntrstand,nacstand);
		}
		
		// ZZ This should be modified for when model.regioneffect == 2
		if(model.regioneffect == 1){
			for(auto th : model.regioneff_param){ 
				if(model.param[th].min != model.param[th].max) area_prop2(samp,burnin,th,L0,areasum,mult,add,paramjumpstand,ntrstand,nacstand);
			}
	
			auto th = model.sigma_param;
			if(model.param[th].min != model.param[th].max) area_prop2(samp,burnin,th,L0,areasum,mult,add,paramjumpstand,ntrstand,nacstand);
		}
	}
	timers.timecovar += clock();
}

void State::area_prop2(unsigned int samp, unsigned int burnin, unsigned int th, double L0, const vector <double> &areasum, const vector < vector <double> >&mult, const vector < vector <double> > &add, vector <float> &paramjumpstand, vector <unsigned int> &ntrstand, 	vector <unsigned int> &nacstand)
{
	auto valst = paramval[th];	
	paramval[th] += normal(0,paramjumpstand[th]);               // Makes a change to a parameter

	double al, Lev_prop=Lev, Pr_prop=Pr;
	vector <double> areafac_prop;
	if(paramval[th] < model.param[th].min || paramval[th] > model.param[th].max) al = 0;
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

	ntrstand[th]++;
	if(ran() < al){
		Lev = Lev_prop;
		Pr = Pr_prop;
		areafac = areafac_prop;
		
		nacstand[th]++;
		if(samp < burnin){ if(samp < 50) paramjumpstand[th] *= 1.05; else paramjumpstand[th] *= 1.01;}
	}
	else{
		paramval[th] = valst;
		if(samp < burnin){ if(samp < 50) paramjumpstand[th] *= 0.975; else paramjumpstand[th] *= 0.995;}
	}
}


	