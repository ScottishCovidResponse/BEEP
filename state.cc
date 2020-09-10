#include <cmath>   
#include <iostream>

using namespace std;

#include "state.hh"


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

/// Calculates the likelihood 
double State::likelihood()
{		
	for(auto c = 0u; c < data.narea; c++){
		for(auto dp = 0u; dp < data.ndemocatpos; dp++){
			auto w = c*data.ndemocatpos + dp;
			popw[w] = data.area[c].ind[dp].size();
		}
	}		
			
	auto L = 0.0;
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
				
				L -= lam[w]*popw[w]*(tmax-t);
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
			L += log(lam[w]);
			if(std::isnan(L)) emsgEC("Chain",47);
			popw[w]--;
			n++;
			
			L += lam[w]*(tmax-t);
		}
		t = tmax;
	} 
	
	return L;
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

		
		


	