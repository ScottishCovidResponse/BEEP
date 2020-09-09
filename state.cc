#include <cmath>   
#include <iostream>

using namespace std;

#include "state.hh"


State::State(const Details &details, const DATA &data, const MODEL &model, const Obsmodel &obsmodel) : details(details), data(data), model(model), obsmodel(obsmodel)
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

/// Checks the observation likelihood and prior are correct
void State::check() const
{
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



	