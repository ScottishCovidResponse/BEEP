//#include <stdio.h>
#include <cmath>   


using namespace std;

#include "state.hh"

State::State(const Details &details, const DATA &data, const MODEL &model) : details(details), data(data), model(model)
{
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
