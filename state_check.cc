/// These are algorithms to check state is working correctly

#include <math.h>  

using namespace std;

#include "state.hh"

/// Checks quantities in the state are correct
void State::check()
{
	if(model.inbounds(paramval) == false) emsg("Parameters not in bounds");
	
	for(auto sett = 0u; sett < details.ndivision; sett++){
		for(auto c = 0u; c < data.narea; c++){
			for(auto tr = 0u; tr < model.trans.size(); tr++){
				auto from = model.trans[tr].from;
				if(model.trans[tr].type == INFECTION_DIST || model.trans[tr].type == EXP_DIST){
					for(auto dp = 0u; dp < data.ndemocatpos; dp++){	 
						auto num = transnum[sett][c][tr][dp];
						if(num < 0) emsgEC("State_check",15);
						
						auto p = pop[sett][c][from][dp];
						if(p <= 0){
							if(num != 0) emsgEC("State_check",74);
						}
						else{
							if(transrate[tr][dp] == 0 ){
								if(num != 0 && model.trans[tr].type != INFECTION_DIST) emsgEC("State_check",75);
							}
						}
					}
				}
			}
		}
	}		
			
	for(auto sett = 0u; sett < details.ndivision-1; sett++){  // Checks pop is correct
		auto popst = pop[sett+1];
		update_pop(sett);
		
		for(auto c = 0u; c < data.narea; c++){
			for(auto co = 0u; co < model.comp.size(); co++){
				for(auto dp = 0u; dp < data.ndemocatpos; dp++){
					if(popst[c][co][dp] != pop[sett+1][c][co][dp]) emsgEC("State_check",1);
				}
			}
		}
	}	
	
	set_Imap(1);  // Checks to make sure Imap is OK
	
	if(details.mode == ABC_SMC || details.mode == ABC_MBP){
		auto EF_temp = -2*obsmodel.calculate(this);
		auto dd = EF - EF_temp;
		if(dd*dd > TINY){
			cout << EF_temp << " " << EF << " obs\n";
			emsgEC("State_check",2);
		}
	}
	
	auto dd = Pr - model.prior(paramval); if(sqrt(dd*dd) > TINY) emsgEC("State_check",59);
	
	auto af = model.create_areafactor(paramval);                      // Checks areafactor
	for(auto c = 0u; c < data.narea; c++){
		if(af[c] != areafactor[c]) emsgEC("State_check",88);
	}
	
	auto sus = model.create_susceptibility(paramval);                 // Checks susceptibility
	for(auto dp = 0u; dp < data.ndemocatpos; dp++){
		if(sus[dp] != susceptibility[dp]) emsgEC("State_check",89);
	}

	for(auto sp = 0u; sp < model.spline.size(); sp++){                // Checks disc_spline
		auto spl = model.create_disc_spline(sp,paramval);  
		for(auto se = 0u; se < details.ndivision; se++){
			if(spl[se] != disc_spline[sp][se]) emsgEC("State_check",90);
		}
	}
	
	auto ctprob = model.create_comptransprob(paramval);               // Checks comptransprob
	for(auto c = 0u; c < model.comp.size(); c++){
		for(auto j = 0u; j < comptransprob[c].prob.size(); j++){
			for(auto k = 0u; k < comptransprob[c].prob[j].size(); k++){
				if(ctprob[c].prob[j][k] != comptransprob[c].prob[j][k]) emsgEC("State_check",91);
			}
		}
	}
} 
