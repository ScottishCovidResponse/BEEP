/// Implements the basic ABC rejection algorithm

#include <iostream>
#include <fstream>
#include <algorithm> 
#include <sstream> 

using namespace std;

#include "abc.hh"
#include "mvn.hh"
#include "inputs.hh"
#include "model.hh"
#include "mpi.hh"

ABC::ABC(const Details &details, const Data &data, const Model &model, Inputs &inputs, const Output &output, const ObservationModel &obsmodel, Mpi &mpi) : state(details,data,model,obsmodel), details(details), model(model), output(output), obsmodel(obsmodel), mpi(mpi)
{
	inputs.find_nrun(nrun);
	inputs.find_nsample_GRmax(Ntot,GRmax,nrun);
	inputs.find_cutoff(cutoff,cutoff_frac);
	if(cutoff_frac != UNSET && GRmax != UNSET) emsgroot("'cutoff_frac' cannot be used with 'GRmax'");
	percentage = UNSET;
}


/// Implements a version of simple ABC rejection algorithm
void ABC::run()
{
	ntr.resize(nrun); nac.resize(nrun); 
	for(auto ru = 0u; ru < nrun; ru++){ ntr[ru] = 0; nac[ru] = 0;}
	
	timer[TIME_ALG].start();
	do{	
		for(auto ru = 0u; ru < nrun; ru++){
			do{
				auto param = model.sample_from_prior();        // Samples parameters from the prior

				state.simulate(param);                         // Simulates a state
			
				ntr[ru]++;
				if(cutoff == UNSET || state.EF < cutoff){      // Stores the state if the error function is below the cutoff    
					particle_store.push_back(state.create_particle(ru));
					nac[ru]++; 
					break;
				}
			}while(true);
		}	
	}while(!terminate());                                // Terminates when sufficient samples are generated
	
	if(cutoff_frac != UNSET) implement_cutoff_frac();    // Implements the acceptance rate to give cut-off

	diagnostic();                                        // Outputs diagnostic information

	output.generate_graphs(particle_store);		           // Outputs the results
	
	timer[TIME_ALG].stop();	
}


/// Determines when to terminate the algorithm
bool ABC::terminate()
{
	bool term = false;

	auto samptot = mpi.sum(long(particle_store.size()));
	if(GRmax != UNSET){
		auto psamp_tot = mpi.gather_psamp(particle_store);
		if(mpi.core == 0){
			auto GR = output.get_Gelman_Rubin_statistic(psamp_tot);
			if(vec_max(GR) < GRmax && samptot >= nrun*20) term = true; 
			cout << "Number of samples: " << samptot << "    Largest GR value:" << vec_max(GR) << "     GRmax: " << GRmax << endl;
		}
	}
	else{
		if(mpi.core == 0){
			if(cutoff != UNSET){ 
				if(samptot >= nrun*Ntot) term = true;
				output.print_percentage(samptot,nrun*Ntot,percentage);
			}
			else{
				if(samptot >= nrun*Ntot/cutoff_frac) term = true;
				output.print_percentage(samptot*cutoff_frac,nrun*Ntot,percentage);
			}
		}
	}
	
	mpi.bcast(term);

	return term;
}


/// Selects EFcut to remove samples at a certain acceptance rate
void ABC::implement_cutoff_frac()
{
	vector <double> EF;                                                     // Calculates the cut-off
	for(auto p = 0u; p < particle_store.size(); p++) EF.push_back(particle_store[p].EF);
	
	auto EFtot = mpi.combine(EF);

	if(mpi.core == 0){
		sort(EFtot.begin(),EFtot.end());
		cutoff = EFtot[nrun*Ntot];
	}
	mpi.bcast(cutoff);
	
	auto p = 0u;                                                            // Removes particles above the cut-off
	while(p < particle_store.size()){
		if(particle_store[p].EF >= cutoff){
			nac[particle_store[p].run]--;
			particle_store.erase(particle_store.begin()+p);
		}
		else p++;
	}
}


/// Diagnostic information
void ABC::diagnostic() const
{
	vector <double> ac_rate(nrun);
	for(auto ru = 0u; ru < nrun; ru++) ac_rate[ru] = mpi.get_acrate(nac[ru],ntr[ru]);     
	
	if(mpi.core == 0){
		auto stat = output.get_statistic(ac_rate);

		cout << endl << "EF cut-off: " << prec(cutoff,2) << "      Fraction accepted: " << per(stat.mean) << endl;	

		vector <double> ME_list; for(auto ru = 0u; ru < nrun; ru++) ME_list.push_back(log(ac_rate[ru]));

		output.final_model_evidence(ME_list,UNSET,cutoff);
	}
}
	