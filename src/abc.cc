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
	inputs.find_nsample(Ntot,1);
	inputs.find_cpu_time(cpu_time); 
	if(Ntot != UNSET && Ntot%mpi.ncore != 0)  emsgroot("'nsample' must be a multiple of the number of cores");
	inputs.find_cutoff(cutoff,cutoff_frac);
	percentage = UNSET;
}


/// Implements a version of simple ABC rejection algorithm
void ABC::run()
{
	ntr = 0; nac = 0;
	timer[TIME_ALG].start();
	do{	
		for(auto loop = 0u; loop < 10; loop++){
			auto param = model.sample_from_prior();          // Samples parameters from the prior

			state.simulate(param);                           // Simulates a state
			ntr++;
			if(cutoff == UNSET || state.EF < cutoff){        // Stores the state if the error function is below the cutoff    
				particle_store.push_back(state.create_particle(UNSET));
				nac++; 
			}
		}
	}while(!terminate());                                // Terminates when sufficient samples are generated
	timer[TIME_ALG].stop();	
		
	if(cutoff_frac != UNSET) implement_cutoff_frac();    // Implements the acceptance rate to give cut-off

	diagnostic();                                        // Outputs diagnostic information

	output.generate_graphs(particle_store,UNSET);		     // Outputs the results
}


/// Determines when to terminate the algorithm
bool ABC::terminate()
{
	bool term = false;

	auto samptot = mpi.sum(long(particle_store.size()));

	if(mpi.core == 0){
		if(cutoff != UNSET){ 
			if(samptot >= Ntot) term = true;
			output.print_percentage(samptot,Ntot,percentage);
		}
		else{
			if(samptot >= Ntot/cutoff_frac) term = true;
			output.print_percentage(samptot*cutoff_frac,Ntot,percentage);
		}
	}
	
	mpi.bcast(term);

	if(cpu_time != UNSET){
		if((clock() - details.time_start)/(60.0*CLOCKS_PER_SEC) > cpu_time){
			emsg("Could not calculate within time limit");
		}
	}
	
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
		cutoff = EFtot[Ntot];
	}
	mpi.bcast(cutoff);
	
	auto p = 0u;                                                            // Removes particles above the cut-off
	while(p < particle_store.size()){
		if(particle_store[p].EF >= cutoff){
			nac--;
			particle_store.erase(particle_store.begin()+p);
		}
		else p++;
	}
}


/// Diagnostic information
void ABC::diagnostic() const
{
	auto ac_rate = mpi.get_acrate(nac,ntr);     
	
	if(mpi.core == 0){
		cout << endl << "EF cut-off: " << prec(cutoff,2) << "      Fraction accepted: " << per(ac_rate) << endl;	

		vector <double> ME_list; ME_list.push_back(log(ac_rate));

		output.final_model_evidence(ME_list,UNSET,cutoff);
	}
}
	