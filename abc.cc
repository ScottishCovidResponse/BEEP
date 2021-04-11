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

ABC::ABC(const Details &details, const Data &data, const Model &model, const Inputs &inputs, const Output &output, const ObservationModel &obsmodel, Mpi &mpi) : state(details,data,model,obsmodel), details(details), model(model), output(output), obsmodel(obsmodel), mpi(mpi)
{
	inputs.find_nsample(Ntot);
	inputs.find_cutoff(cutoff,cutoff_frac);
	percentage = UNSET;
}


/// Implements a version of simple ABC rejection algorithm
void ABC::run()
{
	auto ntr = 0u, nac = 0u;
	
	timer[TIME_ALG].start();
	do{	
		auto param = model.sample_from_prior();            // Samples parameters from the prior

	  state.simulate(param);                             // Simulates a state
		
		ntr++;
		if(cutoff == UNSET || state.EF < cutoff){          // Stores the state if the error function is below the cutoff    
			particle_store.push_back(state.create_particle());
			nac++;
		}
			
		double num = nac; if(cutoff_frac != UNSET) num *= cutoff_frac;
		output.print_percentage(num,Ntot,percentage);      // Prints the percentage to show progress
	}while(!terminate());                                // Terminates when sufficient samples are generated
	
	
	if(cutoff_frac != UNSET) implement_cutoff_frac();    // Implements the acceptance rate to give cut-off
	else cutoff_frac = mpi.get_acrate(nac,ntr);          // Calculates acceptance rate
	
	diagnostic();                                        // Outputs diagnostic information

	output.generate_graphs(particle_store);		           // Outputs the results
	
	timer[TIME_ALG].stop();	
}


/// Determines when to terminate the algorithm
bool ABC::terminate() const
{
	bool term = false;

	auto samptot = mpi.sum(long(particle_store.size()));
	
	if(mpi.core == 0){
		if(cutoff != UNSET){ if(samptot >= Ntot) term = true;}
		else{ if(samptot >= Ntot/cutoff_frac) term = true;}
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
		cutoff = EFtot[Ntot]; 
	}
	mpi.bcast(cutoff);
	
	auto p = 0u;                                                            // Removes particles above the cut-off
	while(p < particle_store.size()){
		if(particle_store[p].EF >= cutoff){
			particle_store.erase(particle_store.begin()+p);
		}
		else p++;
	}
}


/// Diagnostic information
void ABC::diagnostic() const
{
 if(mpi.core == 0){
		cout << endl << "EF cut-off: " << prec(cutoff,2) << "      Fraction accepted: " << per(cutoff_frac) << endl;		
	}
}
	