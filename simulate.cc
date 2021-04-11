// Implements a modified Gillespie algorithm to simulate from the model

#include <iostream>
#include <fstream>
#include <algorithm>
#include <assert.h>
#include <math.h>

using namespace std;

#include "simulate.hh"
#include "data.hh"
#include "model.hh"
#include "output.hh"

/// Initilaises the simulation
Simulate::Simulate(const Details &details, Data &data, const Model &model, const AreaTree &areatree, const Inputs &inputs, const Output &output, const ObservationModel &obsmodel, Mpi &mpi) : state(details,data,model,obsmodel), details(details), data(data), model(model), areatree(areatree), obsmodel(obsmodel), output(output), mpi(mpi)
{	
	if(details.mode == SIM && mpi.ncore != 1) emsgroot("Simulation only requires one core");

	switch(details.mode){
		case SIM: break;
		case MULTISIM:                                                                      // Sets the number of simulations
			nsim = inputs.find_integer("nsimulation",UNSET); 
			if(nsim == UNSET) emsg("'nsimulation' must be specified.");
			break; 
		case COUNTER: inputs.find_nsample(nsample); break;                                  // The number of posterior samples
		case PPC: inputs.find_nsample(nsample); break;
		default: emsgEC("Simulate",14); break;
	}
	
	percentage = UNSET;
}


/// Performs a simulation
void Simulate::run()
{
	auto paramval = model.sample_from_prior();                                   // Sets up the parameters
	
	state.simulate(paramval);                                                    // Simulates the state
	
	auto obs_value = obsmodel.get_obs_value(&state);                             // Copies the observations into the data
	for(auto ob = 0u; ob < data.obs.size(); ob++) data.obs[ob].value = obs_value[ob];
	
	output.simulated_data(obs_value,details.output_directory+"/Simulated_data"); // Outputs the simulated data files
	
	particle_store.push_back(state.create_particle());                           // Generate the pdf output file
	output.generate_graphs(particle_store); 
}


/// Runs multiple simulations
void Simulate::multirun()
{		
	auto nsim_per_core = nsim/mpi.ncore;
	if(nsim_per_core*mpi.ncore != nsim) emsgroot("'nsimulation' must be a multiple of the number of cores");
	
	for(auto s = 0u; s < nsim_per_core; s++){
		auto paramval = model.sample_from_prior();
		
		state.simulate(paramval);
		
		particle_store.push_back(state.create_particle());
		
		output.print_percentage(s+1,nsim,percentage);
	}	

	output.generate_graphs(particle_store);         
}


/// Runs a counter factual analysis or posterior predictive check
void Simulate::counter()
{
	auto nsamp_per_core = nsample/mpi.ncore;

	if(nsamp_per_core*mpi.ncore != nsample) emsgroot("'nsample' must be a multiple of the number of cores");
	
	auto dir = details.output_directory+"/Posterior/samples/";
	for(auto s = 0u; s < nsamp_per_core; s++){
		state.load(dir+"sample"+to_string(nsamp_per_core*mpi.core+s)+".txt");
		
		state.simulate(model.countermod.sett_start,details.ndivision);
		
		particle_store.push_back(state.create_particle());
		
		output.print_percentage(s+1,nsample,percentage);
	}
	
	output.generate_graphs(particle_store);         
}
