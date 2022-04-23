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
Simulate::Simulate(const Details &details, Data &data, const Model &model, Inputs &inputs, const Output &output, const ObservationModel &obsmodel, Mpi &mpi) : state(details,data,model,obsmodel), details(details), data(data), model(model), obsmodel(obsmodel), output(output), mpi(mpi)
{	
	if(details.mode == SIM && mpi.ncore != 1) emsgroot("Simulation only requires one core");

	switch(details.mode){
		case SIM: break;
		case MULTISIM:                                                            // Sets the number of simulations
			nsim = inputs.find_positive_integer("nsimulation",UNSET); 
			if(nsim == UNSET) emsgroot("A value for 'nsimulation' must be set.");
			break; 
		case PREDICTION:                                                          // Sets the simulation per posterior sample
			nsim_per_sample = inputs.find_positive_integer("nsim_per_sample",4); 
			break;
		default: emsgEC("Simulate",1); break;
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

	particle_store.push_back(state.create_particle(0));                          // Generate the pdf output file
	output.generate_graphs(particle_store); 
	
	//state.save(details.output_directory+"/sample");
}


/// Runs multiple simulations
void Simulate::multisim()
{		
	auto nsim_per_core = nsim/mpi.ncore;
	if(nsim_per_core*mpi.ncore != nsim) emsgroot("'nsimulation' must be a multiple of the number of cores");
	
	for(auto s = 0u; s < nsim_per_core; s++){
		auto paramval = model.sample_from_prior();
		
		state.simulate(paramval);
		
		particle_store.push_back(state.create_particle(0));
		
		output.print_percentage((s+1)*mpi.ncore,nsim,percentage);
	}	

	output.generate_graphs(particle_store);         
}


/// Finds info about previously generated posterior samples and check these are consistent with the model
void Simulate::find_posterior_info(const string dir)
{
	nsample_post = UNSET; ndivision_post = UNSET;
	
	auto file = dir+"sample_information.txt";
	ifstream fin(file);
	if(!fin) emsg("Cannot open the file '"+file+"'");
	
	bool flag = false;
	for(auto loop = 0u; loop < 6; loop++){
		string line;
		getline(fin,line);
		auto spl = split(line,'='); if(spl.size() != 2) emsgEC("Simulate",2);
		auto spl2 = split(spl[1],'\''); 
		if(spl2.size() == 3){ if(spl2[2] == "") spl2.pop_back();}    // This removes /r
		if(spl2.size() != 2) emsgEC("Simulate",3);
		
		string quant = spl[0], val = spl2[1]; 
		unsigned int vali;
		if(loop != 1) vali = get_int(val,"Loading posterior samples"); 
		
		switch(loop){
			case 0: 
				if(quant != "# Samples") flag = true;
				nsample_post = vali; 
				break;
			
			case 1:
				if(quant != "TOML file" || val != details.toml_file){
					emsg("The TOML file '"+details.toml_file+"' must be the same as that used for inference '"+val+"'");
				}
				break;
				
			case 2: 
				if(quant != "start" || vali != details.start){
					emsg("The start value '"+details.getdate(details.start+details.start)+"' must be the same as that used for inference '"+details.getdate(details.start+vali)+"'");
				}
				break;
				
			case 3:
				if(quant != "end" || val != to_string(details.end)){
					emsg("The end value '"+details.getdate(details.start+details.end)+"' must be the same as that used for inference '"+details.getdate(details.start+vali)+"'");
				}
				break;
				
			case 4:
				if(quant != "division_per_time" || vali != details.division_per_time){
					emsg("The 'division_per_time' value '"+to_string(details.division_per_time)+"' must be the same as that used for inference '"+val+"'");
				}
				break;
				
			case 5:
				if(quant != "timesteps") flag = true;
				ndivision_post = vali; 
				break;
		}
	}

	if(nsample_post == UNSET) flag = true;
	if(ndivision_post == UNSET) flag = true;
	if(flag == true) emsg("There was a problem loading the posterior samples.");
}


/// Runs a analysis in which the model is modifies
void Simulate::model_modification()
{
	auto dir = details.output_directory+"/Posterior/samples/";
	
	find_posterior_info(dir);

	auto nsamp_per_core = nsample_post/mpi.ncore;

	if(nsamp_per_core*mpi.ncore != nsample_post) emsgroot("'nsample' must be a multiple of the number of cores");

	for(auto s = 0u; s < nsamp_per_core; s++){
		state.load(dir+"sample"+to_string(nsamp_per_core*mpi.core+s),ndivision_post);
		
		for(auto i = 0u; i < nsim_per_sample; i++){
			state.simulate(model.modelmod.pred_start*details.division_per_time,details.ndivision);
		
			particle_store.push_back(state.create_particle(0));
		}
		
		output.print_percentage(s+1,nsamp_per_core,percentage);
	}
	
	output.generate_graphs(particle_store);         
}
