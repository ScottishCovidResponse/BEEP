#ifndef BEEP__SIMULATE_HH
#define BEEP__SIMULATE_HH

#include "struct.hh"
#include "state.hh"

class Simulate
{
	public:	
		Simulate(const Details &details, Data &data, const Model &model, Inputs &inputs, const Output &output, const ObservationModel &obsmodel, Mpi &mpi);	
		void run();
		void multisim();
		void model_modification();
		
	private:
		void find_posterior_info(const string dir);               
	
		unsigned int nsim;                                    // The number of simulations (MULTISIM)
		
		unsigned int nsim_per_sample;                         // The number of simulatation per posterior sample 
		
		unsigned int nsample_post;                            // The number of posterior samples
	
		unsigned int ndivision_post;                          // The number of timesteps loaded in the posterior
		
		vector <Particle> particle_store;                     // Stores particles
		
		State state;                                          // Stores the state
		
		unsigned int percentage;                              // For displaying the percentage complete   
		
		const Details &details;
		Data &data;
		const Model &model;
		const ObservationModel &obsmodel;
		const Output &output;
		Mpi &mpi;
};

#endif
