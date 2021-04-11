#ifndef BEEPMBP__SIMULATE_HH
#define BEEPMBP__SIMULATE_HH

#include "struct.hh"
#include "state.hh"

class Simulate
{
	public:	
		Simulate(const Details &details, Data &data, const Model &model, const AreaTree &areatree, const Inputs &inputs, const Output &output, const ObservationModel &obsmodel, Mpi &mpi);	
		void run();
		void multirun();
		void counter();
		
	private:
		unsigned int nsim;                                    // The number of simulations (MULTISIM)
		
		unsigned int nsample;                                 // The number of simulations (COUNTER)
		
		vector <Particle> particle_store;                     // Stores particles
		
		State state;                                          // Stores the state
		
		unsigned int percentage;                              // For displaying the percentage complete   
		
		const Details &details;
		Data &data;
		const Model &model;
		const AreaTree &areatree;
		const ObservationModel &obsmodel;
		const Output &output;
		Mpi &mpi;
};

#endif
