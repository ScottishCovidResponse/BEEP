#ifndef BEEPMBP__ABC_HH
#define BEEPMBP__ABC_HH

#include "struct.hh"
#include "state.hh"
#include "output.hh"
#include "details.hh"

class ABC
{
public:	
  ABC(const Details &details, const Data &data, const Model &model, Inputs &inputs, const Output &output, const ObservationModel &obsmodel, Mpi &mpi);	
	void run();

private:
	void implement_cutoff_frac();
	bool terminate();
	void diagnostic() const;

	double cutoff;                           // Sets the cut-off used in the rejection sampling

	double cutoff_frac;                      // Sets the acceptance fraction
	
	unsigned int ntr, nac;                   // Gives the number of simulations tried and the number accepted 

	unsigned int Ntot;                       // Sets the total number of samples that need to be generated
	
	double cpu_time;                         // Sets the maximum CPU time for execution
	
	vector <Particle> particle_store;        // Stores the states

	unsigned int percentage;                 // Stores the percentage progress
	
	State state;                             // Stores the state of the system
	
	const Details &details;
	const Model &model;
	const Output &output;
	const ObservationModel &obsmodel;
	Mpi &mpi;
};

#endif
