#ifndef BEEPMBP__ABCSMC_HH
#define BEEPMBP__ABCSMC_HH

#include "struct.hh"
#include "inputs.hh"
#include "mvn.hh"
#include "state.hh"
#include "model.hh"
#include "obsmodel.hh"

class ABCSMC
{
public:	
  ABCSMC(const Details &details, const Data &data, const Model &model, const Inputs &inputs, const Output &output, const ObservationModel &obsmodel, Mpi &mpi);	
	void run();
	
private:
	void calculate_particle_weight(vector <Generation> &generation, const MVN &mvn);
	unsigned int effective_particle_number(const vector <double> &w) const;
	void setup_particle_sampler(const Generation &gen);
	unsigned int particle_sampler() const;
	bool terminate(long nac) const;
	void print_generation(const vector <Generation> &generation, double acrate) const;
	void implement_cutoff_frac(Generation &gen);
	void store_sample(Generation &gen, unsigned int g);
	
	unsigned int G;                          // The total number of generations 
	
	unsigned int Ntot;                       // The total number of particles

	double cutoff_frac;                      // Sets the acceptance fraction
	
	double propsize;                         // The size of the MVN proposals
	
	vector <double> wsum;	                   // Used to sample particles
		
	State state;                             // Stores the state of the system
		
	vector <Particle> particle_store;        // Stores the states in the last generation for output

	const Details &details;
	const Model &model;
	const Output &output;
	const ObservationModel &obsmodel;
	Mpi &mpi;
};

#endif
