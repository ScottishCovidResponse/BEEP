#ifndef BEEP__ABCSMC_HH
#define BEEP__ABCSMC_HH

#include "struct.hh"
#include "inputs.hh"
#include "mvn.hh"
#include "state.hh"
#include "model.hh"
#include "obsmodel.hh"

class ABCSMC
{
public:	
  ABCSMC(const Details &details, const Data &data, const Model &model, Inputs &inputs, const Output &output, const ObservationModel &obsmodel, Mpi &mpi);	
	void run();
	
private:
	double calculate_particle_weight(const vector <double> param_prop, const Generation &gen_last, const MVN &mvn);
	void normalise_particle_weights(Generation &gen);
	
	unsigned int effective_particle_number(const vector <double> &w) const;
	void setup_particle_sampler(const Generation &gen);
	unsigned int particle_sampler() const;
	bool terminate(const Generation &gen) const;
	bool terminate_generation(unsigned int g, double EFcut) const;
	void print_generation(const vector <Generation> &generation, const double acrate) const;
	void implement_cutoff_frac(Generation &gen);
	void remove_small_weights(Generation &gen);
	void store_sample(Generation &gen, const double w);
	void print_model_evidence();
	void results();
	
	vector <Generation> generation;          // Stores information from different generation
	
	unsigned int G;                          // The total number of generations 
	
	double cpu_time;                         // Limit on the cpu time
	
	double cutoff_final;                     // The final error function

	unsigned int Ntot;                       // The total number of particles

	unsigned int Ntot_final;                 // The total number of particles in the final generation
	
	//double GRmax;                          // The maximum value for the Gelman-Rubin statistics

	double cutoff_frac;                      // Sets the acceptance fraction
	double cutoff_frac_init;                 // Sets the acceptance fraction in the initial generation
	
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
