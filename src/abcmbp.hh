#ifndef BEEPMBP__ABCMBP_HH
#define BEEPMBP__ABCMBP_HH

#include "struct.hh"
#include "mbp.hh"
#include "param_prop.hh"

class ABCMBP
{
public:	
  ABCMBP(const Details &details, const Data &data, const Model &model, Inputs &inputs, const Output &output, const ObservationModel &obsmodel, Mpi &mpi);	
	void run();
	
private:
	void EF_cutoff(Generation &gen, vector<Particle> &part, vector <unsigned int> &partcopy);
	bool terminate_generation(unsigned int g, double EFcut) const;
	void store_sample(Generation &gen);
	void model_evidence(vector <Generation> &generation);
	void print_generation(const Generation &gen, const unsigned int g) const;
	void set_time(Generation &gen);

	double cor_max;                          // The maximum correlation  
	double cor_max_last;                     // The maximum correlation in the last generation

	unsigned int G;                          // The total number of generations 

	double cpu_time;                         // Maximum cpu for execution

	double cutoff_final;                     // The final error function

	vector <Generation> generation;          // Stores information about each generation

	unsigned int N;                          // The number of particles per core
	unsigned int Ntot;                       // The total number of particles (across all runs)
	
	unsigned int npart;                      // The number of particles 
	unsigned int nrun;                       // The number of runs
	
	vector <Particle> part;		    	         // These particles store the state
	vector <unsigned int> partcopy;          // Store which particle to copy
	
	unsigned int nproposal;                  // The number of proposals
	unsigned int loop;                       // Stores the number of iterations of MCMC update
	
	State state;                             // Stores the state of the system
		 
	Mbp mbp;                                 // Used for making MBP-MCMC updates
	
	ParamProp paramprop;                     // Stores information about parameter proposals
 	
	//vector <ParamSample> psamp_GR;           // Parameter samples used by Gelman Rubin statistics
	
	const Details &details;
	const Data &data;
	const Model &model;
	const Output &output;
	const ObservationModel &obsmodel;
	Mpi &mpi;
};   

#endif
