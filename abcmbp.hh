#ifndef BEEPMBP__ABCMBP_HH
#define BEEPMBP__ABCMBP_HH

#include "struct.hh"
#include "mbp.hh"
#include "param_prop.hh"

class ABCMBP
{
public:	
  ABCMBP(const Details &details, const Data &data, const Model &model, const AreaTree &areatree, const Inputs &inputs, const Output &output, const ObservationModel &obsmodel, Mpi &mpi);	
	void run();
	
private:
	void EF_cutoff(Generation &gen, vector<Particle> &part, vector <unsigned int> &partcopy);
	void store_sample(Generation &gen);
	void print_generation(const Generation &gen, unsigned int g) const;

	unsigned int G;                          // The total number of generations 

	vector <Generation> generation;          // Stores information about each generation

	unsigned int N;                          // The number of particles per core
	unsigned int Ntot;                       // The total number of particles

	vector <Particle> part;		    	         // These particles store the state
	vector <unsigned int> partcopy;          // Store which particle to copy
	
	unsigned int nproposal;                  // The number of proposals
	
	State state;                             // Stores the state of the system
		 
	Mbp mbp;                                 // Used for making MBP-MCMC updates
	
	ParamProp paramprop;                     // Stores information about parameter proposals
 	
	const Details &details;
	const Data &data;
	const Model &model;
	const AreaTree &areatree;
	const Output &output;
	const ObservationModel &obsmodel;
	Mpi &mpi;
};   

#endif
