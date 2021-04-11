#ifndef BEEPMBP__ABCMBP_GR_HH
#define BEEPMBP__ABCMBP_GR_HH

#include "struct.hh"
#include "mbp.hh"
#include "param_prop.hh"

class ABCMBP_GR
{
public:	
  ABCMBP_GR(const Details &details, const Data &data, const Model &model, const AreaTree &areatree, const Inputs &inputs, const Output &output, const ObservationModel &obsmodel, Mpi &mpi);	
	void run();
	
private:
	void EF_cutoff(Generation &gen, vector<Particle> &part, vector <unsigned int> &partcopy);
	void param_GR_clear();
	bool gelman_rubin_termination(double Nexp);
	void store_sample(Generation &gen);
	void print_generation(const Generation &gen, unsigned int g) const;

	unsigned int G;                          // The total number of generations 

	vector <Generation> generation;          // Stores information about each generation

	unsigned int N;                          // The number of particles per core
	unsigned int Ntot;                       // The total number of particles
	
	unsigned int ngroup;                     // The number of groups the particles are split
	unsigned int groupN;                     // Particles per group
	
	vector <Particle> part;		    	         // These particles store the state
	vector <unsigned int> partcopy;          // Store which particle to copy
	
	unsigned int nproposal;                  // The number of proposals
	
	State state;                             // Stores the state of the system
		 
	Mbp mbp;                                 // Used for making MBP-MCMC updates
	
	ParamProp paramprop;                     // Stores information about parameter proposals
 	
	vector < vector < vector <double> > > param_GR; // Parameter samples used by Gelman Rubin statistics
	
	const Details &details;
	const Data &data;
	const Model &model;
	const AreaTree &areatree;
	const Output &output;
	const ObservationModel &obsmodel;
	Mpi &mpi;
};   

#endif
