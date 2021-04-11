#ifndef BEEPMBP__PAIS_HH
#define BEEPMBP__PAIS_HH

#include "struct.hh"
#include "mbp.hh"
#include "param_prop.hh"

class PAIS
{
public:	
  PAIS(const Details &details, const Data &data, const Model &model, const AreaTree &areatree, const Inputs &inputs, const Output &output, const ObservationModel &obsmodel, Mpi &mpi);	
	void run();
	
private:
	void bootstrap(Generation &gen, vector<Particle> &part, vector <unsigned int> &partcopy, double invT);
	void store_sample(Generation &gen);
	void print_generation(const Generation &gen, unsigned int g) const;

	unsigned int G;                          // The total number of generations 

	double invT_final;                       // The final inverse temperature

	vector <Generation> generation;          // Stores information about each generation

	unsigned int N;                          // The number of particles per core
	unsigned int Ntot;                       // The total number of particles

	vector <Particle> part;		    	         // These particles store the state
	vector <unsigned int> partcopy;          // Store which particle to copy
		
	double quench_factor;                    // Parameter that determines how quickly the PAIC algorithm quenches
		
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
