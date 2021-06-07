#ifndef BEEPMBP__PAIS_HH
#define BEEPMBP__PAIS_HH

#include "struct.hh"
#include "mbp.hh"
#include "param_prop.hh"

class PAIS
{
public:	
  PAIS(const Details &details, const Data &data, const Model &model, Inputs &inputs, const Output &output, const ObservationModel &obsmodel, Mpi &mpi);	
	void run();
	
private:
	void bootstrap(Generation &gen, vector<Particle> &part, vector <unsigned int> &partcopy, const double invT);
	void store_sample(Generation &gen);
	void param_GR_clear();
	bool terminate();
	void model_evidence(vector <Generation> &generation);
	void print_generation(const Generation &gen, const unsigned int g) const;

	unsigned int G;                          // The total number of generations 

	double invT_final;                       // The final inverse temperature

	vector <Generation> generation;          // Stores information about each generation

	unsigned int N;                          // The number of particles per core
	unsigned int Ntot;                       // The total number of particles

	unsigned int npart;                      // The number of particles 
	unsigned int nrun;                       // The number of runs
	
	unsigned int nupdate;                    // Sets the number of updates per generation
	double GRmax;                            // The maximum value for the Gelman-Rubin statistics
	
	vector <Particle> part;		    	         // These particles store the state
	vector <unsigned int> partcopy;          // Store which particle to copy
		
	unsigned int nproposal;                  // The number of proposals
	unsigned int loop;                       // Stores the number of iterations of MCMC update
	
	double quench_factor;                    // Parameter that determines how quickly the PAIC algorithm quenches
		
	State state;                             // Stores the state of the system
		
	Mbp mbp;                                 // Used for making MBP-MCMC updates
	
	ParamProp paramprop;                     // Stores information about parameter proposals
 	
	vector <ParamSample> psamp_GR;           // Parameter samples used by Gelman Rubin statistics
	
	const Details &details;
	const Data &data;
	const Model &model;
	const Output &output;
	const ObservationModel &obsmodel;
	Mpi &mpi;
};   

#endif
