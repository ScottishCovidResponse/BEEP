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
	bool terminate_generation(unsigned int g, double invT) const;
	void model_evidence(vector <Generation> &generation);
	void print_generation(const Generation &gen, const unsigned int g) const;

	double cor_max;                          // The maximum correlation  
	double cor_max_last;                     // The maximum correlation in the last generation

	unsigned int G;                          // The total number of generations 

	double cpu_time;                         // Maximum cpu for execution

	double invT_final;                       // The final inverse temperature

	vector <Generation> generation;          // Stores information about each generation

	unsigned int N;                          // The number of particles per core
	unsigned int Ntot;                       // The total number of particles

	unsigned int npart;                      // The number of particles 
	unsigned int nrun;                       // The number of runs
	
	vector <Particle> part;		    	         // These particles store the state
	vector <unsigned int> partcopy;          // Store which particle to copy
		
	unsigned int nproposal;                  // The number of proposals
	unsigned int loop;                       // Stores the number of iterations of MCMC update
	
	double quench_factor;                    // Parameter that determines how quickly the PAIC algorithm quenches
		
	State state;                             // Stores the state of the system
		
	Mbp mbp;                                 // Used for making MBP-MCMC updates
	
	ParamProp paramprop;                     // Stores information about parameter proposals
 	
	//ofstream invT_vs_EF;                     // Used to generate a plot of invT vs EF
	
	const Details &details;
	const Data &data;
	const Model &model;
	const Output &output;
	const ObservationModel &obsmodel;
	Mpi &mpi;
};   

#endif
