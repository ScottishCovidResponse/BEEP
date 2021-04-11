#ifndef BEEPMBP__MC3_HH
#define BEEPMBP__MC3_HH

#include "struct.hh"
#include "mbp.hh"
#include "param_prop.hh"

struct Chain{                              // Stores information about an MCMC chain (used in MC3)
	Chain(unsigned int num_, const Details &details, const Data &data, const Model &model, const AreaTree &areatree, const Output &output, Mpi &mpi);

	unsigned int num;                        // The number of the chain
	double invT;                             // The inverse temperature
	vector < vector <double> > param_samp;   // Parameter samples

	//vector <double> EF_samp;                 // Likelihood samples generated
	//vector < vector <double> > EF_datatable; // Stores how the error function is divided into datatables
	ParamProp paramprop;                     // Stores information about parameter proposals
	unsigned int nproposal;                  // The number of proposals
	
	unsigned int ntr;                        // The number of times a swap is tried
	unsigned int nac;                        // The number of times a swap is accepted
};

class MC3
{
public:	
	MC3(const Details &details, const Data &data, const Model &model, const AreaTree &areatree, const Inputs &inputs, const Output &output, const ObservationModel &obsmodel, Mpi &mpi);
	
	void run();

private:
	void initialise();
	void update_burnin(unsigned int samp);
	void swap_states();
	void store_sample();
	void diagnostics();
	void set_invT(unsigned int samp);
		
	unsigned int nsample;                    // The number of MCMC samples
	
	unsigned int nburnin;                    // The number of MCMC burnin steps
	
	unsigned int nquench;                    // The number of steps over which quenching is performed
	
	double invT_start;                       // The inverse temperature of the lowest invT (by default zero)
	
	double invT_final;                       // The inverse temperature of the posterior chain
	
	unsigned int thin;                       // Thinning on samples (to avoid memory usage)
	
	ParamUpdate pup;                         // Stores what parameter updates are done during MH
	
	bool burnin;                             // Determines if chain is being burned in
	
	unsigned int N;                          // The number of chains per core
	unsigned int Ntot;                       // The total number of chains

	vector <Particle> part;		    	         // These particles store the state
	
	vector <Particle> part_plot;             // The particles are used for plotting the posterior
	
	vector <Chain> chain;		    	           // These particles store the state

	State state;                             // Stores the state of the system
		
	Mbp mbp;                                 // Used for making MCMC-MBP updates
	
	unsigned int percentage;                 // Stores the percentage progress
	
	const Details &details;
	const Data &data;
	const Model &model;
	const AreaTree &areatree;
	const Output &output;
	const ObservationModel &obsmodel;
	Mpi &mpi;
};   

#endif
