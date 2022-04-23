#ifndef BEEPMBP__MC3_HH
#define BEEPMBP__MC3_HH

#include "struct.hh"
#include "mbp.hh"
#include "param_prop.hh"

class MC3
{
public:	
	MC3(const Details &details, const Data &data, const Model &model, Inputs &inputs, const Output &output, const ObservationModel &obsmodel, Mpi &mpi);
	
	void run();

private:
	void initialise();
	void update_burnin(const unsigned int samp);
	void swap_states();
	void store_param_samp(const unsigned int ch);
	void store_sample();
	void diagnostics();
	void set_invT(const unsigned int samp);
	void model_evidence() const;
	bool terminate(const unsigned int samp);
	void find_EF_range();
	void find_nchain_optimum() const;
		
	unsigned int nsample;                    // The number of MCMC samples
	
	double ESSmin;                            // The minimum effective sample size

	double cpu_time;                         // Maximum cpu for execution

	unsigned int nburnin;                    // The number of MCMC burnin steps
	
	unsigned int nquench;                    // The number of steps over which quenching is performed
	
	double Tpower;                           // The power used to specify the inverse temperature of chains
	 
	double invT_start;                       // The inverse temperature of the lowest invT (by default zero)
	
	double invT_final;                       // The inverse temperature of the posterior chain
	
	unsigned int thin;                       // Thinning on samples (to avoid memory usage)
	
	ParamUpdate pup;                         // Stores what parameter updates are done during MH
	
	bool burnin;                             // Determines if chain is being burned in
	
	unsigned int N;                          // The number of chains per core
	unsigned int Ntot;                       // The total number of chains (across all runs)

	unsigned int nchain;                     // The number of chains per run
	unsigned int nrun;                       // The number of runs
	
	vector <Particle> part;		    	         // These particles store the state
	
	vector <Particle> part_plot;             // The particles are used for plotting the posterior
	
	vector <Chain> chain;		    	           // These particles store the state

	vector <ParamProp> paramprop;            // Stores information about parameter proposals
	
	State state;                             // Stores the state of the system
		
	Mbp mbp;                                 // Used for making MCMC-MBP updates
	
	unsigned int percentage;                 // Stores the percentage progress
	
	const Details &details;
	const Data &data;
	const Model &model;
	const Output &output;
	const ObservationModel &obsmodel;
	Mpi &mpi;
};   

#endif
