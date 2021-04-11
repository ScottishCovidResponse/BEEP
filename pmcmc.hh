#ifndef BEEPMBP__PMCMC_HH
#define BEEPMBP__PMCMC_HH

#include "struct.hh"
#include "param_prop.hh"

class PMCMC
{
public:	
  PMCMC(const Details &details, const Data &data, const Model &model, const AreaTree &areatree, const Inputs &inputs, Output &output, const ObservationModel &obsmodel, Mpi &mpi);	
	
	void run();
	
private:
	double obs_prob(vector <double> &paramv);
	double bootstrap(unsigned int sec, vector <double> &L);
	Particle particle_sample();
	void initialise();
	void initialise_variables();
	void update_burnin(unsigned int samp);
	
	void get_proposals();
	void mcmc_updates();
	double get_al();
	void copy_propose_state();
	void self_proposal(Self &self);
	void mvn_proposal(MVN &mvn);
	void sigma_reff_proposal(MVN &mvn);
	void joint_proposal(Joint &jn);
	void neighbour_proposal(Neighbour &nei);
	void mean_time_proposal(MeanTime& mt);

	double invT;                               // This inverse temperature is used to smooth the likelihood function
	
	double Li;                                 // The unbiased estimator for likelihood
	Particle Pi;                               // A particle for the initial state
	double Pri;                                // The prior for the initial state
		
	double Lp;                                 // The unbiased estimator for the likelihood (proposed state)
	Particle Pp;                               // A particle for the proposed state 
	double Prp;                                // The prior for the proposed state

	bool burnin;                               // Determines if chain is being burned in
	
	ParamUpdate pup;                           // Stores what parameter updates are done during MH
	
	bool invT_dynamic;                         // Determines if invT is dynamically updated
	
	unsigned int core;                         // The number of the core
	
	unsigned int nsample;                      // The number of MCMC samples
	
	unsigned int nburnin;                      // The number of MCMC burnin steps
	
	unsigned int thin;                         // Thinning on samples (to avoid memory usage)
	
	unsigned int N;                            // The number of particles per core
	unsigned int Ntot;                         // The total number of particles

	vector <Proposal> prop_list;               // A list of MCMC proposals
		
	vector <vector <double> > param_samp;      // This is a list of parameter samples to approximate MVN distributions 
	
	vector < vector <unsigned int> > backpart; // How particles are related through the bootstrap step

	vector <State> particle;                   // The vector of states
	
	vector <Particle> particle_store;          // Stores the states for output later
	
	unsigned int buffersize;                   // The size of the buffer for sending / recieving during MPI
	
	unsigned int percentage;                   // Stores the percentage progress
	
	ParamProp paramprop;                       // Stores information about parameter proposals
	
	ofstream trace;                            // Used to output a trace plot
	
	const Details &details;
	const Data &data;
	const Model &model;
	const AreaTree &areatree;
	const Output &output;
	const ObservationModel &obsmodel;
	Mpi &mpi;
};

#endif
