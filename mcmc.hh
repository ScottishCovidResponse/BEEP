#ifndef BEEPMBP__MBP_HH
#define BEEPMBP__MBP_HH

using namespace std;

#include "output.hh"
#include "chain.hh"
#include "model_evidence.hh"

class Data;
class Model;
class AreaTree;
class Output;
class ObservationModel;

enum class proposal_method                                                // Different types of MC3 proposal
{
	allchainsallparams,
	fixednum,
	fixedtime
};

struct ChainInfo                                                          // Stores information about chains (for swapping)
{
	double invT;
	unsigned int ch;
	Jump jump;
};

class Mcmc                                                                // Contains everything related to an MCMC chain
{
	public:	
		Mcmc(const Details &details, const Data &data, const Model &model, const AreaTree &areatree, const Mpi &mpi, 
				 const Inputs &inputs, Output &output, const ObservationModel &obsmodel);
		
		void run();

	private:
		void update();
		void output_measurements(vector <Sample> &opsamp, unsigned int nchain) const;
		void output_parameters(Output &output, vector <ParamSample> &psamp) const;
		void diagnostics() const;
		void chaininfo_set(ChainInfo &chinf, const Chain &cha) const;
		void chain_set(Chain &cha, const ChainInfo &chinf) const;
		void swap_chains();
		void optimise_timeloop();
		unsigned int select_random(const vector <unsigned int> &vec) const;

		vector <Chain> chain; 
		
		unsigned int nsamp;                      // The number of MCMC samples
		unsigned int burnin;                     // The number of burnin samples
		
		unsigned int nchain_total;               // The total number of chains (across all MPI processes);
		unsigned int nchain;                     // The number of chains per core
		
		vector <double> nac_swap;                // Stores the rate at which swaps are accepted
		
		long timeprop, ntimeprop;                // Stores times so that number of proposals can be optimised
		double timeloop;
		
		Model_Evidence model_evidence;           // Stores information used to calculate the model evidence
		
		enum proposal_method propsmethod;        // Stores the type of proposal method
		
		const Details &details;
		const Data &data;
		const Model &model;
		const AreaTree &areatree;
		const Mpi &mpi;
		Output &output;
		const ObservationModel &obsmodel;
};

#endif
