#ifndef BEEPMBP__MBP_HH
#define BEEPMBP__MBP_HH

using namespace std;

#include "output.hh"
#include "chain.hh"
#include "model_evidence.hh"

enum class proposalsmethod
{
	allchainsallparams,
	fixednum,
	fixedtime
};


class DATA;
class MODEL;
class POPTREE;
class Output;
class Obsmodel;

class Mcmc
{
public:	
	Mcmc(const Details &details, DATA &data, const MODEL &model, const POPTREE &poptree, const Mpi &mpi, Inputs &inputs, Output &output, Obsmodel &obsmodel);
	
	void run();

private:
	void output_meas(vector <SAMPLE> &opsamp, unsigned int nchain) const;
	void output_param(Output &output, vector <PARAMSAMP> &psamp) const;
	void diagnostic(vector <double> &nac_swap) const;
	void swap(vector <double> &nac_swap, unsigned int samp, unsigned int nchain);

	vector <Chain> chain; 
	
	unsigned int nsamp;                      // The number of MCMC samples
	unsigned int burnin;                     // The number of burnin samples
	
	unsigned int nchaintot;                  // The total number of chains (across all MPI processes);
	unsigned int nchain;                     // The number of chains per core
	
	vector <double> nac_swap;                // Stores the rate at which swaps are accepted
	
	Model_Evidence model_evidence;           // Stores information used to calculate the model evidence
	
	enum proposalsmethod propsmethod;        // Stores the type of proposal method
	
	const Details &details;
	DATA &data;
	const MODEL &model;
	const POPTREE &poptree;
	const Mpi &mpi;
	Output &output;
	Obsmodel &obsmodel;
};

#endif
