#ifndef BEEPMBP__MBP_HH
#define BEEPMBP__MBP_HH

using namespace std;

#include "output.hh"
#include "chain.hh"

enum class proposalsmethod
{
	allchainsallparams,
	fixednum,
	fixedtime
};

class Mcmc
{
public:	
	Mcmc(Details &details, DATA &data, MODEL &model, POPTREE &poptree, Mpi &mpi, Inputs &inputs, Mode mode, bool verbose);
	
	void run(enum proposalsmethod propmethod);

private:
	void output_meas(vector <SAMPLE> &opsamp) const;
	void output_param(Output &output, vector <vector <double> > &Listore, vector <PARAMSAMP> &psamp) const;
	void diagnostic(vector <vector <double> > &Listore, vector <double> &invTstore, vector <double> &nac_swap) const;
	void swap(vector <double> &nac_swap, unsigned int samp);
	double calcME(vector <vector <double> > &Listore,vector <double> &invTstore) const;

	vector <Chain> chain; 
	
	unsigned int nsamp;                      // The number of MCMC samples
	unsigned int burnin;                     // The number of burnin samples
	unsigned int quench;                     // The number of samples over which the system is quenched
	
	unsigned int nchaintot;                  // The total number of chains (across all MPI processes);
	unsigned int nchain;                     // The number of chains per core
	
	double invTmin, invTmax;                 // The minimum and maximum inverse tenperatures that get run at
	
	Details &details;
	DATA &data;
	MODEL &model;
	POPTREE &poptree;
	Mpi &mpi;
};

#endif
