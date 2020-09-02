#ifndef BEEPMBP__ABC_HH
#define BEEPMBP__ABC_HH

#include "data.hh"
#include "model.hh"

struct PARAMSAMP;
class DATA;
class MODEL;
class POPTREE;
class Output;
class Obsmodel;
class Chain;
struct SAMPLE;

class ABC
{
public:	
  ABC(const Details &details, DATA &data, MODEL &model, const POPTREE &poptree, const Mpi &mpi, Inputs &inputs, Output &output, Obsmodel &obsmodel);	
	void smc();
	void mbp();
	
private:
	void mcmc_updates(Generation &gen, vector <Particle> &part, Chain &chain);
	void cholesky(const vector <vector <double> > &param_samp);
	void cholesky_propose(vector <double> &paramval, double fac);
	double mix(const vector <Particle> &part, unsigned int *partcopy) const;
	vector <double> variance_vector(const vector <vector <double> > &param_samp) const;
	vector <vector <double> > covariance_matrix(const vector <vector <double> > &param_samp) const;
	double mvn_prob(const vector<double> &pend, const vector<double> &pstart, double fac) const;
	void exchange_samples_mpi(Generation &gen);
	double next_generation_mpi(vector<Particle> &part, unsigned int *partcopy);
	void results_mpi(const vector <Generation> &generation, const vector <Particle> &part, Chain &chain) const;
	SAMPLE get_sample(const Particle &part, Chain &chain) const;
	double acceptance(double rate) const;
	void calculate_w(vector <Generation> &generation, double jump);
	unsigned int Neff(vector <double> w);
	
	vector <vector <double> > invert_matrix(const vector <vector <double> > &mat) const;

	unsigned int total_time;                 // The time the algorithm is run for
		
	unsigned int nvar;                       // The number of variables which can change
	vector <unsigned int> param_not_fixed;   // A list of all parameters which actualy change 
	
	vector < vector <double> > Zchol;        // The matrix for Cholesky decompostion
	vector < vector <double> > M; 
	vector < vector <double> > inv; 
		
	vector <double> jumpv;
	
	const Details &details;
	DATA &data;
	MODEL &model;
	const POPTREE &poptree;
	const Mpi &mpi;
	Output &output;
	Obsmodel &obsmodel;
};

#endif
