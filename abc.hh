#ifndef BEEPMBP__ABC_HH
#define BEEPMBP__ABC_HH

#include "data.hh"
#include "model.hh"

struct ParamSample;
class Data;
class Model;
class AreaTree;
class Output;
class ObservationModel;
class Chain;
struct Sample;

class ABC
{
public:	
  ABC(const Details &details, const Data &data, const Model &model, const AreaTree &areatree, const Mpi &mpi, const Inputs &inputs, const Output &output, const ObservationModel &obsmodel);	
	void smc();
	void mbp();
	
private:
	void mcmc_updates(Generation &gen, vector <Particle> &part, Chain &chain);
	void calculate_cholesky_matrix(const vector <vector <double> > &param_samp);
	void mvn_propose(vector <double> &paramval, double fac);
	double calculate_mixing(const vector <Particle> &part, unsigned int *partcopy) const;
	vector <double> variance_vector(const vector <vector <double> > &param_samp) const;
	vector <vector <double> > covariance_matrix(const vector <vector <double> > &param_samp) const;
	double mvn_prob(const vector<double> &pend, const vector<double> &pstart, double fac) const;
	void exchange_samples_mpi(Generation &gen);
	double next_generation_mpi(vector<Particle> &part, unsigned int *partcopy);
	void results_mpi(const vector <Generation> &generation, const vector <Particle> &part, Chain &chain) const;
	Sample get_sample(const Particle &part, Chain &chain) const;
	double acceptance(double rate) const;
	void calculate_particle_weight(vector <Generation> &generation, double jump);
	unsigned int effective_particle_number(vector <double> w);
	
	vector <vector <double> > inv_Mert_matrix(const vector <vector <double> > &mat) const;

	unsigned int total_time;                 // The time the algorithm is run for
	
	unsigned int N;                          // The number of particles per core
	unsigned int Ntot;                       // The total number of particles
	
	unsigned int G;                          // The total number of generations 
	
	unsigned int nvar;                       // The number of variables which can change
	vector <unsigned int> param_not_fixed;   // A list of all parameters which actualy change 
	
	vector < vector <double> > cholesky_matrix;        // The matrix for Cholesky decompostion
	vector < vector <double> > M; 
	vector < vector <double> > inv_M; 
		
	vector <double> jumpv;
	
	const Details &details;
	const Data &data;
	const Model &model;
	const AreaTree &areatree;
	const Mpi &mpi;
	const Output &output;
	const ObservationModel &obsmodel;
};

#endif
