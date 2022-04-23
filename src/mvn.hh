#ifndef BEEPMBP__MVN_HH
#define BEEPMBP__MVN_HH

#include "struct.hh"

class MVN                                      // Multivariate normal
{
public:
	string name;                                 // The name of the multivariate normal
	double ntr, nac, nbo;                        // Keeps track of acceptance of proposals
	double size;                                 // The size of the proposal
	unsigned int number;                         // Used to give the number of proposals for sim-only approach
	ParamType type;                              // The paramter type of the proposal
	MVNType mvntype;                             // Single or multiple variable update
	 
	double num_updates;                          // The number of times performed during an update      
	
 	unsigned int nvar;                           // The number of variables which can change
	vector <unsigned int> var;                   // The variables of the MVN

	vector < vector <double> > prop_vec;         // Stores samples of the successful proposal vector
	vector < vector <double> > after_prop_vec; 	 // Stores samples of the successful proposal vector

	double ac_rate, bo_rate;                     // Acceptance rates (across all MPI processes)
			
	MVN(string name, const vector <unsigned int> &var_, double size_, ParamType type_, MVNType mvntype_);
	void setup(const vector <ParamSample> &param_samp);
	void setup(const vector <ParamSample> &param_samp, const vector <double> &w);

	vector <double> propose(const vector <double> &paramval) const;
	double probability(const vector<double> &pend,const vector<double> &pstart) const;
	vector <double> langevin_shift(const vector <double> &paramv, const Model &model) const;
	Status propose_langevin(vector <double> &param_propose, const vector <double> &paramval, double &probif, const Model &model) const;
	double get_probfi(const vector <double> &param_propose, const vector <double> &paramval, const Model &model);
	Status MH(double al, const ParamUpdate pup);
	Status MH_PMCMC(double al, const double self_ac, const ParamUpdate pup);
	Status sigma_propose(vector <double> &param_propose, const vector <double> &paramval, const Model &model);
	void add_prop_vec(vector <double> &paramv_init, vector <double> &paramv_prop);
	void set_covariance_matrix(const vector < vector <double> > &M);
	vector <double> sample(const vector <double> &mean, const vector < vector <double> > &M);
	void sample_check(const vector <double> &mean, const vector < vector <double> > &_M);
		
private:
	vector < vector <double> > M; 	             // The covariance matrix
	vector < vector <double> > cholesky_matrix;  // The matrix for Cholesky decompostion
	vector < vector <double> > inv_M; 	         // The inverse of the covarience matrixc
	
	void covariance_matrix(const vector <ParamSample> &param_samp);
	void covariance_matrix(const vector <ParamSample> &param_samp, const vector <double> &w);
	void calculate_cholesky_matrix();

	void output_M(const string file);	
};

#endif
