#ifndef BEEPMBP__MVN_HH
#define BEEPMBP__MVN_HH

#include "struct.hh"

class MVN                                      // Multivariate normal
{
public:
	string name;                                 // The name of the multivariate normal
	unsigned int ntr, nac, nbo;                  // Keeps track of acceptance of proposals
	double size;                                 // The size of the proposal
	unsigned int number;                         // Used to give the number of proposals for sim-only approach
	ParamType type;                              // The paramter type of the proposal
	MVNType mvntype;                             // Single or multiple variable update
	 
	double ac_rate, bo_rate;                     // Acceptance rates (across all MPI processes)
			
	MVN(string name, const vector <unsigned int> &var_, double size_, ParamType type_, MVNType mvntype_);
	void setup(const vector <vector <double> > &param_samp);

	vector <double> propose(const vector <double> &paramval) const;
	double probability(const vector<double> &pend,const vector<double> &pstart) const;
	vector <double> langevin_shift(const vector <double> &paramv, const Model &model) const;
	Status propose_langevin(vector <double> &param_propose, const vector <double> &paramval, double &probif, const Model &model) const;
	double get_probfi(vector <double> &param_propose, const vector <double> &paramval, const Model &model);
	Status MH(double al, ParamUpdate pup);
	Status sigma_propose(vector <double> &param_propose, const vector <double> &paramval, const Model &model);
	
private:
//	public:
	unsigned int nvar;                           // The number of variables which can change
	vector <unsigned int> var;                   // The variables of the MVN

	vector < vector <double> > M; 	             // The covariance matrix
	vector < vector <double> > cholesky_matrix;  // The matrix for Cholesky decompostion
	vector < vector <double> > inv_M; 	         // The inverse of the covarience matrixc
	
	void covariance_matrix(const vector <vector <double> > &param_samp);
	void calculate_cholesky_matrix();
	
	void output_M(string file);	
};

#endif
