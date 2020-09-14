#ifndef BEEPMBP__JUMP_HH
#define BEEPMBP__JUMP_HH

#include <vector>
#include <iostream>
#include <algorithm>
#include <fstream> 

using namespace std;

#include "consts.hh"
#include "model.hh"

class Jump                                     // Information about kernal for jumping in parameter space
{
public:
	vector <float> mbp;                          // The size of jumps in parameter space
	vector <unsigned int> mbp_ntr, mbp_nac;      // The number of jumps tried and accepted
	
	vector <float> stand;                        // The size of jumps in parameter space (fixed event sequence)
	vector <unsigned int> stand_ntr, stand_nac;  // The number of jumps tried and accepted

	float naddrem;                             	 // The size of adding and removing events
	unsigned int standev_ntr, standev_nac;    
	
	//
	
	unsigned int nvar;                           // The number of variables which can change
	vector <unsigned int> param_not_fixed;       // A list of all parameters which actualy change 

	vector < vector <double> > cholesky_matrix;        // The matrix for Cholesky decompostion
	vector < vector <double> > M; 
	vector < vector <double> > inv_M; 
		
	vector <double> jumpv;
	
	void init(const vector <double> &paramv, const Model &model);
	vector <double> mbp_prop(const vector <double> &p, unsigned int th);
	
	void setburnin(unsigned int samp, unsigned int _burnin);
	void mbp_accept(unsigned int th);
	void mbp_reject(unsigned int th);
	void stand_accept(unsigned int th);
	void stand_reject(unsigned int th);
	void standev_accept();
	void standev_reject();
	vector <double> variance_vector(const vector <vector <double> > &param_samp) const;
	vector <vector <double> > covariance_matrix(const vector <vector <double> > &param_samp) const;
	void calculate_cholesky_matrix(const vector <vector <double> > &param_samp);
	vector <vector <double> > inv_Mert_matrix(const vector <vector <double> > &mat) const;
	double mvn_prob(const vector<double> &pend,const vector<double> &pstart, double fac) const;
	void mvn_propose(vector <double> &paramval, double fac);
	void output_M(string file);

private:
	Alter alter;
};

#endif
