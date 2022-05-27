/// Store infomation about multivariate normal distributions used to generate proposals 

#include <fstream>
#include <iostream>
#include <sstream>
#include <algorithm>  
#include "stdlib.h"
#include "math.h"
#include "assert.h"

using namespace std;

#include "mvn.hh"
#include "model.hh"
#include "matrix.hh"

MVN::MVN(string name_, const vector <unsigned int> &var_, double size_, ParamType type_, MVNType mvntype_)
{
	name = name_;
	var = var_;
	nvar = var.size();
	size = size_;
	type = type_;
	mvntype =  mvntype_;
	number = 1;
	num_updates = num_updates_min;
}


/// Sets up multivariate normal sampling based on a list of parameter samples
void MVN::setup(const vector <ParamSample> &param_samp)
{
	covariance_matrix(param_samp);
	calculate_cholesky_matrix();
}

	
/// Sets up the multivariate normal sampling based on a list of weighted samples
void MVN::setup(const vector <ParamSample> &param_samp, const vector <double> &w)
{
	covariance_matrix(param_samp,w);
	calculate_cholesky_matrix();
}
	

/// Calculates a covariance matrix from a sets of parameter samples
void MVN::covariance_matrix(const vector <ParamSample> &param_samp)
{
	auto N = param_samp.size();

	M.resize(nvar);
	for(auto i1 = 0u; i1 < nvar; i1++){
		M[i1].resize(nvar);
		auto th1 = var[i1]; 
		for(auto i2 = 0u; i2 < nvar; i2++){
			auto th2 = var[i2];
			
			double av1 = 0, av2 = 0;
			for(auto k = 0u; k < N; k++){
				double val1 = param_samp[k].paramval[th1];
				double val2 = param_samp[k].paramval[th2];
				av1 += val1;
				av2 += val2;
			}
			av1 /= N; av2 /= N; 
			
			double av12 = 0;
			for(auto k = 0u; k < N; k++){
				double val1 = param_samp[k].paramval[th1] - av1;
				double val2 = param_samp[k].paramval[th2] - av2;
				av12 += val1*val2;
			}		
			M[i1][i2] = av12/N;
		
			if(std::isnan(M[i1][i2])) emsgEC("MVN",1);
			if(i1 == i2 && M[i1][i2] < 0) emsgEC("MVN",2);
		}
	}
	
	if(false){
		cout << endl << endl;
		for(auto i1 = 0u; i1 < nvar; i1++){
			for(auto i2 = 0u; i2 < nvar; i2++){
				cout << M[i1][i2] << " ";
			}
			cout << "M covariance matrix" << endl;
		}
	}
}


/// Directly sets the covariance matrix
void MVN::set_covariance_matrix(const vector < vector <double> > &_M)
{
	if(_M.size() != nvar) emsgEC("MVN",58);
	M = _M;
	calculate_cholesky_matrix();
}


/// Calculates a covariance matrix from a sets of parameter samples
void MVN::covariance_matrix(const vector <ParamSample> &param_samp, const vector <double> &w)
{
	auto N = param_samp.size();

	M.resize(nvar);
	for(auto i1 = 0u; i1 < nvar; i1++){
		M[i1].resize(nvar);
		auto th1 = var[i1]; 
		for(auto i2 = 0u; i2 < nvar; i2++){
			auto th2 = var[i2];
			
			double av1 = 0, av2 = 0, wsum = 0;
			for(auto k = 0u; k < N; k++){
				double val1 = param_samp[k].paramval[th1];
				double val2 = param_samp[k].paramval[th2];
				av1 += val1*w[k];
				av2 += val2*w[k];
				wsum += w[k];
			}
			av1 /= wsum; av2 /= wsum; 
			
			double av12 = 0;
			for(auto k = 0u; k < N; k++){
				double val1 = param_samp[k].paramval[th1] - av1;
				double val2 = param_samp[k].paramval[th2] - av2;
				av12 += val1*val2*w[k];
			}		
			M[i1][i2] = av12/wsum;
		
			if(std::isnan(M[i1][i2])) emsgEC("MVN",3);
			if(i1 == i2 && M[i1][i2] < 0) emsgEC("MVN",4);
		}
	}
	
	if(false){
		cout << endl << endl;
		for(auto i1 = 0u; i1 < nvar; i1++){
			for(auto i2 = 0u; i2 < nvar; i2++){
				cout << M[i1][i2] << " ";
			}
			cout << "M" << endl;
		}
	}
}


/// Calculates a lower diagonal matrix used in Cholesky decomposition
void MVN::calculate_cholesky_matrix()
{	
	auto fail = 0u;

	auto fl=0u;
	do{
		if(true){ // Algorithm 1
			cholesky_matrix.resize(nvar);
			for(auto v1 = 0u; v1 < nvar; v1++) cholesky_matrix[v1].resize(nvar);
			
			for(auto i = 0u; i < nvar; i++) {
				for(auto j = 0u; j <= i; j++) {
					auto sum = 0.0; for(auto k = 0u; k < j; k++) sum += cholesky_matrix[i][k]*cholesky_matrix[j][k];

					if (i == j) cholesky_matrix[i][j] = sqrt(M[i][i] - sum);
					else cholesky_matrix[i][j] = (1.0/cholesky_matrix[j][j]*(M[i][j] - sum));
				}
			}
		}
		else{    // Algorithm 2
			vector <vector <double> > A;
			A.resize(nvar);
			cholesky_matrix.resize(nvar);
			for(auto v1 = 0u; v1 < nvar; v1++){
				A[v1].resize(nvar);
				cholesky_matrix[v1].resize(nvar);
				for(auto v2 = 0u; v2 < nvar; v2++){
					A[v1][v2] = M[v1][v2];
					if(v1 == v2) cholesky_matrix[v1][v2] = 1; else cholesky_matrix[v1][v2] = 0;
				}
			}

			vector < vector <double> > Lch, Tch;
			Lch.resize(nvar); Tch.resize(nvar);
			for(auto i = 0u; i < nvar; i++){ Lch[i].resize(nvar); Tch[i].resize(nvar);}
		
			for(auto i = 0u; i < nvar; i++){
				for(auto v1 = 0u; v1 < nvar; v1++){
					for(auto v2 = 0u; v2 < nvar; v2++){
						if(v1 == v2) Lch[v1][v2] = 1; else Lch[v1][v2] = 0;
					}
				}

				double aii = A[i][i];
				Lch[i][i] = sqrt(aii);
				for(auto j = i+1; j < nvar; j++){
					Lch[j][i] = A[j][i]/sqrt(aii);
				}

				for(auto ii = i+1; ii < nvar; ii++){
					for(auto jj = i+1; jj < nvar; jj++){
						A[jj][ii] -= A[ii][i]*A[jj][i]/aii;
					}
				}
				A[i][i] = 1;
				for(auto j = i+1; j < nvar; j++){ A[j][i] = 0; A[i][j] = 0;}

				for(auto v1 = 0u; v1 < nvar; v1++){
					for(auto v2 = 0u; v2 < nvar; v2++){
						double sum = 0u; for(auto ii = 0u; ii < nvar; ii++) sum += cholesky_matrix[v1][ii]*Lch[ii][v2];
					
						Tch[v1][v2] = sum;
					}
				}

				for(auto v1 = 0u; v1 < nvar; v1++){
					for(auto v2 = 0u; v2 < nvar; v2++){
						cholesky_matrix[v1][v2] = Tch[v1][v2];
					}
				}
			}
		}
		
		fl = 0;
		for(auto v1 = 0u; v1 < nvar; v1++){
			for(auto v2 = 0u; v2 < nvar; v2++){
				if(std::isnan(cholesky_matrix[v1][v2])) fl = 1;
			}
		}
		
		if(fl == 1){                                                  // Reduces correlations to allow for convergence
			for(auto v1 = 0u; v1 < nvar; v1++){
				for(auto v2 = 0u; v2 < nvar; v2++){
					if(v1 != v2) M[v1][v2] /= 2;
				}
			}		
			
			int core; MPI_Comm_rank(MPI_COMM_WORLD,&core);
			if(core == 0) cout << "Cholesky convergence" << endl;
			
			
			fail++;
			if(fail == 10 && core == 0){
				for(auto v1 = 0u; v1 < nvar; v1++){
					for(auto v2 = 0u; v2 < nvar; v2++){
						cout << M[v1][v2] << " ";
					}
					cout << "M" << endl;
				}	
				emsg("Could not get Cholesky convergence");
			}
		}
	}while(fl == 1);

	inv_M = invert_matrix(M);
}


/// The log probability of drawing a vector from a MVN distribution
double MVN::probability(const vector<double> &pend, const vector<double> &pstart) const
{
	double sum = 0;
	for(auto v1 = 0u; v1 < nvar; v1++){
		for(auto v2 = 0u; v2 < nvar; v2++){
			double val1 = pend[var[v1]] - pstart[var[v1]];
			double val2 = pend[var[v2]] - pstart[var[v2]];
			sum += (val1/size)*inv_M[v1][v2]*(val2/size);
		}
	}	
	return -0.5*sum;
}


/// Generates a proposed set of parameters from a MVN distribution
Status MVN::propose_langevin(vector <double> &param_prop, const vector <double> &paramval, double &probif, const Model &model) const
{
	auto mean = langevin_shift(paramval,model);
	param_prop = propose(mean);
	probif = probability(param_prop,mean);
	if(model.inbounds(param_prop) == false) return FAIL;
	else return SUCCESS;
}


/// Gets the reverese probability
double MVN::get_probfi(const vector <double> &param_prop, const vector <double> &paramval, const Model &model)
{
	auto mean = langevin_shift(paramval,model);		
	return probability(param_prop,mean);
}	


/// Makes a joint proposal to sigma and regional effects
Status MVN::sigma_propose(vector <double> &param_prop, const vector <double> &paramval, const Model &model)
{
	param_prop = propose(paramval);
	
	auto th = model.region_effect.sigma_param;
	double fac = param_prop[th]/paramval[th]; 
	for(auto th : model.region_effect.param_list) param_prop[th] *= fac;
	
	if(model.inbounds(param_prop) == false) return FAIL;
	else return SUCCESS;
}


/// Perform the Metropolis-Hastings acceptance probability
Status MVN::MH(double al, const ParamUpdate pup)
{
	if(al == -1){ nbo++; al = 0;}
	
	if(al > 1) al = 1;
	if(pup == FAST_UPDATE) update(size,al,eta_fast);
	if(pup == SLOW_UPDATE) update(size,al,eta);
	if(size > sizemax) size = sizemax; 
	
	ntr++; nac += al;
	if(ran() < al) return SUCCESS;
	return FAIL;
}		
	

/// Perform the Metropolis-Hastings acceptance probability
Status MVN::MH_PMCMC(double al, const double self_ac, const ParamUpdate pup)
{
	if(al == -1){ nbo++; al = 0;}
		
	auto target_ac_rate = self_ac/2;
	if(pmcmc_start_param == true) target_ac_rate = 0.3;
	
	if(al > 1) al = 1;
	if(pup == FAST_UPDATE) update(size,al,eta_pmcmc_fast,target_ac_rate);
	if(pup == SLOW_UPDATE) update(size,al,eta_pmcmc,target_ac_rate);
	
	auto sizemin = 0.5*2.4*2.4/nvar;
	if(sizemin > 0.3) sizemin = 0.3;

	if(size > sizemax) size = sizemax;	
	if(pup != FAST_UPDATE){
		if(size < sizemin) size = sizemin;
	}
	
	ntr++; nac += al;

	if(ran() < al) return SUCCESS;
	
	return FAIL;
}		
		
		
/// Generates a proposed set of parameters from a MVN distribution
vector <double> MVN::propose(const vector <double> &paramval) const
{
	auto vec = paramval;

	double norm[nvar];	
	for(auto v = 0u; v < nvar; v++) norm[v] = normal_sample(0,1);

	for(auto v = 0u; v < nvar; v++){
		double dva = 0; for(auto v2 = 0u; v2 <= v; v2++) dva += cholesky_matrix[v][v2]*norm[v2];

		auto th = var[v];
		vec[th] += size*dva;
	}

	return vec;
}


/// Calculates the shift in parameter values to account for gradients in the prior probability
vector <double> MVN::langevin_shift(const vector <double> &paramv, const Model &model) const 
{
	auto vec = paramv;
	
	if(langevin_on == true){			
		vector <double> gr_Pr(nvar);
	
		for(auto v = 0u; v < nvar; v++) gr_Pr[v] = model.dPr_dth(var[v],vec);

		for(auto v = 0u; v < nvar; v++){	
			auto sum = 0.0; for(auto vv = 0u; vv < nvar; vv++) sum += size*size*M[v][vv]*gr_Pr[vv];
	
			vec[var[v]] += 0.5*sum;
		}
		
		if(false){
			for(auto v = 0u; v < nvar; v++) cout << model.param[var[v]].name << " " << paramv[var[v]] << " " << vec[var[v]] << endl;
		}
	}
	
	return vec;
}


/// Adds a vector giving the difference between the initial and proposed states 
void MVN::add_prop_vec(vector <double> &paramv_init, vector <double> &paramv_prop)
{
	vector <double> vec(nvar);
	for(auto v = 0u; v < nvar; v++){
		vec[v] = paramv_prop[var[v]] - paramv_init[var[v]];
	}
	
	prop_vec.push_back(vec);
	
	for(auto v = 0u; v < nvar; v++){
		vec[v] = paramv_prop[var[v]];
	}
	
	after_prop_vec.push_back(vec);
}


/// Samples from a MVN distribution with specified mean and covariance matrix
vector <double> MVN::sample(const vector <double> &mean, const vector < vector <double> > &_M)
{
	M = _M;
	nvar = mean.size();
	var.clear(); for(auto i = 0u; i < nvar; i++) var.push_back(i);
	size = 1;
	calculate_cholesky_matrix();
	
	return propose(mean);
}


/// Checks the sampling function is working correctly
void MVN::sample_check(const vector <double> &mean, const vector < vector <double> > &_M)
{
	auto N = mean.size();
	
	vector <double> mu(N);
	vector < vector <double> > covar;
	
	covar.resize(N);
	for(auto j = 0u; j < N; j++){
		mu[j] = 0;
		covar[j].resize(N); for(auto i = 0u; i < N; i++) covar[j][i] = 0;
	}
	
	const auto loopmax = 1000000u;
	for(auto loop = 0u; loop < loopmax; loop++){ 
		auto p = sample(mean,_M);
		for(auto j = 0u; j < N; j++){
			mu[j] += p[j];
			for(auto i = 0u; i < N; i++) covar[j][i] += p[j]*p[i];
		}
	}
	
	for(auto j = 0u; j < N; j++) mu[j] /= loopmax;
		
	for(auto j = 0u; j < N; j++){
		for(auto i = 0u; i < N; i++) covar[j][i] = covar[j][i]/loopmax - mu[j]*mu[i];
	}
	
	print_matrix("M",M);
	print_matrix("covar",covar);
	emsg("Done");
}
			

/// Outputs the covariance matrix
void MVN::output_M(const string file)
{
	ofstream Mpl(file.c_str());

	for(auto i1 = 0u; i1 < nvar; i1++){
		for(auto i2 = 0u; i2 < nvar; i2++){
			Mpl << M[i1][i2] << " ";
		}
		Mpl << "M" << endl;
	}		
	
	Mpl << "Corelations:" << endl;
	for(auto i1 = 0u; i1 < nvar; i1++){
		for(auto i2 = 0u; i2 < nvar; i2++){
			if(i1 == i2) Mpl  << "--- ";
			else Mpl << M[i1][i2]/sqrt(M[i1][i1]*M[i2][i2]) << " ";
		}
		Mpl << "M" << endl;
	}		
}
