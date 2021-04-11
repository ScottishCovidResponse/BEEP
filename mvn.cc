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

MVN::MVN(string name_, const vector <unsigned int> &var_, double size_, ParamType type_, MVNType mvntype_)
{
	name = name_;
	var = var_;
	nvar = var.size();
	size = size_;
	type = type_;
	mvntype =  mvntype_;
	number = 1;
}


/// Sets up the multivariate normals sampling
void MVN::setup(const vector <vector <double> > &param_samp)
{
	covariance_matrix(param_samp);
	calculate_cholesky_matrix();
}
		

/// Calculates a covariance matrix from a sets of parameter samples
void MVN::covariance_matrix(const vector <vector <double> > &param_samp)
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
				double val1 = param_samp[k][th1];
				double val2 = param_samp[k][th2];
				av1 += val1;
				av2 += val2;
			}
			av1 /= N; av2 /= N; 
			
			double av12 = 0;
			for(auto k = 0u; k < N; k++){
				double val1 = param_samp[k][th1] - av1;
				double val2 = param_samp[k][th2] - av2;
				av12 += val1*val2;
			}		
			M[i1][i2] = av12/N;
		
			if(std::isnan(M[i1][i2])){ cout << N << "N\n"; emsgEC("MVN",10);}
			if(i1 == i2 && M[i1][i2] < 0) emsgEC("MVN",11);
		}
	}
	
	if(false){
		cout << endl << endl;
		for(auto i1 = 0u; i1 < nvar; i1++){
			for(auto i2 = 0u; i2 < nvar; i2++){
				cout << M[i1][i2] << " ";
			}
			cout << "M\n";
		}
	}
}


/// Calculates a lower diagonal matrix used in Cholesky decomposition
void MVN::calculate_cholesky_matrix()
{	
	auto fl=0u;
	do{
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

		double Lch[nvar][nvar];
		double Tch[nvar][nvar];
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
		
		fl = 0;
		for(auto v1 = 0u; v1 < nvar; v1++){
			for(auto v2 = 0u; v2 < nvar; v2++){
				if(std::isnan(cholesky_matrix[v1][v2])) fl = 1;
			}
		}
		
		if(fl == 1){    // Reduces correlations to allow for convergence
			for(auto v1 = 0u; v1 < nvar; v1++){
				for(auto v2 = 0u; v2 < nvar; v2++){
					if(v1 != v2) M[v1][v2] /= 2;
				}
			}		
			
			if(false){
				for(auto v1 = 0u; v1 < nvar; v1++){
					for(auto v2 = 0u; v2 < nvar; v2++){
						cout << M[v1][v2] << " ";
					}
					cout << "M\n";
				}	
			}
			
			int core; MPI_Comm_rank(MPI_COMM_WORLD,&core);
			if(core == 0) cout << "Cholesky converegence" << endl;
		}
	}while(fl == 1);
	
	inv_M = invert_matrix(M);
}


/// The log probability of drawing a vector from a MVN distribution
double MVN::probability(const vector<double> &pend,const vector<double> &pstart) const
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
	
	//if(mvn.type == COVAR_PARAM2) shift_covariate_reff(param_prop);
	
	if(model.inbounds(param_prop) == false) return FAIL;
	else return SUCCESS;
}


/// Gets the reverese probability
double MVN::get_probfi(vector <double> &param_prop, const vector <double> &paramval, const Model &model)
{
	auto mean = langevin_shift(paramval,model);		
	return probability(param_prop,mean);
}	


/// Makes a joint proposal to sigma and regional effects
Status MVN::sigma_propose(vector <double> &param_prop, const vector <double> &paramval, const Model &model)
{
	param_prop = propose(paramval);
	
	double fac = param_prop[model.sigma_param]/paramval[model.sigma_param]; 
	for(auto th : model.region_effect_param) param_prop[th] *= fac;
	
	if(model.inbounds(param_prop) == false) return FAIL;
	else return SUCCESS;
}


/// Perform the Metropolis-Hastings acceptance probability
Status MVN::MH(double al, ParamUpdate pup)
{
	if(al == -1) nbo++;
	ntr++;
	if(ran() < al){
		if(pup == FAST_UPDATE) size *= fac_up_fast;
		if(pup == SLOW_UPDATE) size *= fac_up;
		if(size > sizemax) size = sizemax;
		nac++; 
		return SUCCESS;
	}
	if(pup == FAST_UPDATE) size *= fac_down_fast;
	if(pup == SLOW_UPDATE) size *= fac_down;
	if(size < sizemin) size = sizemin;
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
			for(auto v = 0u; v < nvar; v++) cout << model.param[var[v]].name << " " << paramv[var[v]] << " " << vec[var[v]] <<" \n";
		}
	}
	
	return vec;
}


/// Outputs the covariance matrix
void MVN::output_M(string file)
{
	ofstream Mpl(file.c_str());

	for(auto i1 = 0u; i1 < nvar; i1++){
		for(auto i2 = 0u; i2 < nvar; i2++){
			Mpl << M[i1][i2] << " ";
		}
		Mpl << "M\n";
	}		
	
	Mpl << "Corelations:\n";
	for(auto i1 = 0u; i1 < nvar; i1++){
		for(auto i2 = 0u; i2 < nvar; i2++){
			if(i1 == i2) Mpl  << "--- ";
			else Mpl << M[i1][i2]/sqrt(M[i1][i1]*M[i2][i2]) << " ";
		}
		Mpl << "M\n";
	}		
}
