/// This provides all the functions relating to generating MCMC proposals in parameter space

#include <assert.h>
#include <math.h>

using namespace std;

#include "jump.hh"
#include "utils.hh"

//Jump::Jump(const Details &details, const Model &model) : details(details), model(model){}

/// Quantities in jump are intialised
void Jump::init(const vector <double> &paramv, const Model &model)
{
	auto nparam = paramv.size();
	mbp.resize(nparam); mbp_ntr.resize(nparam); mbp_nac.resize(nparam);     
	stand.resize(nparam); stand_ntr.resize(nparam); stand_nac.resize(nparam);
	for(auto th = 0u; th < nparam; th++){
		mbp[th] = paramv[th]/2; if(mbp[th] == 0) mbp[th] = 0.1;
		mbp_ntr[th] = 0; mbp_nac[th] = 0;
		
		stand[th] = paramv[th]/10; if(stand[th] == 0) stand[th] = 0.1;
		stand_ntr[th] = 0; stand_nac[th] = 0;
	}
	
	naddrem = 20;
	standev_ntr = 0; standev_nac = 0;
	
	param_not_fixed.clear();                                              // Finds the list of model parameteres that change
	for(auto th = 0u; th < model.param.size(); th++){
		if(model.param[th].min != model.param[th].max) param_not_fixed.push_back(th);
	}
	nvar = param_not_fixed.size();
	
	jumpv.resize(nvar); for(auto v = 0u; v < nvar; v++) jumpv[v] = 1;
}

/// Makes a proposal by chaning a given parameter from given parameter set
vector <double> Jump::mbp_prop(const vector <double> &p, unsigned int th)
{
	vector <double> paramv = p;
	paramv[th] += normal_sample(0,mbp[th]);  
	
	return paramv;  
}

void Jump::setburnin(unsigned int samp, unsigned int burnin)
{
	if(samp < burnin){
		if(samp < 50) alter = FAST;
		else alter = SLOW;
	}
	else alter = NONE;
}

/// A MBP is accepted
void Jump::mbp_accept(unsigned int th)
{
	mbp_ntr[th]++;
	mbp_nac[th]++;
	switch(alter){
		case FAST: mbp[th] *= 2; break;
		case SLOW: mbp[th] *= 1.1; break;
		case NONE: break;
	}
}

/// A MBP is rejected
void Jump::mbp_reject(unsigned int th)
{
	mbp_ntr[th]++;
	switch(alter){
		case FAST: mbp[th] *= 0.5; break;
		case SLOW: mbp[th] *= 0.95; break;
		case NONE: break;
	}
}

/// A standand proposal is accepted
void Jump::stand_accept(unsigned int th)
{
	stand_ntr[th]++;
	stand_nac[th]++;
	switch(alter){
		case FAST: stand[th] *= 1.05; break;
		case SLOW: stand[th] *= 1.01; break;
		case NONE: break;
	}
}

/// A standard proposal is rejected
void Jump::stand_reject(unsigned int th)
{
	stand_ntr[th]++;
	switch(alter){
		case FAST: stand[th] *= 0.975; break;
		case SLOW: stand[th] *= 0.995; break;
		case NONE: break;
	}
}

void Jump::standev_accept()
{
	standev_ntr++;
	standev_nac++;
	switch(alter){
		case FAST: naddrem *= 1.05; break;
		case SLOW: naddrem *= 1.05; break;
		case NONE: break;
	}
}

void Jump::standev_reject()
{
	standev_ntr++;
	switch(alter){
		case FAST: naddrem *= 0.95; break;
		case SLOW: naddrem *= 0.95;break;
		case NONE: break;
	}
	if(naddrem < 1) naddrem = 1;
}



/// Generates a covariance matrix from a sets of parameter samples
vector <double> Jump::variance_vector(const vector <vector <double> > &param_samp) const 
{
	vector <double> vec;
	auto N = param_samp.size();                             // Generates the covariance matrix

	vec.resize(nvar);
	for(auto i = 0u; i < nvar; i++){
		auto th = param_not_fixed[i]; 
	
		double av = 0;
		for(auto k = 0u; k < N; k++) av += param_samp[k][th];
		av /= N;
		
		double av2 = 0;
		for(auto k = 0u; k < N; k++){
			double val = param_samp[k][th] - av;
			av2 += val*val;
		}		
		vec[i] = av2/N;
	}
	
	return vec;
}

/// Generates a covariance matrix from a sets of parameter samples
vector <vector <double> > Jump::covariance_matrix(const vector <vector <double> > &param_samp) const 
{
	vector <vector <double> > M;
	
	auto N = param_samp.size();                             // Generates the covariance matrix

	M.resize(nvar);
	for(auto i1 = 0u; i1 < nvar; i1++){
		M[i1].resize(nvar);
		auto th1 = param_not_fixed[i1]; 
		for(auto i2 = 0u; i2 < nvar; i2++){
			auto th2 = param_not_fixed[i2];
			
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
		
			if(std::isnan(M[i1][i2])) emsg("not");
			if(i1 == i2 && M[i1][i2] < 0) emsg("Negative");
		}
	}
	
	return M;
}

/// Calculates a lower diagonal matrix used in Cholesky decomposition
void Jump::calculate_cholesky_matrix(const vector <vector <double> > &param_samp)
{
	M = covariance_matrix(param_samp);
	
	/*
	cout << endl << endl;
	for(auto i1 = 0u; i1 < nvar; i1++){
		for(auto i2 = 0u; i2 < nvar; i2++){
			cout << M[i1][i2] << " ";
		}
		cout << "M\n";
	}
	*/
	//emsg("P");
	
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
			
			/*
			for(auto v1 = 0u; v1 < nvar; v1++){
				for(auto v2 = 0u; v2 < nvar; v2++){
					cout << M[v1][v2] << " ";
				}
				cout << "M\n";
			}	
			*/			
			cout << "Cholesky converegence" << endl;
		}
	}while(fl == 1);
	
	inv_M = inv_Mert_matrix(M);
}

/// Inverts a matrix
vector <vector <double> > Jump::inv_Mert_matrix(const vector <vector <double> > &mat) const    
{
	unsigned int nvar = mat.size();
	vector <vector <double> > inv_M;
	
	double A2[nvar][nvar];

	inv_M.resize(nvar);
  for(auto i = 0u; i < nvar; i++){
		inv_M[i].resize(nvar);
    for(auto j = 0u; j < nvar; j++){
      A2[i][j] = mat[i][j];
      if(i == j) inv_M[i][j] = 1; else inv_M[i][j] = 0;
    }
  }

  for(auto ii = 0u; ii < nvar; ii++){
    double r = A2[ii][ii];
    for(auto i = 0u; i < nvar; i++){
      A2[ii][i] /= r; inv_M[ii][i] /= r; 
    }

    for(auto jj = ii+1; jj < nvar; jj++){
      double r = A2[jj][ii];
      for(auto i = 0u; i < nvar; i++){ 
        A2[jj][i] -= r*A2[ii][i];
        inv_M[jj][i] -= r*inv_M[ii][i];
      }
    }
  }

  for(int ii = nvar-1; ii > 0; ii--){
    for(int jj = ii-1; jj >= 0; jj--){
      double r = A2[jj][ii];
      for(auto i = 0u; i < nvar; i++){ 
        A2[jj][i] -= r*A2[ii][i];
        inv_M[jj][i] -= r*inv_M[ii][i];
      }
    }
  }

	// check inv_Merse
	/*
	for(auto j = 0u; j < nvar; j++){
		for(auto i = 0u; i < nvar; i++){
			double sum = 0; for(auto ii = 0u; ii < nvar; ii++) sum += mat[j][ii]*inv_M[ii][i];
			cout << sum << "\t"; 
		}
		cout << " kk\n";
	}
	*/
	
	return inv_M;
}

/// The probability of drawing a vector from a MVN distribution
double Jump::mvn_prob(const vector<double> &pend,const vector<double> &pstart, double fac) const
{
	double sum = 0;
	for(auto v1 = 0u; v1 < nvar; v1++){
		for(auto v2 = 0u; v2 < nvar; v2++){
			double val1 = pend[param_not_fixed[v1]] - pstart[param_not_fixed[v1]];
			double val2 = pend[param_not_fixed[v2]] - pstart[param_not_fixed[v2]];
			sum += (val1/fac)*inv_M[v1][v2]*(val2/fac);
		}
	}	
	return exp(-0.5*sum);
}

/// Generates a proposed set of parameters from a MVN distribution
void Jump::mvn_propose(vector <double> &paramval, double fac)
{
	double norm[nvar];	
	for(auto v = 0u; v < nvar; v++) norm[v] = normal_sample(0,1);
	
  for(auto v = 0u; v < nvar; v++){
		double dva = 0; for(auto v2 = 0u; v2 <= v; v2++) dva += cholesky_matrix[v][v2]*norm[v2];

		auto th = param_not_fixed[v];
		paramval[th] += fac*dva;
	}
}

void Jump::output_M(string file)
{
	ofstream Mpl(file.c_str());

	for(auto i1 = 0u; i1 < nvar; i1++){
		for(auto i2 = 0u; i2 < nvar; i2++){
			Mpl << M[i1][i2] << " ";
		}
		Mpl << "M\n";
	}		
	
	Mpl << "coorela\n";
	for(auto i1 = 0u; i1 < nvar; i1++){
		for(auto i2 = 0u; i2 < nvar; i2++){
			if(i1 == i2) Mpl  << "--- ";
			else Mpl << M[i1][i2]/sqrt(M[i1][i1]*M[i2][i2]) << " ";
		}
		Mpl << "M\n";
	}		
}

	