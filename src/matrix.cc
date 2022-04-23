// This inclues all functions related to matrices

#include "Eigen/Dense"

using namespace std;

#include "matrix.hh"

unsigned int n_matrix_store;
double **l, **linv, **u, **uinv, **mat2, **matres;
float **lf, **lfinv, **uf, **ufinv, **matf2, **matresf;

/// Performs LU decomposition to return the log of the determinant
vector < vector <double> > invert_determinant_SIMD(const vector < vector <double> > &M, double &det, Accuracy ac) 
{
	timer[TIME_DETERMINANT].start();
	
	auto n = M.size();  

	det = 0.0;
	auto neg = 0u;
	vector < vector <double> > Minv; Minv.resize(n); 
	
	switch(ac){
		case DOUBLE:
			{
				for(auto i = 0u; i < n; i++){
					auto *ui = u[i], *li = l[i];
					
					for(auto j = i; j < n; j++) l[j][i] = M[j][i] - dot_SIMD(l[j],ui,0,i);
					
					u[i][i] = 1;
					for(auto j = i+1; j < n; j++){
						if(li[i] == 0) emsgEC("Matrix",1);
						u[j][i] = (M[i][j] - dot_SIMD(li,u[j],0,i))/li[i];
					}
				}

				for(auto i = 0u; i < n; i++){
					auto valu = u[i][i]; if(valu < 0){ valu = -valu; neg++;}
					auto vall = l[i][i]; if(vall < 0){ vall = -vall; neg++;}
					det += log(valu*vall);
				}
			
				invert_lower_triangular(l,linv,n);
				invert_lower_triangular(u,uinv,n);		
				transpose(linv,n);
				transpose(uinv,n);
			
				for(auto j = 0u; j < n; j++){
					Minv[j].resize(n);
					for(auto i = 0u; i <= j; i++){
						auto kmin = j; if(i > kmin) kmin = i;
						Minv[j][i] = dot_SIMD(uinv[j],linv[i],kmin,n);
						Minv[i][j] = Minv[j][i];
					}
				}	
			}
			break;
			
		case FLOAT:
			{
				for(auto i = 0u; i < n; i++){
					auto *ufi = uf[i], *lfi = lf[i];
					
					for(auto j = i; j < n; j++) lf[j][i] = M[j][i] - dot_SIMD(lf[j],ufi,0,i);
					
					uf[i][i] = 1;
					for(auto j = i+1; j < n; j++){
						if(lfi[i] == 0) emsgEC("Matrix",2);
						uf[j][i] = (M[i][j] - dot_SIMD(lfi,uf[j],0,i))/lfi[i];
					}
				}

				for(auto i = 0u; i < n; i++){
					auto valu = uf[i][i]; if(valu < 0){ valu = -valu; neg++;}
					auto vall = lf[i][i]; if(vall < 0){ vall = -vall; neg++;}
					det += log(valu*vall);
				}
				
				invert_lower_triangular(lf,lfinv,n);
				invert_lower_triangular(uf,ufinv,n);
				transpose(lfinv,n);
				transpose(ufinv,n);
				
				for(auto j = 0u; j < n; j++){
					Minv[j].resize(n);
					for(auto i = 0u; i <= j; i++){
						auto kmin = j; if(i > kmin) kmin = i;
						Minv[j][i] = dot_SIMD(ufinv[j],lfinv[i],kmin,n);
						Minv[i][j] = Minv[j][i];
					}
				}	
			}
			break;
	}

	return Minv;
}


/// Multiplies a matrix and a vector
vector <double> matrix_mult(const vector < vector <double> > &M, const vector <double> &vec)
{
	timer[TIME_MATRIX_MULT].start();
		
	const auto NY = M.size(), NX = M[0].size();
	vector <double> result(NY);
	
	if(vec.size() != NX) emsgEC("Matrix",3);
	
	for(auto j = 0u; j < NY; j++){
		auto sum = 0.0;
		for(auto i = 0u; i < NX; i++) sum += M[j][i]*vec[i];
		result[j] = sum;
	}
	
	timer[TIME_MATRIX_MULT].stop();
	
	return result;
}


/// Multiplies a matrix and a matrix
vector < vector <double> > matrix_mult(const vector < vector <double> > &M1, const vector < vector <double> > &M2)
{
	timer[TIME_MATRIX_MULT].start();
		
	const auto NY = M1.size(), NX = M1[0].size();
	const auto NY2 = M2.size(), NX2 = M2[0].size();
	vector < vector <double> > result;
	
	if(NX != NY2) emsgEC("Matrix",4);
	
	result.resize(NY);
	for(auto j = 0u; j < NY; j++){
		result[j].resize(NX2);
		for(auto k = 0u; k < NX2; k++){
			auto sum = 0.0;
			for(auto i = 0u; i < NX; i++) sum += M1[j][i]*M2[i][k];
			result[j][k] = sum;
		}
	}	
	timer[TIME_MATRIX_MULT].stop();
	
	return result;
}


/// Multiplies a matrix and a matrix (if M1 is sparse this will speed up)
vector < vector <double> > matrix_mult_sparse(const vector < vector <double> > &M1, const vector < vector <double> > &M2)
{
	const auto NY = M1.size(), NX = M1[0].size();
	const auto NY2 = M2.size(), NX2 = M2[0].size();
	vector < vector <double> > result;
	
	if(NX != NY2) emsgEC("Matrix",5);
	
	result.resize(NY);
	for(auto j = 0u; j < NY; j++){
		result[j].resize(NX2);
		for(auto i = 0u; i < NX2; i++) result[j][i] = 0;
	}
	
	for(auto j = 0u; j < NY; j++){
		for(auto i = 0u; i < NX; i++){
			auto val = M1[j][i];
			if(val != 0){
				auto &res1 = result[j];
				const auto &res2 = M2[i];
				for(auto k = 0u; k < NX2; k++) res1[k] += val*res2[k];
			}
		}
	}
	
	return result;
}


/// Sets the matrix (so that is can operated on later
void set_mat(const vector < vector <double> > &M, const Accuracy ac)
{
	auto NY = M.size(), NX = M[0].size();
	for(auto j = 0u; j < NY; j++){
		switch(ac){
			case DOUBLE: for(auto i = 0u; i < NX; i++) mat2[j][i] = M[j][i]; break;
			case FLOAT: for(auto i = 0u; i < NX; i++) matf2[j][i] = M[j][i]; break;
		}
	}
}
	

// Adds the stored matrix to the covariance matrix
void covar_add_mat(vector < vector <double> > &M, const Accuracy ac)
{
	auto NY = M.size(), NX = M[0].size();
	for(auto j = 0u; j < NY; j++){
		switch(ac){
			case DOUBLE: for(auto i = 0u; i < NX; i++) M[j][i] += mat2[j][i] + mat2[i][j]; break;
			case FLOAT: for(auto i = 0u; i < NX; i++) M[j][i] += matf2[j][i] + matf2[i][j]; break;
		}
	}
}
	
/// Multiplies a matrix and a matrix (if M1 is sparse this will speed up)
void covar_matrix_mult_SIMD(const vector < vector <double> > &M1, const unsigned int NX2, const Accuracy ac)
{
	const auto NY = M1.size(), NX = M1[0].size();
	
	if(NX > n_matrix_store || NY > n_matrix_store || NX2 > n_matrix_store) emsgEC("Matrix",6);
		
	switch(ac){
		case DOUBLE:
			{
				for(auto j = 0u; j < NY; j++){
					for(auto i = 0u; i < NX2; i++) matres[j][i] = 0;
				}
				
				for(auto j = 0u; j < NY; j++){
					for(auto i = 0u; i < NX; i++){
						auto val = M1[j][i];
						if(val != 0) sub_row_SIMD(matres[j],mat2[i],val,0,NX2);
					}
				}

				for(auto j = 0u; j < NY; j++){
					for(auto i = 0u; i < NX2; i++) mat2[j][i] = matres[j][i];
				}	
			}
			break;
			
		case FLOAT:
			{
				for(auto j = 0u; j < NY; j++){
					for(auto i = 0u; i < NX2; i++) matresf[j][i] = 0;
				}
				
				for(auto j = 0u; j < NY; j++){
					for(auto i = 0u; i < NX; i++){
						float val = M1[j][i];
						if(val != 0) sub_row_SIMD(matresf[j],matf2[i],val,0,NX2);
					}
				}

				for(auto j = 0u; j < NY; j++){
					for(auto i = 0u; i < NX2; i++) matf2[j][i] = matresf[j][i];
				}	
			}
			break;
	}
}


/// Multiplies two vectors
double vec_mult(const vector <double> &vec1, const vector <double> &vec2)
{
	if(vec1.size() != vec2.size()) emsgEC("Matrix",8);
	
	auto sum = 0.0;
	for(auto i = 0u; i < vec1.size(); i++) sum += vec1[i]*vec2[i];
	
	return sum;
}

/// Adds one matrix onto another
vector <vector <double> > matrix_add(const vector <vector <double> > &M1, const vector <vector <double> > &M2)
{
	const auto NY = M1.size(), NX = M1[0].size();

	if(M2.size() != NY || M2[0].size() != NX) emsgEC("Matrix",9);	
	auto result = M1;
	
	for(auto j = 0u; j < NY; j++){
		for(auto i = 0u; i < NX; i++){
			result[j][i] += M2[j][i];
		}
	}
	
	return  result;
}


/// Calculates the determinant of a matrix
double determinant(const vector < vector <double> > &M) 
{
	timer[TIME_DETERMINANT].start();
		
	auto det = 0.0;
 
	auto n = M.size();
	if(n == 1) return M[0][0];
	if(n == 2) return M[0][0]*M[1][1] - M[0][1]*M[1][0];

	int sign = 1;
	for(auto i = 0u; i < n; i++){
		vector < vector <double> > sub;
		sub.resize(n-1);
		for(auto i = 0u; i < n-1; i++) sub[i].resize(n-1);
		
		auto ii = 0u, jj = 0u;
		for(auto row = 0u; row < n; row++){
			for(auto col = 0u; col < n; col++){
				if(row != 0 && col != i) {
					sub[ii][jj] = M[row][col];
					jj++;
					if(jj == n-1){
						jj = 0;
						ii++;
					}
				}
			}
		}
	 
		det += sign*M[0][i]*determinant(sub);
		sign = -sign;
	}
	
	timer[TIME_DETERMINANT].stop();
		
	return det;
}


/// Performs LU decomposition to return the log of the determinant
double determinant_fast(const vector < vector <double> > &a) 
{
	timer[TIME_DETERMINANT].start();
	
	auto n = a.size();  
	
	vector < vector <double> > l, u; // Performs LU decomposition (converts into a lower and an upper triagular matrix
	
	l.resize(n); u.resize(n);
	for(auto i = 0u; i < n; i++){ l[i].resize(n,0); u[i].resize(n,0);}
	
	auto i = 0u, j = 0u, k = 0u;
	for(i = 0; i < n; i++){
		auto &ui = u[i], &li = l[i];
		
		for(j = i; j < n; j++){
			auto &lj = l[j];
			auto val = a[j][i];
			for(k = 0; k < i; k++) val -= lj[k]*ui[k];
			lj[i] = val;
		}
		
		u[i][i] = 1;
		for(j = i+1; j < n; j++){
			auto &uj = u[j];
			
			if(li[i] == 0) emsgEC("Matrix",10);
			auto fac = 1.0/li[i];
			auto val = fac*a[i][j];
			for(k = 0; k < i; k++) val -= fac*li[k]*uj[k];
			u[j][i] = val;
		}
	}
	
	auto det = 0.0;
	auto neg = 0u;
	for(i = 0; i < n; i++){
		auto valu = u[i][i]; if(valu < 0){ valu = -valu; neg++;}
		auto vall = l[i][i]; if(vall < 0){ vall = -vall; neg++;}
		det += log(valu*vall);
	}
	
	timer[TIME_DETERMINANT].stop();
	
	return det;
}


/// Returns the dot product of two vectors using SIMD to speed up
double dot_SIMD(const double* vec1, const double* vec2, unsigned int n1, unsigned int n2)
{
	const auto si = 4;
	auto val = 0.0;
	auto k = n1;
	
	if(k > 0){
		auto kmin = si*((unsigned int)((k+si-1)/si));
		while(k < kmin && k < n2){ val += vec1[k]*vec2[k]; k++;} 
	}
	
	if(n2-k >= si){
		__m256d v1, v2;
		__m256d summed = _mm256_set1_pd(0.0);
		auto kmax = n2-si+1; 
		while(k < kmax){
			v1 = _mm256_load_pd(vec1+k);
			v2 = _mm256_load_pd(vec2+k);
			summed = _mm256_add_pd( summed, _mm256_mul_pd(v1,v2));
			k += si;
		}
		double* sum = (double*)&summed;
		
		for(auto i = 0u; i < si; i++) val += sum[i];
	}
	
	while(k < n2){ val += vec1[k]*vec2[k]; k++;} 

	return val;
}


/// Returns the dot product of two vectors using SIMD to speed up
float dot_SIMD(const float* vec1, const float* vec2, unsigned int n1, unsigned int n2)
{
	const auto si = 8;
	float val = 0.0;
	auto k = n1;
	
	if(k > 0){
		auto kmin = si*((unsigned int)((k+si-1)/si));
		while(k < kmin && k < n2){ val += vec1[k]*vec2[k]; k++;} 
	}
	
	if(n2-k >= si){
		__m256 v1, v2;
		__m256 summed = _mm256_set1_ps(0.0);
		auto kmax = n2-si+1; 
		while(k < kmax){
			v1 = _mm256_load_ps(vec1+k);
			v2 = _mm256_load_ps(vec2+k);
			summed = _mm256_add_ps( summed, _mm256_mul_ps(v1,v2));
			k += si;
		}
		float* sum = (float*)&summed;
		
		for(auto i = 0u; i < si; i++) val += sum[i];
	}
	
	while(k < n2){ val += vec1[k]*vec2[k]; k++;} 

	return val;
}


/// Performs LU decomposition to return the log of the determinant
double determinant_SIMD(const vector < vector <double> > &a, const Accuracy ac) 
{
	timer[TIME_DETERMINANT].start();
	
	auto n = a.size();  

	auto det = 0.0;
	auto neg = 0u;
	
	if(ac == DOUBLE){
		double **l, **u;
		l = (double**)malloc(n*sizeof(double*));
		u = (double**)malloc(n*sizeof(double*));
		for(auto j = 0u; j < n; j++){
			l[j] = (double*)aligned_alloc(32, n * sizeof(double));
			u[j] = (double*)aligned_alloc(32, n * sizeof(double));
			for(auto i = 0u; i < n; i++){ l[j][i] = 0; u[j][i] = 0;}
		}
		
		for(auto i = 0u; i < n; i++){
			auto *ui = u[i], *li = l[i];
			
			for(auto j = i; j < n; j++){
				l[j][i] = a[j][i] - dot_SIMD(l[j],ui,0,i);
			}
			
			u[i][i] = 1;
			for(auto j = i+1; j < n; j++){
				if(li[i] == 0) emsgEC("Matrix",11);
				u[j][i] = (a[i][j] - dot_SIMD(li,u[j],0,i))/li[i];
			}
		}
	
		for(auto i = 0u; i < n; i++){
			auto valu = u[i][i]; if(valu < 0){ valu = -valu; neg++;}
			auto vall = l[i][i]; if(vall < 0){ vall = -vall; neg++;}
			det += log(valu*vall);
		}
	}
	else{		
	/*
		float **l, **u;
		l = (float**)malloc(n*sizeof(float*));
		u = (float**)malloc(n*sizeof(float*));
		for(auto j = 0u; j < n; j++){
			l[j] = (float*)aligned_alloc(32, n * sizeof(float));
			u[j] = (float*)aligned_alloc(32, n * sizeof(float));
			for(auto i = 0u; i < n; i++){ l[j][i] = 0; u[j][i] = 0;}
		}
		
		for(auto i = 0u; i < n; i++){
			auto *ui = u[i], *li = l[i];
			
			for(auto j = i; j < n; j++){
				l[j][i] = a[j][i] - dot_SIMD(l[j],ui,i);
			}
			
			u[i][i] = 1;
			for(auto j = i+1; j < n; j++){
				if(li[i] == 0) emsgEC("Matrix",12);
				u[j][i] = (a[i][j] - dot_SIMD(li,u[j],i))/li[i];
			}
		}
		
		for(auto i = 0u; i < n; i++){
			auto valu = u[i][i]; if(valu < 0){ valu = -valu; neg++;}
			auto vall = l[i][i]; if(vall < 0){ vall = -vall; neg++;}
			det += log(valu*vall);
		}
		*/
	}
	timer[TIME_DETERMINANT].stop();

	return det;
}

void matrix_operation_initialise(const unsigned int n)
{
	n_matrix_store = n;
	 
	l = (double**)malloc(n*sizeof(double*));
	linv = (double**)malloc(n*sizeof(double*));
	u = (double**)malloc(n*sizeof(double*));
	uinv = (double**)malloc(n*sizeof(double*));
	mat2 = (double**)malloc(n*sizeof(double*));
	matres = (double**)malloc(n*sizeof(double*));
	for(auto j = 0u; j < n; j++){
		l[j] = (double*)aligned_alloc(32, n * sizeof(double));
		linv[j] = (double*)aligned_alloc(32, n * sizeof(double));
		u[j] = (double*)aligned_alloc(32, n * sizeof(double));
		uinv[j] = (double*)aligned_alloc(32, n * sizeof(double));
		mat2[j] = (double*)aligned_alloc(32, n * sizeof(double));
		matres[j] = (double*)aligned_alloc(32, n * sizeof(double));
		for(auto i = 0u; i < n; i++){ l[j][i] = 0; linv[j][i] = 0; u[j][i] = 0; uinv[j][i] = 0;}
	}
	
	lf = (float**)malloc(n*sizeof(float*));
	lfinv = (float**)malloc(n*sizeof(float*));
	uf = (float**)malloc(n*sizeof(float*));
	ufinv = (float**)malloc(n*sizeof(float*));
	matf2 = (float**)malloc(n*sizeof(float*));
	matresf = (float**)malloc(n*sizeof(float*));
	for(auto j = 0u; j < n; j++){
		lf[j] = (float*)aligned_alloc(32, n * sizeof(float));
		lfinv[j] = (float*)aligned_alloc(32, n * sizeof(float));
		uf[j] = (float*)aligned_alloc(32, n * sizeof(float));
		ufinv[j] = (float*)aligned_alloc(32, n * sizeof(float));
		matf2[j] = (float*)aligned_alloc(32, n * sizeof(float));
		matresf[j] = (float*)aligned_alloc(32, n * sizeof(float));
		for(auto i = 0u; i < n; i++){ lf[j][i] = 0; lfinv[j][i] = 0; uf[j][i] = 0; ufinv[j][i] = 0;}
	}
}


/// Inverts a lower triangular
void invert_lower_triangular(double **M, double **Minv, const unsigned int n)
{	
	for(auto j = 0u; j < n; j++){
		auto fac = 1.0/M[j][j];
		for(auto i = 0u; i < n; i++){
			Minv[j][i] = 0;
			M[j][i] *= fac;
		}
		Minv[j][j] = fac;
		M[j][j] = 1;
		
		for(auto i = 0u; i < j; i++){
			sub_row_SIMD(Minv[j],Minv[i], -M[j][i],0,i+1);
		}
	}
}


/// Inverts a lower triangular
void invert_lower_triangular(float **M, float **Minv, const unsigned int n)
{	
	for(auto j = 0u; j < n; j++){
		auto fac = 1.0/M[j][j];
		for(auto i = 0u; i < n; i++){
			Minv[j][i] = 0;
			M[j][i] *= fac;
		}
		Minv[j][j] = fac;
		M[j][j] = 1;
		
		for(auto i = 0u; i < j; i++){
			sub_row_SIMD(Minv[j],Minv[i], -M[j][i],0,i+1);
		}
	}
}


/// Transposes a matrix (in pointer form)
void transpose(double **M, const unsigned int n)
{
	double temp;
	for(auto j = 0u; j < n; j++){
		for(auto i = 0u; i < j; i++){
			temp = M[j][i]; M[j][i] = M[i][j]; M[i][j] = temp;
		}
	}
}


/// Transposes a matrix (in pointer form)
void transpose(float **M, const unsigned int n)
{
	float temp;
	for(auto j = 0u; j < n; j++){
		for(auto i = 0u; i < j; i++){
			temp = M[j][i]; M[j][i] = M[i][j]; M[i][j] = temp;
		}
	}
}


/// Determines if a matrix contains unidentified elements
bool matrix_isnan(const vector < vector <double> > M)
{
	for(auto j = 0u; j < M.size(); j++){
		for(auto i = 0u; i < M[j].size(); i++){
			if(std::isnan(M[j][i])) return true;
		}
	}
	return false;
}


/// Adds one vector onto another
vector <double> vec_add(const vector <double> &vec1, const vector <double> &vec2)
{
	const auto N = vec1.size();
	if(vec2.size() != N) emsgEC("Matrix",13);
	
	auto result = vec1;
	
	for(auto i = 0u; i < N; i++) result[i] += vec2[i];
	
	return result;
}


/// Adds one vector onto another
vector <double> vec_subtract(const vector <double> &vec1, const vector <double> &vec2)
{
	const auto N = vec1.size();
	if(vec2.size() != N) emsgEC("Matrix",15);
	
	auto result = vec1;
	for(auto i = 0u; i < N; i++) result[i] -= vec2[i];
	
	return result;
}


/// Solves a linear set of equations
vector <double> linear_solve(vector <vector <double> > mat, vector <double> vec)   
{
	timer[TIME_LINEAR_EQ].start();
	
	auto nvar = mat.size();

	for(auto ii = 0u; ii < nvar; ii++){
    auto r = mat[ii][ii];
    for(auto i = 0u; i < nvar; i++) mat[ii][i] /= r; 
		vec[ii] /= r;
	
    for(auto jj = ii+1; jj < nvar; jj++){
      auto r = mat[jj][ii];
			if(r != 0){
				for(auto i = 0u; i < nvar; i++) mat[jj][i] -= r*mat[ii][i];
				vec[jj] -= r*vec[ii];
			}
    }
  }
	
	
  for(int ii = nvar-1; ii > 0; ii--){
    for(int jj = ii-1; jj >= 0; jj--){
      auto r = mat[jj][ii];
			if(r != 0){
				for(auto i = 0u; i < nvar; i++) mat[jj][i] -= r*mat[ii][i];
				vec[jj] -= r*vec[ii];
			}
    }
  }
	
	timer[TIME_LINEAR_EQ].stop();
	
	return vec;
}


/// Inverts a matrix
vector < vector <double> > invert_matrix(vector <vector <double> > M)   
{
	timer[TIME_INV_MATRIX].start();
	
	unsigned int nvar = M.size();
	
	vector <vector <double> > inv_M;
	
	inv_M.resize(nvar);
  for(auto i = 0u; i < nvar; i++) inv_M[i].resize(nvar,0);
	for(auto i = 0u; i < nvar; i++) inv_M[i][i] = 1;
   
  for(auto ii = 0u; ii < nvar; ii++){
		auto &Mii = M[ii], &inv_Mii = inv_M[ii];
		
    auto r = Mii[ii];
    for(auto i = ii; i < nvar; i++) Mii[i] /= r; 
	  for(auto i = 0u; i <= ii; i++) inv_Mii[i] /= r; 
	
    for(auto jj = ii+1; jj < nvar; jj++){
      auto r = M[jj][ii];
			if(r != 0){
				auto &Mjj = M[jj], &inv_Mjj = inv_M[jj];
			
				for(auto i = ii; i < nvar; i++) Mjj[i] -= r*Mii[i];
				for(auto i = 0u; i <= ii; i++) inv_Mjj[i] -= r*inv_Mii[i];
			}
    }
  }

  for(int ii = nvar-1; ii > 0; ii--){
		//auto &Mii = M[ii];
		auto &inv_Mii = inv_M[ii];
    for(int jj = ii-1; jj >= 0; jj--){
      auto r = M[jj][ii];
			if(r != 0){
				//auto &Mjj = M[jj];
				auto &inv_Mjj = inv_M[jj];
		
				//for(auto i = 0u; i < nvar; i++) Mjj[i] -= r*Mii[i];
				for(auto i = 0u; i < nvar; i++) inv_Mjj[i] -= r*inv_Mii[i];
			}
    }
  }

	if(false){ // checks inverse
		cout << "Check matrix" << endl;
		for(auto j = 0u; j < nvar; j++){
			for(auto i = 0u; i < nvar; i++){
				double sum = 0; for(auto ii = 0u; ii < nvar; ii++) sum += M[j][ii]*inv_M[ii][i];
				
				cout << sum << " ";
				if(i != j){ if(sum < -TINY || sum > TINY) emsgEC("Matrix",16);}
				else{ if(sum < 1-TINY || sum > 1+TINY) emsgEC("Matrix",17);}		
			}
			cout << endl;
		}
	}
	
	timer[TIME_INV_MATRIX].stop();
	
	return inv_M;
}


/// Subtracts one row from another
void sub_row_SIMD(double *row1, const double *row2, double fac, unsigned int i1, unsigned int i2)
{
	const auto si = 4;
	auto k = i1;
	if(k > 0){
		auto kmin = si*((unsigned int)((k+si-1)/si));
		while(k < kmin){ row1[k] += fac*row2[k]; k++;} 
	}
	
	if(i2-k >= si){
		__m256d v1, v2, res;
		__m256d f = _mm256_set1_pd(fac);
		auto kmax = i2-si+1; 
		while(k < kmax){
			v1 = _mm256_load_pd(row1+k);
			v2 = _mm256_load_pd(row2+k);
			res = _mm256_fmadd_pd(v2,f,v1);	
			_mm256_store_pd(row1+k,res);
			k += si;
		}
	}
	
	while(k < i2){ row1[k] += fac*row2[k]; k++;} 
}


/// Subtracts one row from another
void sub_row_SIMD(float *row1, const float *row2, float fac, unsigned int i1, unsigned int i2)
{
	const auto si = 8;
	auto k = i1;
	
	if(k > 0){
		auto kmin = si*((unsigned int)((k+si-1)/si));
		while(k < kmin && k < i2){ row1[k] += fac*row2[k]; k++;} 
	}
	
	if(i2-k >= si){
		__m256 v1, v2, res;
		__m256 f = _mm256_set1_ps(fac);
		auto kmax = i2-si+1; 
		while(k < kmax){
			v1 = _mm256_load_ps(row1+k);
			v2 = _mm256_load_ps(row2+k);
			res = _mm256_fmadd_ps(v2,f,v1);	
			_mm256_store_ps(row1+k,res);
			k += si; 
		}
	}

	while(k < i2){ row1[k] += fac*row2[k]; k++;} 
}

			
/// Inverts a matrix using SIMD acceleration
vector < vector <double> > invert_matrix_SIMD(const vector <vector <double> > &M_orig, Accuracy ac)   
{
	timer[TIME_INV_MATRIX].start();
	
	unsigned int nvar = M_orig.size();
	
	vector < vector <double> > invM;
	invM.resize(nvar); for(auto j = 0u; j < nvar; j++) invM[j].resize(nvar);
			
	if(ac == DOUBLE){
		double **M, **inv_M;
		M = (double**)malloc(nvar*sizeof(double*));
		inv_M = (double**)malloc(nvar*sizeof(double*));
		for(auto j = 0u; j < nvar; j++){
			M[j] = (double*)aligned_alloc(32, nvar * sizeof(double));
			inv_M[j] = (double*)aligned_alloc(32, nvar * sizeof(double));
			for(auto i = 0u; i < nvar; i++){ M[j][i] = M_orig[j][i]; inv_M[j][i] = 0;}
		}
		
		for(auto i = 0u; i < nvar; i++) inv_M[i][i] = 1;
		 
		for(auto ii = 0u; ii < nvar; ii++){
			auto *Mii = M[ii], *inv_Mii = inv_M[ii];
			
			auto r = Mii[ii];
			for(auto i = ii; i < nvar; i++) Mii[i] /= r; 
			for(auto i = 0u; i <= ii; i++) inv_Mii[i] /= r; 
		
			for(auto jj = ii+1; jj < nvar; jj++){
				auto r = M[jj][ii];
				if(r != 0){
					sub_row_SIMD(M[jj],Mii,-r,ii,nvar);
					sub_row_SIMD(inv_M[jj],inv_Mii,-r,0,ii+1);
				}
			}
		}

		for(int ii = nvar-1; ii > 0; ii--){
			auto *inv_Mii = inv_M[ii];
			for(int jj = ii-1; jj >= 0; jj--){
				auto r = M[jj][ii];
				if(r != 0){
					sub_row_SIMD(inv_M[jj],inv_Mii,-r,0,nvar);
				}
			}
		}

		for(auto j = 0u; j < nvar; j++){	
			for(auto i = 0u; i < nvar; i++) invM[j][i] = inv_M[j][i];
		}
	}
	else{
		float **M, **inv_M;
		M = (float**)malloc(nvar*sizeof(float*));
		inv_M = (float**)malloc(nvar*sizeof(float*));
		for(auto j = 0u; j < nvar; j++){
			M[j] = (float*)aligned_alloc(32, nvar * sizeof(float));
			inv_M[j] = (float*)aligned_alloc(32, nvar * sizeof(float));
			for(auto i = 0u; i < nvar; i++){ M[j][i] = M_orig[j][i]; inv_M[j][i] = 0;}
		}
		
		for(auto i = 0u; i < nvar; i++) inv_M[i][i] = 1;
		 
		for(auto ii = 0u; ii < nvar; ii++){
			auto *Mii = M[ii], *inv_Mii = inv_M[ii];
			
			auto r = Mii[ii];
			for(auto i = ii; i < nvar; i++) Mii[i] /= r; 
			for(auto i = 0u; i <= ii; i++) inv_Mii[i] /= r; 
		
			for(auto jj = ii+1; jj < nvar; jj++){
				auto r = M[jj][ii];
				if(r != 0){
					sub_row_SIMD(M[jj],Mii,-r,ii,nvar);
					sub_row_SIMD(inv_M[jj],inv_Mii,-r,0,ii+1);
				}
			}
		}

		for(int ii = nvar-1; ii > 0; ii--){
			auto *inv_Mii = inv_M[ii];
			for(int jj = ii-1; jj >= 0; jj--){
				auto r = M[jj][ii];
				if(r != 0){
					sub_row_SIMD(inv_M[jj],inv_Mii,-r,0,nvar);
				}
			}
		}
		
		for(auto j = 0u; j < nvar; j++){	
			for(auto i = 0u; i < nvar; i++) invM[j][i] = double(inv_M[j][i]);
		}
		
		//for(auto j = 0u; j < nvar; j++){ free(M[j]); free(inv_M[j]);}
		//free(M); free(inv_M);		
	}
	
	timer[TIME_INV_MATRIX].stop();
	
	return invM;
}


/// Starts by permuting the rows such that the diagonals are non-zero and then does inversion
vector < vector <double> > invert_matrix_permute(const vector <vector <double> > &mat)
{
	auto M = mat;
	auto N = M.size(); 
	
	vector <unsigned int> permute;
	for(auto i = 0u; i < N; i++) permute.push_back(i);
	
	for(auto j = 0u; j < N-1; j++){
		if(M[j][j] != 1){
			auto jj = j+1; while(jj < N && M[jj][j] != 1) jj++;
			if(jj == N) emsgEC("Matrix",18);
			
			for(auto i = 0u; i < N; i++){
				auto temp = M[j][i];
				M[j][i] = M[jj][i];
				M[jj][i] = temp;
			}
			auto temp = permute[j]; permute[j] = permute[jj]; permute[jj] = temp;
		}
	}
	
	auto inv = invert_matrix(M);
	
	vector < vector <double> > answer;
	answer.resize(N);
	for(auto j = 0u; j < N; j++){
		answer[j].resize(N);
		for(auto i = 0u; i < N; i++){
			answer[j][permute[i]] = inv[j][i];
		}
	}
	
	if(true){  // Checks the result
		auto mult = matrix_mult(mat,answer);
		for(auto j = 0u; j < N; j++){
			for(auto i = 0u; i < N; i++){
				auto val = 0.0; if(j == i) val = 1;
				if(mult[j][i] < val-TINY || mult[j][i] > val+TINY){
					print_matrix("diag",mult);
					emsgEC("Matrix",19);
				}
			}
		}
	}
	return answer;
}


/// Inverts a matrix
vector < vector <double> > invert_matrix2(const vector <vector <double> > &mat)   
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
			if(r != 0){
				for(auto i = 0u; i < nvar; i++){ 
					A2[jj][i] -= r*A2[ii][i];
					inv_M[jj][i] -= r*inv_M[ii][i];
				}
			}
    }
  }

  for(int ii = nvar-1; ii > 0; ii--){
    for(int jj = ii-1; jj >= 0; jj--){
      double r = A2[jj][ii];
			if(r != 0){
				for(auto i = 0u; i < nvar; i++){ 
					A2[jj][i] -= r*A2[ii][i];
					inv_M[jj][i] -= r*inv_M[ii][i];
				}
			}
    }
  }

	if(false){ // checks inverse
		for(auto j = 0u; j < nvar; j++){
			for(auto i = 0u; i < nvar; i++){
				double sum = 0; for(auto ii = 0u; ii < nvar; ii++) sum += mat[j][ii]*inv_M[ii][i];
				
				if(i != j){ if(sum < -TINY || sum > TINY) emsgEC("Matrix",20);}
				else{ if(sum < 1-TINY || sum > 1+TINY) emsgEC("Matrix",21);}		
			}
		}
	}
	
	return inv_M;
}


/// Transposes a matrix
vector < vector <double> > transpose(const vector < vector <double> > &M)
{
	vector < vector <double> > T;
	auto Y = M[0].size(), X = M.size();
	
	T.resize(Y);
	for(auto j =0u; j < Y; j++){
		T[j].resize(X);
		for(auto i = 0u; i < X; i++) T[j][i] = M[i][j];
	}
	
	return T;
}


/// Calculates the square root of the inverse matrix using Denman–Beavers iteration
vector < vector <double> > invert_matrix_square_root(const vector < vector <double> > &M)
{
	auto n = M.size();
	
	auto Y = M;
	
	vector < vector <double> > Z;  // Sets the identity matrix
	Z.resize(n);
	for(auto j = 0u; j < n; j++){
		Z[j].resize(n);
		for(auto i = 0u; i < n; i++){
			if(i == j) Z[j][i] = 1;
			else Z[j][i] = 0;
		}
	}

	auto limit = VTINY;
	auto loop = 0u, loopmax = 1000u;
	do{
		auto Yinv = invert_matrix(Y);
		auto Zinv = invert_matrix(Z);
		
		vector < vector <double> > Ynew;
		Ynew.resize(n);
		for(auto j = 0u; j < n; j++){
			Ynew[j].resize(n);
			for(auto i = 0u; i < n; i++){
				Ynew[j][i] = 0.5*(Y[j][i] + Zinv[j][i]);
			}
		}
		
		vector < vector <double> > Znew;
		Znew.resize(n);
		for(auto j = 0u; j < n; j++){
			Znew[j].resize(n);
			for(auto i = 0u; i < n; i++){
				Znew[j][i] = 0.5*(Z[j][i] + Yinv[j][i]);
			}
		}
		
		auto dif = LARGE;
		for(auto j = 0u; j < n; j++){
			for(auto i = 0u; i < n; i++){
				auto d = Znew[j][i]-Z[j][i];
				if(d < 0) d = -d;
				if(d < dif) dif = d;
			}
		}
		
		Y = Ynew; Z = Znew;
		if(dif < limit) break;
		loop++; if(loop%100 == 0) limit *= 10;
	}while(loop < loopmax);
	if(loop == loopmax) emsgEC("Matrix",22);
	
	if(false){
		auto Minv = invert_matrix(M);
		print_matrix("Minv",Minv);
		
		auto N = matrix_mult(Z,Z);
		print_matrix("N",N);
		
		for(auto j = 0u; j < n; j++){
			for(auto i = 0u; i < n; i++){
				auto d = Minv[j][i] - N[j][i];
				if(d < -TINY || d > TINY){
					cout << d << " d\n";
					emsgEC("Matrix",23);
				}
			}
		}		
		emsg("matrix check");
	}
	
	return Z;
}


vector <double> initial_guess;

/// Determines the largest eigenvalue and vector from a matrix
double largest_eigenvalue(const vector < vector <double> > &M, vector <double> &eigenvector)
{
	auto N = M.size();
	
	auto ev = 0.0;

	if(false){
		Eigen::MatrixXd mat(N,N);
		for(auto j = 0u; j < N; j++){
			for(auto i = 0u; i < N; i++) mat(j,i) = M[j][i];
		}
	
		Eigen::EigenSolver<Eigen::MatrixXd> eigensolver(mat);
		auto vec = eigensolver.eigenvalues();
	
		ev=-LARGE;
		auto num=0u;
		for(auto i = 0u; i < N; i++){
			auto evv = vec[i].real();
			if(evv > ev){ ev = evv; num = i;}
		}
		
		auto es = eigensolver.eigenvectors();
		auto evector = es.col(num);
		
		eigenvector.resize(N);
		for(auto i = 0u; i < N; i++) eigenvector[i] = evector[i].real();
		
		auto sum = 0.0; for(auto i = 0u; i < N; i++) sum += eigenvector[i];
		for(auto i = 0u; i < N; i++) eigenvector[i] /= sum;
	}
	else{
		if(initial_guess.size() == 0){
			for(auto i = 0u; i < N; i++) initial_guess.push_back(ran());
		}
		
		vector <double> vec = initial_guess;
		vector <double> vec2(N);
		
		auto loop = 0u;
		do{
			for(auto i = 0u; i < N; i++){
				auto sum = 0.0; for(auto ii = 0u; ii < N; ii++) sum += M[i][ii]*vec[ii];
				vec2[i] = sum;
			}
			
			auto sum = 0.0; for(auto i = 0u; i < N; i++) sum += vec2[i];
			auto ev_new = sum;
			
			for(auto i = 0u; i < N; i++) vec2[i] /= sum;
		
			if(loop%10 == 0){
				auto lim = 100*VTINY;
				if(loop > 100) lim *= 10;
				if(loop > 1000) lim *= 10;
				if(ev_new-ev > -lim && ev_new-ev < lim){
					unsigned int k;
					for(k = 0; k < N; k++){
						auto dif = vec2[k] - vec[k];
						if(dif > lim || dif < -lim) break;
					}
					if(k == N) break;
				}
			}
			
			vec = vec2;	
			ev = ev_new;

			loop++;
			if(loop > 10000){
				print_matrix("M",M);
				emsg("Eigen-vector convergence problem");
			}
		}while(1 == 1);
	
		initial_guess = vec;

		eigenvector = vec;
	}
	
	return ev;
}


/// Prints a matrix
void print_matrix(string name, const vector < vector <double> > &M)
{
	cout << name << ":" << endl;	
	for(auto j = 0u; j < M.size(); j++){
		cout << "   ";
		for(auto i = 0u; i < M[j].size(); i++){
			cout << M[j][i] << " ";
		}
		cout << endl;
	}
}


/// Prints a vector
void print_vector(string name, const vector <double> &vec)
{
	cout << name << ": ";	
	for(auto j = 0u; j < vec.size(); j++){
		cout << vec[j] << ", ";
	}
	cout << " VECTOR" << endl;
}


