#ifndef BEEP__MATRIX_HH
#define BEEP__MATRIX_HH

using namespace std;

#include "struct.hh"

vector < vector <double> > invert_determinant_SIMD(const vector < vector <double> > &a, double &det, Accuracy ac);

vector <double> matrix_mult(const vector < vector <double> > &M, const vector <double> &vec);
vector < vector <double> > matrix_mult(const vector < vector <double> > &M1, const vector < vector <double> > &M2);
vector < vector <double> > matrix_mult_sparse(const vector < vector <double> > &M1, const vector < vector <double> > &M2);
void set_mat(const vector < vector <double> > &M, const Accuracy ac);
void covar_add_mat(vector < vector <double> > &M, const Accuracy ac);
void covar_matrix_mult_SIMD(const vector < vector <double> > &M1, const unsigned int NX2, const Accuracy ac);
double vec_mult(const vector <double> &vec1, const vector <double> &vec2);
double determinant(const vector < vector <double> > &M);
double determinant_fast(const vector < vector <double> > &M);
double dot_SIMD(const double* vec1, const double* vec2, unsigned int n1, unsigned int n2);
float dot_SIMD(const float* vec1, const float* vec2, unsigned int n1, unsigned int n2);
double determinant_SIMD(const vector < vector <double> > &a, const Accuracy ac);
void matrix_operation_initialise(const unsigned int n);

bool matrix_isnan(const vector < vector <double> > M);
vector <vector <double> > matrix_add(const vector <vector <double> > &M1, const vector <vector <double> > &M2);
vector <double> vec_add(const vector <double> &vec1, const vector <double> &vec2);
vector <double> vec_subtract(const vector <double> &vec1, const vector <double> &vec2);
vector <double> linear_solve(vector <vector <double> > mat, vector <double> vec);
vector < vector <double> > invert_matrix_permute(const vector <vector <double> > &mat);
vector <vector <double> > invert_matrix(vector <vector <double> > mat);
void sub_row_SIMD(double *row1, const double *row2, double fac, unsigned int i1, unsigned int i2);
void sub_row_SIMD(float *row1, const float *row2, float fac, unsigned int i1, unsigned int i2);
void invert_lower_triangular(double **M, double **Minv, const unsigned int n);
void invert_lower_triangular(float **M, float **Minv, const unsigned int n);
void transpose(double **M, const unsigned int n);
void transpose(float **M, const unsigned int n);
vector < vector <double> > invert_matrix_SIMD(const vector <vector <double> > &M_orig, Accuracy ac);
vector < vector <double> > invert_matrix2(const vector <vector <double> > &mat);
vector < vector <double> > transpose(const vector < vector <double> > &M);
vector < vector <double> > invert_matrix_square_root(const vector < vector <double> > &M);
double largest_eigenvalue(const vector < vector <double> > &M, vector <double> &eigenvector);
double vec_dot(const vector <double> &vec1, const vector <double> &vec2);
void print_matrix(string name, const vector < vector <double> > &M);
void print_vector(string name, const vector <double> &vec);

#endif
