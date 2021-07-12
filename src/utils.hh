#ifndef BEEPMBP__UTILS_HH
#define BEEPMBP__UTILS_HH

#define USE_MPI 1

#ifdef USE_MPI
// See https://github.com/open-mpi/ompi/issues/5157
#define OMPI_SKIP_MPICXX 1
#include <mpi.h>
#endif

#include <string>
#include <vector>
#include <bits/stdc++.h>

using namespace std;

#include "struct.hh"

/// @cond EMSG
extern bool emsg_throws;
[[ noreturn ]] void emsg(const string &msg);
[[ noreturn ]] void emsgroot(const string &msg);
[[ noreturn ]] void emsgEC(const string &section, unsigned int ec);
/// @endcond
void warning(const string &msg);

double ran();
void sran(const int seed);
double normal_sample(const double mu, const double sd);
double normal_probability(const double x, const double mean, const double var);
double lognormal_sample(const double mean, const double sd);
double lognormal_probability(const double x, const double mean, const double var);
double gamma_sample(const double a, const double b);
double gamma_probability(const double x, const double a, const double b);
double negative_binomial_probability(const unsigned int val, const double m, const double shape);
double exp_sample(const double rate);
double exp_sample_time(const double rate);
int poisson_sample(const double lam);
double poisson_probability(const int i, const double lam);
unsigned int binomial_sample(const double ratio, const unsigned int nn);
double binomial_probability(const double ratio, const unsigned int nn, const unsigned int dn);
void binomial_check();
void strip(string &line);
void rem_pagebreak(string &line);
string toLower(string st);

vector<string> split(const string &s, const char delimiter);
string filebasename(const string &path);
bool stringhasending (std::string const &fullString, std::string const &ending);
unsigned int get_integer(const std::string &st, const unsigned int threshold);
vector <vector <double> > invert_matrix(const vector <vector <double> > &mat);
double largest_eigenvalue(const vector < vector <double> > &M, vector <double> &eigenvector);
unsigned int find_in(const vector <unsigned int> &vec, const unsigned int val);
unsigned int find_in(const vector <string> &vec, const string &val);
void add_vec(vector <unsigned int> &vec, const unsigned int val);
unsigned int find_char(const string st, const string char_in);
double vec_max(const vector <double> &vec);
string prec(const double num, const unsigned int pre);
string per(const string per);
string per(const double per);
bool allow_string(const string st, const string ok_char);
double get_double_with_unset(string st, const string em);
double get_double_with_tobeset(string st, const string em);
double get_double_positive(const string st, const string em);
double get_double(const string st, const string em);
unsigned int get_int(string st, const string em);
unsigned int get_int_error(string st);
double get_data(const string val, const string em, const string &threshold_str, const string &nodata_str);
string replace(const string st, const string s1, const string s2);
ParamSpec ps_one();
ParamSpec ps_zero();
#endif
