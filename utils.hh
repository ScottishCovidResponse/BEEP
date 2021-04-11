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

using namespace std;

/// @cond EMSG
extern bool emsg_throws;
[[ noreturn ]] void emsg(const string& msg);
[[ noreturn ]] void emsgroot(const string& msg);
[[ noreturn ]] void emsgEC(const string& section, unsigned int ec);
/// @endcond

double ran();
void sran(int seed);
double normal_sample(float mu, float sd);
double normal_probability(double x, double mean, double var);
double lognormal_sample(double mean, double sd);
double lognormal_probability(double x, double mean, double var);
double gamma_sample(double a, double b);
double gamma_probability(double x, double a, double b);
double negative_binomial_probability(unsigned int val, double m, double shape);
double exp_sample(double rate);
double exp_sample_time(double rate);
int poisson_sample(double lam);
double poisson_probability(int i, double lam);
unsigned int binomial_sample(double ratio, unsigned int nn);
double binomial_probability(double ratio, unsigned int nn, unsigned int dn);
void binomial_check();
string strip(string line);

vector<string> split(const string& s, char delimiter);
string filebasename(const string &path);
bool stringhasending (std::string const &fullString, std::string const &ending);
unsigned int get_integer(const std::string& st, unsigned int threshold);
vector <vector <double> > invert_matrix(const vector <vector <double> > &mat);
double largest_eigenvalue(const vector < vector <double> > &M, vector <double> &eigenvector);
unsigned int find_in(const vector <unsigned int> vec, const unsigned int val);
unsigned int find_in(const vector <string> vec, const string val);
string prec(double num, unsigned int pre);
string per(string per);
string per(double per);
#endif
