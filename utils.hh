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
double normal(float mu, float sd);
double normalprob(double x, double mean, double var);
double lognormprob(double x, double mean, double var);
double gammasamp(double a, double b);
double gammaprob(double x, double a, double b);
double exp_sample(double rate);

vector<string> split(const string& s, char delimiter);
string filebasename(const string &path);
bool stringhasending (std::string const &fullString, std::string const &ending);
unsigned int get_integer(const std::string& st, unsigned int threshold);

#endif
