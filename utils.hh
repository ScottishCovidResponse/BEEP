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

[[ noreturn ]] void emsg(string msg);
[[ noreturn ]] void emsgroot(string msg);

double ran();
double normal(float mu, float sd);
double normalprob(double x, double mean, double var);
double lognormprob(double x, double mean, double var);
double gammasamp(double a, double b);
double gammaprob(double x, double a, double b);

vector<string> split(const string& s, char delimiter);
void ensuredirectory(const string &path);
#endif
