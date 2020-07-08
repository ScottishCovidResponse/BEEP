#pragma once

#define USE_MPI 1

#ifdef USE_MPI
#include <mpi.h>
#endif

#include <string>
#include <vector>

using namespace std;

[[ noreturn ]] void emsg(string msg);
[[ noreturn ]] void emsgroot(string msg);

double ran();
double normal(float mu, float sd);
double gammasamp(double a, double b);
vector<string> split(const string& s, char delimiter);
