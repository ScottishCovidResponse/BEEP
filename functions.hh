#pragma once

#include <string>

#include "types.hh"

using namespace std;

void emsg(string msg);
double ran();
bool compNEV(NEV lhs, NEV rhs);
double normal(float mu, float sd);
bool compX(long lhs, long rhs);
bool compY(long lhs, long rhs);
double gammasamp(double a, double b);
