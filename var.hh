#pragma once

#include <vector>

using namespace std;

#include "consts.hh"
#include "types.hh"

class PART;

extern long ncase[nregion][tmax/7+1];
extern long timetot, timesim, timeboot;     // Stores the CPU clock times for different parts of the algorithm
