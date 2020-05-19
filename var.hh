#pragma once

#include <vector>

using namespace std;

#include "consts.hh"
#include "types.hh"

class PART;

extern short siminf;                              // Set to 1 for simulation and 0 for inference
extern long ncase[nregion][tmax/7+1];
extern long npart;                                // The number of particles used
extern long timetot, timesim, timeboot;     // Stores the CPU clock times for different parts of the algorithm
