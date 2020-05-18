
#include <vector>

using namespace std;

#include "types.hh"
#include "functions.hh"
#include "consts.hh"

class PART;

short siminf;                              ///< Set to 1 for simulation and 0 for inference
long nsamp;                                ///< The number of PMCMC samples
long ncase[nregion][tmax/7+1];
PART* part[partmax];                       ///< Pointers to each of the particles 
long npart;                                ///< The number of particles used
long timetot=0, timesim=0, timeboot=0;     ///< Stores the CPU clock times for different parts of the algorithm
