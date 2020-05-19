
#include <vector>

using namespace std;

#include "types.hh"
#include "functions.hh"
#include "consts.hh"

class PART;

long ncase[nregion][tmax/7+1];
long timetot=0, timesim=0, timeboot=0;     ///< Stores the CPU clock times for different parts of the algorithm
