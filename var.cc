
#include <vector>

using namespace std;

#include "types.hh"
#include "functions.hh"
#include "consts.hh"

short siminf;                              // Set to 1 for simulation and 0 for inference
long nsamp;                                // The number of PMCMC samples
long areamax;                              // The maximum number of areas 
long ncase[nregion][tmax/7+1];
vector <HOUSE> house;                      // List of all the houses
long **nMval;                              // These are used to store the matrix M of interactions between individuals
float ***Mval;
long ***Mnoderef;
long **naddnoderef;
long ***addnoderef;
vector <LEVEL> lev;
short level;                               // The number of levels of scale in the model
vector <vector <long> > subpop;            // List of all individuals in node on the fine scale
long Cfine;                                // Number of nodes on the fine scale
vector <IND> ind;
PART* part[partmax];                       // Pointers to each of the particles 
long npart;                                // The number of particles used
long timetot=0, timesim=0, timeboot=0;     // Stores the CPU clock times for different parts of the algorithm
double settime[nsettime];
double beta[nsettime];
short nspline;                             // The spline points which are parameters in the model
vector <double> splinet;
vector <PARAM> param;
vector <TRANS> trans;
vector <COMP> comp;	
