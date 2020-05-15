
#include "types.hh"
#include "function_decls.hh"

short siminf;                              // Set to 1 for simulation and 0 for inference

long nsamp;                                // The number of PMCMC samples
 
const short FEV_EV=0, INF_EV=1, SET_EV=2;  // Used to characterise an event type (future event/infection/settime)
	
const double tiny = 0.00000001;            // Used to represent a tiny number
const double large = 10000000;             // Used to represent a big number
const double invT = 1;                     // The inverse temperature (used to relax the observation model)

//const long popsize = 5000000;              // Number of individuals in Scotland (approximately)
//const long nhouse = 1500000;               // Number of houses
const long popsize = 10000;                 // A smaller test case
const long nhouse = 3000;               

long areamax;                              // The maximum number of areas 

const short checkon = 0;                   // Set to one to check algorithm is performing correctly

const double finegridsize = 0.02;          // The range in distance over which the fine grid is used 
const double d0 = 0.05;                    // The minumum distance cut-off for the matrix M
	
const short tmax = 105;                    // The time over which simulation / inference is performed

const short RX = 1, RY = 1;                // When siumlating this gives a hypothetical grid of regions
const short nregion = RX*RY;               // When the real data is analysed these will be Healthboard level regions
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



const long fediv = 1000;                   // The number of divisions into which the global timeline is divided



const long partmax = 10000;                // The maximum number of particles (arbitrarily set)
PART* part[partmax];                       // Pointers to each of the particles 
long npart;                                // The number of particles used

long timetot=0, timesim=0, timeboot=0;     // Stores the CPU clock times for different parts of the algorithm

const short EXP_DIST=0, GAMMA_DIST=1;      // Denotes exponential or gamma distributed 

const long nsettime = 100;                 // The number of time divisions used to represent changes in beta
double settime[nsettime];
double beta[nsettime];

short nspline;                             // The spline points which are parameters in the model
vector <double> splinet;

vector <PARAM> param;

vector <TRANS> trans;

vector <COMP> comp;	

#include "functions.h"
