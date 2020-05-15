
extern short siminf;                              // Set to 1 for simulation and 0 for inference
extern long nsamp;                                // The number of PMCMC samples
extern long areamax;                              // The maximum number of areas 
extern long ncase[nregion][tmax/7+1];
extern vector <HOUSE> house;                      // List of all the houses
extern long **nMval;                              // These are used to store the matrix M of interactions between individuals
extern float ***Mval;
extern long ***Mnoderef;
extern long **naddnoderef;
extern long ***addnoderef;
extern vector <LEVEL> lev;
extern short level;                               // The number of levels of scale in the model
extern vector <vector <long> > subpop;            // List of all individuals in node on the fine scale
extern long Cfine;                                // Number of nodes on the fine scale
extern vector <IND> ind;
extern PART* part[partmax];                       // Pointers to each of the particles 
extern long npart;                                // The number of particles used
extern long timetot, timesim, timeboot;     // Stores the CPU clock times for different parts of the algorithm
extern double settime[nsettime];
extern double beta[nsettime];
extern short nspline;                             // The spline points which are parameters in the model
extern vector <double> splinet;
extern vector <PARAM> param;
extern vector <TRANS> trans;
extern vector <COMP> comp;	
