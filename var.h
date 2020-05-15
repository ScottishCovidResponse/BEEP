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

struct HOUSE {                             // Defines a house
 	double x, y;                             // Position
	vector <long> ind;                       // Individuals which belong to the house
};
vector <HOUSE> house;                      // List of all the houses

long **nMval;                              // These are used to store the matrix M of interactions between individuals
float ***Mval;
long ***Mnoderef;
long **naddnoderef;
long ***addnoderef;

	
struct NODE {                              // Provides information about a node
	vector <long> houseref;                  // References the list of houses within the node
	long parent;                             // The parent node
	vector <long> child;                     // The child nodes
	vector <long> fine;                      // The child nodes on the fine scale
	long popu;                               // The total population in the node
	double x, y;                             // The position of the node (given by the average of all the houses)
	short done;                              // A flag used to determine if this node has been analysed or not
};
 
bool compX(long lhs, long rhs) { return house[lhs].x < house[rhs].x; }
bool compY(long lhs, long rhs) { return house[lhs].y < house[rhs].y; }

struct LEVEL {                             // Stores information about different levels 
	vector <NODE> node;                      // The nodes at a given level
 	vector <long> donelist;                  // List of nodes which have been processed
	vector <double> add;                     // Used when adding up the tree
};
vector <LEVEL> lev;

short level;                               // The number of levels of scale in the model

vector <vector <long> > subpop;            // List of all individuals in node on the fine scale

long Cfine;                                // Number of nodes on the fine scale

struct IND {                               // Provides information about an individual
	long noderef;                            // The node on the finescale to which the individual belongs
	long houseref;                           // The house to which the individual belongs
	short region;                            // The region to which the individual belongs
};
vector <IND> ind;

struct NEV {                               // Information about the immediate next events
  short type; double t;
};
bool compNEV(NEV lhs, NEV rhs) { return lhs.t < rhs.t; }

struct FEV {                               // Stores information about a compartmental transition
  long trans;                              // References the transition type
	long ind;                                // The individual on which the transition happens
	double t;                                // The time of the transition
	short done;                              // Set to 1 if that transition is in the past 
};

const long fediv = 1000;                   // The number of divisions into which the global timeline is divided

class PART                                 // Stores all the things related to a particle 
{
	public:
	long pst;                                // The number of the particle
	
 	double Li;                               // The observation likelihood
 	
	vector <double> ffine;                   // Stores the force of infection on nodes on the fine scale
	vector <vector <long> > indinf;          // Lists all infected individuals  
	vector <vector <long> > pop;             // The total popualtion for nodes on different levels 
	vector <vector <double> > Rtot;          // The total infection rate for nodes on different levels
	vector <vector <double> > addlater;      // A change to the rates Rtot which may be performed later when sampling is performed

	vector < vector <FEV> > fev;             // Stores all compartmental transitions
	
	vector <long> N;                         // The number of individuals in different compartments

	short sett;                              // Index used to track time changes in beta

	long tdnext, tdfnext;                    // Stores when the next future compartmental transition will occur

	public: 
		void gillespie(double ti, double tf);
		void partinit(long p);
		void dofe();
		long nextinfection();
		void addinfc(long c, double t);
		void addfev(double t, long tr, long i);
		vector <long> getnumtrans(string from, string to, short ti, short tf);
		void Lobs(short ti, short tf);
		void copy(long pfrom);
		void simmodel(long i, short enter, double t);
};

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

struct PARAM{                              // Store information about a model parameter
 	string name;                             // Its name
 	double val;                              // The simulation value or starting value for inference
	double sim;                              // The simulation value
	double min;                              // The minimum value (assuming a uniform prior) 
	double max;                              // The maximum value (assuming a uniform prior)
	double jump;                             // The size of proposed changes in PMCMC
	long ntr, nac;                           // Store the number of proposals tried and accepted	
};
vector <PARAM> param;

struct TRANS{                              // Stores information about a compartmental model transition
	short from;                              // Which compartment the individual is coming from
	short to;                                // Which compartment the individual is going to
	short type;                              // The type of transition (exponential or gamma)
	short param1;                            // First characteristic parameter (e.g. rate)
	short param2;                            // Second characteristic parameter (e.g. standard deviation in the case of gamma) 
};
vector <TRANS> trans;

struct COMP{                               // Stores information about a compartment in the model
	string name;                             // Its name
	double infectivity;                      // How infectious that compartment is
	vector <long> trans;                     // The transitions leaving that compartment
};
vector <COMP> comp;	

double ran(){                              // Draws a random number between 0 and 1
	if(RAND_MAX == 32767){
		double v = (double(rand())*32766.0+rand())/(32767.0*RAND_MAX); if(v == 0 || v == 1) return 0.1;
		else return v;
	}
	else return double(0.999999999*rand())/RAND_MAX;
}

// Draws a normally distributed number with mean mu and standard deviation sd
double normal(float mu, float sd){ return mu + sd*sqrt(-2*log(ran()))*cos(2*M_PI*ran());}

// Displays any error messages
void emsg(string msg){ cout << msg << "\n"; exit (EXIT_FAILURE);}

void init();                               // Function declarations
void simulatedata();
void PMCMC();
void readdata();
double bootstrap();
void addcomp(string name, double infectivity);
void addparam(string name, double val, double min, double max);
void addtrans(string from, string to, short type, string param1, string param2);
void betaspline();
double sample();
