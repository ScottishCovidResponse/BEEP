#ifndef BEEPMBP__CONSTS_HH
#define BEEPMBP__CONSTS_HH

using namespace std;

enum Mode { sim, inf, multisim, combinetrace};
	
//const unsigned int MODE_SIM=0, MODE_INF=1, MODE_MULTISIM=2;           // Different modes of operation 
 
const unsigned int FEV_EV=0, INF_EV=1;                           // Event types (future event/infection/settime/external)
const unsigned int SET_EV=2, EXT_EV=3, XIFEV_EV=4, XPFEV_EV=5;

const unsigned int EXP_DIST=1, GAMMA_DIST=2;                     // Denotes exponential or gamma distributed 
const unsigned int LOGNORM_DIST=3, INFECTION=4, TIMEP_DIST=5;

const unsigned int TRANS_DATA=0, POP_DATA=1, MARG_DATA=2;        // Different types of data

const unsigned int TFORM_NUM = 0, TFORM_YMD = 1;                 // Different type of time format

const double vtiny = 0.00000000000000001;                        // Used to represent a very tiny number
const double tiny = 0.00000001;                                  // Used to represent a tiny number
const double large = 10000000;                                   // Used to represent a big number
const unsigned int UNSET = 999999999;                            // A large unsigned integer to represent "Unset"
const unsigned int THRESH = 999999998;                           // Represents a number is under the threshold set
const unsigned int UNKNOWN = 999999997;                          // Represents a number is unknown
const unsigned int ADD = 999999996;                              // Used for summing posterior plots

const double minvar = 5; 																				 // The minimum variance for the observation model
const double Tpower = 4;                                         // The power used for the temerature progression

const unsigned int quenchpl = 0;                                 // Set to 1 if performing a quench plot

const unsigned int checkon = 0;                                  // Set to one to check algorithm is performing correctly
const unsigned int duplicate = 0;                                // Set to one to duplicate chains (this is used as a diagnostic check)

const unsigned int partmax = 10000;                              // The maximum number of particles per core (arbitrarily set)
const unsigned int chainmax = 10000;                             // The maximum number of chains per core (arbitrarily set)

const unsigned int MAX_NUMBERS = 15000000;                       // The maximum buffer size for Send Recv MPI messages

const unsigned int BIN=50;                                       // The number of bins used for plotting probability distributions
	
const unsigned int BOTH=0, PONLY=1, NOT=2;                       // Use to classify particles in MBPs
#endif
