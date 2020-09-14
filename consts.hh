#ifndef BEEPMBP__CONSTS_HH
#define BEEPMBP__CONSTS_HH

using namespace std;

enum Mode { SIM, MULTISIM, MCMCMC, ABC_SMC, ABC_MBP, COMBINE_TRACE};          // Different modes of operation 
	
enum Dist { EXP_DIST, GAMMA_DIST, LOGNORM_DIST, INFECTION_DIST, TIMEP_DIST};  // Different time distributions

enum Alter { FAST, SLOW, NONE };                                              // Different speeds of altering proposals  

enum ParamType { OTHER_PARAM, DISTVAL_PARAM, BRANCHPROB_PARAM};               // Different types of parameters
	
enum DataType { TRANS_DATA, POP_DATA, MARG_DATA};                             // Different types of data

enum Status { SUCCESS, FAIL};                                                 // Determines if successfulful

enum TimeFormat { TIME_FORMAT_NUM, TIME_FORMAT_YMD };                         // Different type of time format

enum IndSus { BOTH_SUSCEPTIBLE, ONLY_PROPOSE_SUSCEPTIBLE         // Classifies if individual is susceptible in initial/propose
						, NOT_SUSCEPTIBLE }; 

const double VTINY = 0.00000000000000001;                        // Used to represent a very tiny number
const double TINY = 0.00000001;                                  // Used to represent a tiny number
const double LARGE = 10000000;                                   // Used to represent a big number
const unsigned int UNSET = 999999999;                            // A large unsigned integer to represent "Unset"
const unsigned int THRESH = 999999998;                           // Represents a number is under the threshold set
const unsigned int UNKNOWN = 999999997;                          // Represents a number is unknown
const unsigned int ADD = 999999996;                              // Used for summing posterior plots

const double MINIMUM_VARIANCE = 5; 													  	 // The minimum variance for the observation model
const double INVT_POWER = 4;                                     // The power used for the temerature progression

const bool checkon = false;                                  // Set to one to check algorithm is performing correctly

const unsigned int BIN=50;                                       // Number of bins used for plotting distributions

const unsigned int smooth_spline = 1;                            // Set to 1 if smoothing of splines is being done
const double smooth = 0.4;                                       // Used for the smoothing priors on splines

#endif
