#ifndef BEEPMBP__CONSTS_HH
#define BEEPMBP__CONSTS_HH

using namespace std;

enum Mode { sim, inf, multisim, abcsmc, abcmbp, combinetrace};                                    // Different modes of operation 
	
enum Dist { exp_dist, gamma_dist, lognorm_dist, infection_dist, timep_dist};              // Different time distributions

enum ParamType { other_paramtype, distval_paramtype, branchprob_paramtype}; 
	
enum DataType { trans_data, pop_data, marg_data};                                         // Different types of data

enum TimeFormat { tform_num, tform_ymd };                                                 // Different type of time format

enum IndSus { both_sus, ponly_sus, not_sus };                                             // Use to classify if individual is susceptible in initial/proposed state

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

const unsigned int MAX_NUMBERS = 15000000;                       // The maximum buffer size for Send Recv MPI messages

const unsigned int BIN=50;                                       // The number of bins used for plotting probability distributions

const unsigned int smooth_spline = 1;                            // Set to 1 if smoothing on splines is being done
const double smooth = 0.4;                                       // Used for the smoothing priors on splines

#endif
