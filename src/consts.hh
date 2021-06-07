#ifndef BEEPMBP__CONSTS_HH
#define BEEPMBP__CONSTS_HH

#include <string>

using namespace std;

const bool sim_only = false;                           // Set if looking at an approach which only uses simualtion

const bool langevin_on = false;                        // Set to true if langevin shift is used when performing MVN proposals

const bool obsmodel_transmean = true;                  // Set to true if the trans mean is used in the observation model

const bool suppress_output = true;                     // Suppresses overly verbose output

const bool diagnotic_output = false;                   // Set to true to display additional diagnostic information

const bool checkon = false;                            // Checks algorithm is performing correctly

const bool all_diagnostics = false;                    // Set to true to show all disgnostic information (for testing) 

const bool plot_sim_param = true;                      // Determines if simulated parameter values are put on plots

enum Mode { SIM, MULTISIM, PREDICTION,                 // Different modes of operation 
            ABC_SIMPLE, ABC_SMC, ABC_MBP, MC3_INF, MCMC_MBP, PAIS_INF, PMCMC_INF};       

enum SimInf { SIMULATE, INFERENCE};                    // Determines if simulation or inference is being performed

enum TransInf { TRANS_INFECTION, TRANS_NOTINFECTION};  // Determines if the transition is an infection or not

enum ErlangPos { FIRST,LAST};                          // Position in Erland distribution

enum ObsModelFunc { NORMAL_OBSMODEL, POISSON_OBSMODEL, NEGBINO_OBSMODEL, SCALE_OBSMODEL}; 
	
enum Dir { X,Y};                                       // Different directions areas sorted by

enum Timers { TIME_TOTAL, TIME_SELF, TIME_MBP, TIME_MBPINIT, TIME_TRANSNUM, TIME_UPDATEPOP, TIME_UPDATEIMAP, TIME_OBSMODEL, TIME_ALG, TIME_MCMCPROP, TIME_WAIT, TIME_GEN, TIME_FIXEDTREE, TIME_SLICETIME, TIME_MEANTIME, TIME_NEIGHBOUR, TIME_JOINT, TIME_SIGMA, TIME_MVN, TIME_RESULTS, TIME_OBSPROB, TIME_PMCMCLIKE, TIME_BOOTSTRAP, TIME_SIMULATE, TIME_PMCMCSWAP, TIME_STATESAMPLE, TIME_SETPARAM, TIME_TRANSMEAN, TIME_INITFROMPART, TIME_SWAP, TIMERMAX};

enum GraphType { GRAPH_TIMESERIES, GRAPH_MARGINAL };
	
enum ParamType { DISTVAL_PARAM, BRANCHPROB_PARAM,               // Different types of parameters 
                   INF_PARAM, RE_PARAM, SIGMA_PARAM, OBS_PARAM,
									 COVAR_PARAM, R_EFOI_PARAM, GEOMIX_PARAM, MODIFY_PARAM, SUSCEPTIBILITY_PARAM, DEMO_SPECIFIC_PARAM, 
									 LEVEL_EFFECT_PARAM, RFACTOR_PARAM, AREA_EFFECT_PARAM, ALL_PARAM, PARAMTYPEMAX};
	
enum Status { SUCCESS, FAIL};                                    // Determines if successful

enum TimeFormat { TIME_FORMAT_NUM, TIME_FORMAT_YMD,             // Different type of time format
									TIME_FORMAT_DMY_SLASH, TIME_FORMAT_DMY_DOT};  

enum Simu_or_mbp {SIMU, MBP};                                    // Whether simulating or doing a MBP

enum MatModType { ROW_COLUMN};
		
enum IndSus { BOTH_SUSCEPTIBLE, ONLY_PROPOSE_SUSCEPTIBLE         // Classifies if individual is susceptible in initial/propose
						, NOT_SUSCEPTIBLE }; 

enum Smooth { NOSMOOTH, LOGSMOOTH, SMOOTH };
						
enum PropType { MVN_PROP, FIXEDTREE_PROP, SLICETIME_PROP,        // Different types of proposal
MEAN_TIME_PROP, NEIGHBOUR_PROP, JOINT_PROP, SELF_PROP};
					
enum MVNType { MULTIPLE, SINGLE };                               // Different types of MVN proposal
					
enum OpType { OUTPUT, DATA};                                     // Different types of datatable
 
enum CovarType { AREA_COVAR, TV_COVAR, AREA_TV_COVAR};           // Different types of covariate

enum Transform { LOG_TRANS, LINEAR_TRANS};                       // Determines how the raw data is transformed

enum DataType { TRANS, POP, POPFRAC, MARGINAL };                 // Different types of data
 
enum NumberType { POPULATION, FRACTION, PERCENT};                // The different types of number

enum ErrorBalance { NOT_BALANCED, BALANCED };                    // Determines if error is balanced across data sources
const ErrorBalance error_balance = NOT_BALANCED; //BALANCED;
 
enum ObsModelMode { CUTOFF, INVT};                               // Determines how observation model is incorporated

enum MbpSimType { ALL_MBP, FIXEDTREE, SLICETIME};                // Different MBP-type proposals
								
enum InfUpdate { INF_UPDATE, INF_DIF_UPDATE};                    // Different ways to update infectivity map								

enum PriorType { FIXED_PRIOR, UNIFORM_PRIOR, EXP_PRIOR,          // Different types of prior
								 NORMAL_PRIOR, DIRICHLET_PRIOR}; 
				
enum JointType { UP_DOWN, SINE};                                 // Different type of joint proposal
	
enum OutputPlotType { OP_GRAPH, OP_GRAPH_MARGINAL, OP_MARGINAL,  // The type of output plot
              OP_PARAMDIST, OP_SPLINE, OP_GENERATION, OP_LOG_GENERATION, OP_CPU, OP_TRACE, OP_ME};

enum LineType { NOLINE, RED_SOLID, RED_DASHED, GREEN_SOLID,      // Different styles of line
                GREEN_DASHED, BLUE_SOLID, BLUE_DASHED, BLACK_SOLID, BLACK_DASHED, 
								RED_DOTTED, RED_DOTDASH, GREEN_DOTTED,     
                GREEN_DOTDASH, BLUE_DOTTED, BLUE_DOTDASH, BLACK_DOTTED, BLACK_DOTDASH,
								YELLOW_SOLID, YELLOW_DASHED, CYAN_SOLID, CYAN_DASHED,MAGENTA_SOLID, MAGENTA_DASHED};

enum LineColour { UNSET_COLOUR, RED, GREEN, BLUE, YELLOW, CYAN, MAGENTA, BLACK};// Different colours of line 
	
enum DistPropType { BINNING, KDE };                              // Different ways to display probability distributions
						 
enum ModificationType { CF_TRANS_RATE, CF_EFOI,                  // DIfferent types of model modification
												CF_SPLINEFAC, CF_SPLINESET, CF_BETA};

enum ParamUpdate { NO_UPDATE, SLOW_UPDATE, FAST_UPDATE};         // Whether to update parameters with MH

const double fac_up_invT = 1.05, fac_down_invT = 0.9;            // These are used to dynamically change invT
const double fac_up_invT_fast = 1.5, fac_down_invT_fast = 0.7;    

const double fac_up = 1.1, fac_down = 0.95;                      // These are used to dynamically change the size of proposals
const double fac_up_fast = 1.5, fac_down_fast = 0.7;  

const double sizemax = 2;                                        // The maximum size of parameter proposals
const double sizemin = 0.2;                                      // The minimum size of parameter proposals

const double VTINY = 0.00000000000000001;                        // Used to represent a very tiny number
const double TINY = 0.00000001;                                  // Used to represent a tiny number
const double SMALL = 0.00001;                                    // Used to represent a small number
const double LARGE = 1000000000;                                 // Used to represent a big number
const unsigned int UNSET = 999999999;                            // A large unsigned integer to represent "Unset"
const unsigned int THRESH = 999999998;                           // Represents a number is under the threshold set
const unsigned int UNKNOWN = 999999997;                          // Represents a number is unknown
const unsigned int TOBESET = 999999996;                          // A large unsigned integer to represent "Unset"

const unsigned int initial_sample_try = 10000;                   // The number of tries to generate initial state before fail
const unsigned int spline_sample_try = 100000;                   // The number of tries to generate parameters before fail
const unsigned int sample_try = 10000;                           // The number of tries to generate spline before fail
const unsigned int initialise_param_samp = 100;                  // Number of random parameter samples to initialise param_samp

const unsigned int graph_step = 2;                               // The step size used when plotting the posterior
#endif
