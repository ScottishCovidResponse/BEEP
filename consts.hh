#ifndef BEEPMBP__CONSTS_HH
#define BEEPMBP__CONSTS_HH

#include <string>

using namespace std;

const bool sim_only = false;                           // Set if looking at an approach which only uses simualtion

const int VALUE_LIMIT = 10;                            // This limits the size of transition take (data is amalgamated)

const bool langevin_on = false;                        // Set to true if langevin shift is used when performing MVN proposals

const bool obsmodel_transmean = true;                  // Set to true if the trans mean is used in the observation model

const bool suppress_output = true;                     // Suppresses overly verbose output

const bool diagnotic_output = false;                   // Set to true to display additional diagnostic information

const bool checkon = false;                            // Checks algorithm is performing correctly

const bool all_diagnostics = false;                    // Set to true to show all disgnostic information (for testing) 

const bool plot_sim_param = true;                      // Determines if simulated parameter values are put on plots

enum Mode { SIM, MULTISIM, COUNTER, PPC, ABC_SIMPLE,   // Different modes of operation 
           ABC_SMC, ABC_MBP, ABC_MBP_GR, MC3_INF, PAIS_INF, PMCMC_INF};       

enum SimInf { SIMULATE, INFERENCE};                    // Determines if simulation or inference is being performed

enum Dist { EXP_DIST, ERLANG_DIST, GAMMA_DIST, LOGNORM_DIST, INFECTION_DIST};  // Different time distributions

enum ObsModelFunc { NORMAL_OBSMODEL, POISSON_OBSMODEL, NEGBINO_OBSMODEL, SCALE_OBSMODEL}; 
	
enum Dir { X,Y};                                       // Different directions areas sorted by

enum Timers { TIME_TOTAL, TIME_SELF, TIME_MBP, TIME_MBPINIT, TIME_TRANSNUM, TIME_UPDATEPOP, TIME_UPDATEIMAP, TIME_OBSMODEL, TIME_ALG, TIME_MCMCPROP, TIME_WAIT, TIME_GEN, TIME_FIXEDTREE, TIME_SLICETIME, TIME_MEANTIME, TIME_NEIGHBOUR, TIME_JOINT, TIME_SIGMA, TIME_MVN, TIME_RESULTS, TIME_OBSPROB, TIME_PMCMCLIKE, TIME_BOOTSTRAP, TIME_SIMULATE, TIME_PMCMCSWAP, TIME_STATESAMPLE, TIME_SETPARAM, TIME_TRANSMEAN, TIME_INITFROMPART, TIME_SWAP, TIMERMAX};

enum GraphType { GRAPH_TIMESERIES, GRAPH_MARGINAL };
	
enum ObsCombineType { COMBINE_VALUE, COMBINE_SECTION};
	
enum ParamType { DISTVAL_PARAM, BRANCHPROB_PARAM, INF_PARAM, RE_PARAM, SIGMA_PARAM, OBS_PARAM, // Different types of parameters
									 COVAR_PARAM, COVAR_PARAM2, R_EFOI_PARAM, GEOMIX_PARAM, MODIFY_PARAM,
									 SUSCEPTIBILITY_PARAM, DEMO_SPECIFIC_PARAM, ALL_PARAM, PARAMTYPEMAX};
	
enum Status { SUCCESS, FAIL};                                    // Determines if successful

enum TimeFormat { TIME_FORMAT_NUM, TIME_FORMAT_YMD };            // Different type of time format

enum Simu_or_mbp {SIMU, MBP};                                    // Whether simulating or doing a MBP

enum MatModType { ROW_COLUMN};
		
enum IndSus { BOTH_SUSCEPTIBLE, ONLY_PROPOSE_SUSCEPTIBLE         // Classifies if individual is susceptible in initial/propose
						, NOT_SUSCEPTIBLE }; 

enum Smooth { NOSMOOTH, LOGSMOOTH, SMOOTH };
						
enum PropType { MVN_PROP, FIXEDTREE_PROP, SLICETIME_PROP,        // Different types of proposal
MEAN_TIME_PROP, NEIGHBOUR_PROP, JOINT_PROP, SELF_PROP};
					
enum MVNType { MULTIPLE, SINGLE };                               // Different types of MVN proposal
					
enum OpType { OUTPUT, DATA};                                     // Different types of datatable
 
enum DataType { TRANS, POP, MARGINAL };                          // Different types of data
 
enum ErrorBalance { NOT_BALANCED, BALANCED };                    // Determines if error is balanced across data sources
const ErrorBalance error_balance = NOT_BALANCED; //BALANCED;
 
enum ModelType { IND_MODEL, POP_MODEL};                          // The type of the model
const ModelType modeltype = POP_MODEL;

enum ObsModelMode { CUTOFF, INVT};                               // Determines how observation model is incorporated

enum MbpSimType { ALL_MBP, FIXEDTREE, SLICETIME};                // Different MBP-type proposals
								
enum InfUpdate { INF_UPDATE, INF_DIF_UPDATE};                    // Different ways to update infectivity map								

enum PriorType { FIXED_PRIOR, UNIFORM_PRIOR, EXP_PRIOR};         // Different types of prior

enum JointType { UP_DOWN, SINE};                                 // Different type of joint proposal
	
enum OutputPlotType { OP_GRAPH, OP_GRAPH_MARGINAL, OP_MARGINAL,  // The type of output plot
              OP_PARAMDIST, OP_SPLINE, OP_GENERATION, OP_LOG_GENERATION, OP_CPU};

enum LineType { NOLINE, RED_SOLID, RED_DASHED, GREEN_SOLID,      // Different styles of line
             GREEN_DASHED, BLUE_SOLID, BLUE_DASHED, BLACK_SOLID,BLACK_DASHED,BLACK_THIN};

enum DistPropType { BINNING, KDE };                              // Different ways to display probability distributions
						 
enum CounterFactType { CF_PPC, CF_TRANS_RATE, CF_EFOI};          // DIfferent types of counterfactual change

enum ParamUpdate { NO_UPDATE, SLOW_UPDATE, FAST_UPDATE};         // Whether to update parameters with MH

const double fac_up_invT = 1.05, fac_down_invT = 0.9;            // These are used to dynamically change invT
const double fac_up_invT_fast = 1.05, fac_down_invT_fast = 0.9;    

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
const unsigned int ADD = 999999996;                              // Used for summing posterior plots

const unsigned int initial_sample_try = 10000;                   // The number of tried to generate initial state before fail
const unsigned int spline_sample_try = 100000;                   // The number of tried to generate spline before fail

const double Tpower = 4;                                         // The power used for the temerature progression

const unsigned int graph_step = 2;                               // The step size used when plotting the posterior
#endif
