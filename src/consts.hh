#ifndef BEEP__CONSTS_HH
#define BEEP__CONSTS_HH

#include <string>

using namespace std;

const bool sim_only = false;                           // Set if looking at an approach which only uses simualtion

const bool langevin_on = false;                        // Set to true if langevin shift is used when performing MVN proposals

const bool obsmodel_transmean = false;                  // Set to true if the trans mean is used in the observation model

const bool suppress_output = true;                     // Suppresses overly verbose output

const bool diagnotic_output = false;                   // Set to true to display additional diagnostic information

const bool checkon = false;                            // Checks algorithm is performing correctly

const bool all_diagnostics = false;                    // Set to true to show all disgnostic information (for testing) 

const bool plot_sim_param_global = true;               // Determines if simulated parameter values are put on plots
 
const double power_obsmodel = 0.7;                     // The power used in the power observation model

const double pmcmc_start_param = false;                // Set to true if pmcmc starts on true parameter value

const double cor_update_num = 0.8;                     // The value of cor used to make changes to num_updates  
const double num_updates_max = 10;                     // The maximum number for num_updates
const double num_updates_min = 1;                      // The minimun number for num_updates

enum Mode { SIM, MULTISIM, PREDICTION,                 // Different modes of operation 
            ABC_SIMPLE, ABC_SMC, ABC_MBP, ABC_DA, ABC_CONT, MC3_INF, MCMC_MBP, PAIS_INF, PMCMC_INF, IMPORTANCE_INF, ML_INF,
						DATAONLY};       

enum SimInf { SIMULATE, INFERENCE, DATAVIEW};          // Determines if simulation or inference is being performed

enum TransInf { TRANS_INFECTION, TRANS_NOTINFECTION};  // Determines if the transition is an infection or not

enum ErlangPos { FIRST, LAST};                          // Position in Erland distribution

enum ObsModelFunc { NORMAL_OBSMODEL, NORMAL_PERCENT_OBSMODEL, POISSON_OBSMODEL, NEGBINO_OBSMODEL}; 
	
enum ObsType { OBS_EXACT, OBS_APPROX };                // Approximate or exact observation (used in likelihood approx)
	
enum Dir { X,Y};                                       // Different directions areas sorted by

enum Timers { TIME_TOTAL, TIME_SELF, TIME_MBP, TIME_MBPINIT, TIME_TRANSNUM, TIME_UPDATEPOP, TIME_UPDATEIMAP, TIME_OBSMODEL, TIME_ALG, TIME_MCMCPROP, TIME_PARAMPROP, TIME_STATEPROP, TIME_UPDATE, TIME_WAIT, TIME_GEN, TIME_FIXEDTREE, TIME_SLICETIME, TIME_MEANTIME, TIME_NEIGHBOUR, TIME_JOINT, TIME_COVAR_AREA, TIME_SIGMA, TIME_MVN, TIME_RESULTS, TIME_OBSPROB, TIME_PMCMCLIKE, TIME_BOOTSTRAP, TIME_SIMULATE, TIME_PMCMCSWAP, TIME_STATESAMPLE, TIME_SETPARAM, TIME_TRANSMEAN, TIME_INITFROMPART, TIME_SWAP, TIME_CREATEN, TIME_BETA_FROM_R, TIME_CUTOFF, TIME_PROP, TIME_MVNSETUP, TIME_MBPUPDATE, TIME_UPDATEPROP, TIME_LIKELIHOOD_APPROX, TIME_OBS_APPROX, TIME_FUTURE_OBS_APPROX, TIME_COVAR, TIME_CORRECT, TIME_GRAD, TIME_POSTERIOR_SAMPLE, TIME_EF_CALCULATE, TIME_TEMP1, TIME_TEMP2, TIME_TEMP3, TIME_TEMP4, TIME_INV_MATRIX, TIME_DETERMINANT, TIME_MATRIX_MULT, TIME_LINEAR_EQ, TIME_ADD_REMOVE_S, TIME_GENERATE_SAMPLES, TIME_CMAES, TIME_SCALE_COVARIANCE, TIMERMAX};

enum Accuracy {DOUBLE, FLOAT};                                  // Sets the level computation

enum GraphType { GRAPH_TIMESERIES, GRAPH_MARGINAL };
	
enum ParamType { DISTVAL_PARAM, BRANCHPROB_PARAM,               // Different types of parameters 
                   INF_PARAM, RE_PARAM, SIGMA_PARAM, OBS_PARAM,
									 COVAR_PARAM, R_EFOI_PARAM, GEOMIX_PARAM, MODIFY_PARAM, SUSCEPTIBILITY_PARAM, DEMO_SPECIFIC_PARAM, 
									 LEVEL_EFFECT_PARAM, RFACTOR_PARAM, AREA_EFFECT_PARAM, ALL_PARAM, DIST_R_JOINT, R_NEIGH, PARAM_PROP, PARAMTYPEMAX};
	
enum Status { SUCCESS, FAIL};                                    // Determines if successful

enum TimeFormat { TIME_FORMAT_NUM, TIME_FORMAT_YMD,             // Different type of time format
									TIME_FORMAT_DMY_SLASH, TIME_FORMAT_DMY_DOT};  

enum Simu_or_mbp {SIMU, MBP};                                    // Whether simulating or doing a MBP

enum MatModType { ROW_COLUMN, PERTURB};
		
enum IndSus { BOTH_SUSCEPTIBLE, ONLY_PROPOSE_SUSCEPTIBLE         // Classifies if individual is susceptible in initial/propose
						, NOT_SUSCEPTIBLE }; 

enum SmoothType { NOSMOOTH, LOGSMOOTH, SMOOTH };                 // The type of smoothing used for a spline
						
enum PropType { MVN_PROP, FIXEDTREE_PROP, SLICETIME_PROP,        // Different types of proposal
MEAN_TIME_PROP, NEIGHBOUR_PROP, JOINT_PROP, COVAR_AREA_PROP, SELF_PROP};
					
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
								 NORMAL_PRIOR, DIRICHLET_PRIOR, DIRICHLET_ALPHA_PRIOR, DIRICHLET_FLAT_PRIOR,
								 MDIRICHLET_PRIOR}; 
				
enum MLAlg { ML_GD, ML_CMAES};                                   // Different types of maximum likelihood algorithm

enum DirType { DIR_NORM, DIR_MODIFIED};                          // The type of Dirichlet distribution

enum MatrixType { DIAG, FULL};

enum JointType { UP_DOWN, SINE};                                 // Different type of joint proposal
	
enum OutputPlotType { OP_GRAPH, OP_GRAPH_MARGINAL, OP_MARGINAL,  // The type of output plot
              OP_PARAMDIST, OP_SPLINE, OP_GENERATION, OP_LOG_GENERATION, OP_CPU, OP_TRACE, OP_ME, 
							OP_ANIM_MAP, OP_AGE_MATRIX, OP_AGE_MATRIX_POST, OP_AGE_MATRIX_DIF, OP_PARAM_TABLE, OP_PRIOR_TABLE, OP_COMP_MODEL, OP_FOI_MODEL,OP_MIXING_WITHIN_ANIM_MAP, OP_MIXING_WITHIN_MAP, OP_MIXING_BETWEEN_ANIM_MAP, OP_MIXING_BETWEEN_MAP, OP_MIXING_POP_ANIM_MAP, OP_MIXING_POP_MAP, OP_DESC,
							OP_AREA_COVAR, OP_TV_COVAR, OP_AREA_TV_COVAR, OP_LEVEL_EFFECT};

enum LineType { NOLINE, RED_SOLID, RED_DASHED, GREEN_SOLID,      // Different styles of line
                GREEN_DASHED, BLUE_SOLID, BLUE_DASHED, BLACK_SOLID, BLACK_DASHED, 
								RED_DOTTED, RED_DOTDASH, GREEN_DOTTED,     
                GREEN_DOTDASH, BLUE_DOTTED, BLUE_DOTDASH, BLACK_DOTTED, BLACK_DOTDASH,
								YELLOW_SOLID, YELLOW_DASHED, CYAN_SOLID, CYAN_DASHED,MAGENTA_SOLID, MAGENTA_DASHED,
								RED_THIN, GREEN_THIN, BLUE_THIN, BLACK_THIN};

enum LineColour { UNSET_COLOUR, RED, GREEN, BLUE, YELLOW, CYAN, MAGENTA, BLACK};// Different colours of line 
	
enum Project { EQUI_PROJ, UNIFROM_PROJ};                         // Different projections for maps

enum FileType{ KML, GEOJSON};                                    // DIfferent types of file
	
enum DistPropType { BINNING, KDE };                              // Different ways to display probability distributions
						 
enum ModificationType { CF_TRANS_RATE, CF_EFOI,                  // DIfferent types of model modification
						CF_SPLINEFAC, CF_SPLINESET, CF_BETA};

enum ParamUpdate { NO_UPDATE, COMBINE_UPDATE,                    // Whether to update parameters with MH
									 COMBINE_DYNAMIC_UPDATE, SLOW_UPDATE, FAST_UPDATE};

enum StateUncertainty{ CI, CURVES};                              // Determines how state uncertainty plotted

const double eta_fast = 1, eta = 0.1;
const double eta_pmcmc_fast = 0.1, eta_pmcmc = 0.01;
const double eta_combine = 4;
const double eta_invT = 0.1;

const double ac_rate = 0.33;                                     // Target acceptance rate
//const double ac_rate_tree = 0.5;                               // Target acceptance rate for trees
const double ac_rate_tree = 0.33;                                // Target acceptance rate for trees
const double ac_rate_slice = 0.33;                               // Target acceptance rate for slicetime
const double ac_rate_self = 0.7;                                 // Target acceptance rate for self proposals (PMCMC);

const double sizemax = 1.5;                                        // The maximum size of parameter proposals

const auto pmcmc_init_samp = 100u;                               // The number of initial PMCMC samples

const double self_cor_PMCMC = 0.99;                              // Gives the correlation in self acceptance when using PMCMC

const double VTINY = 0.000000000000001;                          // Used to represent a very tiny number
const double TINY = 0.00000001;                                  // Used to represent a tiny number
const double SMALL = 0.00001;                                    // Used to represent a small number
const double LARGE = 1000000000;                                 // Used to represent a big number
     
const unsigned int UNSET = 999999999;                            // A large unsigned integer to represent "Unset"
const unsigned int THRESH = 999999998;                           // Represents a number is under the threshold set
const unsigned int UNKNOWN = 999999997;                          // Represents a number is unknown
const unsigned int TOBESET = 999999996;                          // A large unsigned integer to represent "Unset"
const unsigned int ITERATE_GENERATION = 999999995;               // In CMAES iterate generations

const unsigned int initial_sample_try = 10000;                   // The number of tries to generate initial state before fail
const unsigned int spline_sample_try = 100000;                   // The number of tries to generate parameters before fail
const unsigned int sample_try = 10000;                           // The number of tries to generate spline before fail
const unsigned int initialise_param_samp = 100;                  // Number of random parameter samples to initialise param_samp

const unsigned int ML_GENERATION_TERM_COND = 10;                 // The number of generation used in termination (CMAES)

const double map_ratio = 1.22;                                   // The ratio of the map (used when plotting

#endif
