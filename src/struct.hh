/// Stores comonly used data structures

#ifndef BEEP__STRUCT_HH
#define BEEP__STRUCT_HH

#include <vector>

using namespace std;

#include "consts.hh"

class Details;                              // Makes prototypes for the main classes
class Data;
class Model;
class ObservationModel;
class Output;
class Mpi;
class State;
class MVN;
class Inputs;
class Model;
class Mbp;

struct SplineOutput{                       // Used for outputting splines
	string name;                             // The name of the spline
	string desc;                             // The description of the spline
	string fulldesc;                         // A full description of the spline
	string tab, tab2, tab3, tab4;            // Classifiers for menu
	vector <double> splineval;               // The values in the spline
	string spline_param_JSON;                // Stores information about parameters used to inform spline
};

struct DerivedParam{                       // Parameters derived from the model parameters
	string name;                             // The name of the derived parameter
	string desc;                             // The description of the derived parameter
	string file;                             // The output file name
	double value;                            // The value
};

struct RMap{                               // Spatial maps for the reproduction number
	string file;                             // The file of the Rmap
	string fulldesc;                         // The description of the RMap
	string tab, tab2, tab3, tab4;            // Classifiers for menu
	vector < vector <double> > map;
};

struct Sample{                             // Stores information about a state posterior sample
	vector < vector <double> > graph_state;  // Stores the underlying graph state
	vector <SplineOutput> spline_output;     // The splines that need to be output
	vector <DerivedParam> derived_param;     // Derived parameters
	vector <RMap> Rmap;                      // The maps in the reproduction number
	vector < vector < vector <double> > > Nsample;// A sample of the age mixing matrix
};

struct ParamSample{                      	 // Stores information about a parameter sample from the posterior
	unsigned int run;                        // The run from which the samples were generated
	double EF;                               // Stores the error function 
	vector <double> paramval;                // A parameter sample
};

struct ParamSplineRef{                     // References any splines the parameter is within
	unsigned int spline;
	unsigned int i;
};

struct ParamSpec {                         // Gives the specification for a parameter and/or a prior
	string name;
	string value;
	string prior;
	double smooth;
	double factor;
};

struct Param{                            	 // Store information about a model parameter
 	string name;                             // The name
	PriorType priortype;                     // The type of the prior
	double val1;                             // Value for the prior (the minimum value assuming a uniform prior) 
	double val2;                             // Value for the prior (the maximum value assuming a uniform prior)
	unsigned int val1_param;                 // Gives the parameter giving that value
	unsigned int val2_param;                 // Gives the parameter giving that value
	double dir_mean;                         // The mean of a Dirichlet prior
	double dir_sd;                           // The standard deviatiion of a Dirichlet prior
	vector <unsigned int> param_dep;         // Other parameters depenedent on this parameter
	double value;                            // The value used for the simulation
	ParamType type;                          // Set to one if parameter for transition times and 2 for branching probabilities
	ParamSpec ps;                            // The parameter specification (this is needed to generate Erlang)
	vector <unsigned int> greaterthan;       // List of parameters which this one needs to be greater than
	vector <unsigned int> lessthan;          // List of parameters which this one needs to be greater than
	vector <ParamSplineRef> paramsplineref;  // References any splines the parameter is in
	double fixed;                            // This is used to store the fixed value of a parameter (if fixed prior is used)
};

struct Compartment{                        // Stores information about a compartment in the model
	string name;                             // The compartments name
	string name_num;                         // The compartment name and number
	unsigned int num;                        // The position the distribution has in an Erlang distribution
	unsigned int shape;                      // The shape (1 if Exp else Erlang)
 	string mean_dep;                         // The demographic dependency for the mean
	vector <ParamSpec> mean_spec;            // Parameters specificiations for means
	vector <unsigned int> param_mean;        // The parameter for the mean of the distribution (potentially age dependant
	ParamSpec inf_spec;                      // Parameter specification for infectivity
	unsigned int infectivity_param;          // How infectious that compartment is
	vector <unsigned int> trans;             // The transitions leaving that compartment
};

struct CompartmentName{
	string name;                             // The compartments name
	vector <unsigned int> comp;              // List the compartments with this name
};

struct Transition{                         // Stores information about a compartmental model transition
	string name;                             // The name of the transition
	string name_file;                        // The name of the transition used for file (does not have > character)
	unsigned int from;                       // Which compartment the individual is coming from
	unsigned int to;                         // Which compartment the individual is going to
	TransInf inf;                            // Determines if an infection transition
	unsigned int comptrans_ref;              // The reference to comptrans in from
	string prob_dep;                         // The demographic dependency for the mean
	vector <ParamSpec> prob_spec;            // Parameters specificiations for brancing probabilities
	vector <unsigned int> param_prob;        // The parameter for the probability of going down transition (age dependant)
};

struct CompApprox{                         // A unique compartment used in the approximate method
	unsigned int c;                          // The area
	unsigned int co;                         // The compartment
	unsigned int dp;                         // The demographic possibility
};

struct CompCopy{                           // Used to split merged compartments (used for ML approaches)
	unsigned int c;                          // The number of the compartment being copied
	unsigned int c_copy;                     // The new compartment
	unsigned int index;                      // The index for the transition 
};
	
struct TransApprox{                        // A unique compartment used in the approximate method
	unsigned int c;                          // The area
	unsigned int tr;                         // The transition
	unsigned int dp;                         // The demographic possibility
	unsigned int ca_from;                    // The unique compartment from
	unsigned int ca_to;                      // The unique compartment to
};

struct CompProb{                           // Gives compartmental probabilities
	vector <double> value;                   // The probability an indidividuals goes to compartment (used to calculate R)
};

struct SplinePoint{                        // Stores information about a spline point
	double t;                                // The time of the point
	unsigned int param;                      // The parameter which defines the value
	double multfac;                          // A multiplicative factor
	double smooth_value;                     // The value used for smoothing
};

struct Spline{                             // Stores all the information about a spline
	string name;                             // The names of the spline
	string desc;                             // The description of the spline
	SmoothType smoothtype;                   // Determines the type of smoothing
	vector <SplinePoint> p;                  // The points on the spline
	unsigned int param_factor;               // A parameter whihc can multiply the spline
};

struct SplinePriorGrads{                   // Stores the gradients of spline prior
	double dPr_dth;                          // The derivative of the prior
	double d2Pr_dth2;                        // The double derivative of the prior
};

struct SplineInfo{                         // Provides information about R_spline and efoi_splines#
	string name;                             // The name of the spline
	unsigned int spline_ref;                 // A reference to the spline number
	vector <unsigned int> area;              // The areas the spline applies to 
	vector <double> efoi_agedist;            // The age distribution for the force of infection
	vector <double> democatpos_dist;         // The demographic distributions for the selected areas	
	unsigned int total_pop;                  // The total population
};

struct Particle
{
	vector <double> paramval;                // The parameter values for the particle
	vector <vector < vector < vector <double> > > > transnum; // Transition numbers 
	double EF;                               // The value of the error function
	unsigned int run;                        // The run the particle belongs to 

	ParamSample create_param_samp() const
	{
		ParamSample ps; ps.paramval = paramval; ps.run = run; ps.EF = EF;
		return ps;
	}
};

struct ParticleApprox                      // A particle used when generting posterior samples using the approximate method
{
	vector < vector < vector < vector <double> > > > transnum; // Transition numbers 
	vector < vector < vector < vector <double> > > > pop_store;// The populations numbers 
	vector < vector < vector <double> > > pop_end;             // The population at the end
	vector < vector < vector <double> > > pop_end_new;         // The population at the end (new version)
	double L;                                // The likelihood (used for bootstrap)
	double wsum;                             // The sum of weights
};

struct Generation                          // Stores inforamtion about a generation when doing ABC methods
{
	unsigned int num;                        // The number of the generation                 
	vector <ParamSample> param_samp;         // Parameter samples generated
	vector < vector <double> > EF_datatable; // Stores how the error function is divided into datatables
	vector <unsigned int> partcopy;          // Shows which are copied particles
	vector <double> w;                       // The weight of parameter samples (used in ABC-SMC)
	double EFcut;                            // The error function cut-off used (in CUTOFF mode)
	double EFmin, EFmax;                     // Diagnostic information about EF range
	double invT;                             // The inverse temperature used (in INVT mode)
	double time;                             // The execution time
	vector <double> model_evidence;          // Sets the model evidence for each run
};

struct ProbReach                           // The probability of reachin a certain compartment
{
	string name;                             // The name (shown in the output file)
	unsigned int comp;                       // The compartment
};

struct TimePlot {                          // A label can be used for specific dates
	string name;
	unsigned int time;
	string time_str;
};

struct Matrix {                            // Loads a matrix
	unsigned int N;													 // The size of the matrix
	vector <vector <double> > ele;           // The elements of the matrix
	double norm_factor;                      // The factor used when matrix is normalised
};

struct SparseMatrix {                      // Loads a matrix
	unsigned int N;													 // The size of the matrix
	vector <double> diag;                    // The diagonal contribution to the matrix
	vector < vector <unsigned int> > to;     // The value of the element to
	vector < vector <double> > val;          // The the non-diagonal elements of the matrix
};

struct MatrixElement {                     // Identifies an element of a matrix
	unsigned int j;                          // y position in matrix
	unsigned int i;                          // x position in matrix
};

struct MatrixModification {                // Used to modify the mixing matrix
	MatModType type;	                       // The type of modification
	string name;                             // The name of the modification
	string desc;                             // A description of the modification
	vector <ParamSpec> ps;                   // The parameter specifications
	vector <unsigned int> bp;                // Finds the break points for the spline
	vector <unsigned int> ages;              // The group of ages on which the modification acts
	SmoothType smoothtype;                   // The type of smoothing used for the spline
	unsigned int spline_ref;                 // The spline used to generate the modification 
};

struct TreeNode {                         // Provides information about a node
	vector <unsigned int> arearef;          // References the list of areas within the node
	vector <unsigned int> child;            // The child nodes
};

struct GenerateQ {                         // Stores information about age and geographical mixing
	vector <string> N_name;                  // Stores the name of age mixing matrices to be loaded
	vector <MatrixModification> matmod;      // Store potential modifications to the basic age mixing matrix 
	unsigned int nspline;                    // The number of splines used for a particular type
	string M_name;                           // The geographic mixing matrix name
	vector <Matrix> N;                       // Stores the matrices
	SparseMatrix M;                          // The actual loaded geographic mixing matrix 
	vector <double> onlydiag;                // The elements of a completely diagnoal matrix for geomixing
	vector <double> factor;                  // A factor used for geospline
	vector <TreeNode> treenode;              // Stores a tree of nodes which divides area based on matrix M
};

struct Table {                             // Loads a table
	string file; 														 // The file from which the tables was loaded
	unsigned int ncol;                       // The number of columns
	unsigned int nrow;                       // The number if rows
	vector <string> heading;                 // The headings for the columns
	vector <vector <string> > ele;        	 // The elements of the table
};

struct GeographicMap {                     // Geographic maps allow for different combinations of areas 
	string region;                           // The name of the region
	vector <unsigned int> area;              // The areas selected within a given region
};

struct DataFilter {                        // Selects areas and demographic posibilities for a particular column of data
	string desc;                             // The description for the data filter
 	string file;                             // The filename used when saving the file
	string name;                             // The name of the data filter
	string colname;                          // The column name in the file
	vector <unsigned int> area;              // The areas that the data points refer to 
	vector <unsigned int> dp_sel;            // The demographic possibilities that the data points refer to
};

struct DemocatChange{                      // Determines how a demographic category changes with time
	string name;                             // The name of the demographic category
	string file;                             // The name of the file to be loaded
	string democats_filt;                    // The filter used for the demographic categories
	string geo_filt;                         // The geographic filter being used	
	unsigned int shift;                      // A potential shift in the data
	unsigned int d;                          // Stores which democat
	vector <unsigned int> area;              // The areas that the data points refer to 
	vector < vector <double> > frac;         // The fraction in the different dempgraphic groups
	vector < vector <unsigned int> > dp_group;// The dp groups for the particular demographic category
};

struct DataTable {                         // Stores tables of data provided by the user
	OpType optype;                           // Whether it is a true data table of just an output
	DataType type;                           // The type of the data 
	string observation;                      // The quantity being observed 
	string file;                             // The name of the file to be loaded
	int start;                               // The start time for the data (before shifting)
	int end;                                 // The end time for the data (before shifting)
	unsigned int timestep;                   // The time step for transition measurements (e.g. 1=days, 7=weeks) 
	int shift;                               // How much time the data should be shifted before used
	double threshold;                        // Any threshold that is used 
	string democats_filt;                    // The filter used for the demographic categories
	string democats_dep;                     // The demographic dependency
	string geo_filt;                         // The geographic filter being used
	string geo_dep;                          // The geographic dependency
  string factor_spline;                    // The spline factor
	double factor;                           // The factor which multiplies it
	LineColour line_colour;                  // Stores the colour of the line when plotted
	string plot_name;                        // Stores the name of the plot
	string label;                            // Stores the label for the plot (OutputState only)
	vector <string> demolist;                // For marginal distributions stores demographic categories
	
	vector <unsigned int> translist;         // The transitions relating to the data
	vector <unsigned int> complist;          // The compartments relating to the data
	
	vector <unsigned int> graph_ref;         // A list of all the graphs plotted
	
	ObsModelFunc obsmodel;                   // Determines what observaion model is used
	double shape;                            // The shape parameter
	double percent;                          // This stores percentage error
	double sd;                               // This stores sd of error 
	double epsilon_factor;                   // Stores the factor which multiplies epislon
	bool load_sd;                            // Loads the standard deviation from a file
};

struct Observation {                       // A single observation made on the system
	unsigned int factor_spline;              // References the spline number for the factor
	unsigned int datatable;                  // The data table the observation belongs to
	unsigned int graph;                      // The primary graph the observation belongs to
	double value;                            // The value of the observation
	double sd;                               // This is the standard deviation 
	double factor;                           // The factor which multiplies it (used for POPFRAC)
	unsigned int sett_i;                     // The start time
	unsigned int sett_f;                     // The end time
	vector <unsigned int> area;              // The areas involved
	vector <unsigned int> dp_sel;            // The deomgraphic possibilities selected  	
	ObsModelFunc obsmodel;                   // The observation model used
	double shape;                            // Shape paremeter (used for negative binomial).
	double var_approx;                       // The variance used when doing a normal approimation (for approximate methods)
};

struct ObsSlice {                          // Gives all the observations at a particular point in time
	unsigned int sett;                       // The time of the slice
	vector <unsigned int> obs_ref;           // The observations at that time
};

struct GraphPoint {                        // Details of a point on a graph
	double xi;                               // The x position at the start of the observation
	double xf;                               // The x position at the end of the observation
	unsigned int obs;                        // The observations which contribute to the point
};

struct Graph {                             // Details of a graph
	string name;                             // The name of the graph
	string fulldesc;                         // Full descrition
	string fulldesc_data;                    // Full descrition of the data
	string tab, tab2, tab3, tab4;            // Used for menu
	string file;                             // The filename used for the graph     
	string desc;                             // The description of the graph
	string colname;                          // The column the graph relates to
	GraphType type;                          // The type of graph
	vector <GraphPoint> point;               // The points on the graph
	double factor;                           // The factor which multiplies it (used for POPFRAC)
	unsigned int factor_spline;              // References the spline for the factor
	unsigned int datatable;                  // The data table the observation belongs to
	vector <unsigned int> area;              // The areas involved
	vector <unsigned int> dp_sel;            // The deomgraphic possibilities selected  	
};

struct GraphMultiPlot {                    // Allows for multiple graphs to be placed into a single plot
	string plot_name;                        // The name of the plot
	string fulldesc;                         // Full descrition
	string fulldesc_data;                    // Full descrition of the data
	string tab, tab2, tab3, tab4;            // Used for menu
	vector <string> name;                    // The name on the y label
	GraphType type;                          // The type of graph
	double min;                              // The minimum in the graph
	double max;                              // The maximum in the graph
	vector <LineColour> line_colour;         // Stores the colour of the line when plotted
	vector <string> file;                    // The filenames used for the graph     
	vector <string> file_data;               // The filenames for data used for the graph    
	vector <bool> file_data_EB;              // Set to true if data contains an error bar
	vector <string> file_data_thresh;        // The filenames for thresholds used for the graph   
};
	
struct Coord {                             // Stores coordinates
	double x;                                // x
	double y;                                // y
};

struct DemographicCategory {               // Stores demographic categories
	string name;                             // The name of the category
	bool sus_vari;                           // Set to true if there is variation in susceptibility
	vector <string> value;                   // The postential values it can take
	vector <ParamSpec> ps;                   // Information about the parameter specification
	vector <unsigned int> sus_param;         // The parameter which gives the susceptibility
};

struct Strain{                             // Stores information about a strain
	string name;                             // The name of the strain
 	ParamSpec Rfactor_ps;                    // Information about the parameter specification for strain infectivity
	unsigned int Rfactor_param;              // The parameter which gives the strain infectivity	
};

struct LevelEffect {                       // Stores level effect
	bool on;                                 // Set to true if a level effect is on 
	string file;                             // The file used to 
	vector < vector <unsigned int> > param_map;// Gives the name of the parameter for a certain area at a certain time
	vector <double> frac;                    // The fraction of time spend in each of the different states
	vector <ParamSpec> ps;                   // Gives a list of parameter specifications
};


struct AreaEffect {                        // Stores area effect
	bool on;                                 // Set to true if a level effect is on 
	vector <double> frac;                    // The fraction of individuals in different states
	vector <ParamSpec> ps;                   // Gives a list of parameter specifications
};

struct Covariate {                         // Stores the  covariate for the area
	CovarType type;                          // The type of the covariate
	string file;                             // The file for the covariate
	string name;                             // The name of the covariate (i.e. the column in the area data file)
	ParamSpec ps;                            // Gives a parameter specification
	Transform func;                          // The functional transformation
	bool timevary;                           // Set to true is covariate is time-varying
	vector < vector <double> > value;        // The values for the covariates
};

struct Area {                              // Provides information about an area
	string code;                             // The code for the area
  vector < vector <double> > pop_init;     // The initial population in different comparments and demographic categories  
	vector <double> pop_dp;                  // The population within a demographic category
	unsigned int total_pop;                  // The total population in the area
};

struct PartEF                              // Structure used to order particle EFs
{
	unsigned int i;                          // The number of the particle
	double EF;                               // The error function
};

struct Proposal {                          // Stores a proposal to be done    
	PropType type;                           // The type of the proposal
	unsigned int num;                        // The number of the proposal
};

struct EventRefTime{                       // Event reference used for sorting infection events        
	unsigned int ind;                        // Individual
	unsigned int e;	                         // Event number
	double t;	                               // Time
};

struct OutputProp{                         // Properties used for outputting distributions
	DistPropType type;                       // The type of distribution output
	unsigned int nbin;                       // The number of bins used
	double h;                                // The widths used for KDE
};

struct Modification{                       // Used to implement a model modification (NOTE implemented in mpi.copy_data)
	unsigned int start;                      // The start time of a change
	unsigned int end;                        // The end time for a change
	ModificationType type;                   // The type of change
	string spline_name_str;                  // The name of the spline
	string trans_str;                        // The transition
	string strain_str;                       // The strain
	string geo_filt;                         // The filter for geography
	string democats_filt;                    // The filter for demographic categories
	double factor;                           // The factor change
	vector <unsigned int> dp_sel;            // The demographic possibilities which are selected
	vector <unsigned int> area;              // The areas which are selected   
};

struct ModelMod{                           // Specifies modifications made to the model
	unsigned int pred_start;                 // The time at which prediction is started in the model
	vector < vector < vector < vector <double> > > > transmean_mult; // Factor which multiplies transition means
	vector < vector < vector <double> > > beta_mult;                 // Factor which multiplies external force of infection
	vector < vector < vector <double> > > efoi_mult;                 // Factor which multiplies external force of infection
	vector < vector <double> > spline_mult;  // Multiplies a spline by a value
	vector < vector <double> > spline_set;   // Sets a spline to a value
};

struct MCMCUpdate{                           // Determines how mcmc updates are permormed
	bool full_mvn;                           // Determines if mvn update on all parameters is performed
	bool mvn;                                // Determines if mvn updates done on seperate groups of parameters
	bool single;                             // Determines if updates performed of individual parameters
	bool dist_R_joint;                       // Determines if joint updates on distribution values and R performed
	bool demo_spec;                          // Determines if demographic specific updates performed
	bool mean_time;                          // Determines if mean time updates performed
	bool neighbour;                          // Determines if neighbour updates performed 
	bool joint;                              // Determines if joint updates performed 
	bool mvn_multiple;                       // Determines if multiple mvn update performed per update
	double multiple_factor;                  // Factor determining the number of updates
};

struct RegionEffect{
	bool on;                                 // Set to true if the region effect is on
	vector <unsigned int> area_param;        // Parameters for a given area
	vector <unsigned int> param_list;        // List of parameters
	unsigned int sigma_param;                // The parameter for sigma
};

struct DirichletParam                      // A dirichlet parameter
{
	unsigned int th;                         // The fraction than parameter has in the sum
	double frac;	
};

struct Dirichlet{                          // Stores information about a Dirichlet distribution
	DirType dirtype;                         // The type of dirichlet
	vector <DirichletParam> param;           // The parameters which are incorporated into distribution
	unsigned int th_min;                     // Used to check that Dirichlet parameters are unique
};

struct DemoDep{                            // Information about a demographic dependency
	string name;                             // The name of the dependency
	unsigned dp;                             // A demographic possiblity consistent with name
};

struct Chain{                              // Stores information about an MCMC chain (used in MC3)
	Chain(string name_, unsigned int num_);

	string name;                             // The name for the chain
	unsigned int num;                        // The number of the chain
	double invT;                             // The inverse temperature
	vector <ParamSample> param_samp;         // Parameter samples for tuning proposals
	vector <double> EF_samp;                 // Parameter samples for calculating model evidence

	unsigned int nproposal;                  // The number of proposals
	
	unsigned int ntr;                        // The number of times a swap is tried
	unsigned int nac;                        // The number of times a swap is accepted
};

struct AreaPlot {                          // Stores information about plotting areas
	string boundfile;                        // A files which specifies boundaries
	string xcol, ycol;                       // Defines columns which specify the location ofr an area
	Project project;                         // Stores the type of projection
};

struct PointInfo {                         // Provides the likelihood and gradient
	double value;                            // The likelihood value
	vector <double> grad;                    // The gradient in the likelihood
};

struct NormalApprox {                      // A normal approximation
	double mean;                             // The normal mean
	double var;                              // The normal variance
};

struct UsedTomlKey{                        // Used to detemine if a TOML key has been used
	string name;                             // The name of the key
	bool used;	                             // Whether it has been used
};

#include "timers.hh"

#endif
