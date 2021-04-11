/// Stores comonly used data structures

#ifndef BEEPMBP__STRUCT_HH
#define BEEPMBP__STRUCT_HH

#include <vector>

using namespace std;

#include "consts.hh"
#include "utils.hh"

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
class AreaTree;
class Mbp;

struct SplineOutput{                       // Used for outputting splines
	string name;                             // The name of the spline
	string desc;                             // The description of the spline
	vector <double> splineval;               // The values in the spline
};

struct DerivedParam{                       // Parameters derived from the model parameters
	string name;                             // The name of the derived parameter
	string desc;                             // The description of the derived parameter
	string filename;                         // The output file name
	double value;                            // The value
};

struct Sample{                             // Stores information about a state posterior sample
	vector < vector <double> > graph_state;  // Stores the underlying graph state
	vector <SplineOutput> spline_output;     // The splines that need to be output
	vector <DerivedParam> derived_param;     // Derived parameters
};

struct ParamSample{                      	 // Stores information about a parameter sample from the posterior
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

struct Param{                              // Store information about a model parameter
 	string name;                             // The name
	PriorType priortype;                     // The type of the prior
  double val1;                             // Value for the prior (the minimum value assuming a uniform prior) 
	double val2;                             // Value for the prior (the maximum value assuming a uniform prior)
	double value;                            // The value used for the simulation
	ParamType type;                          // Set to one if parameter for transition times and 2 for branching probabilities
	ParamSpec ps;                            // The parameter specification (this is needed to generate Erlang)
	vector <unsigned int> greaterthan;       // List of parameters which this one needs to be greater than
	vector <unsigned int> lessthan;          // List of parameters which this one needs to be greater than
	vector <ParamSplineRef> paramsplineref;  // References any splines the parameter is in
};

struct Compartment{                        // Stores information about a compartment in the model
	string name;                             // The compartments name
	unsigned int infectivity_param;          // How infectious that compartment is
	string inf_trans;                        // A string giving the infection transition on which the infectivity acts
	unsigned int infection_transition;       // The infection transition on which the infectivity acts
	
	vector <unsigned int> trans;             // The transitions leaving that compartment
};

struct CompTransProb                       // If in a compartment this gives probabilities of going down transitions
{
	vector <vector <double> > prob;          // The age-dependent probability of going down transition
	vector <vector <double> > probsum;       // The sum of the age-dependent probability of going down transition
};

struct CompProb{                           // Gives compartmental probabilities
	vector <double> value;                   // The probability an indidividuals goes to compartment (used to calculate R)
};

struct Transition{                         // Stores information about a compartmental model transition
	unsigned int from_origin;                // Which compartment the individual is coming from (if Erlang this is the first one)
	unsigned int from;                       // Which compartment the individual is coming from
	unsigned int to;                         // Which compartment the individual is going to
	unsigned int comptrans_ref;              // The reference to comptrans in from
	double mean_factor;                      // This multiplies the mean (used for Erlang dsitribution)
 	
	unsigned int type;                       // The type of distribution (exponential, gamma or lognormal)
	vector <unsigned int> param_mean;        // The parameter for the mean of the distribution (potentially age dependant)
	vector <unsigned int> param_cv;          // The parameter for the coefficient of variation (if used) (potentially age dependant)
	vector <unsigned int> probparam;         // The parameter for the probability of going down transition (age dependant)
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
	Smooth smooth;                           // Determines the type of smoothing
	vector <SplinePoint> p;                  // The points on the spline
	unsigned int param_factor;               // A parameter whihc can multiply the spline
};

struct SplinePriorGrads{                   // Stores the gradients of spline prior
	double dPr_dth;                          // The derivative of the prior
	double d2Pr_dth2;                        // The double derivative of the prior
};

struct Particle
{
	vector <double> paramval;                // The parameter values for the particle
	vector <vector < vector < vector <int> > > > transnum; // Transition numbers (POP_MODEL)
	double EF;                               // The value of the error function
};

struct Generation                          // Stores inforamtion about a generation when doing ABC methods
{
	vector < vector <double> > param_samp;   // Parameter samples generated
	vector <double> EF_samp;                 // Likelihood samples generated
	vector < vector <double> > EF_datatable; // Stores how the error function is divided into datatables
	vector <unsigned int> partcopy;          // Shows which are copied particles
	vector <double> w;                       // The weight of parameter samples (used in ABC-SMC)
	double EFcut;                            // The error function cut-off used (in CUTOFF mode)
	double EFmin, EFmax;                     // Diagnostic information about EF range
	double invT;                             // The inverse temperature used (in INVT mode)
	long time;                               // The clock time
	
	vector <ParamSample> get_psamp() const;  // Gets a vector of ParamSample from the generation 
};

struct ProbReach                           // The probability of reachin a certain compartment
{
	string name;                             // The name (shown in the output file)
	unsigned int comp;                       // The compartment
	unsigned int inft;                       // The infection transition it relates to 
};

struct TimePlot {                          // A label can be used for specific dates
	string name;
	unsigned int time;
	string time_str;
};

struct Matrix {                            // Loads a matrix
	unsigned int N;													 // The size of the matrix
	vector <vector <double> > ele;           // The elements of the matrix
};

struct SparseMatrix {                      // Loads a matrix
	unsigned int N;													 // The size of the matrix
	vector <double> diag;                    // The diagonal contribution to the matrix
	vector < vector <unsigned int> > to;     // The value of the element to
	vector < vector <double> > val;          // The the non-diagonal elements of the matrix
};

struct MatrixModification {                // Used to modify the mixing matrix
	MatModType type;	                       // The type of modification
	string name;                             // The name of the modification
	string desc;                             // A description of the modification
	vector <ParamSpec> ps;                   // The parameter specifications
	vector <unsigned int> bp;                // Finds the break points for the spline
	vector <unsigned int> ages;              // The group of ages on which the modification acts
	unsigned int spline_ref;                 // The spline used to generate the modification 
};

struct GenerateQ {                         // Stores information about age and geographical mixing
	vector <string> N_name;                  // Stores the name of matrices to be loaded

	vector <MatrixModification> matmod;      // Store potential modifications to the basic mixing matrix 

	unsigned int nspline;                    // The number of splines used for a particular type
	
	string M_name;                           // The geographic mixing matrix name
	
	vector <Matrix> N;                       // Stores the matrices

	SparseMatrix M;                          // The actual loaded geographic mixing matrix 
	
	vector <double> onlydiag;                // The elements of a completely diagnoal matrix for geomixing
	
	vector <double> factor;                  // A factor used for geospline
	
	vector <unsigned int> area_filter_ref;   // When filtering areas provides a reference	
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
	string name;                             // The name of the data filter
	vector <unsigned int> area;              // The areas that the data points refer to 
	vector <unsigned int> dp_sel;            // The demographic possibilities that the data points refer to
};

struct DataTable {                         // Stores tables of data provided by the user
	OpType optype;                           // Whether it is a true data table of just an output
	DataType type;                           // The type of the data 
	string observation;                      // The quantity being observed 
	string geo;                              // The geographic scale the data related to
	string file;                             // The name of the file to be loaded
	unsigned int start;                      // The start time for the data
	unsigned int end;                        // The end time for the data
	
	unsigned int units;                      // The units used (e.g. 1=days, 7=weeks) 
	int shift;
	string democats;
	string fraction_spline;                  // The spline giving the fraction observations made
	
	vector <string> demolist;                // For marginal distributions stores categories
	
	vector <unsigned int> translist;         // The transitions relating to the data
	vector <unsigned int> complist;          // The compartments relating to the data
	
	vector <unsigned int> graph_ref;         // A list of all the graphs plotted
	double weight;                           // The weight attached to the observation model
	
	ObsModelFunc obsmodel;                   // Determines what observaion model is used
	double shape;                            // The shape parameter
	double invT;                             // Scaling done to observation model
};

struct Observation {                       // A single observation made on the system
	unsigned int fraction_spline_ref;        // References the spline number
	unsigned int datatable;                  // The data table the observation belongs to
	unsigned int graph;                      // The primary graph the observation belongs to
	unsigned int value;                      // The value of the observation
	unsigned int sett_i;                     // The start time
	unsigned int sett_f;                     // The end time
	vector <unsigned int> area;              // The areas involved
	vector <unsigned int> dp_sel;            // The deomgraphic possibilities selected  	
	
	double logdif_offset;                    // The offset used when calculating the error
	double w;                                // The weight used when calculating the error
	double invT;                             // Parameter used to scale observation model
	
	ObsModelFunc obsmodel;                   // The observation model used
	double shape;                            // Shape paremeter (used for negative binomial).
};

struct GraphPoint {                        // Details of a point on a graph
	double xi;                               // The x position at the start of the observation
	double xf;                               // The x position at the end of the observation
	vector <unsigned int> obs;               // The observations which contribute to the point
};

struct Graph {                             // Details of a graph
	string name;                             // The name of the graph
	string filename;                         // The filename used for the graph     
	string desc;                             // The description of the graph
	string colname;                          // The column the graph relates to
	GraphType type;
	vector <GraphPoint> point;               // The points on the graph
	
	unsigned int fraction_spline_ref;        // References the spline number
	unsigned int datatable;                  // The data table the observation belongs to
	vector <unsigned int> area;              // The areas involved
	vector <unsigned int> dp_sel;            // The deomgraphic possibilities selected  	
};

struct DemographicCategory {               // Stores demographic categories
	string name;                             // The name of the category
	bool sus_vari;                           // Set to true if there is variation in susceptibility
	vector <string> value;                   // The postential values it can take
	vector <ParamSpec> ps;                   // Information about the parameter specification
};

struct Covariate {                         // Stores the  covariate for the area
	string name;                             // The name of the covariate (i.e. the column in the area data file)
	ParamSpec ps;                            // Gives a parameter specification
	string func;                             // The functional transformation
};

struct Area {                              // Provides information about an area
	string code;                             // The code for the area
	
	unsigned int region_effect;              // If a region has a region_effect
	
 	double x, y;                             // The geographic position
	
  vector <unsigned int> pop;               // The population in different demographic categories          
	vector <double> covar;                   // The covariates for that area

	vector <vector <unsigned int> > ind;     // The individuals in different demographic categories
};

struct Individual {                        // Provides information about an individual
	unsigned int area;                       // The area
	unsigned int dp;                         // The demographic category possibility
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

struct CounterFact{                        // Used to implement a counter factual analysis (NOTE implemented in mpi.copy_data)
	unsigned int start;                      // The start time of a change
	unsigned int end;                        // The end time for a change
	CounterFactType type;                    // The type of change
	string trans_str;                        // The transition
	string geo;                              // The geography
	string democats;                         // The demographic categories
	double factor;                           // The factor change
	vector <unsigned int> dp_sel;            // The demographic possibilities which are selected
	vector <unsigned int> area;              // The areas which are selected   
};

struct CounterMod{                         // Specifies counterfactual modifications made to the model
	unsigned int sett_start;                 // The start time for the first modification
	vector < vector < vector < vector <double> > > > transmean_mult; // Factor which multiplies transition means
	vector < vector < vector <double> > > efoi_mult;                 // Factor which multiplies external force of infection
};

#include "timers.hh"

#endif
