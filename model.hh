#ifndef BEEPMBP__Model_HH
#define BEEPMBP__Model_HH

#include <vector>

using namespace std;

#include "consts.hh"
#include "data.hh"

struct Event {                             // Stores information about a compartmental transition
  unsigned int trans;                      // References the transition type
	unsigned int ind;                        // The individual on which the transition happens
	double t;                                // The time of the transition
	unsigned int timep;                      // The time period in which the transition occurs 
};

struct EventRef {                          // Used to reference an event
	unsigned int ind;                        // The individual
	unsigned int e;	                         // The event number
};

struct DQinfo {                            // Stores information about a change in the Q matrix
	vector <unsigned int> q;                 // Gives which Q tensor is used before and after a transition
	vector <double> fac;                     // Gives the factor which multiplies tensor (accounts for infectivity)
};

struct PriorComp {											   // Stores information about priors on compartmental probabilities
	unsigned int comp;											 // The compartment
	double value;														 // The value of the comparmental probability
	double sd;															 // The standard deviation in the value
};

struct Param{                              // Store information about a model parameter
 	string name;                             // Its name
  double min;                              // The minimum value (assuming a uniform prior) 
	double max;                              // The maximum value (assuming a uniform prior)
	bool used;                               // Determins if the parameter is used in the model.
	unsigned int ntr, nac;                   // Stores the number of proposals tried and accepted	
	double jump;                             // Stores the size of jumps in parameter space
	ParamType type;                          // Set to one if parameter for transition times and 2 for branching probabilities
};

struct Compartment{                        // Stores information about a compartment in the model
	string name;                             // The compartments name
	double infectivity;                      // How infectious that compartment is

	vector <unsigned int> trans;             // The transitions leaving that compartment
	unsigned int transtimep;                 // The transition used to represent a change in timep in that compartment
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
	unsigned int from;                       // Which compartment the individual is coming from
	unsigned int to;                         // Which compartment the individual is going to
	
	unsigned int type;                       // The type of distribution (exponential, gamma or lognormal)
	int param_mean;                          // The parameter for the mean of the distribution
	int param_cv;                            // The parameter for the coefficient of variation (if used)
	
	bool istimep;                            // Set to true if the transition changes the time period
	vector <unsigned int> probparam;         // The parameter for the probability of going down transition (age dependant)
	vector <unsigned int> DQ;                // The change in the Q tensor for going down the transition (age dependant)
};

struct SplinePoint{                        // Stores information about a spline point
	double t;                                // The time of the point
	unsigned int param;                      // The parameter which defines the value
	double multfac;                          // A multiplicative factor
};

struct Particle
{
	vector <double> paramval;                // The parameter values for the particle
	vector <Event> ev;                       // The event sequence for the particle
	double EF;                               // The value of the error function
};

struct Generation                          // Stores inforamtion about a generation when doing ABC methods
{
	vector < vector <double> > param_samp;   // Parameter samples generated
	vector <double> EF_samp;                 // Likelihood samples generated
	vector <unsigned int> partcopy;          // Shows which are copies particles
	vector <double> w;                       // The weight of parameter samples (used in ABC-SMC)
	double EFcut;                            // The error function cut-off used 
	long time;                               // The clock time
};

class Model                                // Stores information about the model
{
	public:
		Model(const Inputs &inputs, const Details &details, Data &data);

		unsigned int ndemocat;                   // The number of demographic categories
		vector <DemographicCategory> democat;    // Stores the demographic categories
		
		vector <vector <SplinePoint> > spline;   // Stores all the splines used in the model (for beta and phi)  
	 
		unsigned int betaspline_ref;             // Denotes which spline refers to variation in beta
		unsigned int phispline_ref;              // Denotes which spline refers to variation in phi

		vector <Param> param;                    // Information about parameters in the model
		vector <PriorComp> priorcomps;           // Priors on compartmental probabilities
		vector <Transition> trans;               // Stores model transitions
		vector <Compartment> comp;	             // Stores model compartments

		unsigned int maximum_infected;           // The maximum number of infected individuals
		
		vector< vector <unsigned int> > suscetibility_param;// The parameters related to fixed effect for susceptibility
		vector <unsigned int> covariate_param;   // The parameters related to covariates for areas
		
		unsigned int region_effect;              // Set to 0 of no effect, 1 if random effect, 2 if fixed effect
		unsigned int sigma_param;                // The standard deviation of the regional effect (if random effect)
		vector <unsigned int> regioneffect_param;// The parameters related to regional effects

		vector <TimePeriod> time_period;         // The timings of changes to Q (e.g. before and after lockdown)
		unsigned int ntime_period;               // The number of time periods
		
		vector <DQinfo> DQ;                      // Keeps track of the change in the Q matrix 
		
		double get_infectivity(const string& name) const;
		void print_to_terminal() const;
		vector <double> sample_from_prior() const;
		vector <double> calculate_R_vs_time(const vector <double> &paramv) const;
		bool do_mbp_events(const vector <double> &parami, const vector <double> &paramp) const;
		void print_events(const string& name, const vector <Event> &ev) const;
		double prior(const vector<double> &paramv) const;
		vector <double> create_disc_spline(unsigned int ref, const vector<double> &paramv) const;
		vector <double> create_susceptibility(const vector<double> &paramv) const;  
		vector <double> create_areafactor(const vector<double> &paramv) const;
		unsigned int create_comptransprob(vector <CompTransProb> &comptransprob, const vector<double> &paramv) const;
		vector <CompProb> create_compprob(const vector <CompTransProb> &comptransprob) const;
		bool inbounds(const vector <double> &paramv) const;
		bool inbounds(double val, unsigned int th) const;
		
	private:
		void add_Q();
		void add_compartment(const string& name, double infectivity);
		void add_parameter(const string& name, double min, double max);
		void add_transition(const string& from, const string& to, const string& prpar,
									unsigned int type, const string& param1, const string& param2);
		unsigned int find_parameter(const string& name);
		void check_data_files() const;
		
		const Details &details;
		Data &data;
};
#endif
