#ifndef BEEPMBP__MODEL_HH
#define BEEPMBP__MODEL_HH

#include <vector>

using namespace std;

#include "consts.hh"
#include "data.hh"


struct FEV {                               // Stores information about a compartmental transition
  unsigned int trans;                      // References the transition type
	unsigned int ind;                        // The individual on which the transition happens
	double t;                                // The time of the transition
	unsigned int timep;                      // The time period in which the transition occurs 
};

struct EVREF {                             // Used to reference an event
	unsigned int ind;                        // The individual
	unsigned int e;	                         // The event number
};

struct DQINFO {                            // Stores information about a change in the Q matrix
	vector <unsigned int> q;
	vector <double> fac;
};

struct PRIORCOMP {											   // Stores information about priors on compartmental probabilities
	unsigned int comp;											 // The compartment
	double value;														 // The value of the comparmental probability
	double sd;															 // The standard deviation in the value
};

struct PARAM{                              // Store information about a model parameter
 	string name;                             // Its name
  double min;                              // The minimum value (assuming a uniform prior) 
	double max;                              // The maximum value (assuming a uniform prior)
	unsigned int used;                       // Determins if the parameter is used in the model.
	unsigned int ntr, nac;                   // Stores the number of proposals tried and accepted	
	double jump;                             // Stores the size of jumps in parameter space
	ParamType type;                          // Set to one if parameter used for transition times and 2 for branching probabilities
};

struct COMP{                               // Stores information about a compartment in the model
	string name;                             // Its name
	double infectivity;                      // How infectious that compartment is

	vector <unsigned int> trans;             // The transitions leaving that compartment
	unsigned int transtimep;                 // The transition used to represent a change in timep in that compartment

	vector <vector <double> > prob, probsum; // The age-dependent probability of going down transition
	vector <vector <double> > probi;         // The age-dependent probability of going down transition for initial state (MBP)
	vector <vector <double> > probp;         // The age-dependent probability of going down transition for proposed state (MBP)
	
	vector <double> probin;                  // The prob ind goes to compartment (used to calculate R0)
	vector <double> infint;                  // The integrated infectivity (used to calculate R0)
};

struct TRANS{                              // Stores information about a compartmental model transition
	unsigned int from;                       // Which compartment the individual is coming from
	unsigned int to;                         // Which compartment the individual is going to
	
	unsigned int type;                       // The type of distribution (exponential or gamma)
	int param_mean;                          // The parameter for the mean of the distribution
	int param_cv;                            // The parameter for the coefficient of variation (if used)
	
	unsigned int istimep;                    // Set to one if the transition is in time period
	vector <unsigned int> num;               // The number of times down transition 
	vector <unsigned int> probparam;         // The parameter for the probability of going down transition (age dependant)
	vector <unsigned int> DQ;                // The change in the Q tensor for going down the transition (age dependant)
	
	                                          // The following are used for making changes to the parameters
	unsigned int numvisittot;                // The number of times the compartment is visited
	double dtsum;                            // Sums up the total time spent
	vector <double> dtlist;                  // Keeps a list of waiting time
};

struct SPLINEP{                            // Stores information about a spline point
	double t;                                // The time of the point
	unsigned int param;                      // The parameter which defines the value
};

struct Particle
{
	vector <double> paramval;                // The parameter values for the particle
	vector <FEV> ev;                         // The event sequence for the particle
	double EF;                               // The value of the error function
};

struct Generation
{
	vector < vector <double> > param_samp;   // Parameter samples generated
	vector <double> EF_samp;                 // Likelihood samples generated
	vector <unsigned int> partcopy;          // Shows which are copies particles
	vector <double> w;                       // The weight of parameter samples (used in ABC-SMC)
	
	double EFcut;                            // The cut off likelihood used 
	
	long time;                               // The clock time
};

class MODEL                                // Stores all the information about the model
{
public:
	MODEL(Inputs &inputs, const Details &details, DATA &data);

	unsigned int ndemocat;                   // The number of demographic categories
	vector <DEMOCAT> democat;                // Stores the demographic categories
	
	vector <double> beta;                    // The value for beta at the various times
	vector <SPLINEP> betaspline;             // The spline used to define beta
	
	vector <double> phi;                     // The value for phi at the various times
	vector <SPLINEP> phispline;              // The spline used to define phi
	
	vector <double> sus;                     // The susceptibility for different demographic categories

	vector <double> areafac;                 // The modification due to area effects
	
	vector <double> betai, betap;            // Under MBPs the values of beta for the initial and proposed states
	vector <double> phii, phip;              // Under MBPs the values of phi for the initial and proposed states
	vector <double> parami, paramp;          // Under MBPs the parameter values for the initial and proposed states
	vector <double> susi, susp;              // Under MBPs the susceptibility for the initial and proposed states
	vector <double> areafaci, areafacp;      // Under MBPs the area factor for the initial and proposed states
				
	vector <PARAM> param;                    // Information about parameters in the model
	vector <PRIORCOMP> priorcomps;           // Priors on compartmental probabilities
	vector <TRANS> trans;                    // Stores model transitions
	vector <COMP> comp;	                     // Stores model compartments

	unsigned int infmax;                     // The maximum number of infected individuals
	
	unsigned int ntr, nac;                   // Gets the base acceptance rate
	
	vector< vector <unsigned int> > sus_param;// The parameters related to fixed effect for susceptibility
	vector <unsigned int> covar_param;       // The parameters related to covariates for areas
	
	unsigned int regioneffect;               // Set to 1 if a regional random effect is put in the force of infection
	unsigned int sigma_param;                // The standard deviation of the regional effect
	vector <unsigned int> regioneff_param;   // The parameters related to regional effects

	vector <TIMEP> timeperiod;               // The timings of changes to Q;
	unsigned int ntimeperiod;
	
	vector <DQINFO> DQ;                      // Keeps track of the change in the Q matrix 
	
	double getinfectivity(const string& name);
	void simmodel(const vector<double> &paramv, vector <FEV> &evlist, unsigned int i, unsigned int c, double t);
	void mbpmodel(vector <FEV> &evlisti, vector <FEV> &evlistp);
	void print_to_terminal() const;
	vector <double> priorsamp();
	unsigned int setup(const vector <double> &paramval);
	void copyi(const vector<double> &paramv);
	void copyp(const vector<double> &paramv);
	vector <double> R0calc(const vector<double> &paramv);
	unsigned int dombpevents();
	void oe(const string& name, vector <FEV> &ev);
	void calcprobin();
	double prior(const vector<double> &paramv);
	void compparam_prop(unsigned int samp, unsigned int burnin, vector <EVREF> &x, vector <vector <FEV> > &indev, vector <double> &paramv,
												   vector <float> &paramjumpxi, vector <unsigned int> &ntrxi,  vector <unsigned int> &nacxi, double &Pri);
													 
private:
	void addQ();
	void checkdata();
	void addcomp(const string& name, double infectivity);
	void addparam(const string& name, double min, double max);
	void addtrans(const string& from, const string& to, const string& prpar,
								unsigned int type, const string& param1, const string& param2);
	void setsus(const vector<double> &paramv); 
	void setarea(const vector<double> &paramv);
	void timevariation(const vector<double> &paramv);
	unsigned int findparam(const string& name);
	unsigned int settransprob(const vector<double> &paramv);
	double likelihood_prob();
	double likelihood_dt(vector <double> &paramv);
	double dlikelihood_dt(vector <double> &paramvi, vector <double> &paramvf);

	const Details &details;
	DATA &data;
};
#endif
