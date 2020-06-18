#pragma once

#include <vector>

using namespace std;

#include "consts.hh"
#include "data.hh"

struct FEV {                               // Stores information about a compartmental transition
  int trans;                               // References the transition type
	int ind;                                 // The individual on which the transition happens
	double t;                                // The time of the transition
	int done;                                // Set to 1 if that transition is in the past 
};

struct PARAM{                              // Store information about a model parameter
 	string name;                             // Its name
 	double valinit;                          // The simulation value or starting value for inference
	double sim;                              // The simulation value
	double min;                              // The minimum value (assuming a uniform prior) 
	double max;                              // The maximum value (assuming a uniform prior)
	int betachange;                          // Set to one if there is a change in the beta spline
	int suschange;                           // Set to one if there us a change in a fixed effect for susceptibility
	int ntr, nac;                            // Store the number of proposals tried and accepted	
	double jump;
};

struct COMP{                               // Stores information about a compartment in the model
	string name;                             // Its name
	double infectivity;                      // How infectious that compartment is

	int type;                                // The type of distribution for waiting in compartment (exponential or gamma)
	int param1;                              // First characteristic parameter (e.g. rate)
	int param2;                              // Second characteristic parameter (e.g. standard deviation in the case of gamma)

	vector <int> trans;                      // The transitions leaving that compartment
	int transtimep;                          // The time period transitions leaving that compartment

	vector <double> prob, probsum;           // The probability of going down transition
	vector <double> probi;                   // The probability of going down transition for initial state (MBP)
	vector <double> probp;                   // The probability of going down transition for proposed state (MBP)
};

struct TRANS{                              // Stores information about a compartmental model transition
	int from;                                // Which compartment the individual is coming from
	int to;                                  // Which compartment the individual is going to
	int probparam;                           // The parameter for the probability of going down transition
	vector <int> DQ;                         // The change in the Q tensor for going down the transition
};

class MODEL                                // Stores all the information about the model
{
public:
	int modelsel;                            // The model used (e.g. irish, old)

	vector <double> settime;                 // The timings at which beta changes
	vector <double> beta;                    // The value for beta at the various times
	vector <double> sus;                     // The susceptibility for different demographic categories
	
	vector <double> betai, betap;            // Under MBPs the values of beta for the initial and proposed states
	vector <double> parami, paramp;          // Under MBPs the parameter values for the initial and proposed states
	vector <double> susi, susp;              // Under MBPs the susceptibility for the initial and proposed states
			
	int phiparam;               					   // Stores which parameters relate to phi and probA
	int nspline;                             // The spline points which are parameters in the model
	vector <double> splinet;                 // The times for the spline points
	vector <PARAM> param;                    // Information about parameters in the model
	vector <double> paramval;                // The values of the parameters
	vector <TRANS> trans;                    // Stores model transitions
	vector <COMP> comp;	                     // Stores model compartments

	int ntr, nac;                            // Gets the base acceptance rate
	
	vector <int> fix_sus_param;              // The parameters related to fixed effect for susceptibility
	
	int ntimeperiod;                         // The number of different time periods (2: before and after lockdown)
	vector <double> timeperiod;              // The timings of changes to Q;
	
	vector <vector <int> > nDQ;               // Stores the changes in the mixing matrix between areas and ages 
	vector <vector< vector <int> > > DQto;
	vector <vector <vector< vector <double> > > > DQval;
	
	double getparam(string name);
	void simmodel(vector <FEV> &evlist, int i, int c, double t);
	void mbpmodel(vector <FEV> &evlisti, vector <FEV> &evlistp);
	void definemodel(DATA &data, int core, double period, int popsize, int mod);
	void addQ(DATA &data);
	void betaspline(double period);
	void priorsamp();
	int settransprob();
	void setsus(DATA &data);           
	void checktransdata(DATA &data);
	
private:
	void addcomp(string name, double infectivity, int type, string param1, string param2);
	void addparam(string name, double val, double min, double max);
	void addtrans(string from, string to, string probparam);
};
