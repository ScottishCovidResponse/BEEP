#pragma once

#include <vector>

using namespace std;

#include "consts.hh"

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
	int infchange;                           // Set to one if there us a change in a fixed effect for infectivity
	int ntr, nac;                            // Store the number of proposals tried and accepted	
	double jump;
};

struct COMP{                               // Stores information about a compartment in the model
	string name;                             // Its name
	double infectivity;                      // How infectious that compartment is
	vector <int> trans;                      // The transitions leaving that compartment
};

struct TRANS{                              // Stores information about a compartmental model transition
	int from;                                // Which compartment the individual is coming from
	int to;                                  // Which compartment the individual is going to
	int type;                                // The type of transition (exponential or gamma)
	int param1;                              // First characteristic parameter (e.g. rate)
	int param2;                              // Second characteristic parameter (e.g. standard deviation in the case of gamma) 
};

class MODEL                                // Stores all the information about the model
{
public:
	int modelsel;                            // The model used (e.g. irish, old)

	vector <double> settime;                 // The timings at which beta changes
	vector <double> beta;                    // The value for beta at the various times
	
	vector <double> betai, betap;            // Under MBPs the values of beta for the initial and proposed states
	vector <double> parami, paramp;          // Under MBPs the parameter values for the initial and proposed states
		
	int phiparam, afracparam, aIparam;       // Stores which parameters relate to phi, afrac and aI
	int nspline;                             // The spline points which are parameters in the model
	vector <double> splinet;                 // The times for the spline points
	vector <PARAM> param;                    // Information about parameters in the model
	vector <double> paramval;                // The values of the parameters
	vector <TRANS> trans;                    // Stores model transitions
	vector <COMP> comp;	                     // Stores model compartments

	int ntr, nac;                            // Gets the base acceptance rate
	
	vector <int> fix_sus_param;              // The parameters related to fixed effect for susceptibility
	vector <int> fix_inf_param;              // The parameters related to fixed effect for infectivity
	
	double getrate(string from, string to);
	void simmodel(vector < vector <FEV> > &fev, int &tdnext, int &tdfnext, int i, double t, double period, vector <int> &N);
	void addfev(vector < vector <FEV> > &fev, int &tdnext, int &tdfnext, double t, int tra, int i, int done, double period, vector <int> &N);
	
	void definemodel(int core, double period, int popsize, int mod);
	void betaspline(double period);

private:
	void addcomp(string name, double infectivity);
	void addparam(string name, double val, double min, double max);
	void addtrans(string from, string to, int type, string param1, string param2);
};
