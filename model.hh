#pragma once

#include <vector>

using namespace std;

#include "consts.hh"

struct PARAM{                              // Store information about a model parameter
 	string name;                             // Its name
 	double val;                              // The simulation value or starting value for inference
	double sim;                              // The simulation value
	double min;                              // The minimum value (assuming a uniform prior) 
	double max;                              // The maximum value (assuming a uniform prior)
	double jump;                             // The size of proposed changes in PMCMC
	short betachange;                        // Set to one if there is a change in the beta spline
	short suschange;                         // Set to one if there us a change in a fixed effect for susceptibility
	short infchange;                         // Set to one if there us a change in a fixed effect for infectivity
	long ntr, nac;                           // Store the number of proposals tried and accepted	
};

struct COMP{                               // Stores information about a compartment in the model
	string name;                             // Its name
	double infectivity;                      // How infectious that compartment is
	vector <long> trans;                     // The transitions leaving that compartment
};

struct TRANS{                              // Stores information about a compartmental model transition
	short from;                              // Which compartment the individual is coming from
	short to;                                // Which compartment the individual is going to
	short type;                              // The type of transition (exponential or gamma)
	short param1;                            // First characteristic parameter (e.g. rate)
	short param2;                            // Second characteristic parameter (e.g. standard deviation in the case of gamma) 
};

class MODEL
{
public:
	void definemodel(short core, double tmax, long popsize);
	void betaspline(double tmax);

	double settime[nsettime];
	double beta[nsettime];
	
	/* MBP*/
	double betai[nsettime], betap[nsettime];
	vector <double> parami, paramp;
	/* MBP END*/
		
	short phiparam, afracparam, aIparam;
	short nspline;                             // The spline points which are parameters in the model
	vector <double> splinet;
	vector <PARAM> param;
	vector <TRANS> trans;
	vector <COMP> comp;	

	long ntr, nac;                             // Gets the base acceptance rate
	
	vector <long> fix_sus_param;               // The parameters related to fixed effect for susceptibility
	vector <long> fix_inf_param;               // The parameters related to fixed effect for infectivity
	
	double getrate(string from, string to);
	
private:
	void addcomp(string name, double infectivity);
	void addparam(string name, double val, double min, double max);
	void addtrans(string from, string to, short type, string param1, string param2);
};
