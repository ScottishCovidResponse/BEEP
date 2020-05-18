#pragma once

struct PARAM{                              // Store information about a model parameter
 	string name;                             // Its name
 	double val;                              // The simulation value or starting value for inference
	double sim;                              // The simulation value
	double min;                              // The minimum value (assuming a uniform prior) 
	double max;                              // The maximum value (assuming a uniform prior)
	double jump;                             // The size of proposed changes in PMCMC
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
	void definemodel();
	void betaspline();

	double settime[nsettime];
	double beta[nsettime];
	short nspline;                             // The spline points which are parameters in the model
	vector <double> splinet;
	vector <PARAM> param;
	vector <TRANS> trans;
	vector <COMP> comp;	

private:
	void addcomp(string name, double infectivity);
	void addparam(string name, double val, double min, double max);
	void addtrans(string from, string to, short type, string param1, string param2);
};
