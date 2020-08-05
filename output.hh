#ifndef BEEPMBP__OUTPUT_HH
#define BEEPMBP__OUTPUT_HH

#include "poptree.hh"

struct SAMPLE{                                        // Stores information about a sample from the posterior
	MEAS meas;                                          // Stores measurements corresponding to the data file
	vector <double> R0;	                                // Time variation in R0
	vector <double> phi;	                              // Time variation in phi
};

struct PARAMSAMP{                                     // Stores information about a sample from the posterior
	vector <double> paramval;                           // A parameter sample
};

struct STAT{                                           // Stores statistical information
	string mean;                                         // The mean
	string CImin, CImax;                                 // The minimum and maximum of the 90% credible interval
	string ESS;                                          // The estimated sample size
};

struct DIST{                                          // Stores a probability distribution
	vector <string> value;
	vector <string> prob;
};

class Obsmodel;

class Output
{
public:
	Output(Details &details, DATA &data, MODEL &model, Obsmodel &obsmodel);
	
	void init();
	void Liinit(unsigned int nchaintot);
	void Li(unsigned int samp, unsigned int nchaintot, double *Litot);
	void plot(string file, vector < vector <FEV> > &xi, double tmin, double period);

	void traceplot(unsigned int samp, double Li, double Pri, unsigned int ninf, vector <double> &paramval);
	void results(vector <PARAMSAMP> &psamp, vector <SAMPLE> &opsamp);
	void eventsample(vector < vector <FEV> > &fev);
	void simulateddata(vector < vector <EVREF> > &trev, vector < vector <FEV> > &indev, string dir);
	void combinedtrace(vector <string> &paramname, vector < vector < vector <double> > > &vals, string file, string distfile, unsigned int burnin);

private:
	STAT getstat(vector <double> &vec); 
	DIST getdist(vector <double> &vec);
	void posterior_plot(vector <SAMPLE> &opsamp, unsigned int d, unsigned int r, unsigned int type);

	ofstream trace, traceLi;
	
	const Details &details;
	const DATA &data;
	MODEL &model;
	Obsmodel &obsmodel;
};
#endif
