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

void outputinit(DATA &data, MODEL &model);
void outputLiinit(DATA &data, unsigned int nchaintot);
void outputLi(unsigned int samp, unsigned int nchaintot, double *Litot);

void outputtrace(DATA &data, MODEL &model, unsigned int samp, double Li, double Pri, unsigned int ninf, vector <double> &paramval);
void outputresults(DATA &data, MODEL &model, vector <PARAMSAMP> &psamp, vector <SAMPLE> &opsamp);
void outputplot(string file, DATA &data, MODEL &model,  vector < vector <FEV> > &xi, double tmin, double period);
void outputeventsample(DATA &data, vector < vector <FEV> > &fev);
void outputsimulateddata(DATA &data, MODEL &model, vector < vector <EVREF> > &trev, vector < vector <FEV> > &indev, string dir);
void outputcombinedtrace(vector <string> &paramname, vector < vector < vector <double> > > &vals, string file, string distfile, unsigned int burnin);
#endif
