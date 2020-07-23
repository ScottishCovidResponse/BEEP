#ifndef BEEPMBP__OUTPUT_HH
#define BEEPMBP__OUTPUT_HH

#include "poptree.hh"

struct SAMPLE{                                        // Stores information about a sample from the posterior
	vector <double> paramval;                           // A parameter sample
	MEAS meas;                                          // Stores measurements corresponding to the data file
	vector <double> R0;	                                // Time variation in R0
};

void outputinit(DATA &data, MODEL &model);
void outputLiinit(DATA &data, unsigned int nchaintot);
void outputLi(unsigned int samp, unsigned int nchaintot, double *Litot);

SAMPLE outputsamp(double invT, unsigned int samp, double Li, double Pri, DATA &data, MODEL &model, POPTREE &poptree, vector <double> &paramval, unsigned int ninf, vector < vector <EVREF> > &trev, vector < vector <FEV> > &indev);
void outputresults(DATA &data, MODEL &model, vector <SAMPLE> &opsamp);
void outputplot(string file, DATA &data, MODEL &model,  vector < vector <FEV> > &xi, double tmin, double period);
void outputeventsample(vector < vector <FEV> > &fev, DATA &data, MODEL &model, POPTREE &poptree);
void outputsimulateddata(DATA &data, MODEL &model, POPTREE &poptree, vector < vector <EVREF> > &trev, vector < vector <FEV> > &indev, string dir);
void outputcombinedtrace(vector <string> &paramname, vector < vector < vector <double> > > &vals, string file);
#endif
