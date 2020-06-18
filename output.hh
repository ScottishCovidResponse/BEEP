#pragma once

struct SAMPLE{                               // Stores information about a sample from the posterior
	vector <double> paramval;                  // A parameter sample
	vector <vector <vector <unsigned int> > > transnum; // Store transition numbers (reflecting the data files)
	vector <double> R0;	                       // Time variation in R0
};

void outputinit(DATA &data, MODEL &model);
void outputLiinit(DATA &data, unsigned int nchaintot);
void outputLi(unsigned int samp, unsigned int nparttot, double *Litot);

SAMPLE outputsamp(double invT, unsigned int samp, double Li, DATA &data, MODEL &model, POPTREE &poptree, vector <double> &paramval, vector < vector <FEV> > &fev);
void outputresults(DATA &data, MODEL &model, vector <SAMPLE> &opsamp);
void outputsimulateddata(DATA &data, MODEL &model, POPTREE &poptree, vector < vector <FEV> > &fev);
void outputplot(string file, DATA &data, MODEL &model,  vector < vector <FEV> > &xi, double tmin, double period);
void outputeventsample(vector < vector <FEV> > &fev, DATA &data, MODEL &model, POPTREE &poptree);
