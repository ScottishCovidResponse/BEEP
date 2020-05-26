#pragma once

struct SAMPLE{                             // Stores information about a sample from the posterior
	vector <double> paramval;                // A parameter sample
	vector <vector <long> > ncase;	         // A case sample
	vector <double> R0;	                     // Time variation in R0
};

void outputinit(MODEL &model, long nparttot);
SAMPLE outputsamp(double invT, long samp, double Li, DATA &data, MODEL &model, POPTREE &poptree, vector < vector <FEV> > &fev);
void outputresults(DATA &data, MODEL &model, vector <SAMPLE> &opsamp, short siminf, long nparttot);
void outputsimulateddata(DATA &data, MODEL &model, POPTREE &poptree, vector < vector <FEV> > &fev);

vector <long> getnumtrans(DATA &data, MODEL &model, POPTREE &poptree, vector < vector <FEV> > &fev, string from, string to, short ti, short tf);
