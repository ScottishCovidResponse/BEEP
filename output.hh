#pragma once

struct SAMPLE{                             // Stores information about a sample from the posterior
	vector <double> paramval;                // A parameter sample
	vector <vector <long> > ncase;	         // A case sample
};

void outputinit(MODEL &model);
SAMPLE outputsamp(long samp, double Li, DATA &data, MODEL &model, POPTREE &poptree, vector < vector <FEV> > &fev);
void outputresults(DATA &data, MODEL &model, vector <SAMPLE> &opsamp);
void outputsimulateddata(DATA &data, MODEL &model, POPTREE &poptree, vector < vector <FEV> > &fev);

vector <long> getnumtrans(DATA &data, MODEL &model, POPTREE &poptree, vector < vector <FEV> > &fev, string from, string to, short ti, short tf);
