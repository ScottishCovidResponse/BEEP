#ifndef BEEPMBP__OBSMODEL_HH
#define BEEPMBP__OBSMODEL_HH

#include <iostream>
#include <vector>
#include <fstream>
#include <sys/stat.h>
#include <sstream>
#include <algorithm>
#include <cmath>

using namespace std;

#include "data.hh"
#include "poptree.hh"

class Obsmodel
{
public:
	Obsmodel(Details &details, DATA &data, MODEL &model);

	vector <unsigned int> getnumtrans(const vector < vector <EVREF> > &trev, const vector < vector <FEV> > &indev, unsigned int tra, unsigned int ti, unsigned int tf, unsigned int d, unsigned int v);

	double Lobs(const vector < vector <EVREF> > &trev, const vector < vector <FEV> > &indev);

	MEAS getmeas(const vector < vector <EVREF> > &trev, const vector < vector <FEV> > &indev);

private:
	double singobs(unsigned int mean, unsigned int val);
	
	Details &details;
	DATA &data;
	MODEL &model;
};

#endif
