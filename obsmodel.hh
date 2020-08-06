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
	Obsmodel(const Details &details, const DATA &data, const MODEL &model);

	vector <unsigned int> getnumtrans(const vector < vector <EVREF> > &trev, const vector < vector <FEV> > &indev, unsigned int tra, unsigned int ti, unsigned int tf, unsigned int d, unsigned int v) const;

	double Lobs(const vector < vector <EVREF> > &trev, const vector < vector <FEV> > &indev) const;

	MEAS getmeas(const vector < vector <EVREF> > &trev, const vector < vector <FEV> > &indev) const ;

private:
	double singobs(unsigned int mean, unsigned int val) const;
	
	const Details &details;
	const DATA &data;
	const MODEL &model;
};

#endif
