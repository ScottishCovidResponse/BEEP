#ifndef BEEPMBP__OBSModel_HH
#define BEEPMBP__OBSModel_HH

#include <iostream>
#include <vector>
#include <fstream>
#include <sys/stat.h>
#include <sstream>
#include <algorithm>
#include <cmath>

using namespace std;

#include "data.hh"
#include "areatree.hh"

class ObservationModel
{
public:
	ObservationModel(const Details &details, const Data &data, const Model &model);

	vector <unsigned int> getnumtrans(const vector < vector <EventRef> > &trev, const vector < vector <Event> > &indev, unsigned int tra, unsigned int ti, unsigned int tf, unsigned int d, unsigned int v) const;

	double Lobs(const vector < vector <EventRef> > &trev, const vector < vector <Event> > &indev) const;

	Measurements getmeas(const vector < vector <EventRef> > &trev, const vector < vector <Event> > &indev) const ;

private:
	double singobs(unsigned int mean, unsigned int val) const;
	
	const Details &details;
	const Data &data;
	const Model &model;
};

#endif
