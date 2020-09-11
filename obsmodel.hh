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

		vector <unsigned int> get_transition_numbers(const vector < vector <EventRef> > &transev, 
																								const vector < vector <Event> > &indev, 
																		unsigned int tra, unsigned int ti, unsigned int tf, unsigned int d, unsigned int v) const;

		double observation_likelihood(const vector < vector <EventRef> > &transev, const vector < vector <Event> > &indev) const;

		Measurements get_measured_quantities(const vector < vector <EventRef> > &transev, const vector < vector <Event> > &indev) const ;

	private:
		double single_observation(unsigned int mean, unsigned int val) const;
		
		const Details &details;
		const Data &data;
		const Model &model;
};

#endif
