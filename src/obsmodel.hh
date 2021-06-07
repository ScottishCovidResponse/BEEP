#ifndef BEEPMBP__OBSModel_HH
#define BEEPMBP__OBSModel_HH

using namespace std;

#include "struct.hh"

class ObservationModel
{
	public:
		ObservationModel(const Details &details, Data &data, const Model &model);
		
		double calculate(const State *state) const;
		double calculate_section(const State *state, unsigned int sec) const;
		vector <double> get_EF_datatable(const State *state) const;
		vector <double> get_obs_value(const State *state) const;
		vector < vector <double> > get_graph_state(const State *state) const;
				
		unsigned int nsection;                                     // The sections the obsevations are split into 
		vector <unsigned int> section_ti;             
		vector <unsigned int> section_tf;
		
	private:
		void initialise_obs_change();
		void split_observations();
		void get_obs_value_section(const State *state,  vector <double> &obs_value, const unsigned int ti, const unsigned int tf) const;
		double obs_prob(double value, const Observation& ob) const;
		
		vector < vector <unsigned int> > section_obs;              // Stores the observations in each section
		
		vector < vector <vector < vector < vector <unsigned int> > > > > obs_trans; // Observation in  transition [sett][c][tr][dp]
		vector < vector <vector < vector < vector <unsigned int> > > > > obs_pop;   // Observation in population [sett][c][co][dp]
				
		const Details &details;
		const Data &data;
		const Model &model;
};

#endif
