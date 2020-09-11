#ifndef BEEPMBP__STATE_HH
#define BEEPMBP__STATE_HH

#include <vector>

using namespace std;

#include "model.hh"
#include "obsmodel.hh"
#include "consts.hh"
#include "jump.hh"

class State
{
public:
	State(const Details &details, const DATA &data, const MODEL &model, const Obsmodel &obsmodel);
		
	double L; 																				// The observation likelihood 
	double EF; 																				// The error function (used in abc methods)
	double Pr; 																		    // The prior probability
	
	vector <double> paramval;                         // The parameter values

	vector < vector <FEV> > indev;                    // The individual event sequences
	vector <EVREF> x;                                 // Ordered list of references to infection events 
	vector < vector <EVREF> > trev;                   // Event references
	
	vector< vector <double> > Qmap;                   // The infectivty map 
	
	// The quantities below are all derived from the parameter values
	vector <double> sus;                              // The susceptibility for different demographic categories
	vector <double> areafac;                          // The modification due to area effects
	vector < vector <double> > disc_spline;           // A discretisation of the splines	
	vector <CompTrans> comptrans; 
	
	// The quantities below are temporary
	double beta, phi;                                 // A temporary store for the values of beta and phi
	vector <double> lam;                              // Total force of infecion for an area
	double Lev; 																	   	// The latent process likelihood 
	
	vector <int> popw;                                    // The population in w

	void simulate_compartmental_transitions(unsigned int i, unsigned int c, double t);
	void set_likelihood();
	void clear();
	void copy(const State &from);
	void check() const;
	Status set_param(const vector <double> &paramv);
	void set_betaphi(unsigned int sett);
	void set_LPr();
	double get_infection_time(unsigned int n) const;
	void add_indev(unsigned int i);
	void sim_model(unsigned int i, unsigned int c, double t);
	void set_Qmap_using_dQ(unsigned int sett, const State &state, const vector <double> &dQmap);
	void set_Qmap(unsigned int check);
	void sort_x();
	void initialise_from_particle(const Particle &part);
	void standard_parameter_prop(Jump &jump);
	void check_LevPr();
	
private:
	void stand_param_betaphi_prop(Jump &jump);	   // These function make up the "standard" proposal
	void stand_param_area_prop(Jump &jump);
	void stand_param_compparam_prop(Jump &jump);
	double likelihood_prob(vector <TransInfo> &transinfo, vector <CompTrans> &comptrans) const;
	double likelihood_dt(vector <TransInfo> &transinfo, vector <double> &paramv) const;
	double dlikelihood_dt(vector <TransInfo> &transinfo, vector <double> &paramvi, vector <double> &paramvf) const;
	
	const vector <COMP> &comp;
	const vector <TRANS> &trans;
	const vector <PARAM> &param;
	
	const Details &details;
	const DATA &data;
	const MODEL &model;
	const Obsmodel &obsmodel;
};

#endif