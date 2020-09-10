#ifndef BEEPMBP__STATE_HH
#define BEEPMBP__STATE_HH

#include <vector>

using namespace std;

#include "model.hh"
#include "obsmodel.hh"
#include "consts.hh"

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

	void setlikelihood();
	void checkLevPr();
	void clear();
	void copy(const State &from);
	void check() const;
	Status setparam(const vector <double> &paramv);
	void setbetaphi(unsigned int sett);
	void setLPr();
	FEV getinfev(unsigned int n) const;
	void addindev(unsigned int i);
	void simmodel(unsigned int i, unsigned int c, double t);
	void setQmapUsingdQ(unsigned int sett, const State &state, const vector <double> &dQmap);
	void setQmap(unsigned int check);
	void sortx();
	void initialise_from_particle(const Particle &part);
	
	// These function make up the "standard" proposal
	void standard_parameter_prop(unsigned int samp, unsigned int burnin, vector <float> &paramjumpstand, vector <unsigned int> &ntrstand, 	vector <unsigned int> &nacstand, float &sigmajump);
	void stand_param_betaphi_prop(unsigned int samp, unsigned int burnin, vector <float> &paramjumpstand, vector <unsigned int> &ntrstand, 	vector <unsigned int> &nacstand);
	
	void stand_param_area_prop(unsigned int samp, unsigned int burnin, vector <float> &paramjumpstand, vector <unsigned int> &ntrstand, 	vector <unsigned int> &nacstand);
	
	void stand_param_compparam_prop(unsigned int samp, unsigned int burnin, vector <float> &paramjumpstand, vector <unsigned int> &ntrstand, vector <unsigned int> &nacstand);
	
	void stand_param_fixarea_prop(unsigned int samp, unsigned int burnin, float &sigmajump);
	
private:
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