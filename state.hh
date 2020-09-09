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

	double likelihood();
	void clear();
	void copy(const State &from);
	void check() const;
	Status setparam(const vector <double> &paramv);
	void setbetaphi(unsigned int sett);
	void setLPr();
	FEV getinfev(unsigned int n) const;
	void addindev(unsigned int i);
	void simmodel(unsigned int i, unsigned int c, double t);
	
	vector <int> popw;                                    // The population in w

private:
	
	const Details &details;
	const DATA &data;
	const MODEL &model;
	const Obsmodel &obsmodel;
};

#endif