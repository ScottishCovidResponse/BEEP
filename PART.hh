#pragma once

#include <vector>

using namespace std;

#include "types.hh"
#include "model.hh"
#include "poptree.hh"

class PART                                 // Stores all the things related to a particle 
{
	public:
	PART(MODEL &model, POPTREE &poptree);

	long pst;                                // The number of the particle
	
 	double Li;                               // The observation likelihood
 	
	vector <double> ffine;                   // Stores the force of infection on nodes on the fine scale
	vector <vector <long> > indinf;          // Lists all infected individuals  
	vector <vector <long> > pop;             // The total popualtion for nodes on different levels 
	vector <vector <double> > Rtot;          // The total infection rate for nodes on different levels
	vector <vector <double> > addlater;      // A change to the rates Rtot which may be performed later when sampling is performed

	vector < vector <FEV> > fev;             // Stores all compartmental transitions
	
	vector <long> N;                         // The number of individuals in different compartments

	short sett;                              // Index used to track time changes in beta

	long tdnext, tdfnext;                    // Stores when the next future compartmental transition will occur

	MODEL  &model;
	vector <TRANS> &trans;
	vector <COMP> &comp;

	POPTREE &poptree;
	vector <LEVEL> &lev;

	public: 
		void gillespie(double ti, double tf, short siminf);
		void partinit(long p);
		void dofe();
		long nextinfection();
		void addinfc(long c, double t);
		void addfev(double t, long tr, long i);
		vector <long> getnumtrans(string from, string to, short ti, short tf);
		void Lobs(short ti, short tf, long ncase[nregion][tmax/7+1]);
		void copy(const PART &other);
		void simmodel(long i, short enter, double t);
};
