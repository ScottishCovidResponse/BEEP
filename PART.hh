#pragma once

#include <vector>

using namespace std;

#include "model.hh"
#include "poptree.hh"
#include "data.hh"

struct FEV {                               // Stores information about a compartmental transition
  long trans;                              // References the transition type
	long ind;                                // The individual on which the transition happens
	double t;                                // The time of the transition
	short done;                              // Set to 1 if that transition is in the past 
};

class PART                                 // Stores all the things related to a particle 
{
	public:
	PART(DATA &data, MODEL &model, POPTREE &poptree);

	long pst;                                // The number of the particle
	
 	double Li;                               // The observation likelihood
 	
	vector <double> ffine;                   // Stores the force of infection on nodes on the fine scale
	vector <vector <long> > indinf;          // Lists all infected individuals 	
		
	vector <vector <double> > sussum;        // The total susceptibility for nodes on different levels 
	//vector <vector <double> > sussumheff;    // The total susceptibility for household effect for nodes on different levels 
	vector <vector <double> > Rtot;          // The total infection rate for nodes on different levels
	vector <vector <double> > addlater;      // A change to the rates Rtot which may be performed later when sampling is performed

	vector < vector <FEV> > fev;             // Stores all compartmental transitions
	
	vector <long> N;                         // The number of individuals in different compartments

	long sett;                               // Index used to track time changes in beta

	long tdnext, tdfnext;                    // Stores when the next future compartmental transition will occur

	MODEL  &model;
	vector <TRANS> &trans;
	vector <COMP> &comp;

	POPTREE &poptree;
	vector <LEVEL> &lev;

	DATA &data;
		
	public: 
		void gillespie(double ti, double tf, short siminf);
		void partinit(long p);
		void dofe();
		long nextinfection(short type);
		void addinfc(long c, double t);
		void addfev(double t, long tr, long i);
		void Lobs(short w, double invT);
		void copy(const PART &other, short fedivmin);
		void simmodel(long i, short enter, double t);
		void check(short num);
		
		void pack(vector <double> &pac, long num);
		void pack(vector <double> &pac, vector <long> &vec);
		void pack(vector <double> &pac, vector <double> &vec);
		void pack(vector <double> &pac, vector< vector <long> > &vec);
		void pack(vector <double> &pac, vector< vector <double> > &vec);
		void pack(vector <double> &pac, vector< vector <FEV> > &vec, short fedivmin, short fedivmax);
		long partpack(double *buffer, short fedivmin);
		long fevpack(double *buffer, short fedivmin, short fedivmax);
		void unpack(long &k, double *buffer, long &num);
		void unpack(long &k, double *buffer, vector <long> &vec);
		void unpack(long &k, double *buffer, vector <double> &vec);
		void unpack(long &k, double *buffer, vector< vector <long> > &vec);
		void unpack(long &k, double *buffer, vector< vector <double> > &vec);
		void unpack(long &k, double *buffer, vector< vector <FEV> > &vec, short fedivmin, short fedivmax);
		void partunpack( double *buffer, int max, short fedivmin);
};
