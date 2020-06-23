#pragma once

#include <vector>

using namespace std;

#include "model.hh"
#include "poptree.hh"
#include "data.hh"

class MBPCHAIN                                          // Stores all the things related to a MBP MCMC chain
{
	public:
	MBPCHAIN(DATA &data, MODEL &model, POPTREE &poptree);
		
	unsigned int ch;                                      // The number of the chain (0=posterior, nchaintot-1=prior)            
	double Li; 																						// The observation likelihood for the current state
	double invT;                                          // The inverse temperature 
	vector <float> paramjump;                             // The size of jumps in parameter space
	vector <unsigned int> ntr, nac;                       // The number of jumps tried and accepted
	long timeprop;                                        // The time for the proposals
	
	vector <EVREF> xi;                                    // Ordered list of references to infection events in init state
	vector <EVREF> xp;                                    // Ordered list of references to infection events in prop state
		                
	vector < vector <EVREF> > trevi;                      // Stores other event reference in initial state
	vector < vector <EVREF> > trevp;                      // Stores other event reference in initial state

 	vector <double> paramval;                             // The values for the parameters

	vector <double> dQmap;                                // The difference in Q between the two states
	vector< vector <double> > Qmapi;                      // The infectivty map for the initial state
	vector< vector <double> > Qmapp;			       		  		// The infectivity map for the proposed state

	vector <unsigned int> dQbuflistv;                     // Used to efficiently calculate dQmap
	vector <unsigned int> dQbuflistq; 
	vector< vector <int> > dQbuf;

	vector <double> lami;                                 // Total force of infecion for an area in the initial state
	vector <double> lamp; 	 															// Total force of infecion for an area in the proposed state
	
	vector < vector <unsigned int> > indbothlist;         // List of individuals in area/demo states which are both susceptible 
	vector <unsigned int> nindbothlist;  
	vector < vector <unsigned int> > indponlylist;        // List of individuals in area/demo states where proposed state is sus
	vector <unsigned int> nindponlylist;       
	vector < vector <unsigned int> > indnotlist;          // List of individuals in area/demo states both 
	vector <unsigned int> nindnotlist;
	vector <unsigned int> indlistref;
	vector <unsigned int> stat;
	
	vector < vector <FEV> > indevi;                       // The individual event sequences for the initial state
	vector < vector <FEV> > indevp;                       // The individual event sequences for the proposed state
	
	vector <unsigned int> indinfi;                        // All the infected individuals in the initial state
	vector <unsigned int> indinfp;                        // All the infected individuals in the initial state
	
	vector <vector <double> > Rtot;                       // Tree giving rate of new infections
	
	vector <int> N;                                       // The number of individuals in different compartments
	
	double betai, phii;                                   // The values of beta and phi in the initial state
	double betap, phip;                                   // The values of beta and phi in the proposed state

	vector <int> popw;                                    // The population in w
	
	DATA &data;
	MODEL &model;
	POPTREE &poptree;
	
	vector <TRANS> &trans;
	vector <COMP> &comp;
	vector <LEVEL> &lev;
	 
	public:
		void init(DATA &data, MODEL &model, POPTREE &poptree, double invTstart, vector < vector <FEV> > &indev, unsigned int chstart);
		void proposal(DATA &data, MODEL &model, POPTREE &poptree, unsigned int th, unsigned int samp, unsigned int burnin);
		void param_prop();
		void setQmapi();
						
	private:
		unsigned int mbp();
		void addindev(unsigned int i, vector <FEV> &indev, vector <EVREF> &x, vector <vector <EVREF> > &trev);
		unsigned int nextinfection();
		void addinfc(unsigned int c, double t);
		void check(unsigned int num, double t, unsigned int sett);
		void updatedQmap(unsigned int sett);
		void setuplists();
		void resetlists();
		void changestat(unsigned int i, unsigned int st);
		void constructRtot(unsigned int sett);
		double likelihood();
};
