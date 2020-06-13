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

	int ch;                                               // The number of the chain (0=posterior, nchaintot-1=prior)            
	double Li; 																						// The observation likelihood for the current state
	double invT;                                          // The inverse temperature 
	vector <float> paramjump;                             // The size of jumps in parameter space
	vector <int> ntr, nac;                                // The number of jumps tried and accepted
	long timeprop;                                        // The time taken for the proporals
	
	int fediv;                                            // The number of divisions the global timeline is divided into
	
 	vector < vector <FEV> > xi;                           // The event timeline in the initial state
	vector < vector <FEV> > xp;                           // The event timeline in the initial state
	
 	vector <double> paramval;                             // The values for the parameters

	vector <double> MIi;                                  // The infectivty map for the initial state
	vector <double> MIp; 																	// The infectivity map for the proposed state

	vector <double> susboth;                              // Total sus. of an area when both initial and proposal states sus.
	vector <double> susp;                                 // Total sus. of an area when only proposed state sus.

	vector <double> lami;                                 // Total force of infecion for an area in the initial state
	vector <double> lamp; 	 															// Total force of infecion for an area in the proposed state
	
	int xitdnext, xitdfnext;                              // Stores the next event in the initial state 
	int xptdnext, xptdfnext;                 					    // Stores the next event in the final state
	
	vector< vector <FEV> > indevi;                        // The individual event sequences for the initial state
	vector< vector <FEV> > indevp;                        // The individual event sequences for the proposed state
	
	int ninftot;                                          // The total number of infected individuals
	int ninftotprop;   
	
	vector <vector <double> > Rtot;                       // Tree giving rate of new infections
	
	vector <int> N;                                       // The number of individuals in different compartments
	
	int sett;                                             // Index used to track time changes in beta
	
	double betai, phii;                                   // The values of beta and phi in the initial state
	double betap, phip;                                   // The values of beta and phi in the proposed state

	DATA &data;
	MODEL &model;
	POPTREE &poptree;
	
	vector <TRANS> &trans;
	vector <COMP> &comp;
	vector <LEVEL> &lev;
	
	public:
		void init(DATA &data, MODEL &model, POPTREE &poptree, double invTstart, vector < vector <FEV> > &xistart, vector < vector <FEV> > &indev, int chstart);
		void proposal(DATA &data, MODEL &model, POPTREE &poptree, int th, int samp, int burnin);
		void mbpinit();
		void mbp();
		int nextinfection();
		void xpdofe(int update);
		void xidofe();
		void addinfc(int c, double t);
		void updateRtot(int c);
		void check(int num);
		void MIupdate(int i, int tra, int upMIi, int upMIp);
		void addxp(FEV fe, double period, double tnow);
};
