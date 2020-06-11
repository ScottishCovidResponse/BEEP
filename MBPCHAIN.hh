#pragma once

#include <vector>

using namespace std;

#include "model.hh"
#include "poptree.hh"
#include "data.hh"

class MBPCHAIN                                 // Stores all the things related to a particle 
{
	public:
	MBPCHAIN(DATA &data, MODEL &model, POPTREE &poptree);

	// Stores all the information swaped between chains
	int ch;                      
	double Li;
	double invT;
	vector <float> paramjump;
	vector <int> ntr, nac;
	long timeprop;
	//
	
	int fediv;
	
	vector < vector <FEV> > xi;
	int ninftot;   //The total number of infected individuals
	vector <double> paramval;

	DATA &data;
	MODEL &model;
	POPTREE &poptree;
	
	vector <TRANS> &trans;
	vector <COMP> &comp;
	vector <LEVEL> &lev;
	
	vector <double> MIi;    
	vector <double> MIp;  
	vector <double> susboth;    
	vector <double> susp; 
	vector <double> lami;    
	vector <double> lamp; 	
	
	vector <vector <double> > Rtot;
	
	int tdnext, tdfnext;                    // Stores when the next future compartmental transition will occur
	int xitdnext, xitdfnext;
	
	vector <int> statilist1;    // Lists all the changes to stati and statp
	vector <int> statplist1, statplist2;
		
	vector <int> stati;
	vector <int> statp;
	
	vector < vector <FEV> > fev;             // Stores all compartmental transitions
	int ninftotprop;   
	
	vector <int> N;                         // The number of individuals in different compartments
	
	int sett;                               // Index used to track time changes in beta
	
	double betai, phii, betap, phip; 

	public:
		void init(DATA &data, MODEL &model, POPTREE &poptree, double invTstart, vector < vector <FEV> > &xistart, int chstart);
		void proposal(DATA &data, MODEL &model, POPTREE &poptree, int th, int samp, int burnin);
		void mbpinit();
		void mbp();
		int nextinfection();
		void dofe(int update);
		void xidofe();
		void addinfc(int c, double t);
		void updateRtot(int c);
		void check(int num);
		void MIupdate(int i, int tra, int upMIi, int upMIp);
};
