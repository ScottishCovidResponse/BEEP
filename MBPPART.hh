#pragma once

#include <vector>

using namespace std;

#include "model.hh"
#include "poptree.hh"
#include "data.hh"

class MBPPART                                 // Stores all the things related to a particle 
{
	public:
	MBPPART(DATA &data, MODEL &model, POPTREE &poptree);

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
	
	long tdnext, tdfnext;                    // Stores when the next future compartmental transition will occur
	long xitdnext, xitdfnext;
	
	vector <long> statilist1;    // Lists all the changes to stati and statp
	vector <long> statplist1, statplist2;
		
	vector <short> stati;
	vector <short> statp;
	
	vector < vector <FEV> > fev;             // Stores all compartmental transitions
	
	vector <long> N;                         // The number of individuals in different compartments

	long sett;                               // Index used to track time changes in beta
	
	double betai, phii, betap, phip; 

	double Li;                               // The observation likelihood
 
	public:
		// copied
		void simmodel(long i, double t);
		void addfev(double t, long tra, long i, short done);
		void Lobs(short w, vector < vector <FEV> > &fev);
	
		void mbpinit(long p, vector < vector <FEV> > &xi);
		void mbp(double ti, double tf, vector < vector <FEV> > &xi);
		void recreate(double ti, double tf, vector < vector <FEV> > &xi, vector < vector <FEV> > &xp);
		long mbpnextinfection();
		void mbpdofe(short update);
		void mbpxidofe(vector < vector <FEV> > &xi);
		void mbpaddinfc(long c, double t);
		void updateRtot(long c);
		void check(short num);
		void MIupdate(long i, long tra, short upMIi, short upMIp);
		void copy(const MBPPART &other, short fedivmin);
		void partpack(short fedivmin);
		void partunpack(short fedivmin);
};
