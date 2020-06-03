// The PMBP algorithm

#include <fstream>
#include <iostream>
#include <sstream>
#include "stdlib.h"
#include "math.h"
#include "assert.h"

using namespace std;

#include "timers.hh"
#include "model.hh"
#include "utils.hh"
#include "PART.hh"
#include "output.hh"
#include "pack.hh"

class MBPPART                                 // Stores all the things related to a particle 
{
	public:
	MBPPART(DATA &data, MODEL &model, POPTREE &poptree);

	DATA &data;
	MODEL  &model;
	POPTREE &poptree;
	
	vector <double> MIi;    
	vector <double> MIp;  
	vector <double> susboth;    
	vector <double> susp; 
	vector <double> lami;    
	vector <double> lamp; 
	long xitdnext, xitdfnext;
	
	vector <short> stati;
	vector <short> statp;
	
	public:
		void mbpinit(long p, vector < vector <FEV> > &xi);
		void mbp(double ti, double tf, vector < vector <FEV> > &xi, long &timegen);
		long mbpnextinfection();
		void mbpdofe();
		void mbpxidofe(vector < vector <FEV> > &xi);
		void mbpaddinfc(long c, double t);
	
}