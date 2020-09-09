#ifndef BEEPMBP__CHAIN_HH
#define BEEPMBP__CHAIN_HH

#include <vector>

using namespace std;

#include "model.hh"
#include "poptree.hh"
#include "data.hh"
#include "simulate.hh"
#include "state.hh"

class Obsmodel;

class Chain                                             // Stores all the things related to an MCMC chain
{
public:
	Chain(const Details &details, const DATA &data, const MODEL &model, const POPTREE &poptree,	Obsmodel &obsmodel, unsigned int chstart);
	void sample_from_prior();
	unsigned int simulate(const vector <double>& paramv);
	void proposal(unsigned int th, unsigned int samp, unsigned int burnin);
	void standard_prop(unsigned int samp, unsigned int burnin, double EFcut=0);
	void setinitialQmap(unsigned int check);
	vector <FEV> event_compress(const vector < vector <FEV> > &indev) const;
	void initialise_from_particle(const Particle &part);
	void generate_particle(Particle &part) const;
	int abcmbp_proposal(const vector <double> &param_propose, double EFcut);
	
	unsigned int ch;                                      // The number of the chain (0=posterior, nchaintot-1=prior)            

	double invT;                                          // The inverse temperature

	long timeprop;                                        // The time for the proposals
	
	vector <float> paramjump;                             // The size of jumps in parameter space
	vector <unsigned int> ntr, nac;                       // The number of jumps tried and accepted

	float numaddrem;                                      // The size of adding and removing events
	unsigned int ntr_addrem, nac_addrem;    
	
	vector <float> paramjumpstand;                           // The size of jumps in parameter space (fixed event sequence)
	vector <unsigned int> ntrstand, nacstand;                   // The number of jumps tried and accepted

	float sigmajump;                                      // Used for jumping in sigma
	
	State initial, propose;                               // The states in the initial and proposed states
	
 	vector <double> paramval;                             // The values for the parameters

private:
	unsigned int mbp();
	void addindev(unsigned int i, vector <FEV> &indev, vector <EVREF> &x, vector <vector <EVREF> > &trev);
	unsigned int nextinfection();
	void addinfc(unsigned int c, double t);
	void check(unsigned int num, double t, unsigned int sett);
	void check_addrem();
	void updatedQmap(vector <EVREF> &trei, vector <EVREF> &trep);
	void setuplists();
	void resetlists();
	void changestat(unsigned int i, unsigned int st, unsigned int updateR);
	void constructRtot(vector <double> &Qmi, vector <double> &Qmp);
	double likelihood(vector < vector<double> > &Qmap, vector <EVREF> &x, vector <vector<FEV> > &indev, vector < vector <double> > &disc_spline, vector <double> &sus, vector <double> &areafac);
	void infsampler(vector< vector<double> > &Qmap);
	void sortx(vector <EVREF> &x, vector <vector <FEV> > &indev);
	void calcproposeQmap();
	void betaphi_prop( unsigned int samp, unsigned int burnin);
	void area_prop(unsigned int samp, unsigned int burnin);
	void area_prop2(unsigned int samp, unsigned int burnin, unsigned int th, double L0, vector <double> &areasum, vector < vector <double> >&mult, vector < vector <double> > &add);
	void fixarea_prop(unsigned int samp, unsigned int burnin);
	void addrem_prop(unsigned int samp, unsigned int burnin, double EFcut=0);
	void proposal_init();
		
	double Levi;         																	// The latent process likelihood
	
	vector < vector <short> > indmap;										  // A map which is used for fast update in updatedQmap 
	    
	vector <double> dQmap;                                // The difference in Q between the two states
	
	vector <unsigned int> dQbuflistv;                     // Used to efficiently calculate dQmap
	vector <unsigned int> dQbuflistq; 
	vector< vector <double> > dQbuf;

	vector <double> lam, lamsum;                          // Used when adding and removing individuals
	
	vector < vector <unsigned int> > indbothlist;         // List of individuals in area/demo states which are both susceptible 
	vector <unsigned int> nindbothlist;  
	vector < vector <unsigned int> > indponlylist;        // List of individuals in area/demo states where proposed state is sus
	vector <unsigned int> nindponlylist;       
	vector < vector <unsigned int> > indnotlist;          // List of individuals in area/demo states both 
	vector <unsigned int> nindnotlist;
	vector <unsigned int> indlistref;
	vector <unsigned int> stat;
	
	vector <vector <double> > Rtot;                       // Tree giving rate of new infections
	
	vector <int> N;                                       // The number of individuals in different compartments
	
	vector <int> popw;                                    // The population in w
	
	const vector <COMP> &comp;
	const vector <LEVEL> &lev;
	const vector <TRANS> &trans;
	
	const Details &details;
	const DATA &data;
	const MODEL &model;
	const POPTREE &poptree;
	Obsmodel &obsmodel;
};
#endif
