#ifndef BEEPMBP__CHAIN_HH
#define BEEPMBP__CHAIN_HH

#include <vector>

using namespace std;

#include "model.hh"
#include "poptree.hh"
#include "data.hh"
#include "simulate.hh"
#include "state.hh"
#include "jump.hh"

class Obsmodel;

class Chain                                             // Stores all the things related to an MCMC chain
{
public:
	Chain(const Details &details, const DATA &data, const MODEL &model, const POPTREE &poptree,	const Obsmodel &obsmodel, unsigned int chstart);
	
	void standard_prop(double EFcut=0); 

	void sample_from_prior();
	Status simulate(const vector <double>& paramv);
	void mbp_proposal(unsigned int th);
	vector <FEV> event_compress(const vector < vector <FEV> > &indev) const;
	void generate_particle(Particle &part) const;
	Status abcmbp_proposal(const vector <double> &param_propose, double EFcut);
	void stand_event_prop(double EFcut);
	
	unsigned int ch;                                      // Indexes the number of the chain 

	double invT;                                          // The inverse temperature

	Jump jump;
	
	State initial, propose;                               // The states in the initial and proposed states
	
private:
	
	Status mbp(const vector<double> &paramv);
	void mbpmodel(vector <FEV> &evlisti, vector <FEV> &evlistp);
	unsigned int nextinfection();
	void addinfc(unsigned int c, double t);
	void check(double t, unsigned int sett) const;
	void updatedQmap(const vector <EVREF> &trei, const vector <EVREF> &trep);
	void setuplists();
	void resetlists();
	void changestat(unsigned int i, unsigned int st, unsigned int updateR);
	void constructRtot(const vector <double> &Qmi, const vector <double> &Qmp);
	void infsampler(const vector< vector<double> > &Qmap);
	void sortx(vector <EVREF> &x, vector <vector <FEV> > &indev) const;
	void calcproposeQmap();
	void area_prop(unsigned int samp, unsigned int burnin);
	void area_prop2(unsigned int samp, unsigned int burnin, unsigned int th, double L0, const vector <double> &areasum, const vector < vector <double> >&mult, const vector < vector <double> > &add);
	void fixarea_prop(unsigned int samp, unsigned int burnin);
	//void proposal_init(const vector <double> &paramv);
		
	double Levi;         																	// The latent process likelihood
	
	vector < vector <short> > indmap;										  // A map which is used for fast update in updatedQmap 
	    
	vector <double> dQmap;                                // The difference in Q between the two states
	
	vector <unsigned int> dQbuflistv;                     // Used to efficiently calculate dQmap
	vector <unsigned int> dQbuflistq; 
	vector< vector <double> > dQbuf;

	vector <double> lam, lamsum;                          // Used when adding and removing individuals
	
	vector < vector <unsigned int> > indbothlist;         // List of individuals which are both susceptible 
	vector <unsigned int> nindbothlist;  
	vector < vector <unsigned int> > indponlylist;        // List of individuals where proposed state is sus
	vector <unsigned int> nindponlylist;       
	vector < vector <unsigned int> > indnotlist;          // List of individuals not sus in either state
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
	const Obsmodel &obsmodel;
};
#endif
