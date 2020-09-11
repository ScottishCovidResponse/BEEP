#ifndef BEEPMBP__CHAIN_HH
#define BEEPMBP__CHAIN_HH

#include <vector>

using namespace std;

#include "model.hh"
#include "areatree.hh"
#include "data.hh"
#include "simulate.hh"
#include "state.hh"
#include "jump.hh"

class ObservationModel;

class Chain                                             // Stores all the things related to an MCMC chain
{
public:
	Chain(const Details &details, const Data &data, const Model &model, const AreaTree &areatree,	const ObservationModel &obsmodel, const Output &output, unsigned int chstart);
	
	void standard_prop(double EFcut=0); 


	void sample_state();
	Status simulate(const vector <double>& paramv);
	void mbp_proposal(unsigned int th);
	vector <Event> event_compress(const vector < vector <Event> > &indev) const;
	void generate_particle(Particle &part) const;
	Status abcmbp_proposal(const vector <double> &param_propose, double EFcut);
	void stand_event_prop(double EFcut);
	
	unsigned int ch;                                      // Indexes the number of the chain 

	double invT;                                          // The inverse temperature

	Jump jump;
	
	State initial, propose;                               // The states in the initial and proposed states
	
private:
	
	void initialise_variables();
	Status mbp(const vector<double> &paramv);
	void copy_event_or_not(unsigned int n);
	void mbp_init();
		
	void mbp_compartmental_transitions(unsigned int i);
	unsigned int area_of_next_infection();
	void add_infection_in_area(unsigned int c, double t);
	void check(double t, unsigned int sett) const;
	void update_dQmap(const vector <EventRef> &trei, const vector <EventRef> &trep);
	void setup_susceptible_lists();
	void reset_susceptible_lists();
	void change_susceptible_status(unsigned int i, unsigned int st, unsigned int updateR);
	void construct_infection_sampler(const vector <double> &Qmi, const vector <double> &Qmp);
	void infsampler(const vector< vector<double> > &Qmap);
	void sortx(vector <EventRef> &x, vector <vector <Event> > &indev) const;
	void calcproposeQmap();
	void area_prop(unsigned int samp, unsigned int burnin);
	void area_prop2(unsigned int samp, unsigned int burnin, unsigned int th, double L0, const vector <double> &areasum, const vector < vector <double> >&mult, const vector < vector <double> > &add);
	void fixarea_prop(unsigned int samp, unsigned int burnin);
	
	vector < vector <short> > indmap;										  // A map which is used for fast update in update_dQmap 
	    
	vector <double> dQmap;                                // The difference in Q between the two states
	
	vector <unsigned int> dQbuflistv;                     // Used to efficiently calculate dQmap
	vector <unsigned int> dQbuflistq; 
	vector< vector <double> > dQbuf;

	vector <double> lam, lamsum;                          // Used when adding and removing individuals
	
	vector < vector <unsigned int> > both_susceptible_list;         // List of individuals which are both susceptible 
	vector <unsigned int> nboth_susceptible_list;  
	vector < vector <unsigned int> > propose_only_susceptible_list; // List of individuals where proposed state is suscptivble
	vector <unsigned int> npropose_only_susceptible_list;       
	vector < vector <unsigned int> > not_susceptible_list;                     // List of individuals not sus in either state
	vector <unsigned int> nnot_susceptible_list;
	vector <unsigned int> susceptible_list_ref;
	vector <unsigned int> susceptible_status;
	
	vector <vector <double> > new_infection_rate;         // Tree giving rate of new infections
	
	vector <int> N;                                       // The number of individuals in different compartments
	
	vector <int> popw;                                    // The population in w
	
	bool do_mbp_event;                                    // Set to true if MBPs on compartmental transitions needed
	
	const vector <Compartment> &comp;
	const vector <Level> &lev;
	const vector <Transition> &trans;
	
	const Details &details;
	const Data &data;
	const Model &model;
	const AreaTree &areatree;
	const ObservationModel &obsmodel;
	const Output &output;
};
#endif
