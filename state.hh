#ifndef BEEPMBP__StatisticsE_HH
#define BEEPMBP__StatisticsE_HH

using namespace std;

#include "struct.hh"
#include "model.hh"
#include "mvn.hh"
#include "obsmodel.hh"

class State                                              // Store information about the state of the system
{
	public:
		State(const Details &details, const Data &data, const Model &model, const ObservationModel &obsmodel);
			
		double EF; 																				   // The error function
		double Pr; 																		       // The prior probability
		
		vector <double> paramval;                            // The parameter values
	
		vector < vector< vector <double> > > Imap;           // The infectivity map coming from other areas
		vector < vector< vector <double> > > Idiag;          // The infectivity coming from within an area
		
		vector < vector < vector < vector <double> > > > transmean;// The mean number of transitions (for Poisson distribution)  
				
		vector < vector < vector < vector <int> > > > transnum;// Realised number of transitions [sett][area][tr][dp]  
		
		vector < vector < vector < vector <int> > > > pop;   // The populations in different compartments [sett][area][comp][dp] 

	 	/* The quantities below are all derived from the model parameter values */
		vector <double> susceptibility;                      // The susceptibility for different demographic categories
		vector <double> areafactor;                          // The modification due to area effects
		vector < vector <double> > disc_spline;              // A discretisation of the splines	
		vector <CompTransProb> comptransprob;                // The probabilities of going down transitions
		vector < vector <double> > transrate;                // Rates for transitions
		vector < vector < vector <double> > > Ntime;         // Time variation in age matrix
		vector < vector <double> > beta_R_ratio;             // The ratio between beta and R [inft][sett]
	
		void set_param(const vector <double> &paramv);
		void set_transmean(unsigned int sett, unsigned int c);
		void set_EF();
		void set_Pr();
		void initialise_from_particle(const Particle &part);
		Particle create_particle() const;
		void obsmodel_prop(MVN &mvn, double EFcut, double invT, ObsModelMode obsmodel_mode);
		void simulate(const vector <double> &paramval);
		void simulate(const unsigned int ti, const unsigned int tf);
		unsigned int get_ninf_total();
		Sample create_sample() const;
		ParamSample create_param_sample() const;
		void save(string file) const;
		void load(string file);
		void check();
		
		// These functions are used for MBPs //
		void set_Imap_sett(unsigned int sett);
		void set_Imap_using_dI(unsigned int sett, const State *state, const vector< vector <double> > &dImap, const vector < vector <double> > &dIdiag);
		void update_I_from_transnum(vector < vector <double> > &Ima, vector < vector <double> > &Idia, const vector < vector < vector <int> > > &dtransnum) const;
		void update_pop(unsigned int sett);
		void pop_init();
		// These functions are used for MBPs //
		
	private:
		void set_Imap(unsigned int check);
		vector <double> get_NMI(unsigned int sett, unsigned int inft, unsigned int c);
		string print_populations(unsigned int sett) const;
		
		const vector <Compartment> &comp;
		const vector <Transition> &trans;
		const vector <Param> &param;
		
		const Details &details;
		const Data &data;
		const Model &model;
		const ObservationModel &obsmodel;
};

#endif