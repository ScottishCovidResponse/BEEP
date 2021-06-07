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
				
		vector < vector < vector < vector <double> > > > transnum;// Realised number of transitions [sett][area][tr][dp]  
		
		vector < vector < vector < vector <double> > > > pop;// The populations in different compartments [sett][area][comp][dp] 

	 	/* The quantities below are all derived from the model parameter values */
		vector <double> susceptibility;                      // The susceptibility for different demographic categories
		vector < vector <double> > areafactor;               // The modification due to area effects at a particular time
		vector < vector <double> > disc_spline;              // A discretisation of the splines	
		vector < vector <double> > transrate;                // Rates for transitions
		vector < vector < vector <double> > > Ntime;         // Time variation in age matrix
		vector < vector < vector <double> > > beta;          // The transmission rate
	
		void set_param(const vector <double> &paramv);
		void set_transmean(const unsigned int sett, const unsigned int c);
		void set_EF();
		void set_Pr();
		void initialise_from_particle(const Particle &part);
		Particle create_particle(const unsigned int run) const;
		void simulate(const vector <double> &paramval);
		void simulate(const unsigned int ti, const unsigned int tf);
		Sample create_sample() const;
		ParamSample create_param_sample(const unsigned int run) const;
		void save(const string file) const;
		void load(const string file, const unsigned int ntimesteps);
		void check();
		
		// START These functions are used for MBPs //
		void set_Imap_sett(const unsigned int sett);
		void set_Imap_using_dI(const unsigned int sett, const State *state, const vector< vector <double> > &dImap, const vector < vector <double> > &dIdiag);
		void update_I_from_transnum(vector < vector <double> > &Ima, vector < vector <double> > &Idia, const vector < vector < vector <double> > > &dtransnum) const;
		void update_pop(const unsigned int sett);
		void pop_init();
		// END //
		
	private:
		void set_Imap(unsigned int check);
		vector <double> get_NMI(const unsigned int sett, const unsigned int inft, const unsigned int c);
		string print_populations(const unsigned int sett) const;
		
		const vector <Compartment> &comp;
		const vector <Transition> &trans;
		const vector <Param> &param;
		
		const Details &details;
		const Data &data;
		const Model &model;
		const ObservationModel &obsmodel;
};

#endif