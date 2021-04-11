#ifndef BEEPMBP__Model_HH
#define BEEPMBP__Model_HH

#include <vector>

using namespace std;

#include "struct.hh"
#include "data.hh"

class Model                                             // Stores information about the model
{
	public:
		Model(const Inputs &inputs, const Details &details, Data &data);
		
		vector <vector <unsigned int> > R_spline_ref;       // Denotes which spline refers to variation in R [inft][area]
		vector <vector <unsigned int> > efoi_spline_ref;    // Denotes which spline refers to variation in extornal force of infection
		
		unsigned int geo_spline_ref;                        // Denotes which spline refers to geographical spread
		
		vector <Param> param;                               // Information about parameters in the model
		vector <Transition> trans;                          // Stores model transitions
		vector <Compartment> comp;	                        // Stores model compartments
		
		unsigned int ninfection_trans;                      // The number of infection transitions
		vector <unsigned int> infection_trans;              // Stores infection transitions in the model
	
		unsigned int start_compartment;                     // This defines the compartment individuals start in
	
		vector <vector <unsigned int> > parameter_type;     // Different type of pararmeters
		
		unsigned int maximum_infected;                      // The maximum number of infected individuals
		
		vector <Spline> spline;                             // Stores the splines used in the model  
	
		vector < vector < vector <double> > > disc_spline_grad;// The gradient in a spline w.r.t. a variable
	
		unsigned int region_effect;                         // Set to 0 of no effect, 1 if random effect, 2 if fixed effect

		vector< vector <unsigned int> > sus_param_list;    	// The parameters related to fixed effect for susceptibility
		vector< vector <double> > sus_param_popfrac;        // Takes into account the population fraction 
		vector< vector <unsigned int> > sus_ref;            // References sus_param_list
		
		vector <unsigned int> covariate_param;              // The parameters related to covariates for areas
		
		unsigned int sigma_param;                           // The standard deviation of the regional effect (if random effect)
		vector <unsigned int> region_effect_param;          // The parameters related to regional effects
		vector < vector <double> > region_effect_shift;     // The shift in region effect to account for change in covariate 

		
		vector <unsigned int> param_not_fixed;              // A list of all parameters which actualy change 

		vector < vector < vector <double> > > efoi_agedist; // The age distribution for the external force of infection
		
		CounterMod countermod;                              // Stores any counerfactual modifications to the model
		
		double get_infectivity_dif(unsigned int inft, unsigned int tr, const vector <double> paramv) const;
		
		double param_prior_sample(unsigned int th) const;
		vector <double> sample_from_prior() const;
		
		vector < vector <double> > calculate_Vinv(const vector < vector <double> > &transrate, const unsigned int inft) const;
		double calculate_area_av(const vector <double> &areafactor) const;
		vector < vector <double> > calculate_F(const vector<double> &paramv, const vector <double> &susceptibility, const vector <vector <double> > &A, const unsigned int inft) const;
		vector < vector <double> > calculate_NGM(const vector < vector<double> > &F, const vector < vector<double> > &Vinv) const;
		double calculate_R_beta_ratio_using_NGM(const vector<double> &paramv, const vector <double> &susceptibility, const vector <vector <double> > &A, const vector < vector <double> > &Vinv, const unsigned int inft) const;
		vector <double> calculate_probreach(const vector<double> &paramv) const;
		vector <double> calculate_external_ninf(const vector<double> &paramv) const;
		vector <SplineOutput> get_spline_output(const vector <double> &paramv) const;
		vector <DerivedParam> get_derived_param(const vector<double> &paramv, const vector <double> &susceptibility, const vector <vector <double> > &A, const vector < vector <double> > &transrate) const;
		double calculate_generation_time(const vector<double> &paramv, const vector <double> &susceptibility, const vector <vector <double> > &A, const vector < vector <double> > &transrate, const unsigned int inft) const;
		bool do_mbp_events(const vector <double> &parami, const vector <double> &paramp) const;
		double prior(const vector<double> &paramv) const;
		vector <double> create_disc_spline(unsigned int ref, const vector<double> &paramv) const;
		vector <double> create_susceptibility(const vector<double> &paramv) const;  
		vector <double> create_areafactor(const vector<double> &paramv) const;
		vector <CompTransProb> create_comptransprob(const vector<double> &paramv) const;
		vector < vector <double> > create_transrate(vector <CompTransProb> &comptransprob, const vector<double> &paramv) const;
		vector <CompProb> create_compprob(const vector <CompTransProb> &comptransprob, unsigned int inft) const;
		void create_Ntime(vector < vector < vector <double> > > &Ntime, const vector < vector<double> > &disc_spline) const;
		bool equal(const vector < vector <double> > &Ntime_1, const vector < vector <double> > &Ntime_2) const;
		vector < vector <double> > create_beta_R_ratio(const vector <double> &susceptibility, const vector <double> &areafactor, const vector <double> &paramv, vector < vector < vector <double> > > &Ntime, const vector < vector <double> > &transrate) const;
		vector < vector <vector <double> > > calculate_R_age(const vector <double> &paramv) const;
		bool inbounds(const vector <double> &paramv) const;
		double dPr_dth(unsigned int th, const vector<double> &paramv) const;
		string print() const;
		void print_parameter_types();
		void print_transition(unsigned int tr) const;
		string trans_str(unsigned int tr) const;
		
	private:
		vector < vector <unsigned int> > inf_state;         // Makes a list of all the infectious states (used to calculate NGM)
		vector < vector <unsigned int> > inf_state_ref;     // Converts from a compartment to inf_state  
	
		vector <ProbReach> prob_reach;                      // Calculates the probability of reaching a certain conpartment
		
		void load_model();
		void add_trans_comps();
		void add_compartment(const string& name, const ParamSpec &inf, const string &inf_trans);
	
		unsigned int add_parameter(const ParamSpec ps, ParamType type);	
		void add_transition(const string& from_origin, const string& from, const string& to, const vector <ParamSpec>& distspec, const vector <ParamSpec>& probspec, unsigned int type, double mean_factor);		
		void add_region_effect();
		void add_splines();
		void add_probreach();
		void set_disc_spline_grad();
		double spline_prior(const vector<double> &paramv) const;
		void spline_sample(vector <double> &paramv) const;
		void complete_datatables();
		void complete_counterfactual();
		void set_susceptibility_param();
		void create_spline(const string name, const string desc, const vector <ParamSpec> &ps, const vector<unsigned int> &bp, const ParamSpec param_fac, Smooth smooth_type, ParamType type);
		void set_infected_uninfected();
		unsigned int get_compartment(string compname);
		unsigned int get_transition(string transname);
		void update_spline_ref(string &name, string &desc, string inft_str, string area_str, const vector <double> efoi_ad, vector < vector <unsigned int> >& spline_ref);
		void counterfactual_modification();
			
		const Details &details;
		Data &data;
		const Inputs &inputs;
};
#endif
