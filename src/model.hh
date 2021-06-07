#ifndef BEEPMBP__Model_HH
#define BEEPMBP__Model_HH

#include <vector>

using namespace std;

#include "utils.hh"
#include "data.hh"

class Model                                             // Stores information about the model
{
	public:
		Model(Inputs &inputs, const Details &details, Data &data, Mpi &mpi);
		
		vector <unsigned int> R_spl_ref;                    // Gives a reference to Rspline_info about Rspline [inf][c]
		vector <SplineInfo> Rspline_info;                   // Stores information about Rsplines (so Reff can be generated) 
	 	
		vector <vector <unsigned int> > efoi_spl_ref;       // Gives a reference to efoi_info about efoi [inf][c]
		vector <SplineInfo> efoispline_info;                // Stores information about efoi splines
		
		unsigned int geo_spline_ref;                        // Denotes which spline refers to geographical spread
		
		vector <Param> param;                               // Information about parameters in the model
		vector <Transition> trans;                          // Stores model transitions
		vector <Compartment> comp;	                        // Stores model compartments
		
		unsigned int infection_trans;                       // Stores infection transitions in the model
	
		unsigned int start_compartment;                     // This defines the compartment individuals start in
	
		vector < vector <DemocatChangeMap> > democat_change_map;// Refers to changes in demography with time
	
		vector <vector <unsigned int> > parameter_type;     // Different type of pararmeters
		
		vector <Spline> spline;                             // Stores the splines used in the model  
	
		vector < vector < vector <double> > > disc_spline_grad;// The gradient in a spline w.r.t. a variable
	
		vector <unsigned int> level_effect_param_list;      // The parameters related to level effects
		
		vector <unsigned int> area_effect_param_list;       // The parameters related to area effects
		
		vector <unsigned int> covariate_param;              // The parameters related to covariates for areas
		
		RegionEffect region_effect;                         // Stores information about regional effects
		
		vector <Dirichlet> dirichlet;                       // Stores information about dirichlet distributions
		
		vector <unsigned int> param_not_fixed;              // A list of all parameters which actualy change 

		ModelMod modelmod;                                  // Stores any modifications to the model
		
		double get_infectivity_dif(const unsigned int tr, const vector <double> &paramv) const;
		double param_prior_sample(const unsigned int th, const vector <double> &paramv) const;
		vector <double> sample_from_prior() const;
		vector < vector <double> > calculate_Vinv(const vector < vector <double> > &transrate, const unsigned int st) const;
		double calculate_area_av(const vector < vector <double> > &areafactor) const;
		vector < vector <double> > calculate_F(const vector<double> &paramv_dir, const vector <double> &susceptibility, const vector <vector <double> > &A, const unsigned int sp, const vector <double> &democatpos_dist) const;
		vector < vector <double> > calculate_NGM(const vector < vector<double> > &F, const vector < vector<double> > &Vinv) const;
		double calculate_R_beta_ratio_using_NGM(const vector<double> &paramv_dir, const vector <double> &susceptibility, const vector <vector <double> > &A, const vector < vector <double> > &Vinv, const unsigned int sp, const vector <double> &democatpos_dist) const;
		vector <double> calculate_probreach(const vector<double> &paramv_dir, const unsigned int st) const;
		vector <double> calculate_external_ninf(const vector<double> &paramv_dir) const;
		vector <SplineOutput> get_spline_output(const vector <double> &paramv_dir, const vector < vector < vector < vector <double> > > > &transnum) const;
		vector <DerivedParam> get_derived_param(const vector<double> &paramv_dir, const vector <double> &susceptibility, const vector <vector <double> > &A, const vector < vector <double> > &transrate) const;
		double calculate_generation_time(const vector<double> &paramv_dir, const vector <double> &susceptibility, const vector <vector <double> > &A, const vector < vector <double> > &transrate, const vector <double> &democatpos_dist, const unsigned int st) const;
		bool do_mbp_events(const vector <double> &parami, const vector <double> &paramp) const;
		double prior(const vector<double> &paramv) const;
		vector < vector <double> > create_disc_spline(const vector<double> &paramv_dir) const;
		vector <double> create_disc_spline(const unsigned int ref, const vector<double> &paramv) const;
		vector <double> create_susceptibility(const vector<double> &paramv_dir) const;
		vector < vector <double> > create_areafactor(const vector<double> &paramv_dir) const;
		vector < vector <double> > create_transrate(const vector<double> &paramv_dir) const;
		vector <CompProb> create_compprob(const vector<double> &paramv_dir, const unsigned int st) const;
		vector < vector < vector <double> > > create_Ntime(const vector < vector<double> > &disc_spline) const;
		bool equal(const vector < vector <double> > &Ntime_1, const vector < vector <double> > &Ntime_2) const;
		vector < vector < vector <double> > > calculate_beta_from_R(const vector <double> &susceptibility, const vector <double> &paramv_dir, const vector < vector < vector <double> > > &Ntime, const vector < vector <double> > &transrate, const vector < vector <double> > &disc_spline) const;
		vector < vector <double> > calculate_R_eff(const vector <double> &paramv_dir, const vector < vector < vector < vector <double> > > > &transnum, const unsigned int st) const;
		vector < vector <double> > calculate_R_age(const vector <double> &paramv_dir) const;
		void setup_democat_change_map();
		bool inbounds(const vector <double> &paramv) const;
		double dPr_dth(const unsigned int th, const vector<double> &paramv) const;
		vector <double> dirichlet_correct(const vector <double> &paramval) const;
		void add_dirichlet(unsigned int th_min, const vector <unsigned int> &th_list, const vector <double> &frac);
		void set_strain();
		string list_dir_param(const vector <DirichletParam> &list) const;
		void set_dirichlet();
		string print() const;
		string print_prior(const unsigned int p) const;
		void print_parameter_types();
		void print_param(const string name, const vector <double> &paramv) const;
		
	private:
		vector <unsigned int> inf_state;                    // Makes a list of all the infectious states (used to calculate NGM)
		vector <unsigned int> inf_state_ref;                // Converts from a compartment to inf_state  
	
		vector <ProbReach> prob_reach;                      // Calculates the probability of reaching a certain conpartment
		
		void load_model();
		void add_comps();
		void add_trans();
		void setup_infectivity();
		void trans_check() const;
		void add_compartment(const string& name, const vector <ParamSpec> &mean_spec, const string &mean_dep, const ParamSpec &inf, const unsigned int num, const unsigned int shape);
		void add_internal_transition(const unsigned int c_from, const unsigned int c_to);
		void add_transition(const string& from, const string& to, const TransInf inf, const vector <ParamSpec>& prob_spec, const string prob_dep);
		double get_val1(const unsigned int th, const vector <double> &paramv) const;
		double get_val2(const unsigned int th, const vector <double> &paramv) const;
		void get_prior_val(const string name, const string st, double &val, unsigned int &valparam);
		unsigned int add_parameter(const ParamSpec ps, const ParamType type);	
		void add_region_effect();
		void add_splines();
		void add_probreach();
		void set_disc_spline_grad();
		double spline_prior(const vector<double> &paramv) const;
		void spline_sample(vector <double> &paramv) const;
		void complete_datatables();
		void setup_modification();
		void create_spline(const string name, const string desc, const vector <ParamSpec> &ps, const vector<unsigned int> &bp, const ParamSpec param_fac, const Smooth smooth_type, const ParamType type);
		void set_infected_uninfected();
		unsigned int get_compartment(const string compname, const ErlangPos pos) const;
		unsigned int get_transition(const string transname);
		unsigned int find_param(string name, string em) const;
		bool prior_order_correct(const vector <double> &paramv) const;
		void prior_order();
		void set_parameter_type();
		void update_spl_ref(string &name, string &desc, string area_str, const vector <double> efoi_ad, vector <unsigned int>& spl_ref, vector <SplineInfo> &spline_info) const;
		void counterfactual_modification();
			
		const Details &details;
		Data &data;
		Inputs &inputs;
		Mpi &mpi;
};
#endif
