#ifndef BEEP__StatisticsE_HH
#define BEEP__StatisticsE_HH

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

		/* Uses in DA and steering */
		vector < vector < vector <double> > > pop_covar;     // Covariance matrix for populations (used in approx)
	
		vector < vector <double> > obs_mean;                 // The mean of the future observation distribution   
		vector < vector < vector <double> > > obs_covar;     // Covariance matrix for future observation model (used in approx)
	
		vector <TransApprox> trans_approx;                    // Stores all unique transitions in the model
		unsigned int ntrans_approx;                          // The total number of transitions
		vector <CompApprox> comp_approx;
		unsigned int ncomp_approx;                           // The total number of compartments
		vector < vector < vector <unsigned int> > > comp_approx_ref; // References compartments to approximate method
		vector < vector < vector <unsigned int> > > trans_approx_ref; // References transition to approximate method
		vector <unsigned int> comp_without_S;                // Lists compartments which are not susceptible
		unsigned int ncomp_without_S;
		vector <unsigned int> comp_with_S;                   // Lists compartments which are susceptible
		vector <vector <unsigned int> > comp_with_S_list;    // Lists compartments which are not susceptible corresponding to sus
		unsigned int ncomp_with_S;
		vector < vector <double> > total_pop;                // The conserved population in areas and demographic categories
		vector < vector <double> > T;                        // Maps transitions to changes in compartment populations
		vector < vector <double> > Tinv_without_S;           // Maps changes in compartment populations with transitions
		vector < vector <double> > dtransmean_dp;            // The gradient in mean number of transitions (steering only)
		vector < vector <unsigned int> > obs_clist;          // Stores the compartments used for observations
		vector < vector <unsigned int> > obs_clist_without_S;// Stores the compartments used for observations (removing susceptible)
		/* End */
	
	 	/* The quantities below are all derived from the model parameter values */
		vector <double> susceptibility;                      // The susceptibility for different demographic categories
		vector < vector <double> > areafactor;               // The modification due to area effects at a particular time
		vector < vector <double> > disc_spline;              // A discretisation of the splines	
		vector < vector <double> > transrate;                // Rates for transitions
		vector < vector < vector <double> > > Ntime;         // Time variation in age matrix
		vector < vector < vector <double> > > beta;          // The transmission rate
	
		double genT;                                         // The generation time
	
		void set_param(const vector <double> &paramv);
		void democat_change_pop_adjust(const unsigned int sett);
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
		
		void check(const unsigned int checknum);
		
		// START These functions are used in steering //
		void set_transmean_gradients(const unsigned int sett);
		void test_transmean_gradients(const unsigned int sett);
		double likelihood_approx(const vector <double> &paramval, vector <ObsSlice> &obs_slice, const double invT, const Accuracy ac);
		void future_obs_approx(vector <ObsSlice> &obs_slice, const double invT, const Accuracy ac);
		Particle posterior_particle_sample(const vector <double> &paramval, vector <ObsSlice> &obs_slice, const unsigned int P, const double invT);
		Particle posterior_sample(const vector <double> &paramval, vector <ObsSlice> &obs_slice, const double invT, const Accuracy ac);
		//void obs_slice_calculate(ObsSlice &os, const double invT);
		void covar_change(vector < vector <double> > &covar_new, const vector < vector <double> > &covar, const vector <double> &mu) const;
		void covar_change_fast(vector < vector <double> > &covar_new, const vector < vector <double> > &covar, const vector <double> &mu, const Accuracy ac, const bool forward) const;
		double apply_observation(const ObsSlice &os, vector < vector <double> > &covar, const double invT);
		double apply_observation_fast(const ObsSlice &os, vector < vector <double> > &covar, const double invT, const Accuracy ac);
		double apply_exact_observation(const ObsSlice &os, vector < vector <double> > &covar);
		void correct_mu_covar(vector <double> &mu, vector < vector <double> > &covar) const;
		vector <double> remove_S(const vector <double> &vec) const;
		vector <double> add_S(const vector <double> &vec) const;
		vector < vector <double> > remove_S(const vector < vector <double> > &M) const;
		vector < vector <double> > add_S(const vector < vector <double> > &M) const;
		void valid_matrix(vector < vector <double> > &M, const vector <double> mean) const;
		double integrated_matrix_prod(const vector <double> &mu1, const vector < vector <double> > &Minv1, const vector <double> &mu2, const vector < vector <double> > &Minv2) const;
		void integrated_matrix_prod_check() const;
		vector <double> set_mu(const unsigned int sett) const;
		vector <double> set_p(const unsigned int sett) const;
		void put_p(const unsigned int sett, const vector <double> &p);
		void initialise_approx();
		// END //
		
		// START These functions are used for MBPs //
		void set_Imap_sett(const unsigned int sett);
		void set_Imap_using_dI(const unsigned int sett, const State *state, const vector< vector <double> > &dImap, const vector < vector <double> > &dIdiag);
		void update_I_from_transnum(vector < vector <double> > &Ima, vector < vector <double> > &Idia, const vector < vector < vector <double> > > &dtransnum) const;
		void set_I_from_pop(const unsigned int sett, const bool check);
		void update_pop(const unsigned int sett);
		void pop_positive(const unsigned int sett);
		void update_pop_rev(const unsigned int sett);
		void pop_init();
		double likelihood() const;
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