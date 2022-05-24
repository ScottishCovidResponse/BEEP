/// This class takes a particle, performs a sequence of MBPs and returns a particle

#ifndef BEEP__MBP_HH
#define BEEP__MBP_HH

using namespace std;

#include "struct.hh"
#include "param_prop.hh"
#include "state.hh"
#include "output.hh"

class ObservationModel;

class Mbp                                                        
{
	public:
		Mbp(ObsModelMode obsmodel_mode_, const Details &details, const Data &data, const Model &model, const ObservationModel &obsmodel, const Output &output, Mpi &mpi);
		
		unsigned int mcmc_updates(vector <Particle> &part, const vector <ParamSample> &param_samp, double EFcut_, double invT_,ParamUpdate pup_, ParamProp &paramprop, const double cor_max);
		unsigned int mcmc_updates_new(vector <Particle> &part, const vector <ParamSample> &param_samp, double EFcut_, double invT_, ParamUpdate pup_, ParamProp &paramprop);
		unsigned int mc3_mcmc_updates(Particle &part, const vector <ParamSample> &param_samp, const double invT_, const ParamUpdate pup_, ParamProp &paramprop);
	
	private:
		void update_particle(Particle &pa, ParamProp &paramprop, unsigned int &nprop);
		//unsigned int mvn_update(vector <Particle> &part, const vector <ParamSample> &param_samp, ParamProp &paramprop, const double cor_max);
		//unsigned int univariate_update(vector <Particle> &part, const vector <ParamSample> &param_samp, ParamProp &paramprop, const double cor_max);
		//unsigned int complex_update(vector <Particle> &part, const vector <ParamSample> &param_samp, ParamProp &paramprop, const double cor_max);
		
		Status mbp(const vector<double> &paramv, const InfUpdate inf_update);	
		Status mbp_rev(const vector<double> &paramv);
		double get_al();
		void initialise_variables();
		void simu_or_mbp_reset();
		void swap_initial_propose_state();
		void mbp_initialise();
		//void update_particle(Particle &pa, const vector <Proposal> &prop_list, ParamProp &paramprop);
		
		void mvn_proposal(MVN &mvn);
		void parameter_proposal(MVN &mvn, double &Li);
		void mvn_proposal_rev(MVN &mvn);
		void dist_R_joint_proposal(MVN &mvn);
		void sigma_reff_proposal(MVN &mvn);
		void mbp_fixedtree(FixedTree &ft);
		void mbp_slicetime(SliceTime &st);
		void mean_time_proposal(MeanTime &mt);
		void neighbour_proposal(Neighbour &rn);
		void joint_proposal(Joint &rn);		
		void covar_area_proposal(CovarArea &ca);
		void copy_initial_to_propse();
	
		void check_MBP();
		
		ObsModelMode obsmodel_mode;                                     // Determines if using invT or EFcut
		
		double invT;                                                    // The inverse temperature (used for MC3 and PAIS)
		double EFcut;                                                   // The error function cutoff (used for ABC methods)
		
		ParamUpdate pup;                                                // Determines how parameters are updated
		
		State state1, state2;                                           // Stores states and swaps references
		State *initial, *propose;                                       // The states in the initial and proposed states

		vector < vector <Simu_or_mbp> > simu_or_mbp;                    // Stores whether doing a simulation or a MBP
				
		vector < vector <double> > dImap;                               // The difference in Imap between the two states
		vector < vector <double> > dIdiag;                          		// The difference in Idiag between the two states
		
		vector < vector < vector <double> > > dtransnum;                // The difference in transnum between state (POP_MODEL)
		
		vector < vector <Simu_or_mbp> > mbp_sim;                        // Used in fixedtree to determine which areas are simulated 
		
		const vector <Compartment> &comp;
		const vector <Transition> &trans;
		
		const Details &details;
		const Data &data;
		const Model &model;
		const ObservationModel &obsmodel;
		const Output &output;
		Mpi &mpi;
};
#endif
