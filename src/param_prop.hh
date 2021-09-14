#ifndef BEEPMBP__PARAMPROP_HH
#define BEEPMBP__PARAMPROP_HH

using namespace std;

#include "struct.hh"

struct FixedTree {        // This stores information about proposals that simulate a given set of areas
	unsigned int n;         // The node on treenode (determines which areas to do MBP and which to simulate)
	double sim_frac;        // The fraction of simulation (vs MBP)
	unsigned int ntr;       // Acceptance rate
	unsigned int nac;
	double ac_rate;
};

struct SliceTime {        // This stores information about proposals that simulate a given time period
	unsigned int sett_i;    // The starting simulation time
	unsigned int sett_f;    // The ending simulation time
	double sim_frac;        // The fraction of simulation (vs MBP)
	unsigned int ntr;       // Acceptance rate
	unsigned int nac;
	double ac_rate;
};

struct Self               // These proposals do not change parameters but resample the state (PMCMC)
{
	unsigned int nbo;
	unsigned int ntr;
	unsigned int nac;
	double bo_rate;
	double ac_rate;
	
	Status MH(double al, double &invT, ParamUpdate pup);
};

struct MeanTime           // These change the time for one transition and do the opposite change for the subsequent transitions
{
	vector <unsigned int> param_mean;
	vector <unsigned int> param_mean_rev;
	
	double size;
	unsigned int nbo;
	unsigned int ntr;
	unsigned int nac;
	double bo_rate;
	double ac_rate;
	
	Status MH(double al, ParamUpdate pup);
	Status propose(vector <double> &param_prop, const vector <double> &paramval, const Model &model);
};

struct Neighbour         // These change a point on a spline with the neighbouring point having the opposite change
{
	unsigned int param1;
	unsigned int param2;

	double size;
	unsigned int nbo;
	unsigned int ntr;
	unsigned int nac;
	double bo_rate;
	double ac_rate;
	
	Status MH(double al, ParamUpdate pup);
	Status propose(vector <double> &param_prop, const vector <double> &paramval, const Model &model);
};


struct Joint            // These make joint changes to all the points on a spline (sinusoidal or all up/down)
{
	vector <unsigned int> var_list;
	JointType type;
	unsigned int sinenum;

	double size;
	unsigned int nbo;
	unsigned int ntr;
	unsigned int nac;
	double bo_rate;
	double ac_rate;
	
	Status MH(double al, ParamUpdate pup);
	Status propose(vector <double> &param_prop, const vector <double> &paramval, const Model &model);
};

	
struct CovarArea        // These make joint proposals on covariates and area effects
{
	unsigned int covar_ref;
	
	double size;
	unsigned int ntr;
	unsigned int nac;
	double ac_rate;
	
	Status MH(double al, ParamUpdate pup);
	Status propose(vector <double> &param_prop, const vector <double> &paramval, const Model &model, const Data &data);
};
	
	
class ParamProp                                  // Information about kernals for propsals in parameter space
{
public:
	ParamProp(const Details &details, const Data &data, const Model &model, const Output &output, Mpi &mpi);

	Self self;                                     // Unchanged parameter proposal (PMCMC)
		
	vector <MVN> mvn;                              // Multivariate proposal distribution
	
	vector <MeanTime> mean_time;                   // Proposals which change compartmental mean transition rates 

	vector <Neighbour> neighbour;                  // Proposals which change neighbouring point on a spline in the opposite way

	vector <Joint> joint;                          // Proposals which change multiple variables

	vector <CovarArea> covar_area;                 // Proposals which simultaneously change covariate and area effects
	
	vector <FixedTree> fixedtree;                  // Proposals which simulate areas distinct areas

	vector <SliceTime> slicetime;                  // Proposals which simulate a time period

	void init(const vector <double> &paramv, const Data &data, const Model &model, const Details &details, unsigned int ncovar, unsigned int nmpp);
	void zero_ntr_nac();
	vector <Proposal> get_proposal_list(const vector <ParamSample> &param_samp);
	void print_prop_list(const vector <Proposal> &prop_list) const;
	void update_proposals();
	void update_sim_proposals();
	void update_fixedtree();
	void update_splicetime();
	void diagnostics() const;
	void set_ac_rate();
	
	vector <double> variance_vector(const vector <ParamSample> &param_samp) const;
	string print_proposal_information(const bool brief) const;
	
private:
	void update(double &val, double acrate);
	void update_high(double &val, double acrate);
	void add_mvn(const string name, const ParamType type, const double size);
	void add_single();
	void add_demographic_specific();
	void randomise(vector <Proposal> &vlist);
	void mean_time_init();
	void neighbour_init();
	void joint_init();
	void covar_area_init();
		
	const Details &details;
	const Data &data;
	const Model &model;
	const Output &output;
	Mpi &mpi;
};

#endif
