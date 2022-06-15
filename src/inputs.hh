#ifndef BEEP__INPUTS_HH
#define BEEP__INPUTS_HH

#include <string>
#include <map>
#include <vector>

using namespace std;

#include "utils.hh"
#include "reader.hh"

#include "fdp/fdp.hxx"

class Inputs
{
	public:
		Inputs(int argc, char** argv);
		~Inputs();
		Inputs(const Inputs&) = delete;
		Inputs(Inputs&&) = delete;
		Inputs& operator=(const Inputs&) = delete;
		Inputs& operator=(Inputs&&) = delete;
		
		unsigned find_positive_integer(const string &key, const int def);
		int find_integer(const string &key, const int def);
		double find_double(const string &key, const double def);
		string find_string(const string &key, const string &def, bool allowcomma=false);	
		unsigned int find_maximum_timestep();
		vector <DataTable> find_datatable(const Details &details);
		LineColour get_line_colour(const string col) const;
		vector <DemocatChange> find_democat_change();
		vector <DemographicCategory> find_demographic_category(vector <Strain> &strain_list);
		LevelEffect find_level_effect();
		vector <Covariate> find_covariates();
		vector <unsigned int> find_agecats(const string ages_str, const vector <string> &age_list, const string name) const;
		void find_genQ(GenerateQ &genQ, const Details &details, const vector <string> &age_list);
		void find_timeplot(vector <TimePlot> &timeplot);
		//unsigned int find_susceptible_compartment();
		SmoothType find_smoothtype(const string st, const string name) const;
		void find_spline(const string name, vector <string> &name_vec, vector < vector <ParamSpec> > &ps_vec, vector <ParamSpec> &factor_param_vec, vector < vector <unsigned int> > &bp_vec, vector <SmoothType> &smooth_type_vec, vector <string> &inft, vector <string> &area, vector < vector <double> > &efoi_agedist_vec, const vector <double> &agedist, const Details &details);
		void find_spline(const string name, vector <ParamSpec> &ps, vector <ParamSpec> &factor_param_vec, vector <unsigned int> &bp, SmoothType &smooth_type, const Details &details);
		void find_probreach(vector <string> &name, vector <string> &comp);
		void find_region_effect(vector <ParamSpec> &ps_vec, ParamSpec &sigma, const vector <Area> &area);
		bool param_not_set(const vector <ParamSpec> &ps_vec) const;
		AreaEffect find_area_effect(const vector <Area> &area);
		void find_generation(unsigned int &G);
		void find_cpu_time(double &cpu_time);
		void find_generation_or_invT_final_or_cpu_time(unsigned int &G, double &invT_final, double &cpu_time);
		void find_generation_or_cutoff_final_or_cpu_time(unsigned int &G, double &cutoff_final, double &cpu_time);
		void find_invT_start_invT_final(double &invT_start, double &invT_final);
		void find_Tpower(double &Tpower);
		void find_invT(double &invT);
		void find_sd(double &sd);
		void find_nrun(unsigned int &nrun);
		void find_algorithm(MLAlg &algorithm, unsigned int &npart, unsigned int &G, double &cpu_time, unsigned int &P, unsigned int &nsample, const unsigned int core, const unsigned int ncore, const unsigned int nvar);
		void find_stateuncer(StateUncertainty &stateuncer);
		void find_nparticle(unsigned int &npart, unsigned int &Ntot, unsigned int &N, const unsigned int nrun, const unsigned int ncore);
		void find_cor_max(double &cor_max, double &cor_max_last);
		void find_prop_max(unsigned &prop_max);
		void find_nparticle_pmcmc(unsigned int &npart, unsigned int &N, const unsigned int ncore);
		void find_nchain(unsigned int &nchain, unsigned int &Ntot, unsigned int &N, const unsigned int nrun, const unsigned int ncore);
		void find_nsample(unsigned int &nsample, unsigned int min);
		void find_nsample_or_ESSmin_or_cpu_time(unsigned int &nsample, double &ESSmin, double &cpu_time);
		void find_nsample_final(unsigned int &nsample_final,unsigned int &nsample);
		void find_nthin(unsigned int &nthin, const unsigned int nsample);
		void find_nburnin(unsigned int &nburnin, const unsigned int &nsample);
		void find_nquench(unsigned int &nquench, const unsigned int &nburnin);
		void find_cutoff_frac(double &cutoff_frac, double &cutoff_frac_init);
		void find_cutoff(double &cutoff, double &cutofffrac);
		void find_quench_factor(double &quench_factor);
		void find_propsize(double &propsize);
		void find_outputprop(OutputProp &prop);
		void find_plot_param_values(bool &plot_param_values);
		void find_modification(const Details &details, vector <Modification> &modification);
		void find_mcmc_update(MCMCUpdate &mcmc_update);
		void find_area_plot(AreaPlot &area_plot);
		vector <vector <unsigned int > > find_demo_ref(const string dep_str, const vector <DemographicCategory> &democat, vector <unsigned int> &dep_order) const;
		Mode mode();
		
		SimInf get_siminf();
		void print_commands_not_used() const;
		
		string inputfilename;                                                 // The name of the TOML file

		static std::string getEnvVar( std::string const & key );

		bool use_datapipeline();
		FairDataPipeline::DataPipeline::sptr datapipeline = nullptr;
		
	private:
		void set_command_line_params(int argc, char *argv[]);  
		void check_for_undefined_parameters(vector<string> allowed, vector<string> given,	const string &context) const;
		
		vector <ParamSpec> find_ct_param_spec(InputNode &it, const string root, const string name, const string value, const string prior, const string dep_str, const vector < vector<unsigned int> > &democatpos, const vector <DemographicCategory> &democat, const string par_root);
		void get_dep(const string root, string &st, string &dep) const;
		string get_dep(const string root, string &name, string &value, string &prior) const;
		string get_dep(InputNode &it, const string root, const string name, const string value, const string prior) const;
		vector <ParamSpec> get_paramspec(InputNode &it, const string root, const string name, const string value, const string prior, const string smooth, const string factor, const unsigned int si, const bool dep_expect);		
		void print_paramspec(const vector <ParamSpec> &ps) const;
		vector <unsigned int> get_breakpoints(string st, const Details &details, const string em);
		void extend_spline_for_prediction(vector <unsigned int> &bp, vector <ParamSpec> &ps_list, const Details &details);
		void check_from_table(string &s);
		vector <string> split_string(string &s, const char delimiter, const unsigned int len);
		vector <double> split_number(string &s, const char delimiter, const unsigned int len);
		bool get_bool(const string st, const string em) const;
		void addused(const string key, vector <UsedTomlKey> &used);
		void check_filename(const string &str) const;

	// Used in inputs_comps
	public:
		void load_compartmental_model(vector <Compartment> &comp, vector <Transition> &trans, vector <CompartmentName> &comp_name, const vector <DemographicCategory> &democat, const vector < vector<unsigned int> > &democatpos, const Mode mode);
		vector <string> find_compartment_names(unsigned int &co_sus, const vector <DemographicCategory> &democat, const vector < vector<unsigned int> > &democatpos);
	
	private:
		void add_comps(vector <Compartment> &comp, const vector <DemographicCategory> &democat, const vector < vector<unsigned int> > &democatpos);
		void add_compartment(vector <Compartment> &comp, const string& name, const vector <ParamSpec> &mean_spec, const string &mean_dep, const ParamSpec &inf, const unsigned int num, const unsigned int shape);
		void add_trans(vector <Transition> &trans, vector <Compartment> &comp, const vector <DemographicCategory> &democat, const vector < vector<unsigned int> > &democatpos);
		void add_internal_transition(vector <Transition> &trans, vector <Compartment> &comp, const unsigned int c_from, const unsigned int c_to);
		void add_transition(vector <Transition> &trans, vector <Compartment> &comp, const string &from, const string &to, const TransInf inf, const vector <ParamSpec> &prob_spec, const string prob_dep);
		void split_merged(vector <Compartment> &comp, vector <Transition> &trans);
		void add_comp_name(vector <Compartment> &comp, vector <CompartmentName> &comp_name);
		void find_compartments(vector <string> &name, vector < vector <ParamSpec> > &mean_spec, vector <string> &mean_dep, vector <unsigned int> &k, vector <ParamSpec> &infectivity, const vector < vector<unsigned int> > &democatpos, const vector <DemographicCategory> &democat);
		void find_transitions(vector <string> &from, vector <string> &to, vector <TransInf> &transinf, vector < vector <ParamSpec> > &prob_spec, vector <string> &prob_dep, const vector < vector<unsigned int> > &democatpos, const vector <DemographicCategory> &democat);
		unsigned int get_compartment(const vector <Compartment> &comp, const string compname, const ErlangPos pos) const;
	
		vector <UsedTomlKey> used;                                        // Stores which keys have been used
	
		map<string,string> cmdlineparams;                                 // A map of all the parameters entered on the command line
		InputData *basedata;
		bool _use_pipeline = false;
}; 

const vector<string> definedparams = {                                // A list of all supported commands
		"baseinputfile",
		"algorithm",//
		"area_covars",//
		"area_effect",//
		"area_plot",
		"area_tv_covars",//
		"age_mixing_matrix",
		"age_mixing_modify",
		"age_mixing_perturb",
		"ages", 
		"areas",
		"cor_max",
		"cor_max_last",
		"comps",
		"cpu_time",
		"cutoff",
		"cutoff_final",
		"cutoff_frac",
		"cutoff_frac_init",
		"datadir", 
		"datapipeline_output_data_product",
		"datapipeline_output_type",
		"data_tables",
		"democats",
		"democat_change",
		"description",
		"dynamics",
		"efoi_factor",
		"efoi_spline",
		"end",
		"ESS_min",
		"geo_mixing_matrix",
		"geo_mixing_modify",
		"init_pop",
		"invT",
		"invT_start",
		"invT_final",
		"invT_power",
		"inputfile",
		"level_effect",
		"mcmc_update",
		"mode",
		"modification",
		"nburnin",
		"nchain",
		"ngeneration",
		"nodata_str",
		"nparticle",
		"nquench",
		"nrun",
		"nsample",
		"nsample_final",
		"nsimulation",
		"nsim_per_sample",
		"nthin",
		"obs_spline",
		"outputdir",
		"output_prop",
		"quench_factor",
		"plot_param_values",
		"posterior_particle",
		"propsize",
		"prediction_end",
		"prediction_start",
		"prior_order",
		"prob_reach",
		"prop_size",
		"prop_max",
		"region_effect",
		"R_spline", 
		"sd",
		"seed",
		"start",
		"state_outputs",
		"state_uncertainty",
		"steps_per_unit_time",
		"strains",
		"threshold_str",
		"trans_combine",
		"time_format",
		"time_labels",
		"trans",//
		"tv_covars",//
		"with_datapipeline"
	};
	
#endif
