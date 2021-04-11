#ifndef BEEPMBP__INPUTS_HH
#define BEEPMBP__INPUTS_HH

#include <string>
#include <map>
#include <vector>

using namespace std;

#include "struct.hh"
#include "reader.hh"

class Inputs
{
	public:
		Inputs(int argc, char** argv);
		~Inputs();
		Inputs(const Inputs&) = delete;
		Inputs(Inputs&&) = delete;
		Inputs& operator=(const Inputs&) = delete;
		Inputs& operator=(Inputs&&) = delete;
		
		int find_integer(const string &key, int def) const;
		double find_double(const string &key, double def) const;
		string find_string(const string &key, const string &def) const;	
		vector <DataTable> find_datatable(const Details &details) const;
		vector <DemographicCategory> find_demographic_category() const;
		vector <Covariate> find_covariate() const;
		void find_genQ(GenerateQ &genQ, const Details &details, const vector <string> &age_list) const;
		void timeplot_add(vector <TimePlot> &timeplot, string time_str, string name) const;
		void find_timeplot(vector <TimePlot> &timeplot) const;
		void find_compartments(vector <string> &name, vector <ParamSpec> &inf, vector <string> &inf_trans) const;
	
		void find_transitions(vector <string> &from, vector <string> &to, vector <int> &type, vector < vector <ParamSpec> > &distspec, vector < vector <ParamSpec> > &probspec,  vector <unsigned int> &k, unsigned int ndemocatpos) const;
		void find_spline(const string name, vector < vector <ParamSpec>>& ps_vec, vector <ParamSpec>& factor_param_vec, vector < vector <unsigned int> >& bp_vec, vector <Smooth>& smooth_type_vec, vector <string>& inft, vector <string>& area, vector < vector <double> >& efoi_agedist_vec, const vector <double>& agedist, const Details &details) const;
		void find_spline(const string name, vector <ParamSpec>& ps, vector <ParamSpec>& factor_param_vec, vector <unsigned int>& bp, Smooth& smooth_type, const Details &details) const;
		void find_probreach(vector <string>& name, vector <string>& comp, vector <string>& inf_trans) const;
		void regional_effects(string& geography, ParamSpec& sigma) const;
		void find_generation(unsigned int &G) const;
		void find_generation_or_invT_final(unsigned int &G, double &invT_final) const;
		void find_invT_start_invT_final(double &invT_start, double &invT_final) const;
		void find_nparticle(unsigned int &Ntot, unsigned int &N, unsigned int ncore) const;
		void find_nchain(unsigned int &Ntot, unsigned int &N, unsigned int ncore) const;
		void find_nsample(unsigned int &nsample) const;
		void find_nthin(unsigned int &nthin, unsigned int nsample) const
		void find_nburnin(unsigned int &nburnin, const unsigned int &nsample) const;
		void find_nquench(unsigned int &nquench, const unsigned int &nburnin) const;
		void find_cutoff_frac(double &cutoff_frac) const;
		void find_cutoff(double &cutoff, double &cutofffrac) const;
		void find_quench_factor(double &quench_factor) const;
		void find_ngroup(unsigned int &ngroup, unsigned int &groupN, unsigned int Ntot) const;
		void find_propsize(double &propsize) const;
		void find_outputprop(OutputProp &prop) const;
		void find_counterfact(const Details &details, vector <CounterFact> &counterfact) const;
		void find_ppc_start(const Details &details, vector <CounterFact> &counterfact) const;

		Mode mode() const;
		SimInf get_siminf() const;
		
	private:
		void set_command_line_params(int argc, char *argv[]);  
		void check_for_undefined_parameters(vector<string> allowed, vector<string> given,	const string &context) const;
		vector <ParamSpec> get_paramspec(string root, unsigned int j, string name, string value, string prior, string smooth, string factor, unsigned int si) const;
		void print_paramspec(const vector <ParamSpec> &ps) const;
		vector <unsigned int> get_breakpoints(string st, const Details &details) const;
		void check_from_table(string& s) const;
		vector<string> split_string(string& s, char delimiter, unsigned int len) const;
		vector<double> split_number(string& s, char delimiter, unsigned int len) const;
		
		map<string,string> cmdlineparams;        // A map of all the parameters entered on the command line
		InputData *basedata;
};

const vector<string> definedparams = {       // A list of all supported parameters (please keep in lexicographic order)
		"baseinputfile",
		"age_mixing_matrix",
		"age_mixing_modify",
		"ages", //*
		"areas", //*
		"comps",//*
		"covars",//
		"counterfactual",//
		"cutoff",
		"cutoff_frac",
		"datadir",//
		"data_tables",//
		"democats",//
		"efoi_factor",
		"efoi_spline",
		"end",//
		"geo_mixing_matrix",//
		"geo_mixing_modify",//
		"geography",//
		"infected_max",
		"invT_final",
		"inputfile",
		"mode",
		"nburnin",
		"nchain",
		"ngeneration",
		"nparticle",
		"nquench",
		"nsample",
		"nthin",
		"outputdir",
		"output_prop",
		"quench_factor",
		"propsize",
		"ppc_start",
		"prob_reach",
		"region_effect",
		"R_spline", // 
		"seed",
		"start",
		"state_outputs",
		"threshold",
		"timeformat",
		"times",
		"trans",
	};
	
#endif
