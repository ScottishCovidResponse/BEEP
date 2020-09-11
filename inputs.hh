#ifndef BEEPMBP__INPUTS_HH
#define BEEPMBP__INPUTS_HH

#include <string>
#include <map>
#include <vector>

using namespace std;

#include "consts.hh"

struct TransitionData;
struct PopulationData;
struct MarginalData;
struct DemographicCategory;
struct Covariate;
struct TimePeriod;
struct GenerateQ;
struct Qtensor;
struct PriorComp;
struct Compartment;
class Details;
class InputData;

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
		vector <TransitionData> find_transition_data(const Details &details) const;
		vector <PopulationData> find_population_data(const Details &details) const;
		vector <MarginalData> find_marginal_data(const Details &details, const vector <DemographicCategory> &democat) const;
		vector <DemographicCategory> find_demographic_category(const Details &details) const;
		vector <Covariate> find_covariate(const Details &details) const;
		vector <TimePeriod> find_time_period(const Details &details) const;
		void find_genQ(GenerateQ &genQ, const Details &details) const;
		void find_Q(vector <Qtensor> &Qvec, const vector <TimePeriod> &time_period, const Details &details) const;
		void find_parameter(vector <string> &name, vector <double> &val) const;
		void find_prior(vector <string> &name, vector <double> &min, vector <double> &max) const;
		void find_compartments(vector <string> &name, vector <double> &infectivty) const;
		void find_transitions(vector <string> &from, vector <string> &to, vector <string> &prpar,
										vector <int> &type, vector <string> &mean, vector <string> &cv) const;
		vector <PriorComp> find_priorcomps(const vector<Compartment> &comp) const;
		void find_spline(const Details &details, string &name, vector <int> &time, vector <string> &param) const;

		Mode mode() const;
		
	private:
		void set_command_line_params(int argc, char *argv[]);  
		void check_for_undefined_parameters(vector<string> allowed, vector<string> given,	const string &context) const;
		
		map<string,string> cmdlineparams;                                                   // A map of all the parameters entered on the command line
		InputData *basedata;
};

const vector<string> definedparams = {       // A list of all supported parameters (please keep in lexicographic order)
		"baseinputfile",
		"Q",
		"agemix",
		"ages",
		"areas",
		"betaspline",
		"burnin",
		"comps",
		"covars",
		"cputime",
		"datadir",
		"democats",
		"distribution",
		"end",
		"genQ",
		"genQoutput",
		"geomix",
		"infmax",
		"invTmin",
		"invTmax",
		"dirs",
		"inputfile",
		"margdata",
		"mode",
		"model",
		"nchain",
		"ngeneration",
		"nparticle",
		"nsamp",
		"output",
		"outputdir",
		"params",
		"phispline",
		"popdata",
		"priorcomps",
		"priors",
		"propsmethod",
		"regions",
		"seed",
		"start",
		"threshold",
		"timeformat",
		"timep",
		"timeunits",
		"trans",
		"transdata",
	};
	
#endif
