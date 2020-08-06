#ifndef BEEPMBP__INPUTS_HH
#define BEEPMBP__INPUTS_HH

#include <string>
#include <map>
#include <vector>

using namespace std;

#include "consts.hh"
#include "toml11/toml.hpp"

struct TRANSDATA;
struct POPDATA;
struct MARGDATA;
struct DEMOCAT;
struct COVAR;
struct TIMEP;
struct GENQ;
struct QTENSOR;
struct PRIORCOMP;
struct COMP;
class Details;

class Inputs
{
public:
	Inputs(int argc, char** argv, bool verbose);

	int find_int(const string &key, int def) const;
	double find_double(const string &key, double def) const;
	string find_string(const string &key, const string &def) const;	
	vector <TRANSDATA> find_transdata(const Details &details) const;
	vector <POPDATA> find_popdata(const Details &details) const;
	vector <MARGDATA> find_margdata(const Details &details, const vector <DEMOCAT> &democat) const;
	vector <DEMOCAT> find_democat(const Details &details) const;
	vector <COVAR> find_covar(const Details &details) const;
	vector <TIMEP> find_timeperiod(const Details &details) const;
	void find_genQ(GENQ &genQ, const Details &details) const;
	void find_Q(vector <QTENSOR> &Qvec, const vector <TIMEP> &timeperiod, const Details &details) const;
	void find_param(vector <string> &name, vector <double> &val) const;
	void find_prior(vector <string> &name, vector <double> &min, vector <double> &max) const;
	void find_comps(vector <string> &name, vector <double> &infectivty) const;
	void find_trans(vector <string> &from, vector <string> &to, vector <string> &prpar,
                 	vector <int> &type, vector <string> &mean, vector <string> &cv) const;
	vector <PRIORCOMP> find_priorcomps(const vector<COMP> &comp) const;
	void find_spline(const Details &details, string &name, vector <int> &time, vector <string> &param) const;

	Mode mode() const;
	
private:
	void set_command_line_params(int argc, char *argv[]);  
	vector<string> get_toml_keys( ) const;
	void check_for_undefined_parameters(vector<string> allowed, vector<string> given,	const string &context) const;
	
	map<string,string> cmdlineparams;                                                   // A map of all the parameters entered on the command line
	toml::basic_value<toml::discard_comments, std::unordered_map, std::vector> tomldata;// Information from the TOML file
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
