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

	Mode mode() const;
	
private:
	void set_command_line_params(int argc, char *argv[]);  
	void read_toml_file(bool verbose); 
	vector<string> get_toml_keys( const toml::basic_value<::toml::discard_comments, std::unordered_map, std::vector> &data) const;
	void check_for_undefined_parameters(vector<string> allowed, vector<string> given,	const string &context) const;
	
	map<string,string> cmdlineparams;                                                   // A map of all the parameters entered on the command line
	toml::basic_value<toml::discard_comments, std::unordered_map, std::vector> tomldata;// Stores information from the TOML file

	vector<string> definedparams;                                                       // Stores all the possible parameters 
};

#endif
