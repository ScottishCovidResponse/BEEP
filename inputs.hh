#ifndef BEEPMBP__INPUTS_HH
#define BEEPMBP__INPUTS_HH

#include <string>
using namespace std;

#include "consts.hh"
#include "toml11/toml.hpp"

class Inputs
{
public:
	Inputs(int argc, char** argv, bool verbose);
	
	Mode mode(bool verbose) const;
	
private:
	map<string,string> get_command_line_params(int argc, char *argv[]);   // Gets the command line parameters
	void read_toml_file(bool verbose);                                                  // Reads the TOML file                    
	string lookup_string_parameter(const map<string,string> &params,
			  												 const toml::basic_value<::toml::discard_comments, std::unordered_map, std::vector> &tomldata,
				  											 const string &key, bool verbose, const string &def="") const;

	vector<string> get_toml_keys( const toml::basic_value<::toml::discard_comments, std::unordered_map, std::vector> &data) const;
	void check_for_undefined_parameters(vector<string> allowed, vector<string> given,	const string &context) const;

	
	map<string,string> cmdlineparams;                                                   // A map of all the parameters entered on the command line

	toml::basic_value<toml::discard_comments, std::unordered_map, std::vector> tomldata;// Stores information from the TOML file

	vector<string> definedparams;                                                       // Stores all the possible parameters 
	
	
};

#endif
