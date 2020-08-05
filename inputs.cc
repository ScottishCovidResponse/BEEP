
#include <iostream>
#include <map>
#include <vector>
#include <string>
#include <sstream>

using namespace std;

#include "inputs.hh"
#include "utils.hh"
#include "consts.hh"
#include "toml11/toml.hpp"

Inputs::Inputs(int argc, char** argv, bool verbose) 
{
	set_command_line_params(argc,argv);                          // Loads up the command line parameters
	
	// A list of all supported parameters (please keep in lexicographic order)
	definedparams = {
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
	
	read_toml_file(verbose);
}

/// Gets the command line parameters and places them into cmdlineparams
void Inputs::set_command_line_params(int argc, char *argv[])
{
	vector <string> commandlist;
	
	for(int op = 1; op < argc; op++){ 
		string str = string(argv[op]);
		unsigned int n = commandlist.size();
		if(n > 0 && (str.substr(0,1) == "=" || commandlist[n-1].substr(commandlist[n-1].length()-1,1) == "=")) commandlist[n-1] += str;
		else{
			commandlist.push_back(str);
		}
	}
	
	// Store the parameters passed on the command line in cmdlineparams
	for(unsigned int op = 0; op < commandlist.size(); op++){ // Goes the various input options
		string str = commandlist[op];
		int j = 0; int jmax = str.length(); while(j < jmax && str.substr(j,1) != "=") j++;
		if(j == jmax){
			stringstream ss; ss << "Cannot understand " << str; 
			emsgroot(ss.str());
		}
		
		string command = str.substr(0,j);
		string value = str.substr(j+1,jmax-(j+1));
		
		if (cmdlineparams.count(command) == 0) {
			cmdlineparams[command] = value;
		} else {
			// Encode repeated parameters as space-separatedd
			cmdlineparams[command] += " " + value;
		}
	}
}

/// Reads TOML parameters and places them into tomldata
void Inputs::read_toml_file(bool verbose)
{
	string inputfilename = "/dev/null";
	if (cmdlineparams.count("inputfile") == 1) {
		inputfilename = cmdlineparams["inputfile"];
	}
	//decltype(toml::parse(inputfilename)) tomldata;
	try {
		tomldata = toml::parse(inputfilename);

		// Allow using values from another TOML file as a base for this one. TODO:
		// make this into functions so you can do this recursively.
		if (tomldata.contains("baseinputfile")) {
			const string basefile = toml::find<string>(tomldata,"baseinputfile");
			// TODO: make the filename relative to the original TOML file
			decltype(toml::parse(basefile)) basetomlddata = toml::parse(basefile);

			for(const auto& p : basetomlddata.as_table())
			{
				if (!tomldata.contains(p.first)) {
					tomldata[p.first] = p.second;
				}
			}
		}

	} catch (const std::exception& e) {
		std::ostringstream oss;
		oss << "toml::parse returns exception\n" << e.what();
		emsg(oss.str());
	}	
	
	vector<string> tomlkeys = get_toml_keys(tomldata);
	check_for_undefined_parameters(definedparams, tomlkeys, "in " + inputfilename);
}

string Inputs::lookup_string_parameter(const map<string,string> &params,
															 const toml::basic_value<::toml::discard_comments, std::unordered_map, std::vector> &tomldata,
															 const string &key, bool verbose, const string &def) const
{
	string val;
	auto val_it = params.find(key);
	if(val_it != params.end()) val = val_it->second;
	else{
		if(tomldata.contains(key)) val = toml::find<string>(tomldata,key);
		else val = def;
	}

	if(verbose) cout << "  " << key << " = " << val << endl;
	
	return val;
}

int Inputs::find_int(const string &key, int def) const
{
	int val;
	auto val_it = cmdlineparams.find(key);
	if (val_it != cmdlineparams.end()) {
		const std::string& valstr = val_it->second;
		try {
			size_t idx;
			val = stoi(valstr,&idx);
			if (idx != valstr.length()) {
				std::ostringstream oss;
				oss << "Should be integer, found '"<< valstr;
				throw std::invalid_argument(oss.str());
			}
		} catch (const std::exception& e) {
			std::ostringstream oss;
			oss << "Bad command-line parameter value for key '"<< key <<"'\n";
			// Add exception description if it's informative
			std::string what = e.what();
			if (what == "stoi") {
				if (valstr == "")
					oss << "Should be integer, found no value\n";
			} else {
				oss << what;
			}
			emsg(oss.str());
		}
	} else {
		if (tomldata.contains(key)) {
			val = toml::find<int>(tomldata,key);
		} else {
			val = def;
		}
	}
	
	return val;
}

double Inputs::find_double(const string &key, double def) const
{
	double val;
	auto val_it = cmdlineparams.find(key);
	if (val_it != cmdlineparams.end()) {
		const std::string& valstr = val_it->second;
		try {
			size_t idx;
			val = stof(valstr,&idx);
			if (idx != valstr.length()) {
				std::ostringstream oss;
				oss << "Should be number, found '"<< valstr;
				throw std::invalid_argument(oss.str());
			}
		} catch (const std::exception& e) {
			std::ostringstream oss;
			oss << "Bad command-line parameter value for key '"<< key <<"'\n";
			// Add exception description if it's informative
			std::string what = e.what();
			if (what == "stoi") {
				if (valstr == "")
					oss << "Should be number, found no value\n";
			} else {
				oss << what;
			}
			emsg(oss.str());
		}
	} else {
		if (tomldata.contains(key)) {
			const auto val_temp = toml::find(tomldata,key);
			if(val_temp.is_floating()) val = val_temp.as_floating(); else val = val_temp.as_integer();	
		} else {
			val = def;
		}
	}
	
	return val;
}

vector<string> Inputs::get_toml_keys( const toml::basic_value<::toml::discard_comments, std::unordered_map, std::vector> &data) const
{
	vector<string> keys;
	for(const auto& p : data.as_table())
	{
		keys.push_back(p.first);
	}
	return keys;
}

void Inputs::check_for_undefined_parameters(vector<string> allowed, vector<string> given,	const string &context) const
{
	vector<string> undefined;

	sort(allowed.begin(), allowed.end());
	sort(given.begin(), given.end());

	set_difference(given.begin(), given.end(),
								 allowed.begin(), allowed.end(),
								 inserter(undefined, undefined.begin()));

	if (undefined.size() != 0) {
		stringstream ss;
		ss << "Unrecognised parameter(s) "+context+":";

		for (const auto &k : undefined) {
			ss << " " << k;
		}
		
		emsgroot(ss.str());
	}
}


Mode Inputs::mode(bool verbose) const
{
	string val = lookup_string_parameter(cmdlineparams, tomldata, "mode", verbose, "UNSET");  
	
	Mode mode;

	map<string,Mode>  modemap{{"sim", sim}, {"inf", inf}, {"multisim", multisim}, {"combinetace", combinetrace}};
	if (modemap.count(val) != 0) {
		mode = modemap[val];
	} else {
		emsgroot("Unrecoginsed value " + val + " for mode parameter");
	}
	return mode;
}
