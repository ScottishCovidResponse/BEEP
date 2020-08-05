
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
#include "data.hh"

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

string Inputs::find_string(const string &key, const string &def) const
{
	string val;
	auto val_it = cmdlineparams.find(key);
	if(val_it != cmdlineparams.end()) val = val_it->second;
	else{
		if(tomldata.contains(key)) val = toml::find<string>(tomldata,key);
		else val = def;
	}

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


Mode Inputs::mode() const
{
	string val = find_string("mode","UNSET");  
	if(val == "UNSET") emsgroot("The 'mode' property must be set");
	
	Mode mode;
	map<string,Mode>  modemap{{"sim", sim}, {"inf", inf}, {"multisim", multisim}, {"combinetace", combinetrace}};
	if (modemap.count(val) != 0) mode = modemap[val];
	else emsgroot("Unrecoginsed value " + val + " for mode parameter");
	
	return mode;
}

/// transdata
vector <TRANSDATA> Inputs::find_transdata(const Details &details) const
{
	vector <TRANSDATA> transdatavec;
	
	if(tomldata.contains("transdata")) {
		const auto tdata = toml::find(tomldata,"transdata");

		for(unsigned int j = 0; j < tdata.size(); j++){
			const auto td = toml::find(tdata,j);
		
			TRANSDATA transdata;
			if(!td.contains("from")) emsgroot("A 'from' property must be specified in 'transdata'.");
			const auto from = toml::find<std::string>(td,"from");
			transdata.fromstr = from;
		
			if(!td.contains("to")) emsgroot("A 'to' property must be specified in 'transdata'.");
			const auto to = toml::find<std::string>(td,"to");
			transdata.tostr = to;
			
			if(!td.contains("area")) emsgroot("An 'area' property must be specified in 'transdata'.");
			const auto area = toml::find<std::string>(td,"area");
			transdata.type = area;
			if(transdata.type != "reg" && transdata.type != "all") emsgroot("Transition data type not recognised"); 
			
			if(!td.contains("file")) emsgroot("A 'file' property must be specified in 'transdata'.");
			const auto file = toml::find<std::string>(td,"file");
			transdata.file = file;

			if(!td.contains("start")) emsgroot("A 'start' property must be specified in 'popdata'.");
			const auto startdata = toml::find<string>(td,"start");
			transdata.start = details.gettime(startdata)-details.start;
			
			if(!td.contains("units")) emsgroot("A 'units' property must be specified in 'transdata'.");
			const auto units = toml::find<std::string>(td,"units");
			if(units == "days") transdata.units = 1;
			else{
				if(units == "weeks") transdata.units = 7;
				else emsgroot("Units in 'transdata' not recognised");
			}
			
			if(details.mode != inf){
				transdata.rows = (unsigned int)((details.period - transdata.start)/transdata.units);
				if(transdata.rows == 0) emsgroot("Transition data '"+file+"' cannot be generated because the time period is not sufficiently long.");
			}
			
			transdatavec.push_back(transdata);
		}
	}
	return transdatavec;
}

/// popdata
vector <POPDATA> Inputs::find_popdata(const Details &details) const
{
	vector <POPDATA> popdatavec;

	if(tomldata.contains("popdata")) {
		const auto pdata = toml::find(tomldata,"popdata");

		POPDATA popdata;
		for(unsigned int j = 0; j < pdata.size(); j++){
			const auto pd = toml::find(pdata,j);
			
			if(!pd.contains("comp")) emsgroot("A 'comp' property must be specified in 'popdata'.");
			const auto comp = toml::find<std::string>(pd,"comp");
			popdata.compstr = comp;
			
			if(!pd.contains("area")) emsgroot("An 'area' property must be specified in 'popdata'.");
			const auto area = toml::find<std::string>(pd,"area");
			popdata.type = area;
			if(popdata.type != "reg" && popdata.type != "all") emsgroot("popition data type not recognised"); 
			
			if(!pd.contains("file")) emsgroot("A 'file' property must be specified in 'popdata'.");
			const auto file = toml::find<std::string>(pd,"file");
			popdata.file = file;

			if(!pd.contains("start")) emsgroot("A 'start' property must be specified in 'popdata'.");
			const auto startdata = toml::find<string>(pd,"start");
			popdata.start = details.gettime(startdata)-details.start;
			
			if(!pd.contains("units")) emsgroot("A 'units' property must be specified in 'popdata'.");
			const auto units = toml::find<std::string>(pd,"units");
			
			if(units == "days") popdata.units = 1;
			else{
				if(units == "weeks") popdata.units = 7;
				else emsgroot("Units in 'popdata' not recognised");
			}
			
			if(details.mode != inf){
				popdata.rows = (unsigned int)((details.period - popdata.start)/popdata.units);
				if(popdata.rows == 0) emsgroot("popition data '"+file+"' cannot be generated because the time period is not sufficiently long.");
			}
			
			popdatavec.push_back(popdata);
		}
	}
	
	return popdatavec;
}

/// margdata
vector <MARGDATA> Inputs::find_margdata(const Details &details, const vector <DEMOCAT> &democat) const
{
	vector <MARGDATA> margdatavec;
	
	if(tomldata.contains("margdata")) {
		const auto mdata = toml::find(tomldata,"margdata");

		for(unsigned int j = 0; j < mdata.size(); j++){
			const auto md = toml::find(mdata,j);
			
			MARGDATA margdata;
			if(!md.contains("from")) emsgroot("A 'from' property must be specified in 'margdata'.");
			const auto from = toml::find<std::string>(md,"from");
			margdata.fromstr = from;
		
			if(!md.contains("to")) emsgroot("A 'to' property must be specified in 'margdata'.");
			const auto to = toml::find<std::string>(md,"to");
			margdata.tostr = to;
			
			if(!md.contains("area")) emsgroot("An 'area' property must be specified in 'margdata'.");
			const auto area = toml::find<std::string>(md,"area");
			margdata.type = area;
			if(margdata.type != "reg" && margdata.type != "all") emsgroot("Marginal data type not recognised"); 
			
			if(!md.contains("type")) emsgroot("An 'type' property must be specified in 'margdata'.");
			const auto type = toml::find<std::string>(md,"type");
			unsigned int d;
			for(d = 0; d < democat.size(); d++) if(type == democat[d].name) break;
			if(d == democat.size()) emsg("The 'type' property must be 'age' or a demographic property.");
			margdata.democat = d;
				
			if(!md.contains("file")) emsgroot("A 'file' property must be specified in 'margdata'.");
			const auto file = toml::find<std::string>(md,"file");
			margdata.file = file;
			
			margdatavec.push_back(margdata);
		}
	}
	
	return margdatavec;
}

/// democats
vector <DEMOCAT> Inputs::find_democat(const Details &details) const
{
	vector <DEMOCAT> democatvec;
	
	if(tomldata.contains("ages")){                           // Age categories
		const auto ages = toml::find(tomldata,"ages");
		
		DEMOCAT democat;
		democat.name = "age";
		for(unsigned int j = 0; j < ages.size(); j++){
			const auto ag = toml::find(ages,j);
			
			if(!ag.contains("range")) emsgroot("A 'range' must be specified in 'ages'.");
			const auto range = toml::find<std::string>(ag,"range");
			democat.value.push_back(range);
			
			if(!ag.contains("sus")) emsgroot("A 'sus' must be specified in 'ages'.");
			const auto sus = toml::find<std::string>(ag,"sus");
			democat.param.push_back(sus);
		}
		democatvec.push_back(democat);
	}
	else emsgroot("The 'ages' parameter must be set.");
	
	if(tomldata.contains("democats")){                        // Other demographic possibilities
		const auto democats = toml::find(tomldata,"democats");
	
		for(unsigned int k = 0; k < democats.size(); k++){
			const auto democ = toml::find(democats,k);
			
			DEMOCAT democat;
			democat.name="";
			for(unsigned int j = 0; j < democ.size(); j++){
				const auto demoval = toml::find(democ,j);
				
				if(!demoval.contains("value")) emsgroot("A 'value' must be specified in 'democats'.");
				const auto value = toml::find<std::string>(demoval,"value");
				democat.value.push_back(value);
				
				if(!demoval.contains("sus")) emsgroot("The property 'sus' must be specified in 'democats'.");
				const auto sus = toml::find<std::string>(demoval,"sus");
				democat.param.push_back(sus);
			}
			
			democatvec.push_back(democat);
		}
	}
	
	return democatvec;
}

void DATA::addcovar(string name, string param, string func)
{
	COVAR cov;
	cov.name = name;
	cov.param = param;
	cov.func = func;
	cov.col = UNSET;
	
	covar.push_back(cov);
	ncovar = covar.size();
}

// area covariates
vector <COVAR> Inputs::find_covar(const Details &details) const
{
	vector <COVAR> covarvec;
	
	if(tomldata.contains("covars")){
		const auto covars = toml::find(tomldata,"covars");
		
		COVAR cov;
		for(unsigned int j = 0; j < covars.size(); j++){
			const auto covar = toml::find(covars,j);
			
			if(!covar.contains("name")) emsgroot("A 'name' must be specified in 'covars'.");
			cov.name = toml::find<std::string>(covar,"name");
	
			if(!covar.contains("param")) emsgroot("A 'param' must be specified in 'covars'.");
			cov.param = toml::find<std::string>(covar,"param");

			if(!covar.contains("func")) emsgroot("A 'func' must be specified in 'covars'.");
			cov.func = toml::find<std::string>(covar,"func");
			cov.col = UNSET;
			
			covarvec.push_back(cov);
		}
	}
	
	return covarvec;
}
