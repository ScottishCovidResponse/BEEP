/*
Load mpi: module load mpi/openmpi-x86_64

Compile using: make

Simulation:  ./beepmbp inputfile="examples/sim.toml"        

Inference:    mpirun -n 20 ./beepmbp inputfile="examples/inf.toml" nchain=20
*/

#include <iostream>
#include <sstream>
#include <math.h>
#include <map>
#include <algorithm>
#include <vector>
#include <iterator>

#include "toml11/toml.hpp"

#include "stdlib.h"
#include "time.h"

#include "inputs.hh"
#include "details.hh"
#include "data.hh"
#include "poptree.hh"
#include "model.hh"
#include "output.hh"
#include "obsmodel.hh"

#include "utils.hh"
#include "timers.hh"

#include "simulate.hh"
#include "mcmc.hh"
#include "consts.hh"

#include "combinetrace.hh"

#ifdef USE_DATA_PIPELINE
#include "pybind11/embed.h"
#include "datapipeline.hh"
#endif

using namespace std;

string lookup_string_parameter(const map<string,string> &params,
															 const toml::basic_value<::toml::discard_comments, std::unordered_map, std::vector> &tomldata,
															 const string &key, bool verbose, const string &def="")                        // zz
{
	string val;
	auto val_it = params.find(key);
	if (val_it != params.end()) {
		val = val_it->second;
	} else {
		if (tomldata.contains(key)) {
			val = toml::find<string>(tomldata,key);
		} else {
			val = def;
			// emsgroot("ERROR: Parameter \'"+key+"\' must be supplied");
		}
	}
	if (verbose)
		cout << "  " << key << " = " << val << endl;
	return val;
}


int lookup_int_parameter(const map<string,string> &params,
												 const toml::basic_value<::toml::discard_comments, std::unordered_map, std::vector> &tomldata,
												 const string &key, bool verbose, int def=-1)
{
	int val;
	auto val_it = params.find(key);
	if (val_it != params.end()) {
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
			// emsgroot("ERROR: Parameter \'"+key+"\' must be supplied");
		}
	}
	if (verbose)
		cout << "  " << key << " = " << val << endl;
	return val;
}

double lookup_double_parameter(const map<string,string> &params,
												 const toml::basic_value<::toml::discard_comments, std::unordered_map, std::vector> &tomldata,
												 const string &key, bool verbose, int def=-1) //zz
{
	double val;
	auto val_it = params.find(key);
	if (val_it != params.end()) {
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
			// emsgroot("ERROR: Parameter \'"+key+"\' must be supplied");
		}
	}
	if (verbose)
		cout << "  " << key << " = " << val << endl;
	return val;
}


vector<string> string_split(const string &s)
{
	std::stringstream ss(s);
	std::istream_iterator<std::string> begin(ss);
	std::istream_iterator<std::string> end;
	std::vector<std::string> vstrings(begin, end);
	// std::copy(vstrings.begin(), vstrings.end(), std::ostream_iterator<std::string>(std::cout, "\n"));
	return vstrings;
}

vector<string> lookup_stringlist_parameter(
	const map<string,string> &params,
	const toml::basic_value<::toml::discard_comments, std::unordered_map, std::vector> &tomldata,
	const string &key, bool verbose)
{
	vector<string> val;
	auto val_it = params.find(key);
	if (val_it != params.end()) {
		val = string_split(val_it->second);
	} else {
		if (tomldata.contains(key)) {
			val = toml::find<vector<string>>(tomldata,key);
		} else {
			emsgroot("ERROR: Parameter \'"+key+"\' must be supplied");
		}
	}

	if (verbose) {
		cout << "  " << key << " = ";

		for (auto it = val.begin(); it != val.end(); it++) {
			cout << *it << " ";
		}
		cout << endl;
	}
	return val;
}

vector<string> get_toml_keys( // zz
	const toml::basic_value<::toml::discard_comments, std::unordered_map, std::vector> &data)
{
	vector<string> keys;
	for(const auto& p : data.as_table())
	{
		keys.push_back(p.first);
	}
	return keys;
}

void check_for_undefined_parameters(vector<string> allowed, vector<string> given,
																		const string &context) // zz
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

map<string,string> get_command_line_params(int argc, char *argv[])    // zz
{
	map<string,string> cmdlineparams;
	
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
	return cmdlineparams;
}

vector<string> keys_of_map(map<string,string> m)
{
	vector<string> keys;
	for (const auto &p : m) {
		keys.push_back(p.first);
	}
	return keys;
}

string gitversion();

int main(int argc, char** argv)
{
	#ifdef USE_MPI
  MPI_Init(&argc, &argv);
  #endif

	Mpi mpi;                                 // Stores mpi information (core and ncore)
	
	bool verbose = (mpi.core == 0);            // Paramater which ensure that only core 0 outputs results
		
	if(verbose) cout << "BEEPmbp version " << gitversion() << endl;

	Inputs inputs(argc,argv,verbose);                 // Loads up the command line arguments
	
	
	Details details(inputs);
	
	Mode mode = details.mode;
	//unsigned int nsamp=UNSET;                 // The number of samples for inference
	//unsigned int nchain=UNSET;                // The number of chains per core
	
	//unsigned int mode=UNSET;                  // Sets the mode of operation (sim/inf)

	
	bool param_verbose = (mpi.core == 0);         // Paramater which ensure that only core 0 outputs results

	unsigned int j, k, d;
	//TRANSDATA transdata;
	//POPDATA popdata;
	//MARGDATA margdata;
	string str, command, value;//, startstr, endstr;

	
#ifdef USE_DATA_PIPELINE
	pybind11::scoped_interpreter guard{};

	using namespace pybind11::literals;

  pybind11::module::import("logging").attr("basicConfig")("level"_a="DEBUG", "format"_a="%(asctime)s %(filename)s:%(lineno)s %(levelname)s - %(message)s");

	// TODO: get the path to the config file from the options
 	DataPipeline *dp = new DataPipeline(
		"dpconfig.yaml", "https://github.com/ScottishCovidResponse/BEEPmbp",
		gitversion());

	DATA data(inputs,details,mpi,dp);    // The following file names will need to be read in by the interface:
#else
	DATA data(inputs,details,mpi);    // The following file names will need to be read in by the interface:
#endif


	if(mode == combinetrace){
		combine_trace(data,inputs);
		return 0;
	}
	
	MODEL model(inputs,details,data);
		
	data.read_data_files(inputs,model,mpi);
	
	
	//data.threshold=UNSET;

	// A list of all supported parameters (please keep in lexicographic order)   // zz
	vector<string>  definedparams {
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
	
	// Read command line parameters
	map<string,string> cmdlineparams = get_command_line_params(argc, argv); //zz
	check_for_undefined_parameters(definedparams, keys_of_map(cmdlineparams), "on command line");//zz
	
	
	// Read TOML parameters    zz
	string inputfilename = "/dev/null";
	if (cmdlineparams.count("inputfile") == 1) {
		inputfilename = cmdlineparams["inputfile"];
	}
	decltype(toml::parse(inputfilename)) tomldata;
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

	// The code could be simplified by reading the TOML parameters into a
	// map<string,string> and merging it with cmdlineparams. However, this would
	// require casting all the TOML values to strings, and we would lose the
	// vectors and types in the TOML file, and have to parse them out of
	// strings. An alternative would be to interpret the command line as a TOML
	// fragment, and work with just TOML values. That might be a better approach.

	/*********************************************************************************
	/ Process parameters
  **********************************************************************************/



	if (param_verbose)
		cout << endl << "Settings:" << endl;

	string key,val;

	// DETAILS
	
	// mode
	/*
	val = lookup_string_parameter(cmdlineparams, tomldata, "mode", param_verbose, "UNSET");     // zz
	map<string,int>  modemap{{"sim", MODE_SIM}, {"inf", MODE_INF}, {"multisim", MODE_MULTISIM}};
	if (modemap.count(val) != 0) {
		mode = modemap[val];
	} else {
		emsgroot("Unrecoginsed value " + val + " for mode parameter");
	}
*/
	
	/*
	if(mode == multisim || mode == inf){
		// nsamp
		nsamp = lookup_int_parameter(cmdlineparams, tomldata, "nsamp", param_verbose);
	}
	*/
	
	//model.infmax = large;
	
	
	// start
	/*
	startstr = lookup_string_parameter(cmdlineparams, tomldata, "start", param_verbose);
	details.start = details.gettime(startstr);

	// end
	endstr = lookup_string_parameter(cmdlineparams, tomldata, "end", param_verbose);
	data.end = details.gettime(endstr);
	
	// sets the period
	details.period = data.end - details.start;
*/

	
	/*
	// seed
	if(tomldata.contains("seed")){
		seed = lookup_int_parameter(cmdlineparams, tomldata, "seed", param_verbose);
	}
	*/
	
	


	// outputdir
	//data.outputdir = lookup_string_parameter(cmdlineparams, tomldata, "outputdir", param_verbose, data.outputdir);

	// data directory
	//data.datadir = lookup_string_parameter(cmdlineparams, tomldata, "datadir", param_verbose, "UNSET");
	
	// region data
	//data.regiondatafile = lookup_string_parameter(cmdlineparams, tomldata, "regions", param_verbose, "UNSET");

	// area data
	//data.areadatafile = lookup_string_parameter(cmdlineparams, tomldata, "areas", param_verbose, "UNSET");

	
	// THE MODEL
	
	
	if(verbose) data.print_to_terminal();
	
	MPI_Barrier(MPI_COMM_WORLD);

	//data.read_data_files(inputs,mpi);

	POPTREE poptree;
	poptree.init(data,mpi.core);	

	model.definemodel(mpi.core,details.period,data.popsize,tomldata);
	model.addQ();
	model.checkdata();

	if(verbose) cout << "Running...." << endl;

	timersinit();
	timers.timetot = -clock();
	
	MPI_Barrier(MPI_COMM_WORLD);

	unsigned int seed = inputs.find_int("seed",0);
	sran(mpi.core*10000+seed);
	
	if(mode != inf && mpi.ncore != 1) emsgroot("Simulation only requires one core");
	
	Obsmodel obsmodel(details,data,model);
	Output output(details,data,model,obsmodel);
	
	switch(mode){
	case sim:
		{		
			Simulate simu(details,data,model,poptree,mpi,inputs,output,obsmodel,mode,verbose);
			simu.run();
		}
		break;
	
	case multisim:
		{
			Simulate simu(details,data,model,poptree,mpi,inputs,output,obsmodel,mode,verbose);
			simu.multirun();
		}
		break;
			
	case inf: 
		{
			Mcmc mcmc(details,data,model,poptree,mpi,inputs,output,obsmodel,mode,verbose);
			mcmc.run();
		}
		break;

	default: emsgroot("Mode not recognised"); break;
 	}
	
	timers.timetot += clock();
	
	if(verbose) cout << double(timers.timetot)/CLOCKS_PER_SEC << " Total time (seconds)" << endl;
	
#ifdef USE_DATA_PIPELINE
	delete dp;
#endif

#ifdef USE_MPI
	MPI_Finalize();
#endif
}
