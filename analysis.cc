/*
Load mpi: module load mpi/openmpi-x86_64
Compile using: make

INPUTS:

Simulation:         
mpirun -n 1 ./run mode=sim model=irish simtype=smallsim seed=0 period=16 transdata=I,H,reg,cases.txt transdata=H,D,all,deaths.txt

MBP Inference:    
mpirun -n 20 ./run mode=mbp model=irish simtype=smallsim nchain=20 nsamp=1001 period=16 transdata=I,H,reg,cases.txt transdata=H,D,all,deaths.txt

mpirun -n 1 ./run mode=mbp model=irish simtype=smallsim nchain=1 nsamp=10001 period=16 transdata=I,H,reg,cases.txt transdata=H,D,all,deaths.txt

PMCMC Inference:  
mpirun -n 20 ./run mode=pmcmc model=irish simtype=smallsim npart=20 nsamp=1000 period=16 transdata=I,H,reg,cases.txt transdata=H,D,all,deaths.txt

The flag -n 1 sets the number of cores (set to 1 for simulation or more for inference)

Here is a description of the various inputs:

mode - Defines how the code operates:
		"sim" generates simulated data.
		"pmcmc" performs inference using particle MCMC.
		"mbp" performs inference using multi-temperature MBP MCMC.

model - This defines the compartmental model being used:
		"irish" defines a SEAPIRHD model for asymtopmatic / presymptomatic / symptomatic individuals.
	
area - Determines	the number of areas into which the houses are divided (should be a power of 4)

npart - The total number of particles used when performing PMCMC (should be a multiple of the number of cores).

nchain - The total number of chains used when performing multi-temperature MBP MCMC (should be a multiple of the number of cores).

nsamp - The number of samples used for inference (note, burnin is assumed to be a quater this value).

seed - Sets the random seed when performing inference (this is set to zero by default)

outputdir - Gives the name of the output directory (optional).

transdata - Transition data. Gives the "from" then "to" compartments, the type ("reg" means regional data and "all" means global) and then the file name. More than one set of transition data can be added for an analysis.

housedata - House data. Gives the file name for a file giving the positions of houses.

period - The period of time over which simulation/inference is performed (e.g. measured in weeks).

simtype - Determines the system on which simulation is performed (simulation mode only):
		"smallsim" a small system with 1024 houses, 10000 individuals and a 2x2 grid of data regions (used for testing)
		"scotsim" a Scotland-like system with 1.5 million houses, 5.5 million individuals and a 4x4 grid of data regions (used for testing)
		"uksim" a UK-like system with 20 million houses, 68 million individuals and a 10x10 grid of data regions (used for testing)
	
 
OUTPUTS:

Simulation - This creates the specified 'transdata' and 'housedata' files along with output directory containing:
1) Plots for the transitions corresponding to the 'transdata' files.
2) "R0.txt", which gives time variation in R0.
3) "parameter.txt", which gives the parameter values used in the simulation.

Inference - The output directory contains postior information (with means and 90% credible intervals) for:
1) Plots for the transitions corresponding to the 'transdata' files.
2) "R0.txt", which gives posterior plotstime variation in R0.
3) "parameter.txt", which gives information about parameters.
4) "trace.txt", which gives trace plots for different models.
5) "traceLi.txt", which gives trace plots for the likelihoods on different chains (MBPs only).
6) "MCMCdiagnostic.txt", which gives diagnostic information on the MCMC algorthm.
*/

/*
 Commands for running on DiRAC:
 ./submit-csd3 --groups 65536 --seed 1 --samples 10 --nprocs 32 --nnode 1  --walltime 1:00:00 --dir ~/rds/rds-dirac-dc003/dc-pool1/test
 dos2unix submit-csd3
 squeue -u dc-pool1
 scancel <jobid>
 mybalance
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

#include "utils.hh"
#include "timers.hh"
#include "poptree.hh"
#include "model.hh"
#include "data.hh"

#include "simulate.hh"
#include "PMCMC.hh"
#include "gitversion.hh"
#include "MBP.hh"
#include "consts.hh"

using namespace std;
 


string lookup_string_parameter(const map<string,string> &params,
															 const toml::basic_value<::toml::discard_comments, std::unordered_map, std::vector> &tomldata,
															 const string &key, bool verbose, const string &def="")
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
			// emsg("ERROR: Parameter \'"+key+"\' must be supplied");
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
		val = stoi(val_it->second);
	} else {
		if (tomldata.contains(key)) {
			val = toml::find<int>(tomldata,key);
		} else {
			val = def;
			// emsg("ERROR: Parameter \'"+key+"\' must be supplied");
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
			emsg("ERROR: Parameter \'"+key+"\' must be supplied");
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

vector<string> get_toml_keys(
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
																		const string &context)
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
		
		emsg(ss.str());
	}
}

map<string,string> get_command_line_params(int argc, char *argv[])
{
	map<string,string> cmdlineparams;

	// Store the parameters passed on the command line in cmdlineparams
	for(int op = 1; op < argc; op++){ // Goes the various input options
		string str = string(argv[op]);
		int j = 0; int jmax = str.length(); while(j < jmax && str.substr(j,1) != "=") j++;
		if(j == jmax){
			stringstream ss; ss << "Cannot understand " << str; 
			emsg(ss.str());
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


int main(int argc, char** argv)
{
	unsigned int ncore, core;                 // Stores the number of cores and the core of the current process
	unsigned int nsamp=UNSET;                 // The number of samples for inference
	unsigned int nchain=UNSET;                // The number of chains per core (MBP only)
	unsigned int npart=UNSET;                 // The number of particles per core (PMCMC only)
	
	unsigned int mode=UNSET;                  // Sets the mode of operation (sim/PMCMC/MBP)
  unsigned int modelsel=UNSET;              // Sets the model used
	unsigned int period=UNSET;                // The period of time for simulation / inference
	
	int seed=0;                               // Sets the random seed for simulation
	
	unsigned int j, jst, jmax, flag;
	int op;
	TRANSDATA transdata;
	string str, command, value;

	#ifdef USE_MPI
  MPI_Init(&argc, &argv);
	int num;
	MPI_Comm_size(MPI_COMM_WORLD,&num); ncore =num;
  MPI_Comm_rank(MPI_COMM_WORLD,&num); core =num;
  #endif

	#ifndef USE_MPI
	ncore = 1;
	core = 0;
	#endif

	if (core == 0)
		cout << "CoronaPMCMC version " << GIT_VERSION << endl;
	
	DATA data;    // The following file names will need to be read in by the interface:
	
	data.democatfile = "Data_small/democat.txt";
	data.regiondatafile = "Data_small/regiondata.txt";  
	data.areadatafile = "Data_small/areadata.txt";  
	data.Mdatafile = "Data_small/Mdata.txt";
	data.Ndatafile = "Data_small/Ndata.txt";   
	
/*
	//data.democatfile = "Data_scotland/democat.txt";
	data.democatfile = "Data_scotland/democat_noage.txt";
	data.regiondatafile = "Data_scotland/regiondata.txt";  
	//data.areadatafile = "Data_scotland/areadata.txt";  
	data.areadatafile = "Data_scotland/areadata_noage.txt";  
	data.Mdatafile = "Data_scotland/Mdata.txt";
	data.Ndatafile = "Data_scotland/Ndata.txt";   
	*/
	
	data.outputdir="Output";                // The default output directory
		
	// A list of all supported parameters
	vector<string>  definedparams {"mode", "model", "simtype", "npart", "nchain", "nsamp",
																 "period", "seed", "transdata", "outputdir", "inputfile"};

	// Read command line parameters
	map<string,string> cmdlineparams = get_command_line_params(argc, argv);
	check_for_undefined_parameters(definedparams, keys_of_map(cmdlineparams), "on command line");
	
	// Read TOML parameters
	string inputfilename = "/dev/null";
	if (cmdlineparams.count("inputfile") == 1) {
		inputfilename = cmdlineparams["inputfile"];
	}
  auto tomldata = toml::parse(inputfilename);
	vector<string>  tomlkeys = get_toml_keys(tomldata);
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

	bool param_verbose = (core == 0);

	if (param_verbose)
		cout << endl << "Parameters:" << endl;

	string key,val;

	// mode
	val = lookup_string_parameter(cmdlineparams, tomldata, "mode", param_verbose, "UNSET");
	map<string,int>  modemap{{"sim", MODE_SIM}, {"pmcmc", MODE_PMCMC}, {"mbp", MODE_MBP}};
	if (modemap.count(val) != 0) {
		mode = modemap[val];
	} else {
		emsg("Unrecoginsed value " + val + " for mode parameter");
	}

	// model
	val = lookup_string_parameter(cmdlineparams, tomldata, "model", param_verbose, "UNSET");
	if (val == "irish") {
		modelsel = MOD_IRISH;
	} else {
		emsg("Unrecognised value "+val+" for model parameter");
	}

	// simtype
	data.simtype = lookup_string_parameter(cmdlineparams, tomldata, "simtype", param_verbose, "UNSET");
	if (!(data.simtype == "smallsim" || data.simtype == "scotsim" ||
				data.simtype == "uksim" || data.simtype == "")) {
		emsg("Unrecognised value \'"+data.simtype+"\' for simtype parameter");
	}

	// npart
	int nparttot = lookup_int_parameter(cmdlineparams, tomldata, "npart", param_verbose, npart);
	if (npart != UNSET) {
		if(nparttot%ncore != 0) emsg("The number of particles must be a multiple of the number of cores");
		npart = nparttot/ncore;
	}

	// nchain
	int nchaintot = lookup_int_parameter(cmdlineparams, tomldata, "nchain", param_verbose, UNSET);
	if (nchaintot != UNSET) {
		if(nchaintot%ncore != 0) emsg("The number of chains must be a multiple of the number of cores");
		nchain = nchaintot/ncore;
	}

	// nsamp
	nsamp = lookup_int_parameter(cmdlineparams, tomldata, "nsamp", param_verbose);

	// period
	period = lookup_int_parameter(cmdlineparams, tomldata, "period", param_verbose);

	// seed
	seed = lookup_int_parameter(cmdlineparams, tomldata, "seed", param_verbose);

	// transdata
	vector<string> slval = lookup_stringlist_parameter(cmdlineparams, tomldata, "transdata",
																										 param_verbose);

	for (auto it = slval.begin(); it != slval.end(); it++) {

		string value = *it;

		j = 0; jmax = value.length();
		while(j < jmax && value.substr(j,1) != ",") j++;
		if(j == jmax) emsg("Problem with transition data");
		transdata.from = value.substr(0,j);
		j++;
			
		jst = j; while(j < jmax && value.substr(j,1) != ",") j++;
		if(j == jmax) emsg("Problem with transition data");
		transdata.to = value.substr(jst,j-jst);
		j++;
			
		jst = j; while(j < jmax && value.substr(j,1) != ",") j++;
		if(j == jmax) emsg("Problem with transition data");
		transdata.type = value.substr(jst,j-jst);
		if(transdata.type != "reg" && transdata.type != "all") emsg("Transition data type not recognised"); 
		j++;
			
		transdata.file = value.substr(j,jmax-j);
		data.transdata.push_back(transdata);
	}	

	// outputdir
	data.outputdir = lookup_string_parameter(cmdlineparams, tomldata, "outputdir", param_verbose,
																					 data.outputdir);

	// End of parameters
	if (param_verbose)
		cout << endl;
	
	if(mode == MODE_SIM && ncore != 1) emsg("Simulation only requires one core");
	
	if(core == 0) cout << "Initialising...." << endl;
 
	if(mode == UNSET) emsg("The property 'mode' must be set"); 
 
	MPI_Barrier(MPI_COMM_WORLD);

	if(period == UNSET) emsg("The property 'period' must be set");
	//if(data.housefile=="") emsg("The property 'housedata' must be set");

	data.readdata(core,ncore,mode,period);

	//if(narea==-1) emsg("The property 'narea' must be set"); 
	
	POPTREE poptree;
	poptree.init(data,core);	
		
	MODEL model;
	model.definemodel(data,core,data.period,data.popsize,modelsel);
	model.addQ(data);
	
	model.checktransdata(data);
	
	model.setsus(data);

	if(core == 0) cout << "Running...." << endl;

	timersinit();
	timers.timetot = -clock();
	
	MPI_Barrier(MPI_COMM_WORLD);

	switch(mode){
	case MODE_SIM:
		srand(seed); 
		simulatedata(data,model,poptree);
		break;
		
	case MODE_PMCMC:
		if(nsamp == UNSET) emsg("The number of samples must be set");
		if(npart == UNSET) emsg("The number of particles must be set");

		PMCMC(data,model,poptree,nsamp,core,ncore,npart);
		break;

	case MODE_MBP: 
		if(nsamp == UNSET) emsg("The number of samples must be set");
		if(nchain == UNSET) emsg("The number of chains must be set");
		MBP(data,model,poptree,nsamp,core,ncore,nchain);
		break;

	default: emsg("Mode not recognised"); break;
 	}
	
	timers.timetot += clock();
	
	if(core == 0){
		cout << double(timers.timetot)/CLOCKS_PER_SEC << " Total time (seconds)" << endl;
		if(mode != MODE_MBP) cout << double(timers.timesim)/CLOCKS_PER_SEC << " Simulation time (seconds)" << endl;
		if(mode == MODE_PMCMC) cout << double(timers.timeboot)/CLOCKS_PER_SEC << " Bootstrap time (seconds)" << endl;
		if(mode == MODE_MBP) cout << double(timers.timewait)/CLOCKS_PER_SEC << " MBP waiting time (seconds)" << endl;
		if(mode == MODE_MBP) cout << double(timers.timembp)/CLOCKS_PER_SEC << " MBP time (seconds)" << endl;
		if(mode == MODE_MBP) cout << double(timers.timembpinit)/CLOCKS_PER_SEC << " MBP init (seconds)" << endl;
		if(mode == MODE_MBP) cout << double(timers.timembpQmap)/CLOCKS_PER_SEC << " MBP Qmap (seconds)" << endl;
		if(mode == MODE_MBP) cout << double(timers.timembpprop)/CLOCKS_PER_SEC << " MBP prop (seconds)" << endl;
	}
	
	#ifdef USE_MPI
	MPI_Finalize();
	#endif
}
