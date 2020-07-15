/*
Load mpi: module load mpi/openmpi-x86_64
Compile using: make

Simulation:  
 ./beepmbp inputfile="examples/sim.toml"        
  ./beepmbp inputfile="examples/infMSOA_noage.toml" nchain=1
 ./beepmbp inputfile="examples/simMSOA_noage.toml" mode="sim"
  mode="inf"  nsamp=10000
Inference:    
mpirun -n 2 ./beepmbp inputfile="examples/inf.toml" nchain=2
mpirun -n 20 ./beepmbp inputfile="examples/inf.toml" nchain=20

mpirun -n 20 ./beepmbp inputfile="examples/simMSOA_noage.toml" mode="inf" nchain=20 nsamp=10000
*/

// nohup mpirun -n 20 ./beepmbp inputfile="examples/simMSOA_noage.toml" nchain=20  mode="inf"  nsamp=10000 &

// nohup mpirun -n 20 gdb --batch --quiet -ex "run" -ex "bt" -ex "quit" --args ./beepmbp inputfile="examples/simMSOA_noage.toml" nchain=20  mode="inf"  nsamp=10000 &

// mpirun -n 2 gdb --batch --quiet -ex "run" -ex "bt" -ex "quit" --args ./beepmbp inputfile="examples/simMSOA.toml" nchain=2  mode="inf"  nsamp=10000 
// make -j 10
/*
 Commands for running on DiRAC:
 ./submit-csd3 --groups 65536 --seed 1 --samples 10 --nprocs 32 --nnode 1  --walltime 1:00:00 --dir ~/rds/rds-dirac-dc003/dc-pool1/test
 dos2unix submit-csd3
 squeue -u dc-pool1
 scancel <jobid>
 mybalance
*/
// sed -i 's/\r$//g' gitversion.sh    gets rid of \r
// sed -i 's/\r$//g' convert_sim.sh    gets rid of \r


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
#include "generateQ.hh"

#include "simulate.hh"
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
		
		emsgroot(ss.str());
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

int main(int argc, char** argv)
{
	unsigned int ncore, core;                 // Stores the number of cores and the core of the current process
	unsigned int nsamp=UNSET;                 // The number of samples for inference
	unsigned int nchain=UNSET;                // The number of chains per core
	
	unsigned int mode=UNSET;                  // Sets the mode of operation (sim/inf)
	
	int seed=0;                               // Sets the random seed for simulation
	
	unsigned int j, k, d;
	TRANSDATA transdata;
	POPDATA popdata;
	MARGDATA margdata;
	string str, command, value, tformat, startstr, endstr;

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
		cout << "BEEPmbp version " << GIT_VERSION << endl;
	
	DATA data;    // The following file names will need to be read in by the interface:
	MODEL model(data);
		
	data.outputdir="Output";                // The default output directory
	data.threshold=UNSET;

	// A list of all supported parameters
	vector<string>  definedparams {"threshold", "covars", "infmax", "datadir", "trans", "betaspline", "phispline", "comps", "params", "priorcomps", "priors", "democats", "ages", "regions", "areas", "genQ", "agemix", "geomix", "genQoutput", "Q", "timep", "mode", "model", "nchain", "nsamp", "timeformat", "start", "end", "timeunits", "seed", "margdata", "transdata", "popdata", "outputdir", "inputfile", "propsmethod"};

	// Read command line parameters
	map<string,string> cmdlineparams = get_command_line_params(argc, argv);
	check_for_undefined_parameters(definedparams, keys_of_map(cmdlineparams), "on command line");
	
	// Read TOML parameters
	string inputfilename = "/dev/null";
	if (cmdlineparams.count("inputfile") == 1) {
		inputfilename = cmdlineparams["inputfile"];
	}
	decltype(toml::parse(inputfilename)) tomldata;
	try {
		tomldata = toml::parse(inputfilename);
	} catch (const std::exception& e) {
		std::ostringstream oss;
		oss << "toml::parse returns exception\n" << e.what();
		emsg(oss.str());
	}
		
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
		cout << endl << "Settings:" << endl;

	string key,val;

	// DETAILS
	
	// mode
	val = lookup_string_parameter(cmdlineparams, tomldata, "mode", param_verbose, "UNSET");
	map<string,int>  modemap{{"sim", MODE_SIM}, {"inf", MODE_INF}};
	if (modemap.count(val) != 0) {
		mode = modemap[val];
	} else {
		emsgroot("Unrecoginsed value " + val + " for mode parameter");
	}

	if(mode != MODE_SIM){
		// nchain
		int nchaintot = lookup_int_parameter(cmdlineparams, tomldata, "nchain", param_verbose, UNSET);
		if (nchaintot != UNSET) {
			if(nchaintot%ncore != 0) emsgroot("The number of chains must be a multiple of the number of cores");
			nchain = nchaintot/ncore;
		}
	}
	
	if(mode != MODE_SIM){
		// nsamp
		nsamp = lookup_int_parameter(cmdlineparams, tomldata, "nsamp", param_verbose);
		
		// infmax
		if(tomldata.contains("infmax")) {
			model.infmax = lookup_int_parameter(cmdlineparams, tomldata, "infmax", param_verbose);
		}
		else emsgroot("Input file must contain a limit on the maximum number of individuals through 'infmax'.");
	}
	
	tformat = lookup_string_parameter(cmdlineparams, tomldata, "timeformat", param_verbose);
	if(tformat == "number"){ data.tform = TFORM_NUM; data.tformat = "time";}
	else{ 
		data.tformat = "date";
		if(tformat == "year-month-day") data.tform = TFORM_YMD; 
		else emsgroot("Do not recognise time format '"+tformat+"'.");
	}
	
	// start
	startstr = lookup_string_parameter(cmdlineparams, tomldata, "start", param_verbose);
	data.start = data.gettime(startstr);

	// end
	endstr = lookup_string_parameter(cmdlineparams, tomldata, "end", param_verbose);
	data.end = data.gettime(endstr);
	
	// sets the period
	data.period = data.end - data.start;

	// seed
	if(mode == MODE_SIM){
		seed = lookup_int_parameter(cmdlineparams, tomldata, "seed", param_verbose);
	}
		
	// proposals method
	string propsmethod_str = lookup_string_parameter(cmdlineparams, tomldata, "propsmethod", param_verbose, "fixedtime");

	enum proposalsmethod propsmethod;

	if (propsmethod_str == "allchainsallparams") {
		propsmethod = proposalsmethod::allchainsallparams;
	} else if (propsmethod_str == "fixednum") {
		propsmethod = proposalsmethod::fixednum;
	} else if (propsmethod_str == "fixedtime") {
		propsmethod = proposalsmethod::fixedtime;
	} else {
		emsg("Parameter propsmethod set to an unrecognised value \""+propsmethod_str+"\"");
	}

	// outputdir
	data.outputdir = lookup_string_parameter(cmdlineparams, tomldata, "outputdir", param_verbose, data.outputdir);

	// data directory
	data.datadir = lookup_string_parameter(cmdlineparams, tomldata, "datadir", param_verbose, "UNSET");
	
	// region data
	data.regiondatafile = lookup_string_parameter(cmdlineparams, tomldata, "regions", param_verbose, "UNSET");

	// area data
	data.areadatafile = lookup_string_parameter(cmdlineparams, tomldata, "areas", param_verbose, "UNSET");

	data.genQ.onoff = lookup_string_parameter(cmdlineparams, tomldata, "genQ", param_verbose, "UNSET");

	// threshold
	if(tomldata.contains("threshold")){
		data.threshold = lookup_int_parameter(cmdlineparams, tomldata, "threshold", param_verbose);
		data.thres_h = log(1.0/(data.threshold + 0.5*sqrt(2*M_PI*minvar*varfac)));
	}

	// End of parameters
	if (param_verbose)
		cout << endl;
	
	// THE MODEL
	
	// age categories
	if(tomldata.contains("ages")) {
		const auto ages = toml::find(tomldata,"ages");
		vector<string> pos, possus;
		for(j = 0; j < ages.size(); j++){
			const auto ag = toml::find(ages,j);
			
			if(!ag.contains("range")) emsgroot("A 'range' must be specified in 'ages'.");
			const auto range = toml::find<std::string>(ag,"range");
			pos.push_back(range);
			
			if(!ag.contains("sus")) emsgroot("A 'sus' must be specified in 'ages'.");
			const auto sus = toml::find<std::string>(ag,"sus");
			possus.push_back(sus);
		}
		data.adddemocat("age",pos,possus);
	}
	else emsgroot("The 'ages' parameter must be set.");
	
	if(core == 0){	
		cout << "Age categories: " << endl << "  ";
		for(j = 0; j < data.democat[0].value.size(); j++){
			if(j != 0) cout << ", ";
			cout << data.democat[0].value[j] << " sus='" <<  data.democat[0].param[j] << "'";
		}
		cout << endl << endl;
	}
	
	// democats
	if(tomldata.contains("democats")) {
		const auto democats = toml::find(tomldata,"democats");
	
		for(k = 0; k < democats.size(); k++){
			vector<string> pos, possus;
			
			const auto democat = toml::find(democats,k);
			
			for(j = 0; j < democat.size(); j++){
				const auto demoval = toml::find(democat,j);
				
				if(!demoval.contains("value")) emsgroot("A 'value' must be specified in 'democats'.");
				const auto value = toml::find<std::string>(demoval,"value");
				pos.push_back(value);
				
				if(!demoval.contains("sus")) emsgroot("The property 'sus' must be specified in 'democats'.");
				const auto sus = toml::find<std::string>(demoval,"sus");
				vector<string> pos;
				possus.push_back(sus);
			}
			
			data.adddemocat("",pos,possus);
		}
	}
	
	if(core == 0 && data.democat.size() > 1){
		cout << "Demographic categories: " << endl;
		for(k = 1; k < data.democat.size(); k++){
			cout << "  ";
			for(j = 0; j < data.democat[k].value.size(); j++){
				if(j != 0) cout << ", ";
				cout << data.democat[k].value[j] << " sus='" <<  data.democat[k].param[j] << "'";
			}	
			cout << endl;
		}
		cout << endl;
	}
	
	// area covariates
	if(tomldata.contains("covars")) {
		const auto covars = toml::find(tomldata,"covars");
		for(j = 0; j < covars.size(); j++){
			const auto covar = toml::find(covars,j);
			
			if(!covar.contains("name")) emsgroot("A 'name' must be specified in 'covars'.");
			const auto name = toml::find<std::string>(covar,"name");
	
			if(!covar.contains("param")) emsgroot("A 'param' must be specified in 'covars'.");
			const auto par = toml::find<std::string>(covar,"param");

			if(!covar.contains("func")) emsgroot("A 'func' must be specified in 'covars'.");
			const auto func = toml::find<std::string>(covar,"func");

			data.addcovar(name,par,func);
		}
		
		if(core == 0){
			cout << "Area covariates: " << endl;
			cout << "  ";
			for(j = 0; j < data.covar.size(); j++) cout << data.covar[j].name << "   param='" << data.covar[j].param << "'" << endl; 
			cout << endl;
		}
	}

	// timep
	if(tomldata.contains("timep")) {
		const auto timep = toml::find(tomldata,"timep");
		for(j = 0; j < timep.size(); j++){
			const auto tim = toml::find(timep,j);
			
			if(!tim.contains("name")) emsgroot("A 'name' must be specified in 'timep'.");
			const auto name = toml::find<std::string>(tim,"name");
			
			if(!tim.contains("tend")) emsgroot("'tend' must be specified in 'timep'.");
			auto tendstr = toml::find<string>(tim,"tend");
			int tend = data.gettime(tendstr) - data.start;
			
			if(tend < 0 || tend > (int)data.period) emsg("Time '"+tendstr+"' is out of range."); 
			if(j > 0){
				if(tend < data.timeperiod[j-1].tend) emsg("'timep' is not time ordered.");
			}
			
			if(j == timep.size()-1){
				if(tend != (int)data.period) emsg("'tend' in 'timep' must finish with the end time.");
			}
			data.addtimep(name,tend);
		}
	}
	else emsgroot("Property 'timep' defining time periods must be set.");

	if(core == 0){
		cout << "Time periods defined:" << endl;
		for(j = 0; j < data.timeperiod.size(); j++){
			cout << "  ";
			cout << data.timeperiod[j].name << ": ";
			if(j == 0) cout << "0"; else cout << data.timeperiod[j-1].tend;
			cout << " - " <<  data.timeperiod[j].tend << endl;
		}
		cout << endl;
	}

	// genQ
	if(data.genQ.onoff == "on"){
		if(tomldata.contains("agemix")) {
			const auto agemix = toml::find(tomldata,"agemix");
			
			if(!agemix.contains("Nall")) emsgroot("'Nall' must be specified in 'agemix'.");
			const auto Nall = toml::find<std::string>(agemix,"Nall");
			data.genQ.Nall = Nall;
			
			if(!agemix.contains("Nhome")) emsgroot("'Nhome' must be specified in 'agemix'.");
			const auto Nhome = toml::find<std::string>(agemix,"Nhome");
			data.genQ.Nhome = Nhome;
			
			if(!agemix.contains("Nother")) emsgroot("'Nother' must be specified in 'agemix'.");
			const auto Nother = toml::find<std::string>(agemix,"Nother");
			data.genQ.Nother = Nother;
			
			if(!agemix.contains("Nschool")) emsgroot("'Nschool' must be specified in 'agemix'.");
			const auto Nschool = toml::find<std::string>(agemix,"Nschool");
			data.genQ.Nschool = Nschool;
			
			if(!agemix.contains("Nwork")) emsgroot("'Nwork' must be specified in 'agemix'.");
			const auto Nwork = toml::find<std::string>(agemix,"Nwork");
			data.genQ.Nwork = Nwork;
		}
		else emsgroot("'agemix' must be specified.");
	
		if(tomldata.contains("geomix")) {
			const auto geomix = toml::find(tomldata,"geomix");
			
			if(!geomix.contains("M")) emsgroot("'M' must be specified in 'geomix'.");
			const auto M = toml::find<std::string>(geomix,"M");
			data.genQ.M = M;
		}
		else emsgroot("'geomix' must be specified.");
		
		if(tomldata.contains("genQoutput")) {
			const auto qout = toml::find(tomldata,"genQoutput");
			
			if(!qout.contains("localhome")) emsgroot("'localhome' must be specified in 'genQoutput'.");
			const auto localhome = toml::find<std::string>(qout,"localhome");
			data.genQ.localhome = localhome;
	
			if(!qout.contains("flowall")) emsgroot("'flowall' must be specified in 'genQoutput'.");
			const auto flowall = toml::find<std::string>(qout,"flowall");
			data.genQ.flowall = flowall;
		}
		else emsgroot("'genQoutput' must be specified.");
		
		data.genQ.nage = data.democat[0].value.size();
		data.genQ.datadir = data.datadir;
		data.genQ.outputdir = data.outputdir;
		data.genQ.areadata = data.areadatafile;
		if(core == 0) generateQ(data.genQ); 
	}

	// Q
	if(tomldata.contains("Q")) {
		const auto Qlist = toml::find(tomldata,"Q");
		for(j = 0; j < Qlist.size(); j++){
			const auto Q = toml::find(Qlist,j);
			
			if(!Q.contains("timep")) emsgroot("A 'timep' must be specified in 'Q'.");
			const auto timep = toml::find<std::string>(Q,"timep");
			
			if(!Q.contains("comp")) emsgroot("'comp' must be specified in 'Q'.");
			const auto comp = toml::find<std::string>(Q,"comp");
			
			if(!Q.contains("file")) emsgroot("'file' must be specified in 'Q'.");
			const auto file = toml::find<std::string>(Q,"file");
			
			data.addQtensor(timep,comp,file);
		}
	}
	else emsgroot("Property 'timep' defining time periods must be set.");
	
	if(core == 0){
		cout << "Q tensors loaded:" << endl;
		for(j = 0; j < data.Q.size(); j++){
			cout << "  ";
			cout << "timep: " << data.timeperiod[data.Q[j].timep].name << "  ";
			cout << "compartment: " << data.Q[j].comp << "  ";
			cout << "file: " << data.Q[j].file << "  ";
			cout << endl;
		}
		cout << endl;
	}

	// THE DATA
	
	// transdata
	if(tomldata.contains("transdata")) {
		const auto tdata = toml::find(tomldata,"transdata");

		for(j = 0; j < tdata.size(); j++){
			const auto td = toml::find(tdata,j);
			
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
			transdata.start = data.gettime(startdata)-data.start;
			
			if(!td.contains("units")) emsgroot("A 'units' property must be specified in 'transdata'.");
			const auto units = toml::find<std::string>(td,"units");
			if(units == "days") transdata.units = 1;
			else{
				if(units == "weeks") transdata.units = 7;
				else emsgroot("Units in 'transdata' not recognised");
			}
			
			if(mode == MODE_SIM){
				transdata.rows = (unsigned int)((data.period - transdata.start)/transdata.units);
				if(transdata.rows == 0) emsgroot("Transition data '"+file+"' cannot be generated because the time period is not sufficiently long.");
			}
			
			data.transdata.push_back(transdata);
		}
	}

	// popdata
	if(tomldata.contains("popdata")) {
		const auto pdata = toml::find(tomldata,"popdata");

		for(j = 0; j < pdata.size(); j++){
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
			popdata.start = data.gettime(startdata)-data.start;
			
			if(!pd.contains("units")) emsgroot("A 'units' property must be specified in 'popdata'.");
			const auto units = toml::find<std::string>(pd,"units");
			
			if(units == "days") popdata.units = 1;
			else{
				if(units == "weeks") popdata.units = 7;
				else emsgroot("Units in 'popdata' not recognised");
			}
			
			if(mode == MODE_SIM){
				popdata.rows = (unsigned int)((data.period - popdata.start)/popdata.units);
				if(popdata.rows == 0) emsgroot("popition data '"+file+"' cannot be generated because the time period is not sufficiently long.");
			}
			
			data.popdata.push_back(popdata);
		}
	}

	if(!tomldata.contains("transdata") && !tomldata.contains("popdata"))  emsgroot("'transdata' and/or 'popdata' must be set.");
			
	// margdata
	
	if(tomldata.contains("margdata")) {
		const auto mdata = toml::find(tomldata,"margdata");

		for(j = 0; j < mdata.size(); j++){
			const auto md = toml::find(mdata,j);
			
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
			for(d = 0; d < data.ndemocat; d++) if(type == data.democat[d].name) break;
			if(d == data.ndemocat) emsg("The 'type' property must be 'age' or a demographic property.");
			margdata.democat = d;
				
			if(!md.contains("file")) emsgroot("A 'file' property must be specified in 'margdata'.");
			const auto file = toml::find<std::string>(md,"file");
			margdata.file = file;
			
			data.margdata.push_back(margdata);
		}
	}

	if(mode == MODE_SIM && ncore != 1) emsgroot("Simulation only requires one core");
	
	if(mode == UNSET) emsgroot("The property 'mode' must be set"); 
 
	MPI_Barrier(MPI_COMM_WORLD);

	data.readdata(core,ncore,mode);

	POPTREE poptree;
	poptree.init(data,core);	

	model.definemodel(core,data.period,data.popsize,tomldata);
	model.addQ();
	model.checkdata();

	if(core == 0) cout << "Running...." << endl;

	timersinit();
	timers.timetot = -clock();
	
	MPI_Barrier(MPI_COMM_WORLD);

	srand(core*10000+seed);
	
	switch(mode){
	case MODE_SIM:
		simulatedata(data,model,poptree);
		break;

	case MODE_INF: 
		if(nsamp == UNSET) emsgroot("The number of samples must be set");
		if(nchain == UNSET) emsgroot("The number of chains must be set");
		MBP(data,model,poptree,nsamp,core,ncore,nchain,propsmethod);
		break;

	default: emsgroot("Mode not recognised"); break;
 	}
	
	timers.timetot += clock();
	
	if(core == 0){
		cout << double(timers.timetot)/CLOCKS_PER_SEC << " Total time (seconds)" << endl;
		if(mode == MODE_SIM){
			cout << double(timers.timesim)/CLOCKS_PER_SEC << " Simulation time (seconds)" << endl;
		}
		else{
			cout << double(timers.timewait)/CLOCKS_PER_SEC << " MBP waiting time (seconds)" << endl;
			cout << double(timers.timembp)/CLOCKS_PER_SEC << " MBP time (seconds)" << endl;
			cout << double(timers.timembpinit)/CLOCKS_PER_SEC << " MBP init (seconds)" << endl;
			cout << double(timers.timembpQmap)/CLOCKS_PER_SEC << " MBP Qmap (seconds)" << endl;
			cout << double(timers.timembpprop)/CLOCKS_PER_SEC << " MBP prop (seconds)" << endl;
		}
	}
	
	#ifdef USE_MPI
	MPI_Finalize();
	#endif
}
