
#include <algorithm>
#include <iostream>
#include <map>
#include <vector>
#include <string>
#include <sstream>

using namespace std;

#include "inputs.hh"
#include "reader.hh"
#include "utils.hh"
#include "consts.hh"

#include "data.hh"
#include "model.hh"

Inputs::~Inputs()
{
	delete basedata;
}

/// /// Reads TOML and command line parameters
Inputs::Inputs(int argc, char** argv, bool /* verbose */) 
{
	set_command_line_params(argc,argv);                          // Loads up the command line parameters
	
	string inputfilename = "/dev/null";
	if (cmdlineparams.count("inputfile") == 1) {
		inputfilename = cmdlineparams["inputfile"];
	}

	try {
		basedata = new InputData{inputfilename};

	} catch (const std::exception& e) {
		std::ostringstream oss;
		oss << "toml::parse returns exception\n" << e.what();
		emsgroot(oss.str());
	}	
	
	vector<string> tomlkeys = basedata->keys();
	check_for_undefined_parameters(definedparams, tomlkeys, "in " + inputfilename);
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

/// Finds a string from the TOML file (or uses the default 'def' value)
string Inputs::find_string(const string &key, const string &def) const
{
	string val;
	auto val_it = cmdlineparams.find(key);
	if(val_it != cmdlineparams.end()) val = val_it->second;
	else{
		if(basedata->contains(key))
			val = basedata->data.stringfield_unchecked(key);
		else val = def;
	}

	return val;
}

/// Finds an integer from the TOML file (or uses the default 'def' value)
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
			emsgroot(oss.str());
		}
	} else {
		if (basedata->contains(key)) {
			val = basedata->data.intfield_unchecked(key);
		} else {
			val = def;
		}
	}
	
	return val;
}

/// Finds a double from the TOML file (or uses the default 'def' value)
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
			if (what == "stof") {
				if (valstr == "")
					oss << "Should be number, found no value\n";
			} else {
				oss << what;
			}
			emsgroot(oss.str());
		}
	} else {
		if (basedata->contains(key)) {
			val = basedata->data.numberfield_unchecked(key);
		} else {
			val = def;
		}
	}
	
	return val;
}

/// Checks for unrecognised parameters
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

/// Returns the mode of operation
Mode Inputs::mode() const
{
	string val = find_string("mode","UNSET");  
	if(val == "UNSET") emsgroot("The 'mode' property must be set");
	
	Mode mode;
	map<string,Mode>  modemap{{"sim", sim}, {"inf", inf}, {"multisim", multisim}, {"abcsmc", abcsmc}, {"abcmbp", abcmbp}, {"combinetrace", combinetrace}};
	if (modemap.count(val) != 0) mode = modemap[val];
	else emsgroot("Unrecoginsed value " + val + " for mode parameter");
	
	return mode;
}

/// Finds and returns 'transdata'
vector <TRANSDATA> Inputs::find_transdata(const Details &details) const
{
	vector <TRANSDATA> transdatavec;
	
	if(basedata->contains("transdata")) {
		const auto tdata = basedata->open("transdata");

		for(unsigned int j = 0; j < tdata.size(); j++){
			const auto td = tdata[j];
		
			TRANSDATA transdata;
			transdata.trans = UNSET;
			transdata.fromstr = td.stringfield("from");
			transdata.tostr = td.stringfield("to");
			transdata.type = td.stringfield("area");

			if(transdata.type != "reg" && transdata.type != "all") emsgroot("Transition data type not recognised"); 
			
			transdata.file = td.stringfield("file");

			const auto startdata = td.stringfield("start");
			transdata.start = details.gettime(startdata)-details.start;
			
			const auto units = td.stringfield("units");
			if(units == "days") transdata.units = 1;
			else{
				if(units == "weeks") transdata.units = 7;
				else emsgroot("Units in 'transdata' not recognised");
			}
			
			if(details.mode != inf && details.mode != abcsmc && details.mode != abcmbp){
				transdata.rows = (unsigned int)((details.period - transdata.start)/transdata.units);
				if(transdata.rows == 0) emsgroot("Transition data '"+transdata.file+"' cannot be generated because the time period is not sufficiently long.");
			} else {
				transdata.rows = UNSET;
			}
			
			transdatavec.push_back(transdata);
		}
	}
	return transdatavec;
}

/// Finds and returns 'popdata'
vector <POPDATA> Inputs::find_popdata(const Details &details) const
{
	vector <POPDATA> popdatavec;

	if(basedata->contains("popdata")) {
		const auto pdata = basedata->open("popdata");

		for(unsigned int j = 0; j < pdata.size(); j++){
			const auto pd = pdata[j];

			POPDATA popdata;
			popdata.comp = UNSET;
			popdata.compstr = pd.stringfield("comp");
			popdata.type = pd.stringfield("area");

			if(popdata.type != "reg" && popdata.type != "all") emsgroot("popition data type not recognised"); 
			
			popdata.file = pd.stringfield("file");

			const auto startdata = pd.stringfield("start");
			popdata.start = details.gettime(startdata)-details.start;
			
			const auto units = pd.stringfield("units");
			if(units == "days") popdata.units = 1;
			else{
				if(units == "weeks") popdata.units = 7;
				else emsgroot("Units in 'popdata' not recognised");
			}
			
			if(details.mode != inf && details.mode != abcsmc && details.mode != abcmbp){
				popdata.rows = (unsigned int)((details.period - popdata.start)/popdata.units);
				if(popdata.rows == 0) emsgroot("popition data '"+popdata.file+"' cannot be generated because the time period is not sufficiently long.");
			} else {
				popdata.rows = UNSET;
			}
			
			popdatavec.push_back(popdata);
		}
	}
	
	return popdatavec;
}

/// Finds and returns 'margdata'
vector <MARGDATA> Inputs::find_margdata(const Details & /*details */, const vector <DEMOCAT> &democat) const
{
	vector <MARGDATA> margdatavec;
	
	if(basedata->contains("margdata")) {
		const auto mdata = basedata->open("margdata");

		for(unsigned int j = 0; j < mdata.size(); j++){
			const auto md = mdata[j];
			
			MARGDATA margdata;
			margdata.trans = UNSET;
			margdata.fromstr = md.stringfield("from");
			margdata.tostr = md.stringfield("to");
			
			margdata.type = md.stringfield("area");
			if(margdata.type != "reg" && margdata.type != "all") emsgroot("Marginal data type not recognised"); 
			
			const auto type = md.stringfield("type");
			unsigned int d;
			for(d = 0; d < democat.size(); d++) if(type == democat[d].name) break;
			if(d == democat.size()) emsgroot("The 'type' property must be 'age' or a demographic property.");
			margdata.democat = d;
				
			margdata.file = md.stringfield("file");
			
			margdatavec.push_back(margdata);
		}
	}
	
	return margdatavec;
}

/// Finds and returns 'democats'
vector <DEMOCAT> Inputs::find_democat(const Details & /* details */) const
{
	vector <DEMOCAT> democatvec;
	
	if(basedata->contains("ages")){                           // Age categories
		const auto ages = basedata->open("ages");
		
		DEMOCAT democat;
		democat.name = "age";
		for(unsigned int j = 0; j < ages.size(); j++){
			const auto ag = ages[j];
			
			const auto range = ag.stringfield("range");
			democat.value.push_back(range);
			
			const auto sus = ag.stringfield("sus");
			democat.param.push_back(sus);
		}
		democatvec.push_back(democat);
	}
	else emsgroot("The 'ages' parameter must be set.");
	
	if(basedata->contains("democats")){                        // Other demographic possibilities
		const auto democats = basedata->open("democats");
	
		for(unsigned int k = 0; k < democats.size(); k++){
			const auto democ = democats[k];
			
			DEMOCAT democat;
			democat.name="";
			for(unsigned int j = 0; j < democ.size(); j++){
				const auto demoval = democ[j];
				
				const auto value = demoval.stringfield("value");
				democat.value.push_back(value);
				
				const auto sus = demoval.stringfield("sus");
				democat.param.push_back(sus);
			}
			
			democatvec.push_back(democat);
		}
	}
	
	return democatvec;
}

/// Finds and returns 'covar'
vector <COVAR> Inputs::find_covar(const Details &/*details*/) const
{
	vector <COVAR> covarvec;
	
	if(basedata->contains("covars")){
		const auto covars = basedata->open("covars");
		
		COVAR cov;
		for(unsigned int j = 0; j < covars.size(); j++){
			const auto covar = covars[j];
			
			cov.name = covar.stringfield("name");
			cov.param = covar.stringfield("param");
			cov.func = covar.stringfield("func");
			cov.col = UNSET;
			
			covarvec.push_back(cov);
		}
	}
	
	return covarvec;
}

/// Finds and returns 'timep'
vector <TIMEP> Inputs::find_timeperiod(const Details &details) const
{
	vector <TIMEP> timeperiodvec;
	
	if(basedata->contains("timep")) {
		const auto timep = basedata->open("timep");
		for(unsigned int j = 0; j < timep.size(); j++){
			const auto tim = timep[j];
			
			TIMEP timeperiod;
			timeperiod.name = tim.stringfield("name");
			
			auto tendstr = tim.stringfield("tend");
			timeperiod.tend = details.gettime(tendstr) - details.start;
			
			if(timeperiod.tend < 0 || timeperiod.tend > (int)details.period) emsgroot("Time '"+tendstr+"' is out of range."); 
			if(j > 0){
				if(timeperiod.tend < timeperiodvec[j-1].tend) emsgroot("'timep' is not time ordered.");
			}
			
			if(j == timep.size()-1){
				if(timeperiod.tend != (int)details.period) emsgroot("'tend' in 'timep' must finish with the end time.");
			}
			timeperiodvec.push_back(timeperiod);
		}
	}
	else emsgroot("Property 'timep' defining time periods must be set.");

	return timeperiodvec;
}

/// Finds properties of 'genQ'
void Inputs::find_genQ(GENQ &genQ, const Details &/*details*/) const
{
	if(basedata->contains("agemix")) {
		const auto agemix = basedata->open("agemix");
		
		genQ.Nall = agemix.stringfield("Nall");
		genQ.Nhome = agemix.stringfield("Nhome");
		genQ.Nother = agemix.stringfield("Nother");
		genQ.Nschool = agemix.stringfield("Nschool");
		genQ.Nwork = agemix.stringfield("Nwork");
	}
	else emsgroot("'agemix' must be specified.");

	if(basedata->contains("geomix")) {
		const auto geomix = basedata->open("geomix");		
		genQ.M = geomix.stringfield("M");
	}
	else emsgroot("'geomix' must be specified.");
	
	if(basedata->contains("genQoutput")) {
		const auto qout = basedata->open("genQoutput");
		
		genQ.localhome = qout.stringfield("localhome");
		genQ.flowall = qout.stringfield("flowall");
	}
	else emsgroot("'genQoutput' must be specified.");
}

/// Sets up Q
void Inputs::find_Q(vector <QTENSOR> &Qvec, const vector <TIMEP> &timeperiod, const Details & /* details */) const
{
	if(basedata->contains("Q")) {
		const auto Qlist = basedata->open("Q");
		for(unsigned int j = 0; j < Qlist.size(); j++){
			QTENSOR qten;
		
			const auto Q = Qlist[j];
			
			const auto timep = Q.stringfield("timep");
			unsigned int tp = 0;
			while(tp < timeperiod.size() && timeperiod[tp].name != timep)
				tp++;
			if(tp == timeperiod.size()) emsgroot("Cannot find '"+timep+"' as a time period defined using the 'timep' command in the input TOML file.");
			qten.timep = tp;
	
			qten.comp = Q.stringfield("comp");			
			qten.name = Q.stringfield("name");
		
		  qten.Qtenref = UNSET;
			
			Qvec.push_back(qten);
		}
	}
	else emsgroot("Property 'timep' defining time periods must be set.");
}

/// Finds 'params'
void Inputs::find_param(vector <string> &name, vector <double> &val) const
{
	if(basedata->contains("params")){
		const auto paramsin = basedata->open("params");
		for(unsigned int j = 0; j < paramsin.size(); j++){
			const auto params = paramsin[j];
			string nam = params.stringfield("name");
			
			double value = params.numberfield("value");
			
			name.push_back(nam);
			val.push_back(value);
		}
	}
	else{ emsgroot("The input file must contain parameter values through 'params'.");}
}

/// Finds 'priors'
void Inputs::find_prior(vector <string> &name, vector <double> &min, vector <double> &max) const
{
	if(basedata->contains("priors")){
		const auto paramsin = basedata->open("priors");
		for(unsigned int j = 0; j < paramsin.size(); j++){
			const auto params = paramsin[j];
			string nam = params.stringfield("name");

			double mi, ma;
			if(params.contains("value")){
				double value = params.numberfield("value");
				mi = value; ma = value;
			}
			else{
				if(!params.contains("type")) emsgroot("The prior '"+nam+"' must have a 'value' or a 'type'");
				
				string type = params.stringfield_unchecked("type");
				if(type == "uniform"){
					mi = params.numberfield("min");
					ma = params.numberfield("max");
				}
				else emsgroot("In 'priors', the prior type '"+type+"' is not recognised.");
			}
			
			name.push_back(nam); min.push_back(mi); max.push_back(ma);
		}
	}
	else{ emsgroot("The input file must contain quantity 'priors'.");}
}

/// Finds 'comps'
void Inputs::find_comps(vector <string> &name, vector <double> &infectivity) const
{
	if(basedata->contains("comps")) {
		const auto compsin = basedata->open("comps");
		for(unsigned int j = 0; j < compsin.size(); j++){
			const auto comps = compsin[j];

			string nam = comps.stringfield("name");

			double infect = comps.numberfield("inf");
	
			name.push_back(nam); infectivity.push_back(infect);
		}
	}
	else{ emsgroot("The input file must contain compartment definitions through 'comps'");}
}

/// Finds 'trans'
void Inputs::find_trans(vector <string> &from, vector <string> &to, vector <string> &prpar, vector <int> &type, vector <string> &mean, vector <string> &cv) const
{
	if(basedata->contains("trans")){
		const auto transin = basedata->open("trans");
		for(unsigned int j = 0; j < transin.size(); j++){
			const auto trans = transin[j];
			
			string fr_temp = trans.stringfield("from");
			
			string to_temp = trans.stringfield("to");
		
			string name = fr_temp+"→"+to_temp;
			if(!trans.contains("dist")) emsgroot("For the '"+name+"' transition the 'dist' distribution must be set.");
			
			string dist = trans.stringfield_unchecked("dist");
			
			unsigned int distval = UNSET;
			string mean_temp="", cv_temp="";
			
			if(dist == "infection"){
				distval = infection_dist;
			}
			
			if(dist == "exp"){
				distval = exp_dist;
				mean_temp = trans.stringfield("mean");				
			}
			
			if(dist == "lognorm"){
				distval = lognorm_dist;
				mean_temp = trans.stringfield("mean");				
				cv_temp = trans.stringfield("cv");
			}
			
			if(dist == "gamma"){
				distval = gamma_dist;
				mean_temp = trans.stringfield("mean");
				cv_temp = trans.stringfield("cv");
			}
			
			if(distval == UNSET) emsgroot("For the '"+name+"' transition the distribution '"+dist+"' is not recognised.");
			
			string prob="";
			if(trans.contains("prob"))
				prob = trans.stringfield_unchecked("prob");

			from.push_back(fr_temp); to.push_back(to_temp); prpar.push_back(prob);
			type.push_back(distval); mean.push_back(mean_temp); cv.push_back(cv_temp);
		}
	}
	else{ emsgroot("The input file must contain transition definitions through the 'trans' quantity.");}
}

/// Finds 'priorcomps'
vector <PRIORCOMP> Inputs::find_priorcomps(const vector<COMP> &comp) const
{
	vector <PRIORCOMP> priorcompvec;
	
	if(basedata->contains("priorcomps")){
		const auto prcomps = basedata->open("priorcomps");
		for(unsigned int j = 0; j < prcomps.size(); j++){
			const auto prcomp = prcomps[j];
			
			PRIORCOMP pricomp;
			string co = prcomp.stringfield("comp");
			unsigned int c = 0; while(c < comp.size() && comp[c].name != co) c++;
			if(c == comp.size()) emsgroot("Cannot find '"+co+"' in 'priorcomps'");
			pricomp.comp = c;
			
			pricomp.value = prcomp.numberfield("inf");
			pricomp.sd = prcomp.numberfield("sd");
	
			priorcompvec.push_back(pricomp);
		}
	}
	
	return priorcompvec;
}

/// Finds 'betaspline' and 'phispline'
void Inputs::find_spline(const Details &details, string &name, vector <int> &time, vector <string> &param) const
{
	time.clear(); param.clear();
	if(basedata->contains(name)) {
		const auto bespin = basedata->open(name);
		for(unsigned int j = 0; j < bespin.size(); j++){
			const auto besp = bespin[j];

			const auto nam = besp.stringfield("param");
			const auto timstr = besp.stringfield("time");
			
			int tim = details.gettime(timstr) - details.start;
			
			if(j == 0 && tim != 0) emsgroot("The first point in '"+name+"' must be at the 'start' time.");
			if(j == bespin.size()-1 && tim != (int)details.period) emsgroot("The last '"+name+"' point must be at the 'end' time.");
			if(tim < 0 || tim > (int)details.period) emsgroot("The '"+name+"' points must be within the time period set for the simulation/inference.");
			
			time.push_back(tim);
			param.push_back(nam);
		}
	}
	else emsgroot("'"+name+"' must be specified");
}
	
