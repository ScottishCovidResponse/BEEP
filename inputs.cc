
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
#include "model.hh"

class InputData {
public:
	typedef toml::basic_value<toml::discard_comments,
														std::unordered_map, std::vector> Node;
	InputData(const std::string& inputfilename);
	InputData(const InputData& data) = delete;
	InputData(InputData&& data) = delete;
	InputData& operator=(const InputData& data) = delete;
	InputData& operator=(InputData&& data) = delete;
	~InputData() {}
	bool contains(const std::string& name) const
		{
			return tomldata.contains(name);
		}
	Node open(const std::string& name)
		{
			return toml::find(tomldata, name);
		}
	vector<string> get_keys() const;
	Node tomldata;// Information from the TOML file
};

InputData::InputData(const std::string& inputfilename) :
	tomldata(toml::parse(inputfilename))
{
	// Allow using values from another TOML file as a base for this one. TODO:
	// make this into functions so you can do this recursively.
	if (contains("baseinputfile")) {
		const string basefile = toml::find<string>(tomldata,"baseinputfile");
		// TODO: make the filename relative to the original TOML file
		decltype(toml::parse(basefile)) basetomlddata = toml::parse(basefile);
		
		for(const auto& p : basetomlddata.as_table())
		{
			if (!contains(p.first)) {
				tomldata[p.first] = p.second;
			}
		}
	}
}

Inputs::~Inputs()
{
	delete basedata;
}

std::string stringfield(
	const InputData::Node& td,
	const char *title,
	const char *name)
{
	if(!td.contains(name)) {
		ostringstream oss;
		oss << "A '" << name <<
			"' property must be specified in '" << title << "'.";
		emsgroot(oss.str().c_str());
	}
	return toml::find<std::string>(td,name);
}

const InputData::Node& opennamedtable(
	const InputData::Node& t, const char* name)
{
	return toml::find(t,name);
}
const InputData::Node& openindexedtable(
	const InputData::Node& t, unsigned int index)
{
	return toml::find(t,index);
}

/// /// Reads TOML and command line parameters
Inputs::Inputs(int argc, char** argv, bool verbose) 
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
	
	vector<string> tomlkeys = basedata->get_keys();
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
		if(basedata->contains(key)) val = toml::find<string>(basedata->tomldata,key);
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
			val = toml::find<int>(basedata->tomldata,key);
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
			if (what == "stoi") {
				if (valstr == "")
					oss << "Should be number, found no value\n";
			} else {
				oss << what;
			}
			emsgroot(oss.str());
		}
	} else {
		if (basedata->contains(key)) {
			const auto val_temp = toml::find(basedata->tomldata,key);
			if(val_temp.is_floating()) val = val_temp.as_floating(); else val = val_temp.as_integer();	
		} else {
			val = def;
		}
	}
	
	return val;
}

/// Gets a list of all the keys
vector<string> InputData::get_keys( ) const
{
	vector<string> keys;
	for(const auto& p : tomldata.as_table())
	{
		keys.push_back(p.first);
	}
	return keys;
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
	map<string,Mode>  modemap{{"sim", sim}, {"inf", inf}, {"multisim", multisim}, {"combinetace", combinetrace}};
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
			const auto td = openindexedtable(tdata,j);
		
			TRANSDATA transdata;
			transdata.fromstr = stringfield(td,"transdata","from");
			transdata.tostr = stringfield(td,"transdata","to");
			transdata.type = stringfield(td,"transdata","area");

			if(transdata.type != "reg" && transdata.type != "all") emsgroot("Transition data type not recognised"); 
			
			transdata.file = stringfield(td,"transdata","file");

			const auto startdata = stringfield(td,"transdata","start");
			transdata.start = details.gettime(startdata)-details.start;
			
			const auto units = stringfield(td,"transdata","units");
			if(units == "days") transdata.units = 1;
			else{
				if(units == "weeks") transdata.units = 7;
				else emsgroot("Units in 'transdata' not recognised");
			}
			
			if(details.mode != inf){
				transdata.rows = (unsigned int)((details.period - transdata.start)/transdata.units);
				if(transdata.rows == 0) emsgroot("Transition data '"+transdata.file+"' cannot be generated because the time period is not sufficiently long.");
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

		POPDATA popdata;
		for(unsigned int j = 0; j < pdata.size(); j++){
			const auto pd = toml::find(pdata,j);
			
			popdata.compstr = stringfield(pd,"popdata","comp");
			popdata.type = stringfield(pd,"popdata","area");

			if(popdata.type != "reg" && popdata.type != "all") emsgroot("popition data type not recognised"); 
			
			popdata.file = stringfield(pd,"popdata","file");

			const auto startdata = stringfield(pd,"popdata","start");
			popdata.start = details.gettime(startdata)-details.start;
			
			const auto units = stringfield(pd,"popdata","units");;			
			if(units == "days") popdata.units = 1;
			else{
				if(units == "weeks") popdata.units = 7;
				else emsgroot("Units in 'popdata' not recognised");
			}
			
			if(details.mode != inf){
				popdata.rows = (unsigned int)((details.period - popdata.start)/popdata.units);
				if(popdata.rows == 0) emsgroot("popition data '"+popdata.file+"' cannot be generated because the time period is not sufficiently long.");
			}
			
			popdatavec.push_back(popdata);
		}
	}
	
	return popdatavec;
}

/// Finds and returns 'margdata'
vector <MARGDATA> Inputs::find_margdata(const Details &details, const vector <DEMOCAT> &democat) const
{
	vector <MARGDATA> margdatavec;
	
	if(basedata->contains("margdata")) {
		const auto mdata = basedata->open("margdata");

		for(unsigned int j = 0; j < mdata.size(); j++){
			const auto md = toml::find(mdata,j);
			
			MARGDATA margdata;
			margdata.fromstr = stringfield(md,"margdata","from");
		
			margdata.tostr = stringfield(md,"margdata","to");
			
			margdata.type = stringfield(md,"margdata","area");
			if(margdata.type != "reg" && margdata.type != "all") emsgroot("Marginal data type not recognised"); 
			
			const auto type = stringfield(md,"margdata","type");
			unsigned int d;
			for(d = 0; d < democat.size(); d++) if(type == democat[d].name) break;
			if(d == democat.size()) emsgroot("The 'type' property must be 'age' or a demographic property.");
			margdata.democat = d;
				
			margdata.file = stringfield(md,"margdata","file");
			
			margdatavec.push_back(margdata);
		}
	}
	
	return margdatavec;
}

/// Finds and returns 'democats'
vector <DEMOCAT> Inputs::find_democat(const Details &details) const
{
	vector <DEMOCAT> democatvec;
	
	if(basedata->contains("ages")){                           // Age categories
		const auto ages = basedata->open("ages");
		
		DEMOCAT democat;
		democat.name = "age";
		for(unsigned int j = 0; j < ages.size(); j++){
			const auto ag = toml::find(ages,j);
			
			const auto range = stringfield(ag,"ages","range");
			democat.value.push_back(range);
			
			const auto sus = stringfield(ag,"ages","sus");
			democat.param.push_back(sus);
		}
		democatvec.push_back(democat);
	}
	else emsgroot("The 'ages' parameter must be set.");
	
	if(basedata->contains("democats")){                        // Other demographic possibilities
		const auto democats = basedata->open("democats");
	
		for(unsigned int k = 0; k < democats.size(); k++){
			const auto democ = toml::find(democats,k);
			
			DEMOCAT democat;
			democat.name="";
			for(unsigned int j = 0; j < democ.size(); j++){
				const auto demoval = toml::find(democ,j);
				
				const auto value = stringfield(demoval,"democats","value");
				democat.value.push_back(value);
				
				const auto sus = stringfield(demoval,"democats","sus");
				democat.param.push_back(sus);
			}
			
			democatvec.push_back(democat);
		}
	}
	
	return democatvec;
}

/// Finds and returns 'covar'
vector <COVAR> Inputs::find_covar(const Details &details) const
{
	vector <COVAR> covarvec;
	
	if(basedata->contains("covars")){
		const auto covars = basedata->open("covars");
		
		COVAR cov;
		for(unsigned int j = 0; j < covars.size(); j++){
			const auto covar = toml::find(covars,j);
			
			cov.name = stringfield(covar,"covars","name");
			cov.param = stringfield(covar,"covars","param");
			cov.func = stringfield(covar,"covars","func");
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
			const auto tim = toml::find(timep,j);
			
			TIMEP timeperiod;
			timeperiod.name = stringfield(tim,"timep","name");
			
			auto tendstr = stringfield(tim,"timep","tend");
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
void Inputs::find_genQ(GENQ &genQ, const Details &details) const
{
	if(basedata->contains("agemix")) {
		const auto agemix = basedata->open("agemix");
		
		genQ.Nall = stringfield(agemix,"agemix","Nall");
		genQ.Nhome = stringfield(agemix,"agemix","Nhome");
		genQ.Nother = stringfield(agemix,"agemix","Nother");
		genQ.Nschool = stringfield(agemix,"agemix","Nschool");
		genQ.Nwork = stringfield(agemix,"agemix","Nwork");
	}
	else emsgroot("'agemix' must be specified.");

	if(basedata->contains("geomix")) {
		const auto geomix = basedata->open("geomix");		
		genQ.M = stringfield(geomix,"geomix","M");
	}
	else emsgroot("'geomix' must be specified.");
	
	if(basedata->contains("genQoutput")) {
		const auto qout = basedata->open("genQoutput");
		
		genQ.localhome = stringfield(qout,"genQoutput","localhome");
		genQ.flowall = stringfield(qout,"genQoutput","flowall");
	}
	else emsgroot("'genQoutput' must be specified.");
}

/// Sets up Q
void Inputs::find_Q(vector <QTENSOR> &Qvec, const vector <TIMEP> &timeperiod, const Details &details) const
{
	if(basedata->contains("Q")) {
		const auto Qlist = basedata->open("Q");
		for(unsigned int j = 0; j < Qlist.size(); j++){
			QTENSOR qten;
		
			const auto Q = toml::find(Qlist,j);
			
			const auto timep = stringfield(Q,"Q","timep");
			unsigned int tp = 0;
			while(tp < timeperiod.size() && timeperiod[tp].name != timep)
				tp++;
			if(tp == timeperiod.size()) emsgroot("Cannot find '"+timep+"' as a time period defined using the 'timep' command in the input TOML file.");
			qten.timep = tp;
	
			qten.comp = stringfield(Q,"Q","comp");			
			qten.name = stringfield(Q,"Q","name");
		
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
			const auto params = toml::find(paramsin,j);
			if(!params.contains("name")) emsgroot("The quantity 'params' must contain a 'name' definition.");
			string nam = stringfield(params,"params","name");
			
			if(!params.contains("value")) emsgroot("The quantity 'params' must contain a 'value' definition.");
			const auto value_temp = toml::find(params,"value");
			double value;
			if(value_temp.is_floating()) value = value_temp.as_floating(); else value = value_temp.as_integer();
			
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
			const auto params = toml::find(paramsin,j);
			string nam = stringfield(params,"priors","name");

			double mi, ma;
			if(params.contains("value")){
				const auto value_temp = toml::find(params,"value");
				double value;
				if(value_temp.is_floating()) value = value_temp.as_floating(); else value = value_temp.as_integer();	
				mi = value; ma = value;
			}
			else{
				if(!params.contains("type")) emsgroot("The prior '"+nam+"' must have a 'value' or a 'type'");
				
				string type = toml::find<std::string>(params,"type");
				if(type == "uniform"){
					if(!params.contains("min")) emsgroot("For the prior '"+nam+"', the uniform distribution must contain a 'min' definition.");
					const auto min_temp = toml::find(params,"min");
					if(min_temp.is_floating()) mi = min_temp.as_floating(); else mi = min_temp.as_integer();	
		
					if(!params.contains("max")) emsgroot("For the prior '"+nam+"', the uniform distribution must contain a 'max' definition.");
					const auto max_temp = toml::find(params,"max");
					if(max_temp.is_floating()) ma = max_temp.as_floating(); else ma = max_temp.as_integer();	
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
			const toml::value comps = toml::find(compsin,j);

			string nam = stringfield(comps,"comps","name");

			if(!comps.contains("inf")) emsgroot("The compartment 'name' must contain an 'inf' definition.");
			
			const auto inf_temp = toml::find(comps,"inf");
			double infect;
			if(inf_temp.is_floating()) infect = inf_temp.as_floating(); else infect = inf_temp.as_integer();	
	
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
			const auto trans = toml::find(transin,j);
			
			string fr_temp = stringfield(trans, "trans", "from");
			
			string to_temp = stringfield(trans, "trans", "to");
		
			string name = fr_temp+"→"+to_temp;
			if(!trans.contains("dist")) emsgroot("For the '"+name+"' transition the 'dist' distribution must be set.");
			
			string dist = toml::find<std::string>(trans, "dist");
			
			unsigned int distval = UNSET;
			string mean_temp="", cv_temp="";
			
			if(dist == "infection"){
				distval = infection_dist;
			}
			
			if(dist == "exp"){
				distval = exp_dist;
				mean_temp = stringfield(trans,name.c_str(),"mean");				
			}
			
			if(dist == "lognorm"){
				distval = lognorm_dist;
				mean_temp = stringfield(trans,name.c_str(),"mean");				
				cv_temp = stringfield(trans,name.c_str(),"cv");
			}
			
			if(dist == "gamma"){
				distval = gamma_dist;
				mean_temp = stringfield(trans,name.c_str(),"mean");
				cv_temp = stringfield(trans,name.c_str(),"cv");
			}
			
			if(distval == UNSET) emsgroot("For the '"+name+"' transition the distribution '"+dist+"' is not recognised.");
			
			string prob="";
			if(trans.contains("prob")) prob = toml::find<std::string>(trans, "prob");

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
			const auto prcomp = toml::find(prcomps,j);
			
			PRIORCOMP pricomp;
			string co = stringfield(prcomp,"priorcomps","comp");
			unsigned int c = 0; while(c < comp.size() && comp[c].name != co) c++;
			if(c == comp.size()) emsgroot("Cannot find '"+co+"' in 'priorcomps'");
			pricomp.comp = c;
			
			if(!prcomp.contains("value")) emsgroot("'priorcomps' must contain a 'value' definition.");
			double val;
			const auto val_temp = toml::find(prcomp,"inf");
			if(val_temp.is_floating()) val = val_temp.as_floating(); else val = val_temp.as_integer();	
	
			pricomp.value = val;
			
			if(!prcomp.contains("sd")) emsgroot("'priorcomps' must contain a 'sd' standard deviation definition.");
			double sd;
			const auto sd_temp = toml::find(prcomp,"sd");
			if(sd_temp.is_floating()) sd = sd_temp.as_floating(); else sd = sd_temp.as_integer();	
			pricomp.sd = sd;
	
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
			const auto besp = toml::find(bespin,j);

			const auto nam = stringfield(besp,name.c_str(),"param");
			const auto timstr = stringfield(besp,name.c_str(),"time");
			
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
	
