// This file deals with loading quantities from the TOML input file.

#include <algorithm>
#include <iostream>
#include <fstream> 
#include <map>
#include <vector>
#include <string>
#include <sstream>

using namespace std;

#include "inputs.hh"
#include "details.hh"

Inputs::~Inputs()
{
	delete basedata;
}

/// Reads TOML and command line parameters
Inputs::Inputs(int argc, char** argv) 
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
		auto n = commandlist.size();
		if(n > 0 && (str.substr(0,1) == "=" || commandlist[n-1].substr(commandlist[n-1].length()-1,1) == "=")) commandlist[n-1] += str;
		else{
			commandlist.push_back(str);
		}
	}
	
	// Store the parameters passed on the command line in cmdlineparams
	for(const auto& str : commandlist ){ // Goes the various input options
		auto j = 0u; auto jmax = str.length(); while(j < jmax && str.substr(j,1) != "=") j++;
		if(j == jmax){
			stringstream ss; ss << "Cannot understand " << str; 
			emsgroot(ss.str());
		}
		
		string command = str.substr(0,j);
		string value = str.substr(j+1,jmax-(j+1));
		
		if (cmdlineparams.count(command) == 0) cmdlineparams[command] = value;
		else cmdlineparams[command] += " " + value;
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
int Inputs::find_integer(const string &key, int def) const
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
				ostringstream oss;
				oss << "Should be number, found '" << valstr;
				throw std::invalid_argument(oss.str());
			}
		} catch (const std::exception& e) {
			ostringstream oss;
			oss << "Bad command-line parameter value for key '"<< key <<"'\n";
			// Add exception description if it's informative
			std::string what = e.what();
			if(what == "stof"){
				if (valstr == "") oss << "Should be number, found no value\n";
			} 
			else oss << what;
			emsgroot(oss.str());
		}
	} 
	else {
		if(basedata->contains(key)) val = basedata->data.numberfield_unchecked(key);
		else val = def;
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
	map<string,Mode>  modemap{{"sim", SIM}, {"multisim", MULTISIM}, {"counter", COUNTER}, {"ppc", PPC}, {"abc", ABC_SIMPLE}, {"abcsmc", ABC_SMC}, {"abcmbp", ABC_MBP}, {"abcmbp_gr", ABC_MBP_GR}, {"mc3", MC3_INF}, {"pais", PAIS_INF}, {"pmcmc", PMCMC_INF}};
	if (modemap.count(val) != 0) mode = modemap[val];
	else emsgroot("Unrecoginsed value " + val + " for mode parameter");
	
	return mode;
}


/// Returns whether simulation or inference is being performed
SimInf Inputs::get_siminf() const
{
	switch(mode()){                                                                      // Determines if inference or simulation
	case SIM: case MULTISIM: case COUNTER: case PPC: return SIMULATE;
	default: return INFERENCE;
	}
	return INFERENCE;
}
	
	
/// Finds information about datatables
vector <DataTable> Inputs::find_datatable(const Details &details) const
{
	vector <DataTable> datatable;
	
	if(basedata->contains("data_tables")) {
		const auto datatables = basedata->open("data_tables");

		for(auto j = 0u; j < datatables.size(); j++){
			const auto td = datatables[j];
		
			DataTable datatab;
			datatab.optype = DATA;
			
			auto type =  td.stringfield("type");
			if(type == "population") datatab.type = POP;
			else{
				if(type == "transition") datatab.type = TRANS;
				else{
					if(type == "marginal") datatab.type = MARGINAL;
					else emsg("Data type '"+type+"' not recognised");
				}
			}
			
			datatab.observation = td.stringfield("observation");
			
			datatab.start = 0;
			if(td.contains("start")){
				const auto startdata = td.stringfield("start");
				datatab.start = details.gettime(startdata)-details.start;
			}
			
			datatab.fraction_spline = "";
			if(td.contains("fraction")) datatab.fraction_spline = td.stringfield("fraction");
			
			datatab.democats = "all"; 
			if(td.contains("democats")) datatab.democats = td.stringfield("democats");
			
			datatab.file = td.stringfield("file");

			datatab.weight = 1;
			if(td.contains("weight")) datatab.weight = td.numberfield("weight");
			
			datatab.file = td.stringfield("file");
			
			switch(datatab.type){
				case POP: case TRANS:
					{
						datatab.geo = "all";
						if(td.contains("geo")) datatab.geo = td.stringfield("geo");
						
						datatab.shift = 0;
						if(td.contains("shift")) datatab.shift = atoi(td.stringfield("shift").c_str());
						datatab.start += datatab.shift;
						
						datatab.units = 1;
						if(td.contains("units")){
							const auto units = td.stringfield("units");
							if(units == "days") datatab.units = 1;
							else{
								if(units == "weeks") datatab.units = 7;
								else emsgroot("Units in 'datatables' not recognised");
							}
						}
					}
					break;
					
				case MARGINAL:	
					{				
						datatab.units = UNSET;
						datatab.shift = UNSET;
				
						datatab.geo = "all";
						if(td.contains("geo")) datatab.geo = td.stringfield("geo");
						
						const auto enddata = td.stringfield("end");
						datatab.end = details.gettime(enddata)-details.start;
					}
					break;
			}
	
			datatab.invT = 1;	
			if(td.contains("invT")) datatab.invT = atof(td.stringfield("invT").c_str());
	
			datatab.obsmodel = SCALE_OBSMODEL;
			if(td.contains("obsmodel")){
				auto obsmodel = td.stringfield("obsmodel");
				if(obsmodel == "normal") datatab.obsmodel = NORMAL_OBSMODEL;
				else{
					if(obsmodel == "poisson") datatab.obsmodel = POISSON_OBSMODEL;
					else{
						if(obsmodel == "negbin"){ 
							datatab.obsmodel = NEGBINO_OBSMODEL; 
							datatab.shape = atof(td.stringfield("shape").c_str());
						}
						else{
							if(obsmodel == "scale") datatab.obsmodel = SCALE_OBSMODEL;
							else emsg("obsmodel='"+obsmodel+"' is not recognised.");
						}
					}
				}			
			}
			
			datatable.push_back(datatab);
		}
	}
	
	if(basedata->contains("state_outputs")) {
		const auto datatables = basedata->open("state_outputs");

		for(auto j = 0u; j < datatables.size(); j++){
			const auto td = datatables[j];
		
			DataTable datatab;
			datatab.optype = OUTPUT;
				
			auto type =  td.stringfield("type");
			if(type == "population") datatab.type = POP;
			else{
				if(type == "transition") datatab.type = TRANS;
				else{
					if(type == "marginal") datatab.type = MARGINAL;
					else emsg("Data type '"+type+"' not recognised");
				}
			}
			
			datatab.observation = td.stringfield("observation");
			
			if(td.contains("start")) emsg("'state_outputs' should not set a start.");
			datatab.start = 0;
			
			if(td.contains("fraction")) emsg("'state_outputs' should not set a fraction.");
			datatab.fraction_spline = "";
			
			datatab.democats = "all"; 
			if(td.contains("democats")) datatab.democats = td.stringfield("democats");
			
			if(td.contains("file")) emsg("'state_outputs' should not set a file.");
			datatab.file = "";

			if(td.contains("weight")) emsg("'state_outputs' should not set a weight.");
			datatab.weight = 1;
		
			switch(datatab.type){
				case POP: case TRANS:
					{
						datatab.geo = td.stringfield("geo");
						
						datatab.shift = 0;
						if(td.contains("shift")) emsg("'state_outputs' should not set a shift.");
						
						if(td.contains("units")) emsg("'state_outputs' should not set units.");
						datatab.units = 1;
					}
					break;
					
				case MARGINAL:	
					{				
						datatab.units = UNSET;
						datatab.shift = UNSET;
				
						datatab.geo = "all";
						if(td.contains("geo")) datatab.geo = td.stringfield("geo");
						
						if(td.contains("units")) emsg("'state_outputs' should not set an end.");

						datatab.end = details.end-details.start;
					}
					break;
			}
	
			datatable.push_back(datatab);
		}
	}
	
	return datatable;
}


/// Finds and returns 'democats'
vector <DemographicCategory> Inputs::find_demographic_category() const
{
	vector <DemographicCategory> democatvec;
	
	DemographicCategory democat;
	democat.name = "age";
	democat.sus_vari= false;
		
	if(basedata->contains("ages")){                           // Age categories
		const auto ages = basedata->open("ages");
		
		for(auto j = 0u; j < ages.size(); j++){
			const auto ag = ages[j];
			
			const auto range = ag.stringfield("range");
			democat.value.push_back(range);
			
			if(ag.contains("sus")){
				democat.sus_vari = true;
				const auto sus = ag.stringfield("sus");
			
				ParamSpec ps;
				ps.name = sus; ps.prior = "Exp(1)";
				if(ag.contains("sus_value")) ps.value = ag.stringfield("sus_value");
				democat.ps.push_back(ps);
			}
			else{
				if(democat.sus_vari == true) emsg("A susceptibility parameter needs to be set for each age range.");
				democat.sus_vari = false;
			}
		}
	
	}
	else{
		democat.value.push_back("all");
	}
	democatvec.push_back(democat);
	
	if(basedata->contains("democats")){                        // Other demographic possibilities
		const auto democats = basedata->open("democats");
	
		for(auto k = 0u; k < democats.size(); k++){
			const auto democ = democats[k];
			
			DemographicCategory democat;
			democat.name="";
			democat.sus_vari = false;
			
			for(auto j = 0u; j < democ.size(); j++){
				const auto demoval = democ[j];
				
				if(demoval.contains("name")){
					const auto name = demoval.stringfield("name");
					if(democat.name == "") democat.name = name;
					else{
						if(democat.name != name) emsg("The name '"+democat.name+"' is inconsistent");
					}
				}
				
				if(demoval.contains("sus")){
					democat.sus_vari = true;
					const auto sus = demoval.stringfield("sus");
			
					const auto value = demoval.stringfield("sus_value");
					democat.value.push_back(value);
					
					ParamSpec ps;
					ps.name = sus; ps.prior = "Exp(1)";
					if(demoval.contains("sus_value")) ps.value = demoval.stringfield("sus_value");
					democat.ps.push_back(ps);
				}
				else{
					if(democat.sus_vari == true) emsg("A susceptibility parameter needs to be set for each age range.");
					democat.sus_vari = false;
				}
			}
			
			democatvec.push_back(democat);
		}
	}
	
	return democatvec;
}


/// Finds and returns 'covar'
vector <Covariate> Inputs::find_covariate() const
{
	vector <Covariate> covarvec;
	
	if(basedata->contains("covars")){
		const auto covars = basedata->open("covars");
		
		Covariate cov;
		for(auto j = 0u; j < covars.size(); j++){
			const auto covar = covars[j];
			
			cov.name = covar.stringfield("name");
			
			auto vec = get_paramspec("covars",j,"param","value","prior","","",1);
			 
			cov.ps = vec[0];
			
			cov.func = covar.stringfield("func");
		
			covarvec.push_back(cov);
		}
	}
	
	return covarvec;
}


/// Finds the break points from a string 
vector <unsigned int> Inputs::get_breakpoints(string st, const Details &details) const
{
	check_from_table(st);

	auto bp = split(st,'|');
	
	vector <unsigned int> vec;
	for(auto j = 0u; j < bp.size(); j++) vec.push_back(details.gettime(bp[j])-details.start); 

	return vec;	
}


/// Finds properties of 'genQ'
void Inputs::find_genQ(GenerateQ &genQ, const Details &details, const vector <string> &age_list) const
{
	genQ.nspline = 0;
	
	if(basedata->contains("age_mixing_matrix")) {
		auto st = find_string("age_mixing_matrix","");
		genQ.N_name = split(st,'|');
	}
	else{
		if(age_list.size() != 1) emsgroot("'age_mixing_matrix' must be specified.");
		else genQ.N_name.push_back("UNIT MATRIX");
	}

	if(basedata->contains("age_mixing_modify")) {
		const auto matrix_mod = basedata->open("age_mixing_modify");

		MatrixModification matmod;
		for(auto j = 0u; j < matrix_mod.size(); j++){
			const auto mm = matrix_mod[j];
		
			string type = mm.stringfield("type");
			if(type == "row_column"){
				matmod.type = ROW_COLUMN;
				
				auto nage = age_list.size();
				auto ages = split(mm.stringfield("ages"),'|');
				for(auto st : ages){
					auto a = 0u; while(a < nage && age_list[a] != st) a++;
					if(a == nage) emsg("Cannot find '"+st+"'");
					matmod.ages.push_back(a);
				}
				
				matmod.desc = "Factor multipling the "+mm.stringfield("ages")+" rows and columns of the contact matrix";
				matmod.name = mm.stringfield("ages")+" factor";
				
				matmod.bp = get_breakpoints(mm.stringfield("bp"),details);
			
				matmod.ps = get_paramspec("age_mixing_modify",j,"param","value","prior","smooth","factor",matmod.bp.size());
				
				genQ.nspline++;
			}
			else{
				emsg("The 'age_mixing_modify' type is not recognised");
			}
			
			genQ.matmod.push_back(matmod);
		}
	}

	if(basedata->contains("geo_mixing_matrix")) {
		genQ.M_name = find_string("geo_mixing_matrix","");
	}
	else{
		genQ.M_name = "";
	}
}


/// Adds a new label onto timeplot
void Inputs::timeplot_add(vector <TimePlot> &timeplot, string time_str, string name) const 
{
	auto i = 0u; while(i < timeplot.size() && timeplot[i].time_str != time_str) i++;
	if(i == timeplot.size()){
		TimePlot tp;	
		tp.name = name;
		tp.time_str = time_str;
		timeplot.push_back(tp);
	}
}


/// Finds times (used for plots)
void Inputs::find_timeplot(vector <TimePlot> &timeplot) const
{
	if(basedata->contains("times")){
		const auto times = basedata->open("times");
		for(auto j = 0u; j < times.size(); j++){
			const auto time = times[j];
			TimePlot tp;
			tp.name = time.stringfield("name");
			tp.time_str = time.stringfield("time");
			timeplot.push_back(tp);
		}
	}
	
	if(mode() == COUNTER){
		if(basedata->contains("counterfactual")){
			const auto pr = basedata->open("counterfactual");
		
			for(auto j = 0u; j < pr.size(); j++){
				const auto cf = pr[j];
		
				if(cf.contains("start")) timeplot_add(timeplot,cf.stringfield("start"),"Counterfactual");
				if(cf.contains("end")) timeplot_add(timeplot,cf.stringfield("end"),"Counterfactual");
			}
		}
	}
	
	if(mode() == PPC){
		auto start_str = find_string("ppc_start","UNSET"); 
		if(start_str != "UNSET") timeplot_add(timeplot,start_str,"PPC");
	}	
}


/// Finds information about regional effects
void Inputs::regional_effects(string& geography, ParamSpec& sigma) const
{
	if(basedata->contains("region_effect")){
		const auto re = basedata->open("region_effect");

		geography = re.stringfield("geography");
		sigma.name = re.stringfield("sigma");
		if(re.contains("sigma_value")) sigma.value = re.stringfield("sigma_value");
		if(re.contains("sigma_prior")) sigma.value = re.stringfield("sigma_prior");
	}
}	


/// Gets a parameter specification from various elements
vector <ParamSpec> Inputs::get_paramspec(string root, unsigned int j, string name, string value, string prior, string smooth, string factor, unsigned int si) const
{
	const auto ro = basedata->open(root);
	const auto it = ro[j];
	
	string name_string = it.stringfield(name.c_str());
	string sim_string = "", prior_string = "", smooth_string = "", factor_string = "1";
	
	switch(get_siminf()){
		case SIMULATE:
			if(it.contains(value.c_str())){
				sim_string = it.stringfield(value.c_str());
			}			
			else{
				emsg("The parameter '"+name+"' must have a value set.");
			}
			break;
		
		case INFERENCE:
			if(it.contains(value.c_str())){
				sim_string = it.stringfield(value.c_str());
			}			
			
			if(it.contains(prior.c_str())){
				prior_string = it.stringfield(prior.c_str());
			}
			else{
				if(!it.contains(value.c_str())){   
					emsg("The prior for '"+name+"' must be set.");
				}
			}
			break;
	}
	
	if(it.contains(smooth.c_str())){   
		smooth_string = it.stringfield(smooth.c_str());  
	}
			
	if(it.contains(factor.c_str())){   
		factor_string = it.stringfield(factor.c_str());  
	}

	auto na = split_string(name_string,'|',si);
	if(na.size() == 0) emsg("'"+root+"' does not have a name.");
	
	unsigned int num = 1;
	for(auto i = 0u; i < na.size(); i++){   // This fills in '.' when performing auto numbering
		if(na[i].substr(na[i].length()-1,1) == "."){
			stringstream ss; ss << na[i].substr(0,na[i].length()-1) << num;
			na[i] = ss.str(); num++;
		}
	}
	
	auto val = split_string(sim_string,'|',si);
	
	if(get_siminf() == INFERENCE && prior_string == "" && sim_string != ""){  // Creates fixed  prior if value is set
		for(auto i = 0u; i < si; i++){
			if(i != 0) prior_string += "|";
			prior_string += "Fixed("+val[i]+")";
		}
	}
	
	auto pri = split_string(prior_string,'|',si);
	auto smo = split_number(smooth_string,'|',si);
	auto fac = split_number(factor_string,'|',si);
	
	vector <ParamSpec> vec; 
	for(auto i=0u; i < si; i++){
		ParamSpec ps; ps.name = na[i]; ps.value = val[i]; ps.prior = pri[i]; ps.smooth = smo[i]; ps.factor = fac[i];
		vec.push_back(ps);
	}
	
	return vec;
}


/// Prints off a vector of parameter specifications
void Inputs::print_paramspec(const vector <ParamSpec> &ps) const
{
	for(auto i = 0u; i < ps.size(); i++){
		cout << i << " Name:" << ps[i].name << "  Value:" <<  ps[i].value;
		cout << " Prior:" <<  ps[i].prior << " Smooth:" <<  ps[i].smooth << " ps\n";
	}
}


/// Finds 'comps'
void Inputs::find_compartments(vector <string> &name, vector <ParamSpec> &infectivity, vector <string> &inf_trans) const
{
	if(basedata->contains("comps")) {
		const auto compsin = basedata->open("comps");
		for(auto j = 0u; j < compsin.size(); j++){
			const auto comps = compsin[j];

			string nam = comps.stringfield("name");

			ParamSpec inf; inf.name = "zero"; inf.value = "0"; inf.prior = "Fixed(0)";
			string inf_tr = "";
			
			if(comps.contains("inf")){		
				auto infvec = get_paramspec("comps",j,"inf","inf_value","inf_prior","","",1);
				inf = infvec[0];
			}
	
			if(comps.contains("inf_trans")) inf_tr = comps.stringfield("inf_trans");
			
			name.push_back(nam); 
			infectivity.push_back(inf);
			inf_trans.push_back(inf_tr);
		}
	}
	else{ emsgroot("The input file must contain compartment definitions through 'comps'");}
}


/// Finds 'trans'
void Inputs::find_transitions(vector <string> &from, vector <string> &to, vector <int> &type, vector < vector <ParamSpec> > &distspec, vector < vector <ParamSpec> > &probspec,  vector <unsigned int> &k, unsigned int ndemocatpos) const
{
	if(basedata->contains("trans")){
		const auto transin = basedata->open("trans");
		for(auto j = 0u; j < transin.size(); j++){
			
			vector <ParamSpec> distsp, probsp;
			
			const auto trans = transin[j];
			
			string fr_temp = trans.stringfield("from");
			
			string to_temp = trans.stringfield("to");
		
			string name = fr_temp+"→"+to_temp;
			if(!trans.contains("dist")) emsgroot("For the '"+name+"' transition the 'dist' distribution must be set.");
			
			string dist = trans.stringfield_unchecked("dist");
			
			unsigned int distval = UNSET;
			unsigned int k_temp = UNSET;
			
			if(dist == "Infection"){
				distval = INFECTION_DIST;
			}
			else{
				distsp = get_paramspec("trans",j,"mean","mean_value","mean_prior","","",ndemocatpos);
				if(trans.contains("prob")){
					probsp = get_paramspec("trans",j,"prob","prob_value","prob_prior","","",ndemocatpos);
				}
			}
		
			if(dist == "Exp"){
				distval = EXP_DIST;
			}
			
			if(dist == "Lognorm"){
				distval = LOGNORM_DIST;
				emsg("Lognorm not working yet");
			}
			
			if(dist == "Gamma"){
				distval = GAMMA_DIST;
				emsg("Gamma not working yet");
			}
			
			if(dist == "Erlang"){
				distval = ERLANG_DIST;
				k_temp = trans.numberfield("k");
			}
			
			if(distval == UNSET) emsgroot("For the '"+name+"' transition the distribution '"+dist+"' is not recognised.");
			
			string prob="";
			if(trans.contains("prob"))
				prob = trans.stringfield_unchecked("prob");

			from.push_back(fr_temp); to.push_back(to_temp);
			type.push_back(distval); 
			
			distspec.push_back(distsp);
			probspec.push_back(probsp);
			
			k.push_back(k_temp);
		}
	}
	else{ emsgroot("The input file must contain transition definitions through the 'trans' quantity.");}
}


/// Find information about probreach
void Inputs::find_probreach(vector <string>& name, vector <string>& comp, vector <string>& inf_trans) const
{
	if(basedata->contains("prob_reach")){
		const auto pr = basedata->open("prob_reach");
		
		for(auto num = 0u; num < pr.size(); num++){   // Allows for different spline for different areas / infectious transitions
			auto probr = pr[num];
			name.push_back(probr.stringfield("name"));
			
			comp.push_back(probr.stringfield("comp"));
			
			if(probr.contains("inf_trans")) inf_trans.push_back(probr.stringfield("inf_trans"));
			else inf_trans.push_back(probr.stringfield(""));
		}
	}
}
 
 
/// Finds spline information
void Inputs::find_spline(const string name, vector < vector <ParamSpec>>& ps_vec, vector <ParamSpec>& factor_param_vec, vector < vector <unsigned int> >& bp_vec, vector <Smooth>& smooth_type_vec, vector <string>& inft, vector <string>& area, vector < vector <double> >& efoi_agedist_vec, const vector <double>& agedist, const Details &details) const
{
	ps_vec.clear(); factor_param_vec.clear(); 
	bp_vec.clear(); smooth_type_vec.clear(); inft.clear(); area.clear(); efoi_agedist_vec.clear();
	
	if(basedata->contains(name)) {
		const auto bespin = basedata->open(name);
		
		for(auto num = 0u; num < bespin.size(); num++){   // Allows for different spline for different areas / infectious transitions
			auto besp = bespin[num];
			
			vector <unsigned int> bp;
			if(besp.contains("bp")) bp = get_breakpoints(besp.stringfield("bp"),details);
			else bp = get_breakpoints("start|end",details);

			bp_vec.push_back(bp);
		
			auto ps_list = get_paramspec(name,num,"param","value","prior","smooth","factor",bp.size());
			if(name == "efoi_spline"){ for(auto& ps : ps_list) ps.factor /= details.efoi_factor;}
			ps_vec.push_back(ps_list);
		
			if(besp.contains("factor_param")){
				auto fpps = get_paramspec(name,num,"factor_param","factor_value","factor_prior","","",1);
				factor_param_vec.push_back(fpps[0]);
			}
			else{
				ParamSpec one; one.name = "one"; one.value = "1"; one.prior = "Fixed(1)";
				factor_param_vec.push_back(one);
			}
		
			Smooth st = LOGSMOOTH;
			if(besp.contains("smooth_type")){
				string val = besp.stringfield("smooth_type");
				if(val == "log_smooth") st = LOGSMOOTH;
				else{
					if(val == "smooth") st = SMOOTH;
					else{
						if(val == "no_smooth") st = NOSMOOTH;
						else{
							emsg("The 'smooth_type' value '"+val+"' is not recognised.");
						}
					}
				}
			}
			smooth_type_vec.push_back(st);
			
			if(besp.contains("inf_trans")) inft.push_back(besp.stringfield("inf_trans")); else inft.push_back("");
	
			if(besp.contains("area")) area.push_back(besp.stringfield("area")); else area.push_back("");
		
			unsigned int nage = agedist.size();
			if(besp.contains("age_dist")){
				if(nage == 0) emsg("'age_dist' is not expected.");
			
				string aged = besp.stringfield("age_dist");
				auto ad = split_number(aged,'|',nage);
				
				auto sum = 0.0; for(auto a = 0u; a < nage; a++) sum += ad[a];
				for(auto a = 0u; a < nage; a++) ad[a] /= (sum*agedist[a]);

				if(false){
					for(auto a = 0u; a < nage; a++) cout << ad[a] << "  age dist\n"; 
					emsg("Age distribution");
				}
	
				efoi_agedist_vec.push_back(ad);    
			}
			else{
				vector <double> ad(nage);
				for(auto a = 0u; a < nage; a++) ad[a] = 1;
				efoi_agedist_vec.push_back(ad);
			}
		}
	}
}

/*
/// Finds unit spline (which is just set to 1)
void Inputs::find_unit_spline(const string name, vector < vector <ParamSpec>>& ps_vec, vector <ParamSpec>& factor_param_vec, vector < vector <unsigned int> >& bp_vec, vector <Smooth>& smooth_type_vec, vector <string>& inft, vector <string>& area, vector < vector <double> >& efoi_agedist_vec, const vector <double>& agedist, const Details &details) const
{
	ps_vec.clear(); factor_param_vec.clear(); 
	bp_vec.clear(); smooth_type_vec.clear(); inft.clear(); area.clear(); efoi_agedist_vec.clear();
				
	auto bp = get_breakpoints("start|end",details);
	bp_vec.push_back(bp);
		
			auto ps_list = get_paramspec(name,num,"param","value","prior","smooth","factor",bp.size());
			if(name == "efoi_spline"){ for(auto& ps : ps_list) ps.factor /= details.efoi_factor;}
			ps_vec.push_back(ps_list);
		
			if(besp.contains("factor_param")){
				auto fpps = get_paramspec(name,num,"factor_param","factor_value","factor_prior","","",1);
				factor_param_vec.push_back(fpps[0]);
			}
			else{
				ParamSpec one; one.name = "one"; one.value = "1"; one.prior = "Fixed(1)";
				factor_param_vec.push_back(one);
			}
		
			Smooth st = LOGSMOOTH;
			if(besp.contains("smooth_type")){
				string val = besp.stringfield("smooth_type");
				if(val == "log_smooth") st = LOGSMOOTH;
				else{
					if(val == "smooth") st = SMOOTH;
					else{
						if(val == "no_smooth") st = NOSMOOTH;
						else{
							emsg("The 'smooth_type' value '"+val+"' is not recognised.");
						}
					}
				}
			}
			smooth_type_vec.push_back(st);
			
			if(besp.contains("inf_trans")) inft.push_back(besp.stringfield("inf_trans")); else inft.push_back("");
	
			if(besp.contains("area")) area.push_back(besp.stringfield("area")); else area.push_back("");
		
			unsigned int nage = agedist.size();
			if(besp.contains("age_dist")){
				if(nage == 0) emsg("'age_dist' is not expected.");
			
				string aged = besp.stringfield("age_dist");
				auto ad = split_number(aged,'|',nage);
				
				auto sum = 0.0; for(auto a = 0u; a < nage; a++) sum += ad[a];
				for(auto a = 0u; a < nage; a++) ad[a] /= (sum*agedist[a]);

				if(false){
					for(auto a = 0u; a < nage; a++) cout << ad[a] << "  age dist\n"; 
					emsg("Age distribution");
				}
	
				efoi_agedist_vec.push_back(ad);    
			}
			else{
				vector <double> ad(nage);
				for(auto a = 0u; a < nage; a++) ad[a] = 1;
				efoi_agedist_vec.push_back(ad);
			}
		}
	}
}
*/

/// Checks if a string needs to be read from a table
void Inputs::check_from_table(string& s) const
{
	if(s.length() > 4){
		if(s.substr(s.length()-4,4) == ".txt"){
			auto spl = split(s,'|');          
			if(spl.size() != 2) emsg("'"+s+"' does not use the format 'Column|file.txt'");
			
			auto head = spl[0];
			auto file = spl[1];
			
			vector < vector <string> > matrix;
			
			auto datadir = find_string("datadir","UNSET");
			if(datadir == "UNSET") emsg("Data directory is not set");
			
			auto file_full = datadir+"/"+file;
			ifstream in(file_full);
			if(!in) emsg("Cannot open the file '"+file_full+"'.");

			unsigned int ncol = UNSET;
			do{
				vector <string> vec;
				string line;
				getline(in,line);
				if(in.eof()) break;
				
				stringstream ss(line);
				do{
					string st;
					getline(ss,st,'\t'); st = strip(st);
					vec.push_back(st);
					if(ss.eof()) break;
				}while(true);
				
				if(ncol == UNSET) ncol = vec.size();
				else{
					if(vec.size() != ncol) emsg("Rows in file '"+spl[1]+"' do not all share the same number of columns.");
				}
		
				matrix.push_back(vec);
			}while(true);
			
			if(matrix.size() == 0) emsg("The file "+file+" does not contain any information.");
			auto c=0u; while(c < ncol && head != matrix[0][c]) c++;

			if(c == ncol) emsg("The heading '"+head+"' could not be found in '"+file+"'.");
		
			s = "";
			for(auto r = 1u; r < matrix.size(); r++){ 
				if(r != 1) s += "|";
				s += matrix[r][c];
			}
		}
	}
}
	
	
/// Split up a string at a specified delimiter
vector<string> Inputs::split_string(string& s, char delimiter, unsigned int len) const                                
{                  
	check_from_table(s);
	auto spl = split(s, delimiter);

	if(spl.size() == 0){
		for(auto i = 0u; i < len; i++) spl.push_back("");
	}
	else{
		if(spl.size() == 1 && len > 1){
			for(auto i = 0u; i < len-1; i++) spl.push_back(spl[0]);
		}
		
		if(spl.size() != len){
			stringstream ss; ss << "The string '" << s << "' has the wrong number of elements. It is " << spl.size() << " instead of " << len << ","; 
			emsg(ss.str());
		}
	}
	
	return spl;                                           
}


/// Split up a string at a specified delimiter
vector<double> Inputs::split_number(string& s, char delimiter, unsigned int len) const                                                            
{                   
	check_from_table(s);
	auto spl =  split(s, delimiter);

	vector <double> num;
	if(spl.size() == 0){
		for(auto i = 0u; i < len; i++) num.push_back(UNSET);
	}
	else{
		for(auto i = 0u; i < spl.size(); i++) num.push_back(atof(spl[i].c_str()));
		
		if(num.size() == 1 && len > 1){
			for(auto i = 0u; i < len-1; i++) num.push_back(num[0]);
		}
	
		if(num.size() != len){
			stringstream ss; ss << "The string '" << s << "' has the wrong number of elements. It is " << spl.size() << " instead of " << len << ","; 
			emsg(ss.str());
		}
	}
	
  return num;                                           
}


/// Finds the number of generations used 
void Inputs::find_generation(unsigned int &G) const
{
	G = find_integer("ngeneration",UNSET); 
	if(G == UNSET) emsgroot("'ngeneration' must be set");
}


/// Finds the number of gnerations used or final invT (PAIS algorithm) 
void Inputs::find_generation_or_invT_final(unsigned int &G, double &invT_final) const
{
	G = find_integer("ngeneration",UNSET);
	invT_final = find_double("invT_final",UNSET);  
	if(G == UNSET && invT_final == UNSET) emsgroot("'ngeneration' or 'invT_final' must be set");
	if(G != UNSET && invT_final != UNSET) emsgroot("'ngeneration' and 'invT_final' cannot both be set");
}


/// Finds the number of gnerations used or final invT (PAIS algorithm) 
void Inputs::find_invT_start_invT_final(double &invT_start, double &invT_final) const
{
	invT_start = find_double("invT_start",0);  
	invT_final = find_double("invT_final",UNSET);  
	if(invT_final == UNSET) emsgroot("'invT_final' must be set");
}

	
/// Finds the number of particles
void Inputs::find_nparticle(unsigned int &Ntot, unsigned int &N, unsigned int ncore) const
{
	Ntot = find_integer("nparticle",UNSET);                 // Sets the total number of mcmc particles
	if(Ntot == UNSET) emsgroot("The number of particles must be set");
	if(Ntot%ncore != 0) emsgroot("The number of particles must be a multiple of the number of cores");
	if(Ntot%2 != 0) emsg("The number of particles must be even.");
		
	N = Ntot/ncore;                                            // The number of particles per core
	if(N == 0) emsgroot("'nparticle' must be non-zero");
}


/// Finds the number of chains (MC3)
void Inputs::find_nchain(unsigned int &Ntot, unsigned int &N, unsigned int ncore) const
{
	Ntot = find_integer("nchain",UNSET);                 // Sets the total number of mcmc particles
	if(Ntot == UNSET) emsgroot("'nchain' must be set");
	if(Ntot%ncore != 0) emsgroot("'nchain' must be a multiple of the number of cores");
	
	N = Ntot/ncore;                                            // The number of particles per core
	if(N == 0) emsgroot("'nchain' must be non-zero");
}


/// Finds the number of samples
void Inputs::find_nsample(unsigned int &nsample) const
{
	nsample = find_integer("nsample",UNSET);
	if(nsample == UNSET) emsgroot("'nsample' must be set");
	if(nsample < 10) emsg("'nsample' must be at least 10.");
}


/// Finds the number of burnin steps (PMCMC)
void Inputs::find_nburnin(unsigned int &nburnin, const unsigned int &nsample) const
{
	nburnin = find_integer("nburnin",nsample/4);                 // Sets the total number of mcmc particles
	if(nburnin >= nsample) emsg("'nburnin' must be less than 'nsample'");
}


/// Finds the thinning of samples when using PMCMC or MC3
void Inputs::find_nthin(unsigned int &nthin) const
{
	nthin = find_integer("nthin",1);
	if(nthin >= nsample/10) emsg("'nthin' is too large");
	if(nthin < 1) emsg("'nthin' must be 1 or above");
}


/// Finds the number of burnin steps (PMCMC)
void Inputs::find_nquench(unsigned int &nquench, const unsigned int &nburnin) const
{
	nquench = find_integer("nquench",nburnin/2);                 // Sets the period of quenching (MC3)
	if(nquench >= nburnin) emsg("'nquench' must be less than 'nburnin'");
}


/// Finds the cutoff in EF or the fraction of samples accepted (for ABC)
void Inputs::find_cutoff(double &cutoff, double &cutoff_frac) const
{
	cutoff = find_double("cutoff",UNSET);
	cutoff_frac = find_double("cutoff_frac",UNSET);
	
	if(cutoff == UNSET && cutoff_frac == UNSET) emsgroot("Either 'cutoff' or 'cutoff_frac' need to be set.");
	if(cutoff != UNSET && cutoff_frac != UNSET) emsgroot("'cutoff' and 'cutoff_frac' cannot both be set.");
	
	if(cutoff_frac != UNSET){
		if(cutoff_frac <= 0) emsg("'cutoff_frac' must be positive");
		if(cutoff_frac > 1) emsg("'cutoff_frac' must be less than one.");
	}
	else{
		if(cutoff <= 0) emsg("The error function cutoff must be positive.");
	}
}


/// Finds the proposal size when performing ABCSMC 
void Inputs::find_propsize(double &propsize) const
{
	propsize = find_double("propsize",1);
	
	if(propsize <= 0) emsg("'propsize' must be positive.");
}


/// Finds the fraction of samples accepted per generation (for ABCSMC)
void Inputs::find_cutoff_frac(double &cutoff_frac) const
{
	cutoff_frac = find_double("cutoff_frac",UNSET);
	if(cutoff_frac == UNSET) emsgroot("'cutoff_frac' need to be set.");
		
	if(cutoff_frac <= 0) emsg("'cutoff_frac' must be positive");
	if(cutoff_frac > 1) emsg("'cutoff_frac' must be less than one.");
}


/// Finds the quenching factor for the PAIC algorithm
void Inputs::find_quench_factor(double &quench_factor) const 
{
	quench_factor = find_double("quench_factor",0.5);
}


/// Finds the number of groups particles are divided into (ABC_MBP_GR algorithm)
void Inputs::find_ngroup(unsigned int &ngroup, unsigned int &groupN, unsigned int Ntot) const 
{
	ngroup = find_integer("ngroup",4);
	if(Ntot%ngroup != 0) emsg("'nparticle' goes not split equally into groups");
	if(ngroup <= 1) emsg("'ngroup' must be more than one.");
	groupN = Ntot/ngroup;
}


/// Finds the way in which the output is displayed
void Inputs::find_outputprop(OutputProp &prop) const
{
	prop.type = KDE; prop.nbin = 200; prop.h = 10; 
	
	if(basedata->contains("output_prop")){
		const auto pr = basedata->open("output_prop");
		
		if(pr.contains("probdist")){
			auto pd = pr.stringfield("probdist");
			if(pd == "kde") prop.type = KDE;
			else{
				if(pd == "bin") prop.type = BINNING;
			}
		}
		
		if(pr.contains("nbin")){
			prop.nbin = atoi(pr.stringfield("nbin").c_str());
		}
		
		if(pr.contains("h")){
			prop.h = atof(pr.stringfield("h").c_str());
		}
	}
}


/// Finds any counterfactuals which are being analysed
void Inputs::find_counterfact(const Details &details, vector <CounterFact> &counterfact) const
{
	if(basedata->contains("counterfactual")){
		const auto pr = basedata->open("counterfactual");
		
		for(auto j = 0u; j < pr.size(); j++){
			CounterFact cofa;
			const auto cf = pr[j];
		
			cofa.start = 0;
			if(cf.contains("start")){
				auto start = cf.stringfield("start");
				cofa.start = details.gettime(start)-details.start;
			}
			
			cofa.end = details.end - details.start;
			if(cf.contains("end")){
				auto end = cf.stringfield("end");
				cofa.end = details.gettime(end)-details.start;
			}
		
			if(cf.contains("geo")) cofa.geo = cf.stringfield("geo"); else cofa.geo = "all";
			
			if(cf.contains("democats")) cofa.democats = cf.stringfield("democats"); else cofa.democats = "all";
			
			if(cf.contains("factor")) cofa.factor = atof(cf.stringfield("factor").c_str()); 
			else emsg("'counterfactual' must contain a 'factor' definition");
			
			if(cf.contains("type")){
				auto type = cf.stringfield("type");
				if(type == "trans_rate_fac"){
					cofa.type = CF_TRANS_RATE;
					if(cf.contains("trans")) cofa.trans_str = cf.stringfield("trans");
					else emsg("'counterfactual' must contain a 'trans' definition");
				}
				else{
					if(type == "efoi_fac"){
						cofa.type = CF_EFOI;
						if(cf.contains("trans")) cofa.trans_str = cf.stringfield("trans");
					}
					else emsg("'"+type+"' in 'counterfactual' is not recognised");
				}
			}
				
			counterfact.push_back(cofa);
		}
	}
	else emsg("The input file must contain a 'counterfact' defninition.");
}


/// Finds posterior predictive check start time
void Inputs::find_ppc_start(const Details &details, vector <CounterFact> &counterfact) const
{
	CounterFact cofa;
	cofa.type = CF_PPC;
	
	auto start_str = find_string("ppc_start","UNSET"); 
	if(start_str == "UNSET") cofa.start = 0;
	else cofa.start = details.gettime(start_str)-details.start;
		
	counterfact.push_back(cofa);
}
