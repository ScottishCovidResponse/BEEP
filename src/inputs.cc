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
	set_command_line_params(argc,argv);                               // Loads up the command line parameters
	
	inputfilename = "/dev/null";
	if(cmdlineparams.count("inputfile") == 1) {
		inputfilename = cmdlineparams["inputfile"];
	}

	ifstream ip(inputfilename);
	if(!ip) emsgroot("The file '"+inputfilename+"' could not be opened");
	
	try {
		basedata = new InputData{inputfilename};
	} catch (const std::exception& e) {
		std::ostringstream oss;
		oss << "Error reading TOML file" << endl << e.what();
		emsgroot(oss.str());
	}	
	
	
	auto tomlkeys = basedata->keys();
	check_for_undefined_parameters(definedparams, tomlkeys, "in input TOML file");
	
	for(const auto &cmd : cmdlineparams){
		if(find_in(tomlkeys,cmd.first)) tomlkeys.push_back(cmd.first);
		if(find_in(definedparams,cmd.first) == UNSET) emsgroot("'"+cmd.first+"' is not a valid command");
	}	
	
	for(auto i = 0u; i < tomlkeys.size(); i++){
		UsedTomlKey utk; utk.name = tomlkeys[i]; utk.used = false;
		used.push_back(utk);
	}
	
	addused("inputfile",used);
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
	for(const auto &str : commandlist ){ // Goes the various input options
		auto j = 0u; auto jmax = str.length(); while(j < jmax && str.substr(j,1) != "=") j++;
		if(j == jmax){
			stringstream ss; ss << "Cannot understand the command line option '" << str << "'"; 
			emsgroot(ss.str());
		}
		
		string command = str.substr(0,j);
		string value = str.substr(j+1,jmax-(j+1));
		
		if (cmdlineparams.count(command) == 0) cmdlineparams[command] = value;
		else cmdlineparams[command] += " " + value;
	}
}


/// Finds a string from the TOML file (or uses the default 'def' value)
string Inputs::find_string(const string &key, const string &def, bool allowcomma)
{
	string val;
	auto val_it = cmdlineparams.find(key);
	if(val_it != cmdlineparams.end()){
		val = val_it->second;
		addused(key,used);
	}
	else{
		if(basedata->contains(key)){
			if(allowcomma == false) val = basedata->data.stringfield_unchecked(key);
			else val = basedata->data.stringfield_allowcomma(key,"");
			addused(key,used);
		}
		else val = def;
	}

	return val;
}

/// Returns a positive integer
unsigned int Inputs::find_positive_integer(const string &key, const int def)
{
	auto val = find_integer(key,def);
	if(val <= 0) emsgroot("The value of '"+key+"' should be positive.");
	return val;
}


/// Finds an integer from the TOML file (or uses the default 'def' value)
int Inputs::find_integer(const string &key, const int def)
{
	int val;
	auto val_it = cmdlineparams.find(key);
	if (val_it != cmdlineparams.end()) {
		const string &valstr = val_it->second;
		addused(key,used);
		
		try {
			size_t idx;
			val = stoi(valstr,&idx);
			if (idx != valstr.length()) {
				ostringstream oss;
				oss << "Should be integer, found '"<< valstr;
				throw std::invalid_argument(oss.str());
			}
		} catch (const std::exception &e) {
			std::ostringstream oss;
			oss << "Bad command-line parameter value for key '"<< key << "'" << endl;
		
			std::string what = e.what();
			if(what == "stoi"){
				if(valstr == "") oss << "Should be integer" << endl;
			} 
			else oss << what;
	
			emsgroot(oss.str());
		}
	} 
	else{
		if(basedata->contains(key)){
			val = basedata->data.intfield_unchecked(key);
			addused(key,used);
		} 
		else{
			val = def;
		}
	}
	
	return val;
}


/// Finds a double from the TOML file (or uses the default 'def' value)
double Inputs::find_double(const string &key, const double def)
{
	double val;
	auto val_it = cmdlineparams.find(key);
	if (val_it != cmdlineparams.end()) {
		const std::string &valstr = val_it->second;
		addused(key,used);
		try {
			size_t idx;
			val = stof(valstr,&idx);
			if (idx != valstr.length()) {
				ostringstream oss;
				oss << "Should be number, found '" << valstr;
				throw std::invalid_argument(oss.str());
			}
		} catch (const std::exception &e) {
			ostringstream oss;
			oss << "Bad command-line parameter value for key '"<< key <<"'" << endl;

			std::string what = e.what();
			if(what == "stof"){
				if(valstr == "") oss << "Should be number, found no value" << endl;
			} 
			else oss << what;
			emsgroot(oss.str());
		}
	} 
	else{
		if(basedata->contains(key)){
			val = basedata->data.numberfield_unchecked(key);
			addused(key,used);
		}
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

	set_difference(given.begin(),given.end(),allowed.begin(),allowed.end(),inserter(undefined,undefined.begin()));

	if (undefined.size() != 0) {
		stringstream ss;
		ss << "Unrecognised parameter(s) "+context+":";
		for (const auto &k : undefined){ ss << " " << k;}
		
		emsgroot(ss.str());
	}
}


/// Returns the mode of operation
Mode Inputs::mode()
{
	string val = find_string("mode","UNSET");  
	if(val == "UNSET") emsgroot("The 'mode' property must be set");
	
	Mode mode;
	map<string,Mode>  modemap{{"sim", SIM}, {"multisim", MULTISIM}, {"prediction", PREDICTION}, {"data", DATAONLY}, {"abc", ABC_SIMPLE}, {"abcsmc", ABC_SMC}, {"abcmbp", ABC_MBP}, {"abcda", ABC_DA}, {"abccont", ABC_CONT}, {"mc3", MC3_INF}, {"mcmcmbp", MCMC_MBP}, {"pas", PAS_INF}, {"pmcmc", PMCMC_INF}, {"importance", IMPORTANCE_INF}, {"map", ML_INF}};
	if (modemap.count(val) != 0) mode = modemap[val];
	else emsgroot("Unrecoginsed value '" + val + "' for 'mode'");
	
	return mode;
}


/// Returns whether simulation or inference is being performed
SimInf Inputs::get_siminf()
{
	switch(mode()){
	case SIM: case MULTISIM: case PREDICTION: return SIMULATE;
	case DATAONLY: return DATAVIEW;
	default: return INFERENCE;
	}
	return INFERENCE;
}
	
	
/// Returns the maximum timestep used in data_tables (used in PMCMC to split up observations
unsigned int Inputs::find_maximum_timestep()
{
	auto max = 1u;
	if(basedata->contains("data_tables")) {
		const auto datatables = basedata->open("data_tables",used);
		
		for(auto j = 0u; j < datatables.size(); j++){
			auto td = datatables[j];
			auto timestep = td.stringfield("timestep","");
			if(timestep != ""){
				auto val = get_int(timestep,"In 'data_tables'");
				if(val > max) max = val;
			}
		}
	}
	
	return max;
}

/// Finds information about datatables
vector <DataTable> Inputs::find_datatable(const Details &details)
{
	vector <DataTable> datatable;
	
	if(basedata->contains("data_tables")) {
		const auto datatables = basedata->open("data_tables",used);

		for(auto j = 0u; j < datatables.size(); j++){
			auto td = datatables[j]; td.set_used();
		
			DataTable datatab;
			datatab.optype = DATA;
			
			auto type = td.stringfield("type","In 'data_tables'");
			if(type == "population") datatab.type = POP;
			else{
				if(type == "population_fraction") datatab.type = POPFRAC;
				else{
					if(type == "transition") datatab.type = TRANS;
					else{
						if(type == "marginal") datatab.type = MARGINAL;
						else emsgroot("In 'data_tables' the value of 'type="+type+"' is not recognised (it should be 'transition', 'population', 'population_fraction', or 'marginal')");
					}
				}
			}
			
			datatab.file = td.stringfield("file","In 'data_tables'");
			check_filename(datatab.file);
			
			datatab.threshold = UNSET;
			auto thresh = td.stringfield("threshold","");
			if(thresh != "") datatab.threshold = get_double_positive(thresh,"In 'data_tables' and for 'threshold'");
			
			datatab.factor = 1;
			auto factor = td.stringfield("factor","");
			if(factor != "") datatab.factor = get_double_positive(factor,"In 'data_tables' and for 'factor'");
			
			datatab.observation = td.stringfield_allowcomma("observation","In 'data_tables'");
			
			datatab.shift = 0;
			switch(datatab.type){
				case POP: case POPFRAC: case TRANS:
					{
						auto shift = td.stringfield("shift","");
						if(shift != "") datatab.shift = get_pos_neg_int(shift,"In 'data_tables' and 'shift'");
					}
					break;
					
				case MARGINAL:	
					{				
						auto shift = td.stringfield("shift","");
						if(shift != "") emsgroot("In 'data_tables' for 'marginal' data the property 'shift' cannot be used.");
					}
					break;
			}
	
			datatab.start = UNSET;
			const auto startdata = td.stringfield("start","");
			if(startdata == ""){
				if(details.mode == SIM) datatab.start = -datatab.shift;	
				else{
					if(datatab.type == MARGINAL) datatab.start = 0;
				}
			}
			else{
				datatab.start = details.gettime(startdata,"In 'data_tables' with 'start'") - details.start;
				
				if(datatab.start+datatab.shift < 0){
					emsgroot("In 'data_tables' the value for 'start' (after shifting) must be after 'start' specified for the analysis");
				}
				
				if(datatab.start+datatab.shift >= (int)details.period){
					emsgroot("In 'data_tables' the value for 'start' (after shifting) must be before 'end' specified for the analysis");
				}
			}
			
			datatab.timestep = UNSET;
			auto timestep = td.stringfield("timestep","");
			if(timestep == ""){
				if(details.mode == SIM || datatab.type != TRANS) datatab.timestep = 1;	
			}
			else{
				if(datatab.type == MARGINAL) emsgroot("In 'data_tables' the property 'timestep' cannot be set for 'type=\"marginal\"'");
				datatab.timestep = get_int(timestep,"In 'data_tables'");
			}
			
			datatab.end = UNSET;
			auto enddata = td.stringfield("end","");
			if(enddata == ""){
				if(details.mode == SIM){
					datatab.end = datatab.start + details.period;
					if(datatab.timestep != UNSET){
						datatab.end -= (datatab.end - datatab.start)%datatab.timestep;
						datatab.end -= datatab.timestep;
						if(datatab.end == datatab.start) emsgroot("In 'data_table' the value for 'timestep' is too large");
					}
				}
				else{
					if(datatab.type == MARGINAL) datatab.end = datatab.start + details.period;
				}
			}
			else{
				datatab.end = details.gettime(enddata,"In 'data_tables' for 'end'") - details.start;
				if(datatab.timestep != UNSET && (datatab.end - datatab.start)%datatab.timestep != 0){
					emsgroot("In 'data_tables' the difference between 'start' and 'end' must be a multiple of 'timestep'");
				}
				
				if(datatab.end+datatab.shift < 0){
					emsgroot("In 'data_tables' the value for 'end' must be after 'start' specified for the analysis");
				}
				
				if(datatab.end+datatab.shift >= (int)details.period){
					emsgroot("In 'data_tables' the value for 'end' (after shifting) must be before 'end' specified for the analysis");
				}
			}
		
			datatab.factor_spline = td.stringfield("factor_spline","");
			
			datatab.democats_filt = td.stringfield_allowcomma("democats_filt",""); 
			
			datatab.democats_dep = td.stringfield("democats_dep",""); 
			
			datatab.geo_filt = td.stringfield_allowcomma("geo_filt","");
						
			datatab.geo_dep = td.stringfield("geo_dep","");
						
			if(datatab.democats_dep != "" && datatab.geo_dep != ""){
				emsgroot("In 'data_tables' values for 'democats_dep' and 'geo_dep' cannot both be set");
			}

			datatab.plot_name = td.stringfield("plot_name","");
					
			datatab.line_colour = get_line_colour(td.stringfield("line_colour",""));
	
			datatab.percent = UNSET;
			datatab.epsilon_factor = 0.05;
			auto epsilon_factor = td.stringfield("epsilon_factor","");
			if(epsilon_factor != "") datatab.epsilon_factor = get_double_positive(epsilon_factor,"In 'data_tables' for 'epsilon_factor'");
	
			datatab.load_sd = false;
			auto obsmodel = td.stringfield("obsmodel","");
			if(obsmodel != ""){
				auto spl = split(obsmodel,' ');
				if(spl.size() > 2){
					emsgroot("In 'data_tables' the value '"+obsmodel+"' for 'obsmodel' is invalid.");
				}
				
				auto type = toLower(spl[0]);
				if(type == "normal"){
					datatab.obsmodel = NORMAL_OBSMODEL;
					datatab.percent = UNSET;
				
					if(spl.size() == 2){
						auto value = toLower(spl[1]);
						
						if(value == "load"){
							datatab.load_sd = true;
						}
						else{
							if(value.substr(value.length()-1,1) == "%"){
								datatab.obsmodel = NORMAL_PERCENT_OBSMODEL;
								auto num = value.substr(0,value.length()-1);
								if(num != ""){
									datatab.percent = get_double(num,"In 'obsmodel'");
								}
							}
							else{
								if(value != ""){
									datatab.sd = get_double(value,"In 'obsmodel'");
									cout << datatab.sd << " sd\n";
								}
							}
						}
					}
				}					
				else{
					if(obsmodel == "poisson"){
						datatab.obsmodel = POISSON_OBSMODEL;
						if(spl.size() > 1){
							emsgroot("In 'data_tables' the value '"+obsmodel+"' for 'obsmodel' is invalid.");
						}
					}
					else{
						if(obsmodel == "negbin"){ 
							datatab.obsmodel = NEGBINO_OBSMODEL; 
							
							if(spl.size() != 2){
								emsgroot("In 'data_tables' the value '"+obsmodel+"' for 'obsmodel' must contain a definition for the shape parameter");	
							}
							
							auto shape = spl[1];
						
							datatab.shape = get_double_positive(shape,"In 'data_tables' for 'shape'");
						}
						else{
							emsgroot("In 'data_tables' the value 'obsmodel="+obsmodel+"' is not recognised (it should be '.%', 'normal', 'poisson', 'negbin' or 'scale').");
						}
					}
				}
			}
			else{
				emsgroot("In 'data_tables' a value for 'obsmodel' must be set.");
			}
			
			td.check_used("datatable");
			datatable.push_back(datatab);
		}
	}
	
	if(basedata->contains("state_outputs")) {
		const auto datatables = basedata->open("state_outputs",used);

		for(auto j = 0u; j < datatables.size(); j++){
			auto td = datatables[j]; td.set_used();
		
			DataTable datatab;
			datatab.optype = OUTPUT;
				
			auto type = td.stringfield("type","In 'state_outputs'");
			if(type == "population") datatab.type = POP;
			else{
				if(type == "population_fraction") datatab.type = POPFRAC;
				else{
					if(type == "transition") datatab.type = TRANS;
					else{
						if(type == "marginal") datatab.type = MARGINAL;
						else emsgroot("In 'state_outputs' the value 'type="+type+"' is not recognised (it should be 'population', 'transition' or 'marginal'");
					}
				}
			}
			
			datatab.observation = td.stringfield_allowcomma("observation","In 'state_outputs'");
			
			datatab.label = td.stringfield("label","");
			
			if(td.contains("start")) emsgroot("In 'state_outputs' the value 'start' should not be set.");
			datatab.start = 0;
			datatab.end = details.end-details.start;
			
			if(td.contains("factor_spline")) emsgroot("In 'state_outputs' the value 'factor_spline' should not be set.");
			
			datatab.factor = 1;
			datatab.factor_spline = "";
			
			datatab.democats_filt = td.stringfield_allowcomma("democats_filt",""); 
			
			datatab.democats_dep = ""; 
			
			datatab.geo_filt = td.stringfield_allowcomma("geo_filt","");
					
			datatab.geo_dep = "";
			
			datatab.plot_name = td.stringfield("plot_name","");
					
			datatab.line_colour = get_line_colour(td.stringfield("line_colour",""));
			
			if(td.contains("file")) emsgroot("In 'state_outputs' the value 'file' should not be set.");
			datatab.file = "";
		
			if(td.contains("timestep")) emsgroot("In 'state_outputs' a value for 'timestep' should not be set.");
			if(td.contains("shift")) emsgroot("In 'state_outputs' the value 'shift' should not be set.");
						
			switch(datatab.type){
				case POP: case POPFRAC: case TRANS:
					datatab.shift = 0;
					datatab.timestep = 1;
					break;
					
				case MARGINAL:	
					datatab.shift = UNSET;
					datatab.timestep = UNSET;
					break;
			}
	
			td.check_used("state_outputs");
			
			datatable.push_back(datatab);
		}
	}
	
	return datatable;
}


/// Gets a line colour from a string
LineColour Inputs::get_line_colour(const string col) const
{ 
	if(col == "red") return RED;
	else{
		if(col == "green") return GREEN;
		else{
			if(col == "blue") return BLUE;
			else{
				if(col == "yellow") return YELLOW;
				else{
					if(col == "cyan") return CYAN;
					else{
						if(col == "magenta") return MAGENTA;
						else{
							if(col == "black") return BLACK;
							else{
								if(col != "") emsgroot("In 'line_colour' do not recognise the colour '"+col+"'");
							}
						}
					}
				}
			}
		}
	}
	return UNSET_COLOUR;
}

				
/// Finds information about changes in demographic catergory
vector <DemocatChange> Inputs::find_democat_change()
{
	vector <DemocatChange> democatchange;
	
	if(basedata->contains("democat_change")) {
		const auto democh = basedata->open("democat_change",used);

		for(auto j = 0u; j < democh.size(); j++){
			auto td = democh[j]; td.set_used();
		
			DemocatChange dcc;
			dcc.name = td.stringfield("name","In 'data_tables'");
			dcc.file = td.stringfield("file","In 'data_tables'");
			check_filename(dcc.file);
			
			dcc.democats_filt = td.stringfield_allowcomma("democats_filt",""); 
			dcc.geo_filt = td.stringfield_allowcomma("geo_filt","");
			dcc.shift = 0;
			auto shift = td.stringfield("shift","");
			if(shift != "") dcc.shift = get_int(shift,"In 'democat_change' and 'shift'");
						
			td.check_used("democat_change");
			democatchange.push_back(dcc);
		}
	}
	
	return democatchange;
}


/// Finds and returns information about demographic categories
vector <DemographicCategory> Inputs::find_demographic_category(vector <Strain> &strain_list)
{
	vector <DemographicCategory> democatvec;
	
	DemographicCategory democat;
	democat.name = "age";
	democat.sus_vari= false;
		
	if(basedata->contains("ages")){                                  // Age categories
		auto ages = basedata->open("ages",used); ages.set_used();
			
		auto cats = ages.stringfield("cats","In 'ages'");
		democat.value = split(cats,'|');
			
		if(ages.contains("sus_param") || ages.contains("sus_value") || ages.contains("sus_prior")){
			democat.sus_vari = true;
			democat.ps = get_paramspec(ages,"ages","sus_param","sus_value","sus_prior","","",democat.value.size(),false);
		}
		ages.check_used("ages");
	}
	else{
		democat.value.push_back("population");
	}
	democatvec.push_back(democat);
	
	if(basedata->contains("democats")){                              // Other demographic possibilities
		const auto democats = basedata->open("democats",used);
	
		for(auto k = 0u; k < democats.size(); k++){
			auto democ = democats[k]; democ.set_used();
		
			DemographicCategory democat;
			democat.name = democ.stringfield("name","In 'democats'");
					
			auto cats = democ.stringfield("cats","In 'democats'");
			democat.value = split(cats,'|');
		
			democat.sus_vari = false;
			
			if(democ.contains("sus_param") || democ.contains("sus_value") || democ.contains("sus_prior")){
				democat.sus_vari = true;
				democat.ps = get_paramspec(democ,"democats","sus_param","sus_value","sus_prior","","",democat.value.size(),false);
			}
			
			democ.check_used("democats");
			democatvec.push_back(democat);
		}
	}
	
	DemographicCategory strain;                                         // Different strains
	strain.name = "strain";
	strain.sus_vari = false;
	if(basedata->contains("strains")){ 
		auto strains = basedata->open("strains",used); strains.set_used();
		
		auto cats = strains.stringfield("cats","In 'strains'");
		auto strain_name = split(cats,'|');
		strain.value = strain_name;
			
		if(strains.contains("sus_param") || strains.contains("sus_value") || strains.contains("sus_prior")){
			strain.sus_vari = true;
			strain.ps = get_paramspec(strains,"strains","sus_param","sus_value","sus_prior","","",strain.value.size(),false);
		}
		
		for(auto i = 0u; i < strain_name.size(); i++){
			Strain str; str.name = strain_name[i]; str.Rfactor_ps = ps_one();
			strain_list.push_back(str);
		}
		
		if(strains.contains("Rfactor_param") || strains.contains("Rfactor_value") || strains.contains("Rfactor_prior")){
			auto ps_vec = get_paramspec(strains,"strains","Rfactor_param","Rfactor_value","Rfactor_prior","","",strain_name.size(),false);
			for(auto i = 0u; i < strain_name.size(); i++){
				strain_list[i].Rfactor_ps = ps_vec[i];
				if(strain_list[i].Rfactor_ps.name == "") strain_list[i].Rfactor_ps.name = "Rfactor_"+strain_name[i];
			}
		}
		
		strains.check_used("strains");
	}
	else{
		strain.value.push_back("all");
		
		Strain str; str.name = "all"; str.Rfactor_ps = ps_one();
		strain_list.push_back(str);
	}
	democatvec.push_back(strain);
		
	for(auto d = 0u; d < democatvec.size(); d++){                      // Checks names are not repeated
		for(auto dd = 0u; dd < democatvec.size(); dd++){
			for(auto i = 0u; i < democatvec[d].value.size(); i++){
				for(auto ii = 0u; ii < democatvec[dd].value.size(); ii++){
					if(!(d == dd && i == ii) && democatvec[d].value[i] == democatvec[dd].value[ii]){
						string type;
						if(d == 0) type = "'ages'";
						else{ 
							if(d == democatvec.size()-1) type = "strains"; 
							else type = "'democats'";
						}
						
						if(d != dd){
							if(dd == 0) type += "and 'ages'";
							else{ 
								if(dd == democatvec.size()-1) type += "and 'strains'"; 
								else type += "and 'democats'";
							}
						}
						
						emsgroot("For "+type+" in 'cats' the expression '"+democatvec[d].value[i]+"' cannot be used more than once.");
					}
				}
			}
		}
	}
		
	for(auto &democat : democatvec){                                  // Sets names of susceptibility parameters (if not set)
		if(democat.sus_vari == true){
			if(param_not_set(democat.ps)){
				for(auto i = 0u; i < democat.ps.size(); i++){
					democat.ps[i].name = democat.value[i]+" susceptibility";
				}
			}
		}
	}
	
	return democatvec;
}


/// Finds and returns information about level effects
LevelEffect Inputs::find_level_effect()
{
	LevelEffect le;
	
	le.on = false;
	if(basedata->contains("level_effect")){
		le.on = true;
		
		auto lin = basedata->open("level_effect",used); lin.set_used();
	
		le.file = lin.stringfield("file","In 'level_effect'");
		check_filename(le.file);
	
		auto param = lin.stringfield("param","In 'level_effect'");
		auto spl = split(param,'|'); 
		auto si = spl.size();
	
		le.ps = get_paramspec(lin,"level_effect","param","value","prior","","",si,false);
	
		lin.check_used("level_effect");
	}	
	
	return le;
}


/// Finds and returns information about area covariates
vector <Covariate> Inputs::find_covariates()
{
	vector <Covariate> covarvec;
	
	if(basedata->contains("area_covars")){
		const auto covars = basedata->open("area_covars",used);
		
		Covariate cov; cov.type = AREA_COVAR; cov.file = find_string("areas","UNSET");
		check_filename(cov.file);
		for(auto j = 0u; j < covars.size(); j++){
			auto covar = covars[j]; covar.set_used();
			
			cov.name = covar.stringfield("name","In 'area_covars'");
			cov.timevary = false;
			 
			auto vec = get_paramspec(covar,"area_covars","param","value","prior","","",1,false);
			 
			if(param_not_set(vec)) vec[0].name = "Covariate: "+cov.name;   // Sets parameter name (if uninitialised)
	
			cov.ps = vec[0];
			
			auto funct = covar.stringfield("func","In 'area_covars'");
			if(funct == "log") cov.func = LOG_TRANS; 
			else{ 
				cov.func = LINEAR_TRANS; 
				if(funct != "linear") emsgroot("In 'area_covars' the value 'func="+funct+"' is not recognised."); 
			}
			
			covar.check_used("area_covars");
			covarvec.push_back(cov);
		}
	}
	
	if(basedata->contains("tv_covars")){
		const auto covars = basedata->open("tv_covars",used);
		
		Covariate cov; cov.type = TV_COVAR;
		for(auto j = 0u; j < covars.size(); j++){
			auto covar = covars[j]; covar.set_used();
			
			cov.file = covar.stringfield("file","In 'tv_covars'");
			check_filename(cov.file);
			
			cov.name = covar.stringfield("name","In 'tv_covars'");
			cov.timevary = true;
			
			auto vec = get_paramspec(covar,"tv_covars","param","value","prior","","",1,false);
			 
			if(param_not_set(vec)) vec[0].name = "TV Covariate: "+cov.name;   // Sets parameter name (if uninitialised)
	
			cov.ps = vec[0];
			
			auto funct = covar.stringfield("func","In 'tv_covars'");
			if(funct == "log") cov.func = LOG_TRANS; 
			else{ 
				cov.func = LINEAR_TRANS; 
				if(funct != "linear") emsgroot("In 'tv_covars' the value 'func="+funct+"' is not recognised."); 
			}
			
			covar.check_used("tv_covars");
			covarvec.push_back(cov);
		}
	}
	
	if(basedata->contains("area_tv_covars")){
		const auto covars = basedata->open("area_tv_covars",used);
		
		Covariate cov; cov.type = AREA_TV_COVAR;
		for(auto j = 0u; j < covars.size(); j++){
			auto covar = covars[j]; covar.set_used();
			
			cov.file = covar.stringfield("file","In 'area_tv_covars'");
			check_filename(cov.file);
			
			cov.name = covar.stringfield("name","In 'area_tv_covars'");
			
			auto vec = get_paramspec(covar,"area_tv_covars","param","value","prior","","",1,false);
			 
			if(param_not_set(vec)) vec[0].name = "Area TV Covariate: "+cov.name;   // Sets parameter name (if uninitialised)
	
			cov.ps = vec[0];
			
			auto funct = covar.stringfield("func","In 'area_tv_covars'");
			if(funct == "log") cov.func = LOG_TRANS; 
			else{ 
				cov.func = LINEAR_TRANS; 
				if(funct != "linear") emsgroot("In 'area_tv_covars' the value 'func="+funct+"' is not recognised."); 
			}
			
			covar.check_used("area_tv_covars");
			covarvec.push_back(cov);
		}
	}
	
	return covarvec;
}


/// Finds the break points from a string 
vector <unsigned int> Inputs::get_breakpoints(string st, const Details &details, const string em)
{
	check_from_table(st);

	auto bp = split(st,'|');
	
	vector <unsigned int> vec;
	for(auto j = 0u; j < bp.size(); j++){
	
		auto t = details.gettime(bp[j],em+" for bp=\""+st+"\"");
		
		if(j == 0 && t != details.start){
			emsgroot(em+" 'bp' must start at 'start'.");
		}
		
		if(j == bp.size()-1){
			if(details.mode == PREDICTION){
				if(t > details.pred_end){
					emsgroot(em+" 'bp' must not exceed time for prediction.");
				}
			}
			else{
				if(t-details.start != details.period){
					emsgroot(em+" 'bp' must end at 'end'.");
				}
			}
		}
		
		if(t > details.end) emsgroot(em+" '"+bp[j]+"' in 'bp' is outside of the analysis time period.");
		if(j > 0){
			if(t-details.start < vec[j-1]) emsgroot(em+" times in 'bp' are not in the right order.");
		}
		vec.push_back(t-details.start);
	}		

	return vec;	
}


/// Finds the age catergories from a string
vector <unsigned int> Inputs::find_agecats(const string ages_str, const vector <string> &age_list, const string name) const
{
	vector <unsigned int> vec;
	
	auto nage = age_list.size();
	auto ages = split(ages_str,'|');
	for(auto st : ages){
		auto a = 0u; while(a < nage && age_list[a] != st) a++;
		if(a == nage){
			emsgroot("In '"+name+"' the value '"+st+"' in 'agecat' does not specifiy a valid age classification");
		}
		vec.push_back(a);
	}
	
	return vec;
}


/// Finds properties of 'genQ'
void Inputs::find_genQ(GenerateQ &genQ, const Details &details, const vector <string> &age_list)
{
	genQ.nspline = 0;
	
	if(basedata->contains("age_mixing_matrix")) {
		auto st = find_string("age_mixing_matrix","");
		genQ.N_name = split(st,'|');
	}
	else{
		if(age_list.size() != 1) emsgroot("When more than one age group exists a value for 'age_mixing_matrix' must be set.");
		else genQ.N_name.push_back("UNIT MATRIX");
	}

	if(basedata->contains("age_mixing_modify")) {
		const auto matrix_mod = basedata->open("age_mixing_modify",used);

		string em = "In 'age_mixing_modify'";
		
		for(auto j = 0u; j < matrix_mod.size(); j++){
			auto mm = matrix_mod[j]; mm.set_used();
		
			string type = mm.stringfield("type",em);
			if(type == "row_column"){
				MatrixModification matmod;
				matmod.type = ROW_COLUMN;
				
				auto ages_str = mm.stringfield("agecat",em);
				matmod.ages = find_agecats(ages_str,age_list,"age_mixing_modify");
				
				matmod.desc = "Factor multipling the "+ages_str+" rows and columns of the contact matrix";
				matmod.name = "MatMod "+ages_str+" factor";
				
				matmod.bp = get_breakpoints(mm.stringfield("bp",em),details,em);
			
				matmod.ps = get_paramspec(mm,"age_mixing_modify","param","value","prior","smooth","factor",matmod.bp.size(),false);
				
				if(param_not_set(matmod.ps)){                                  // Sets parameter name (if uninitialised)
					for(auto i = 0u; i < matmod.bp.size(); i++) matmod.ps[i].name = matmod.name+" t="+to_string(matmod.bp[i]);
				}
	
				genQ.nspline++;
				
				genQ.matmod.push_back(matmod);
			}
			else{
				emsgroot("In 'age_mixing_modify' the value of 'type="+type+"' is not recognised");
			}
			
			mm.check_used("age_mixing_modify");
		}
	}

	if(basedata->contains("age_mixing_perturb")) {
		const auto matrix_mod = basedata->open("age_mixing_perturb",used);
		
		string em = "In 'age_mixing_perturb";
		
		for(auto j = 0u; j < matrix_mod.size(); j++){
			auto mm = matrix_mod[j]; mm.set_used();

			MatrixModification matmod;
			matmod.type = PERTURB;
			
			auto ages_str = mm.stringfield("agecat",em);
			matmod.ages = find_agecats(ages_str,age_list,"age_mixing_perturb");
		
			matmod.smoothtype = find_smoothtype(mm.stringfield("smooth_type",""),"age_mixing_perturb");
			matmod.desc = "Factor multipling the "+ages_str+" rows and columns of the contact matrix";
			matmod.name = "MatMod "+ages_str+" factor";
			matmod.bp = get_breakpoints(mm.stringfield("bp",em),details,em);

			matmod.ps = get_paramspec(mm,"age_mixing_perturb","param","value","prior","smooth","factor",matmod.bp.size(),false);
	
			if(param_not_set(matmod.ps)){                                  // Sets parameter name (if uninitialised)
				for(auto i = 0u; i < matmod.bp.size(); i++) matmod.ps[i].name = matmod.name+" t="+to_string(matmod.bp[i]);
			}
			
			genQ.nspline++;
			genQ.matmod.push_back(matmod);
			
			mm.check_used("age_mixing_perturb");
		}
	}
	
	if(basedata->contains("geo_mixing_matrix")) {
		genQ.M_name = find_string("geo_mixing_matrix","");
	}
	else{
		genQ.M_name = "";
	}
}


/// Finds times (used for plots)
void Inputs::find_timeplot(vector <TimePlot> &timeplot)
{
	if(basedata->contains("time_labels")){
		const auto times = basedata->open("time_labels",used);
		for(auto j = 0u; j < times.size(); j++){
			auto time = times[j]; time.set_used();
			TimePlot tp;
			tp.name = time.stringfield("name","In 'time_labels'");
			tp.time_str = time.stringfield("time","In 'time_labels'");
			timeplot.push_back(tp);
			time.check_used("time_labels");
		}
	}
}


/// Finds information about regional effects
void Inputs::find_region_effect(vector <ParamSpec> &ps_vec, ParamSpec &sigma, const vector <Area> &area)
{
	if(basedata->contains("region_effect")){
		auto re = basedata->open("region_effect",used);
		
		auto vec = get_paramspec(re,"region_effect","sigma_param","sigma_value","sigma_prior","","",1,false);
		sigma = vec[0];
		 
		string par = "param";
		if(!re.contains(par)){
			par = ""; for(auto c = 0u; c < area.size(); c++){ if(c != 0) par += "|"; par += "RE_"+area[c].code;}
		}
		
		string pri = "prior";
		if(!re.contains(pri) && sigma.name != "") pri = "Normal(0,"+sigma.name+")";
		
		ps_vec = get_paramspec(re,"region_effect",par,"value",pri,"","",area.size(),false);
	}
}	


/// Determines if the parameters have not been set
bool Inputs::param_not_set(const vector <ParamSpec> &ps_vec) const
{
	auto nnot_set = 0u;
	for(const auto &ps : ps_vec){ if(ps.name == "") nnot_set++;}
	if(nnot_set == ps_vec.size()) return true;
	if(nnot_set != 0) emsgroot("Not all parameter values are set");
	return false;
}


/// Finds information about area effects
AreaEffect Inputs::find_area_effect(const vector <Area> &area)
{
	AreaEffect ae;
	
	ae.on = false;
	if(basedata->contains("area_effect")){
		auto areaeff = basedata->open("area_effect",used); areaeff.set_used();
		
		ae.on = true;
		
		ae.ps = get_paramspec(areaeff,"area_effect","param","value","prior","","",area.size(),false);
		
		if(param_not_set(ae.ps)){                                        // Sets parameter name (if uninitialised)
			for(auto c = 0u; c < area.size(); c++) ae.ps[c].name = "Area effect "+area[c].code;
		}
		
		ae.frac.resize(area.size());
		auto tot_pop = 0.0;
		for(auto c = 0u; c < area.size(); c++){
			ae.frac[c] = 0; for(const auto &vec : area[c].pop_init){ for(auto val : vec) ae.frac[c] += val;}
			tot_pop += ae.frac[c];
		}
		for(auto c = 0u; c < area.size(); c++) ae.frac[c] /= tot_pop;
		
		 areaeff.check_used("area effect");
	}	
	
	return ae;
}	


/// Gets the dependency from a string
void Inputs::get_dep(const string root, string &st, string &dep) const 
{
	if(st != ""){
		auto splbr = split(st,'[');
		auto spl = split(splbr[0],':');
		if(spl.size() > 2) emsgroot("In '"+root+"' there was a problem with '"+st+"'");
		if(spl.size() == 2){
			st = spl[1]; 
			if(dep == "") dep = spl[0];
			else{ if(dep != spl[0]) emsgroot("In '"+root+"' cannot have conflicting dependencies");}
		}
	}
}


/// Gets the dependency from the strings
string Inputs::get_dep(const string root, string &name, string &value, string &prior) const
{
	string dep="";
	get_dep(root,name,dep);
	get_dep(root,value,dep);
	get_dep(root,prior,dep);
	
	return dep;
}


/// Gets the dependency from InputNode
string Inputs::get_dep(InputNode &it, const string root, const string name, const string value, const string prior) const
{	
	auto name_string = it.stringfield(name,"");	
	auto value_string = it.stringfield(value,"");
	auto prior_string = it.stringfield_allowcomma(prior,"");
	
	return get_dep(root,name_string,value_string,prior_string);
}
		
		
/// Gets a parameter specification from a particular key and an index j
vector <ParamSpec> Inputs::get_paramspec(InputNode &it, const string root, const string name, const string value, const string prior, const string smooth, const string factor, const unsigned int si, const bool dep_expect)
{
	auto name_string = it.stringfield(name,"");	
	auto value_string = it.stringfield(value,"");
	auto prior_string = it.stringfield_allowcomma(prior,"");
	
	auto dep = get_dep(root,name_string,value_string,prior_string);
	if(dep != "" && dep_expect == false) emsgroot("In '"+root+"' the dependency '"+dep+"' cannot be used");
	
	auto smooth_string = it.stringfield(smooth,"");
	auto factor_string = it.stringfield(factor,""); if(factor_string == "") factor_string = "1";

	if(root == "region_effect"){
		if(name_string == "" && name != "sigma_param") name_string = name;
		if(prior_string == "" && prior != "sigma_prior") prior_string = prior;
	}
	
	switch(get_siminf()){
		case SIMULATE:
			if(value_string == ""){
				if(root == "trans"){
					emsgroot("In '"+root+"' a value for '"+value+"' must be set for the transition "+it.stringfield("from","")+" -> "+it.stringfield("to",""));	
				}
				else{
					emsgroot("In '"+root+"' a value for '"+value+"' must be set.");
				}
			}
			break;
		
		case INFERENCE: case DATAVIEW:
			if(prior_string == "" && value_string == "") emsgroot("In '"+root+"' either '"+prior+"' or '"+value+"' must be set");
			break;
	}
	
	
	auto na = split_string(name_string,'|',si);
	auto val = split_string(value_string,'|',si);
	
	// Creates fixed  prior if value is set
	if((get_siminf() == INFERENCE || get_siminf() == DATAVIEW) && prior_string == "" && value_string != ""){ 
		for(auto i = 0u; i < si; i++){
			if(i != 0) prior_string += "|";
			prior_string += "Fixed("+val[i]+")";
		}
	}
	
	auto pri = split_string(prior_string,'|',si);
	auto smo = split_number(smooth_string,'|',si);
	auto fac = split_number(factor_string,'|',si);

	vector <ParamSpec> vec; 
	for(auto i = 0u; i < si; i++){
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
		cout << " Prior:" <<  ps[i].prior << " Smooth:" <<  ps[i].smooth << " ps" << endl;
	}
}


/// Given a string this works out a list of states consistent with that
vector <vector <unsigned int > > Inputs::find_demo_ref(const string dep_str, const vector <DemographicCategory> &democat, vector <unsigned int> &dep) const
{
	vector <vector <unsigned int > > demo_ref;
	
	auto ndemocat = democat.size();
	
	vector <unsigned int> index(ndemocat);
	
	if(dep_str == ""){
		for(auto i = 0u; i < ndemocat; i++) index[i] = UNSET;
		demo_ref.push_back(index);
	}
	else{
		auto spl = split(dep_str,',');

		auto ndep = spl.size();
		
		dep.resize(ndep);
		for(auto i = 0u; i < ndep; i++){
			auto k = 0u; while(k < ndemocat && democat[k].name != spl[i]) k++;
			if(k == ndemocat){
				emsgroot("The dependency '"+spl[i]+"' cannot be found");
			}
			
			if(democat[k].value.size() == 1){
				emsgroot("The dependency '"+dep_str+":' cannot be used because '"+democat[k].name+"' only contains one value.");
			}
			
			dep[i] = k; 
			for(auto jj = 0u; jj < i; jj++){
				if(dep[jj] == k) emsgroot("Ihe value '"+dep_str+"' contains multiple dependencies");
			}
		}
		
		for(auto i = 0u; i < ndemocat; i++) index[i] = UNSET;
		
		for(auto i = 0u; i < ndep; i++) index[dep[i]] = 0;
		
		unsigned int i;
		do{
			demo_ref.push_back(index);
			 
			i = 0; 
			bool flag;
			do{
				flag = false;
				index[dep[i]]++; if(index[dep[i]] == democat[dep[i]].value.size()){ index[dep[i]] = 0; flag = true; i++;}
			}while(flag == true && i < ndep);
		}while(i < ndep);
	}
	
	if(false){
		for(auto k = 0u; k < demo_ref.size(); k++){
			cout << k << " ";
			for(auto i = 0u; i < ndemocat; i++){ 
				if(demo_ref[k][i] == UNSET) cout << "unset ";
				else cout << democat[i].value[demo_ref[k][i]] << " ";
			}
			cout << "demo" << endl;
		} 
	}
	
	return demo_ref;
}
		
		
/// Finds parameter specifications for transitions / branching probabilities
vector <ParamSpec> Inputs::find_ct_param_spec(InputNode &it, const string root, const string name, const string value, const string prior, const string dep_str, const vector < vector<unsigned int> > &democatpos, const vector <DemographicCategory> &democat, const string par_root)
{
	vector <ParamSpec> vec;
	
	auto ndemocatpos = democatpos.size();
	auto ndemocat = democat.size();

	if(it.contains(name) || it.contains(value) || it.contains(prior)){ 
		vector <unsigned int> dep_order;
		auto demo_ref = find_demo_ref(dep_str,democat,dep_order);
		
		auto param_spec = get_paramspec(it,root,name,value,prior,"","",demo_ref.size(),true);
		
		if(param_not_set(param_spec)){                                  // Sets the parameter name (if not intialised)
			for(auto ref = 0u; ref < demo_ref.size(); ref++){
				string par_name;
				for(auto i = 0u; i < dep_order.size(); i++){
					auto c = dep_order[i];
					if(i > 0) par_name += ",";
					par_name += democat[c].value[demo_ref[ref][c]];
				}
				if(dep_order.size() != 0) par_name += " ";
				par_name += par_root;
				param_spec[ref].name = par_name;
			}
		}
		
		for(auto d = 0u; d < ndemocatpos; d++){
			vector <unsigned int> list;
			for(auto ref = 0u; ref < demo_ref.size(); ref++){
				auto i = 0u; while(i < ndemocat && (democatpos[d][i] == demo_ref[ref][i] || demo_ref[ref][i] == UNSET)) i++;
				if(i == ndemocat) list.push_back(ref);
			}
			if(list.size() != 1) emsgEC("Input",1);
		
			vec.push_back(param_spec[list[0]]);
		}
	}
	
	return vec;
}


/// Find information about probreach
void Inputs::find_probreach(vector <string> &name, vector <string> &comp)
{
	if(basedata->contains("prob_reach")){
		const auto pr = basedata->open("prob_reach",used);
		
		for(auto num = 0u; num < pr.size(); num++){   // Allows for different spline for different areas / infectious transitions
			auto probr = pr[num];
			name.push_back(probr.stringfield("name","In 'prob_reach'"));
			comp.push_back(probr.stringfield("comp","In 'prob_reach'"));
		}
	}
}


// Converts from a string to a smooth type 
SmoothType Inputs::find_smoothtype(const string st, const string name) const
{
	if(st != ""){
		if(st == "log_smooth") return LOGSMOOTH;
		else{
			if(st == "smooth") return SMOOTH;
			else{
				if(st == "no_smooth") return NOSMOOTH;
				else{
					emsgroot("In '"+name+"' the value of 'smooth_type="+st+"' is not recognised.");
				}
			}
		}
	}
	
	return NOSMOOTH;
}

 
/// Finds spline information
void Inputs::find_spline(const string name, vector <string> &name_vec, vector < vector <ParamSpec> > &ps_vec, vector <ParamSpec> &factor_param_vec, vector < vector <unsigned int> > &bp_vec, vector <SmoothType> &smooth_type_vec, vector <string> &strain, vector <string> &area, vector < vector <double> > &efoi_agedist_vec, const vector <double> &agedist, const Details &details)
{
	name_vec.clear(); ps_vec.clear(); factor_param_vec.clear(); 
	bp_vec.clear(); smooth_type_vec.clear(); strain.clear(); area.clear(); efoi_agedist_vec.clear();
	
	if(basedata->contains(name)) {
		const auto bespin = basedata->open(name,used);
		
		for(auto num = 0u; num < bespin.size(); num++){   // Allows for different spline for different areas / infectious transitions
			auto besp = bespin[num];
			
			name_vec.push_back(besp.stringfield("name",""));
			
			auto strain_name = besp.stringfield("strain","");
			strain.push_back(strain_name);
		
			auto geo_filt = besp.stringfield_allowcomma("geo_filt","");
			area.push_back(geo_filt);
	
			auto bp_str = besp.stringfield("bp",""); if(bp_str == "") bp_str = "start|end";
			
			
			auto bp = get_breakpoints(bp_str,details,"In '"+name+"'");
			
			auto ps_list = get_paramspec(besp,name,"param","value","prior","smooth","factor",bp.size(),false);
			
			if(details.mode == PREDICTION) extend_spline_for_prediction(bp,ps_list,details);
		
			bp_vec.push_back(bp);
		
			if(param_not_set(ps_list)){                                  // Sets the parameter name (if not intialised)
				for(auto i = 0u; i < bp.size(); i++){
					auto par_name = replace(name,"_","-");
				
					if(name == "obs_spline") par_name = name_vec[name_vec.size()-1];
					
					if(geo_filt != "") par_name += " "+geo_filt;
					if(strain_name != "") par_name += " "+strain_name;
					
					if(bp.size() > 2) par_name += " t="+to_string(bp[i]);
					ps_list[i].name = par_name;
				}
			}
	
			ps_vec.push_back(ps_list);
		
			if(besp.contains("factor_param") || besp.contains("factor_value") || besp.contains("factor_prior")){
				auto fpps = get_paramspec(besp,name,"factor_param","factor_value","factor_prior","","",1,false);
				
				if(param_not_set(fpps)) fpps[0].name = name+" factor";              // Sets the parameter name (of not intialised)

				factor_param_vec.push_back(fpps[0]);
			}
			else{
				factor_param_vec.push_back(ps_one());
			}
		
			auto st = find_smoothtype(besp.stringfield("smooth_type",""),name);
			smooth_type_vec.push_back(st);
			
			unsigned int nage = agedist.size();
			
			auto aged = besp.stringfield("age_dist","");
			if(aged != ""){
				if(nage <= 1) emsgroot("In '"+name+"' the value of 'age_dist' is not required.");
			
				auto ad = split_number(aged,'|',nage);
				
				auto sum = 0.0; for(auto a = 0u; a < nage; a++) sum += ad[a];
				for(auto a = 0u; a < nage; a++) ad[a] /= (sum*agedist[a]);

				if(false){
					for(auto a = 0u; a < nage; a++) cout << ad[a] << "  age dist" << endl; 
					emsgroot("Age distribution");
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
	else{
		if(name == "efoi_spline"){          // If no external foi specified then set to zero
			name_vec.push_back("");
			
			vector <ParamSpec> ps_list; ps_list.push_back(ps_zero()); ps_list.push_back(ps_zero());
			ps_vec.push_back(ps_list);
			factor_param_vec.push_back(ps_one());
			
			auto bp_str = "start|end";
			auto bp = get_breakpoints(bp_str,details,"In '"+name+"'");
			bp_vec.push_back(bp);
			
			auto st = find_smoothtype("",name);
			smooth_type_vec.push_back(st);
			
			strain.push_back("");
			area.push_back("");
			
			unsigned int nage = agedist.size();
			vector <double> ad(nage);
			for(auto a = 0u; a < nage; a++) ad[a] = 1;
			efoi_agedist_vec.push_back(ad);
		}
	}
}


/// If in prediction mode then extends splines to match up with the predition time
void Inputs::extend_spline_for_prediction(vector <unsigned int> &bp, vector <ParamSpec> &ps_list, const Details &details)
{
	if(bp[bp.size()-1] < details.pred_end){
		auto flag = false;
		if(bp.size() >= 2){
			const auto &ps1 = ps_list[ps_list.size()-2];
			const auto &ps2 = ps_list[ps_list.size()-1];
			if(ps1.name == ps2.name && ps1.value == ps2.value && ps1.prior == ps2.prior
				&& ps1.smooth == ps2.smooth && ps1.factor == ps2.factor){
				flag = true;
			}
		}
		
		if(flag == true) bp[bp.size()-1] = details.pred_end;
		else{
			bp.push_back(details.pred_end);
			ps_list.push_back(ps_list[ps_list.size()-1]);
		}
	}
}


/// Checks if a string needs to be read from a table
void Inputs::check_from_table(string &s)
{
	strip(s);
	if(s.length() > 6 && s.substr(0,1) == "[" && s.substr(s.length()-1,1) == "]"){
		s = s.substr(1,s.length()-2);
		auto spl = split(s,':');
		if(spl.size() != 2) emsgroot("The expression '"+s+"' does not use the format '[column:file.csv]'");
		auto head = spl[0];
		auto file = spl[1];
		check_filename(file);
		
		auto end = file.substr(file.length()-4,4);
		if(end != ".txt" && end != ".csv") emsgroot("The expression '"+s+"' does not use the format '[column:file.csv]'");
			
		vector < vector <string> > matrix;
		
		auto datadir = find_string("datadir","UNSET");
		if(datadir == "UNSET") emsgroot("'datadir' is not set");
		
		auto file_full = datadir+"/"+file;
		ifstream in(file_full);
		if(!in) emsgroot("Cannot open the file '"+file_full+"'.");

		char del; if(end == ".txt") del = '\t'; else del = ',';
		
		unsigned int ncol = UNSET;
		do{
			string line;
			getline(in,line);
			if(in.eof()) break;
			auto vec = split(line,del);
			
			if(ncol == UNSET) ncol = vec.size();
			else{
				if(vec.size() != ncol) emsgroot("The rows in file '"+spl[1]+"' do not all share the same number of columns.");
			}
	
			matrix.push_back(vec);
		}while(true);
		
		if(matrix.size() == 0) emsgroot("The file '"+file+"' does not contain any information.");
		auto c=0u; while(c < ncol && head != matrix[0][c]) c++;

		if(c == ncol) emsgroot("The heading '"+head+"' in file '"+file+"' could not be found.");
	
		s = "";
		for(auto r = 1u; r < matrix.size(); r++){ 
			if(r != 1) s += "|";
			s += matrix[r][c];
		}
	}
}
	
	
/// Split up a string at a specified delimiter
vector<string> Inputs::split_string(string &s, const char delimiter, const unsigned int len)                        
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
			stringstream ss;
			ss << "The string '" << s << "' has the wrong number of elements. It is ";
			ss << spl.size() << " instead of " << len << ".";
			if(len == 2) ss << " Perhaps 'bp' needs to be specified.";
			emsgroot(ss.str());
		}
	}
	
	return spl;                                           
}


/// Split up a string into numbers at a specified delimiter
vector<double> Inputs::split_number(string &s, const char delimiter, const unsigned int len)                                                            
{                  
	check_from_table(s);
	auto spl = split(s,delimiter);

	vector <double> num;
	if(spl.size() == 0){
		for(auto i = 0u; i < len; i++) num.push_back(UNSET);
	}
	else{
		for(auto i = 0u; i < spl.size(); i++){
			num.push_back(get_double_with_tobeset(spl[i],"In '"+s+"'"));
		}
		
		if(num.size() == 1 && len > 1){
			for(auto i = 0u; i < len-1; i++) num.push_back(num[0]);
		}
	
		if(num.size() != len){
			stringstream ss; ss << "The string '" << s << "' has the wrong number of elements. It is " << spl.size() << " instead of " << len << ","; 
			emsgroot(ss.str());
		}
	}
	
  return num;                                           
}


/// Finds the number of generations used 
void Inputs::find_generation(unsigned int &G)
{
	G = find_positive_integer("ngeneration",UNSET)+1; 
	if(G == UNSET) emsgroot("A value for 'ngeneration' must be set");
}


/// Finds the number of gnerations used or final invT (PAS algorithm) 
void Inputs::find_generation_or_invT_final_or_cpu_time(unsigned int &G, double &invT_final, double &cpu_time)
{
	G = find_positive_integer("ngeneration",UNSET);
	invT_final = find_double("invT_final",UNSET);  
	cpu_time = find_double("cpu_time",UNSET); 

	auto num = 0u;
	if(G != UNSET) num++;
	if(invT_final != UNSET) num++;
	if(cpu_time != UNSET) num++;
	
	if(num == 0) emsgroot("A value for 'ngeneration', 'invT_final', or 'cpu_time' must be set");
	if(num > 1) emsgroot("A value for either 'ngeneration', 'invT_final', or 'cpu_time' must be set");
	if(G != UNSET) G++;
}


/// Finds the cpu time limit
void Inputs::find_cpu_time(double &cpu_time)
{
	cpu_time = find_double("cpu_time",UNSET); 
}


/// Finds the number of gnerations used or final invT (ABC-MBP algorithm) 
void Inputs::find_generation_or_cutoff_final_or_cpu_time(unsigned int &G, double &cutoff_final, double &cpu_time)
{
	G = find_positive_integer("ngeneration",UNSET);
	cutoff_final = find_double("cutoff_final",UNSET); 
	cpu_time = find_double("cpu_time",UNSET); 

	auto num = 0u;
	if(G != UNSET) num++;
	if(cutoff_final != UNSET) num++;
	if(cpu_time != UNSET) num++;
	
	if(num == 0) emsgroot("A value for 'ngeneration', 'cutoff_final', or 'cpu_time' must be set");
	if(num > 1) emsgroot("A value for either 'ngeneration', 'cutoff_final', or 'cpu_time' must be set");
	if(G != UNSET) G++;
}


/// Finds the number of gnerations used or final invT (PAS algorithm) 
void Inputs::find_invT_start_invT_final(double &invT_start, double &invT_final)
{
	invT_start = find_double("invT_start",0);  
	invT_final = find_double("invT_final",UNSET);  
	if(invT_final == UNSET) emsgroot("A value for 'invT_final' must be set");
}


/// Finds the number of gnerations used or final invT (PAS algorithm) 
void Inputs::find_Tpower(double &Tpower)
{
	Tpower = find_double("invT_power",4);  
	if(Tpower < 0) emsgroot("The value for 'invT_power' must be positive");
}


/// Finds the invT used (MCMCMBP algorithm) 
void Inputs::find_invT(double &invT)
{
	invT = find_double("invT",UNSET);  
	if(invT == UNSET) emsgroot("A value for 'invT' must be set");
}


/// Finds the number of runs
void Inputs::find_nrun(unsigned int &nrun)
{
	nrun = find_positive_integer("nrun",1);                    // Sets the number of runs
	if(nrun < 1) emsgroot("'nrun' must be positive");
	if(nrun > 100) emsgroot("'nrun' is too large");
}


/// Checks if a filename is valid or not 
void Inputs::check_filename(const string &str) const
{
	string invalid = "#£%&{}\\/$!'\":+`|=<>*?@";
	
	for(auto i = 0u; i < str.length(); i++){
		for(auto j = 0u; j < invalid.length(); j++){
			if(str.substr(i,1) == invalid.substr(j,1)){
				emsgroot("The filename '"+str+"' is not valid because of the '"+invalid.substr(j,1)+"' character.");
			}
		}
	}
}


/// Finds the algorithm type for maximum likelihood approaches
void Inputs::find_algorithm(MLAlg &algorithm, unsigned int &npart, unsigned int &G, double &cpu_time, unsigned int &P, unsigned int &nsample_final, const unsigned int core, const unsigned int ncore, const unsigned int nvar)
{
	string val = toLower(find_string("algorithm","cmaes"));

	if(val == "gd") algorithm = ML_GD;
	else{
		if(val == "cmaes") algorithm = ML_CMAES;
		else emsgroot("The value for 'algorithm' is not recognised");
	}
	
	npart = UNSET;
	G = UNSET;
	cpu_time = UNSET;
	P = UNSET;
	if(algorithm == ML_CMAES){
		G = find_positive_integer("ngeneration",UNSET);
		cpu_time = find_double("cpu_time",UNSET); 
		if(G != UNSET && cpu_time != UNSET) emsgroot("For 'cmaes' the quantities 'ngeneration' and 'cpu_time' cannot both be set");
		
		
		if(G == UNSET && cpu_time == UNSET){
			G = ITERATE_GENERATION; 
			if(core == 0) cout << "By default generations are iterated until convergence." << endl;
		}
		
		npart = find_positive_integer("nparticle",UNSET);
		if(npart == UNSET){
			npart = 4+int(5*log(nvar)); 
			npart = int((npart+ncore-1)/ncore)*ncore;
			if(core == 0) cout << "By default 'nparticle' set to " << npart << "." << endl;
		}
	
		P = find_positive_integer("posterior_particle",UNSET);
		if(P == UNSET){
			P = 100;
			//WP = 20;
			if(core == 0) cout << "By default 'posterior_particle' is set to " << P << "." << endl;
			//emsgroot("'posterior_particle' must be set to a positive integer. This determines the number of particles used when generating the final posterior samples (this would typically be over 100).");
		}
		
		nsample_final = find_positive_integer("nsample_final",UNSET);
		if(nsample_final == UNSET){
			nsample_final = (int((400+ncore-1)/(ncore)))*ncore;
			if(core == 0) cout << "By default 'nsample_final' is set to " << nsample_final << "." << endl;
		}
		else{
			if(nsample_final%ncore != 0) emsgroot("'nsample_final' must be a multiple of the number of cores");
		}
	}		
}


/// Finds the standard deviation (used for testing PMCMC)
void Inputs::find_sd(double &sd)
{
	sd = find_double("sd",0.5);
}


/// Finds how uncertainty in 
void Inputs::find_stateuncer(StateUncertainty &stateuncer)
{
	string val = toLower(find_string("state_uncertainty","CI"));  
	
	if(val == "ci") stateuncer = CI;
	else{
		if(val == "curves") stateuncer = CURVES;
		else emsgroot("The value for 'state_uncertainty' is not recognised");
	}
}


/// Finds the maximum correlations (used in ABCMBP / PAS / ABCDA)
void Inputs::find_cor_max(double &cor_max, double &cor_max_last)
{
	cor_max = find_double("cor_max",0.8);
	cor_max_last = find_double("cor_max_last",cor_max);
}


/// Finds the the maximum number of proposals before an error is generated (used in ABCMBP / PAS)
void Inputs::find_prop_max(unsigned &prop_max)
{
	prop_max = find_positive_integer("prop_max",10000);
}


/// Finds the number of particles
void Inputs::find_nparticle(unsigned int &npart, unsigned int &Ntot, unsigned int &N, const unsigned int nrun, const unsigned int ncore)
{
	npart = find_positive_integer("nparticle",UNSET);          // Sets the total number of mcmc particles
	if(npart == UNSET) emsgroot("A value for 'nparticle' must be set");
	if(npart%2 != 0) emsgroot("'nparticle' must be an even number.");
		
	Ntot = npart*nrun;                                         // The total number of particles (across all runs)
	if(Ntot%ncore != 0) emsgroot("The total number of particles ('nparticle' times 'nrun') must be a multiple of the number of cores");
	N = Ntot/ncore;                                            // The number of particles per core
	if(N == 0) emsgroot("'nparticle' must be non-zero");
}


/// Finds the number of particles
void Inputs::find_nparticle_pmcmc(unsigned int &npart, unsigned int &N, const unsigned int ncore)
{
	npart = find_positive_integer("nparticle",UNSET);          // Sets the total number of mcmc particles
	if(npart == UNSET) emsgroot("'nparticle' must be set");

	if(npart%ncore != 0) emsgroot("'nparticle' must be a multiple of the number of cores");
	N = npart/ncore;                                           // The number of particles per core
	if(N == 0) emsgroot("'nparticle' must be non-zero");
}


/// Finds the number of chains (MC3)
void Inputs::find_nchain(unsigned int &nchain, unsigned int &Ntot, unsigned int &N, const unsigned int nrun, const unsigned int ncore)
{
	nchain = find_positive_integer("nchain",UNSET);            // Sets the number of chains per run
	if(nchain == UNSET) emsgroot("A value for 'nchain' must be set");
	if(nchain < 2) emsgroot("'nchain' must be at least 2");
	Ntot = nchain*nrun;                                        // Sets the total number of chains (across all runs)
	if(Ntot%ncore != 0) emsgroot("'nchain' must be a multiple of the number of cores");
	N = Ntot/ncore;                                            // The number of particles per core
	if(N == 0) emsgroot("'nchain' must be non-zero");
}


/// Finds the number of samples
void Inputs::find_nsample(unsigned int &nsample, unsigned int min)
{
	nsample = find_positive_integer("nsample",UNSET);
	if(nsample == UNSET) emsgroot("A value for 'nsample' must be set");
	if(nsample < min) emsgroot("'nsample' must be at least 10.");
}


/// Finds the number of samples
void Inputs::find_nsample_final(unsigned int &nsample_final,unsigned int &nsample)
{
	nsample_final = find_positive_integer("nsample_final",nsample);
}


/// Finds the number of samples
void Inputs::find_nsample_or_ESSmin_or_cpu_time(unsigned int &nsample, double &ESSmin, double &cpu_time)
{
	nsample = find_positive_integer("nsample",UNSET);
	ESSmin = find_double("ESS_min",UNSET);
	cpu_time = find_double("cpu_time",UNSET);
			
	auto num = 0u;
	if(nsample != UNSET) num++;
	if(ESSmin != UNSET) num++;
	if(cpu_time != UNSET) num++;
	
	if(num == 0) emsgroot("A value for 'nsample', 'ESS_min', or 'cpu_time' must be set");
	if(num > 1) emsgroot("Either 'nsample', 'ESS_min', or 'cpu_time' must be set");
}


/// Finds the number of burnin steps (PMCMC)
void Inputs::find_nburnin(unsigned int &nburnin, const unsigned int &nsample)
{
	nburnin = find_positive_integer("nburnin",UNSET);           // Sets the total number of mcmc particles
	//if(nburnin == UNSET && nsample != UNSET) nburnin = nsample/4;
	if(nburnin == UNSET && nsample != UNSET) nburnin = nsample/2; // TO DO
	if(nburnin >= nsample && nsample != UNSET) emsgroot("'nburnin' must be less than 'nsample'");
	if(nburnin == UNSET) emsgroot("'nburnin' must be set");
}


/// Finds the thinning of samples when using PMCMC or MC3
void Inputs::find_nthin(unsigned int &nthin, const unsigned int nsample)
{
	auto def = 1u;
	if(nsample != UNSET){
		def = (unsigned int)(nsample/1000.0);
		if(def == 0) def = 1;
	}
		
	nthin = find_positive_integer("nthin",def);
	if(nthin >= nsample/10) emsgroot("'nthin' is too large");
	if(nthin < 1) emsgroot("'nthin' must be 1 or above");
}


/// Finds the number of burnin steps (PMCMC)
void Inputs::find_nquench(unsigned int &nquench, const unsigned int &nburnin)
{
	nquench = find_positive_integer("nquench",nburnin/2);         // Sets the period of quenching (MC3)
	if(nquench >= nburnin) emsgroot("'nquench' must be less than 'nburnin'");
}


/// Finds the cutoff in EF or the fraction of samples accepted (for ABC)
void Inputs::find_cutoff(double &cutoff, double &cutoff_frac)
{
	cutoff = find_double("cutoff",UNSET);
	cutoff_frac = find_double("cutoff_frac",UNSET);
	
	if(cutoff == UNSET && cutoff_frac == UNSET) emsgroot("Either 'cutoff' or 'cutoff_frac' needs to be set.");
	if(cutoff != UNSET && cutoff_frac != UNSET) emsgroot("'cutoff' and 'cutoff_frac' cannot both be set.");
	
	if(cutoff_frac != UNSET){
		if(cutoff_frac <= 0) emsgroot("'cutoff_frac' must be positive");
		if(cutoff_frac > 1) emsgroot("'cutoff_frac' must be less than one.");
	}
	else{
		if(cutoff <= 0) emsgroot("The error function cutoff must be positive.");
	}
}


/// Finds the proposal size when performing ABCSMC 
void Inputs::find_propsize(double &propsize)
{
	propsize = find_double("prop_size",1);
	
	if(propsize <= 0) emsgroot("'prop_size' must be positive.");
}


/// Finds the fraction of samples accepted per generation (for ABCSMC)
void Inputs::find_cutoff_frac(double &cutoff_frac, double &cutoff_frac_init)
{
	cutoff_frac = find_double("cutoff_frac",0.3);	
	if(cutoff_frac <= 0) emsgroot("'cutoff_frac' must be positive");
	if(cutoff_frac > 1) emsgroot("'cutoff_frac' must be less than one.");
	cutoff_frac_init = find_double("cutoff_frac_init",cutoff_frac);	
}


/// Finds the quenching factor for the PAIC algorithm
void Inputs::find_quench_factor(double &quench_factor)
{
	quench_factor = find_double("quench_factor",1.0);
}


/// Finds the way in which the output is displayed
void Inputs::find_outputprop(OutputProp &prop)
{
	prop.type = KDE; prop.nbin = 200; prop.h = 10; 
	
	if(basedata->contains("output_prop")){
		auto pr = basedata->open("output_prop",used); pr.set_used();
		
		auto pd = pr.stringfield("probdist","");
		if(pd != ""){
			if(pd == "kde") prop.type = KDE;
			else{
				if(pd == "bin") prop.type = BINNING;
			}
		}
		
		auto nbin = pr.stringfield("nbin","");
		if(nbin != "") prop.nbin = get_int(nbin,"In 'output_prop' for 'nbin'");
		
		auto h = pr.stringfield("h","");
		if(h != "") prop.h = get_double_positive(h,"In 'output_prop' for 'h'");
		pr.check_used("output_prop");
	}
}

/// Finds if parameter values are plotted onto graphs
void Inputs::find_plot_param_values(bool &plot_param_values)
{
	plot_param_values = false; 
	auto plot_param_values_str = find_string("plot_param_values","false");
	if(plot_param_values_str == "true") plot_param_values = true;
	else{
		if(plot_param_values_str != "false") emsgroot("'plot_param_values' must be 'true' or 'false'");
	}
	
	if(get_siminf() == SIMULATE) plot_param_values = true;
}

/// Finds any modifications made to the model
void Inputs::find_modification(const Details &details, vector <Modification> &modification)
{
	if(basedata->contains("modification")){
		const auto pr = basedata->open("modification",used);
		
		for(auto j = 0u; j < pr.size(); j++){
			Modification cofa;
			auto cf = pr[j]; cf.set_used();
		
			cofa.start = 0;
			auto start = cf.stringfield("start","");
			if(start != "") cofa.start = details.gettime(start,"In 'modification' for 'start'")-details.start;
			
			cofa.end = details.period;
			auto end = cf.stringfield("end", "");
			if(end != "") cofa.end = details.gettime(end,"In 'modification' for 'end'")-details.start;
		
			cofa.strain_str = cf.stringfield("strain","");
		
			cofa.geo_filt = cf.stringfield_allowcomma("geo_filt",""); 
			
			cofa.democats_filt = cf.stringfield_allowcomma("democats_filt",""); 
		
			auto type = cf.stringfield("type","");
			if(type != ""){
				if(type == "trans_rate_fac"){
					cofa.type = CF_TRANS_RATE;
					cofa.trans_str = cf.stringfield("trans","In 'modification'");
					
					auto fac = cf.stringfield("factor","In 'modification'");
					cofa.factor = get_double_positive(fac,"In 'modification' for 'factor'"); 
				}
				else{
					if(type == "beta_fac"){
						cofa.type = CF_BETA;
						
						auto fac = cf.stringfield("factor","In 'modification'");
						cofa.factor = get_double_positive(fac,"In 'modification' for 'factor'"); 
					}
					else{
						if(type == "efoi_fac"){
							cofa.type = CF_EFOI;
							
							auto fac = cf.stringfield("factor","In 'modification'");
							cofa.factor = get_double_positive(fac,"In 'modification' for 'factor'"); 
						}
						else{
							if(type == "spline_fac"){
								cofa.type = CF_SPLINEFAC;
								cofa.spline_name_str = cf.stringfield("name","In 'modification'");
								
								auto fac = cf.stringfield("factor","In 'modification'");
								cofa.factor = get_double_positive(fac,"In 'modification' for 'factor'");
							}
							else{
								if(type == "spline_set"){
									cofa.type = CF_SPLINESET;
									cofa.spline_name_str = cf.stringfield("name","In 'modification'");
								
									auto fac = cf.stringfield("value","In 'modification'");
									cofa.factor = get_double_positive(fac,"In 'modification' for 'factor'");
								}
								else{
									emsgroot("'"+type+"' in 'modification' is not recognised (it should be 'trans_rate_fac' or 'efoi_fac').");
								}
							}
						}
					}
				}
			}
			
			cf.check_used("modification");
			modification.push_back(cofa);
		}
	}
}


/// Gets a boolean from a string
bool Inputs::get_bool(const string st, const string em) const
{
	if(st == "on") return true;
	if(st != "off" && st != "") emsgroot(em+" '"+st+"' must either be 'on' or 'off'");
	return false;
}


/// Find properties of the mcmc update
void Inputs::find_mcmc_update(MCMCUpdate &mcmc_update)
{
	mcmc_update.full_mvn = false;
	mcmc_update.mvn = false;
	mcmc_update.single = true;
	mcmc_update.dist_R_joint = true;
	mcmc_update.demo_spec = true;
	mcmc_update.mean_time = true;
	mcmc_update.neighbour = true;
	mcmc_update.joint = true;
	mcmc_update.mvn_multiple = true;
	mcmc_update.multiple_factor = 2;
	
	if(basedata->contains("mcmc_update")) {
		auto mup = basedata->open("mcmc_update",used); mup.set_used();
		mcmc_update.full_mvn = get_bool(mup.stringfield("full_mvn",""),"In 'mcmc_update'");
		mcmc_update.mvn = get_bool(mup.stringfield("mvn",""),"In 'mcmc_update'");
		mcmc_update.single = get_bool(mup.stringfield("single",""),"In 'mcmc_update'");
		mcmc_update.dist_R_joint = get_bool(mup.stringfield("dist_R_joint",""),"In 'mcmc_update'");
		mcmc_update.demo_spec = get_bool(mup.stringfield("demo_spec",""),"In 'mcmc_update'");
		mcmc_update.mean_time = get_bool(mup.stringfield("mean_time",""),"In 'mcmc_update'");
		mcmc_update.neighbour = get_bool(mup.stringfield("neighbour",""),"In 'mcmc_update'");
		mcmc_update.joint = get_bool(mup.stringfield("joint",""),"In 'mcmc_update'");
		mcmc_update.mvn_multiple = get_bool(mup.stringfield("mvn_multiple",""),"In 'mcmc_update'");
		auto num = mup.stringfield("multiple_factor",""); 
		if(num != ""){
			mcmc_update.multiple_factor = get_int(num,"In 'mcmc_update' and 'multiple_factor'");
		}
		mup.check_used("mcmc_update");
	}
}


void Inputs::find_area_plot(AreaPlot &area_plot)
{
	if(basedata->contains("area_plot")){
		auto mup = basedata->open("area_plot",used); mup.set_used();
	
		area_plot.boundfile = mup.stringfield("boundary", "");
		area_plot.xcol = mup.stringfield("x_column","");
		area_plot.ycol = mup.stringfield("y_column","");
		if(area_plot.boundfile != "" && (area_plot.xcol != "" || area_plot.ycol != "")){
			emsgroot("When specifying 'area_plot' cannot set 'boundary' and 'x_column'/'y_column'");
		}		
		
		if(area_plot.boundfile == ""){
			if(area_plot.xcol == "" || area_plot.ycol == ""){
				emsgroot("When specifying 'area_plot' either 'boundary' or 'x_column'/'y_column' must be set");
			}
		}
		
		area_plot.project = UNIFROM_PROJ;
		
		auto pro = mup.stringfield("projection", ""); 
		if(pro != ""){
			if(pro == "equirectangular") area_plot.project = EQUI_PROJ;
			else{
				if(pro == "uniform") area_plot.project = UNIFROM_PROJ;
				else emsgroot("In 'area_plot' the value '"+pro+"' is not recognised");
			}
		}
		mup.check_used("area_plot");
	}
}


/// Updates which TOML keys have been used
void Inputs::addused(const string key, vector <UsedTomlKey> &used)
{
	for(auto &us : used){ if(us.name == key) us.used = true;}
}
		

/// Prints all the commands which are not used
void Inputs::print_commands_not_used() const
{
	vector <string> not_used; 
	for(const auto &us : used){ if(us.used == false) not_used.push_back(us.name);}
	
	if(not_used.size() > 0){
		cout << "WARNING: The TOML ";
		
		switch(not_used.size()){
			case 1: cout << "command '" << not_used[0] << "' is not used."; break;
			default:
				cout << "commands ";
				for(auto i = 0u; i < not_used.size(); i++){
					if(i != 0){
						if(i < not_used.size()-1) cout << ", "; else cout << " and ";
					}
					cout << "'" << not_used[i] << "'";		
				}
				cout << " are not used." << endl;
				break;
		}
		cout << endl << endl;
	}
}
