/// Functions related to making a MVN approximiation to the populations

//#include <math.h>  

using namespace std;

#include "inputs.hh"

/// Loads up the compartmental model and provides a list of compartments and transitions
void Inputs::load_compartmental_model(vector <Compartment> &comp, vector <Transition> &trans, vector <CompartmentName> &comp_name, const vector <DemographicCategory> &democat, const vector < vector<unsigned int> > &democatpos, const Mode mode)
{
	add_comps(comp,democat,democatpos);                                  // Adds compartments to the model

	add_trans(trans,comp,democat,democatpos);                            // Adds transitions to the model

	if(mode == ML_INF) split_merged(comp,trans);                         // Expands model to ensure no merging 
	 
	add_comp_name(comp,comp_name);                                       // Generates unique compartment names
}
	

/// Adds comaprtments to the model
void Inputs::add_comps(vector <Compartment> &comp, const vector <DemographicCategory> &democat, const vector < vector<unsigned int> > &democatpos)
{
	vector <string> name; 
	vector <ParamSpec> inf;
	vector < vector <ParamSpec> > mean_spec;
	vector <unsigned int>	shape;
	vector <string> mean_dep;
	
	find_compartments(name,mean_spec,mean_dep,shape,inf,democatpos,democat); 

	for(auto c = 0u; c < name.size(); c++){
		auto j = 0u; while(j < comp.size() && comp[j].name != name[c]) j++;
		if(j < comp.size()) emsgroot("The compartment '"+name[c]+"' cannot be defined multiple times.");
		
		auto k = shape[c]; 
		if(k == UNSET){
			add_compartment(comp,name[c],mean_spec[c],mean_dep[c],inf[c],UNSET,UNSET);	
		}
		else{
			for(auto i = 0u; i < k; i++){
				add_compartment(comp,name[c],mean_spec[c],mean_dep[c],inf[c],i,k);	
			}
		}
	}
}

	
/// Adds a compartment to the model
void Inputs::add_compartment(vector <Compartment> &comp, const string& name, const vector <ParamSpec> &mean_spec, const string &mean_dep, const ParamSpec &inf, const unsigned int num, const unsigned int shape)
{
	Compartment co;
	co.name = name;
	co.num = num;
	co.name_num = name; if(num > 0) co.name_num += to_string(num);
	co.shape = shape;
	co.mean_dep = mean_dep;
	co.inf_spec = inf;

	co.mean_spec = mean_spec;
	
	comp.push_back(co);	
}


/// Adds transitions to the model
void Inputs::add_trans(vector <Transition> &trans, vector <Compartment> &comp, const vector <DemographicCategory> &democat, const vector < vector<unsigned int> > &democatpos)
{
	vector <string> from, to;                            
	vector <vector <ParamSpec> > prob_spec;
	vector <string> prob_dep;
	vector <TransInf> transinf;
	
	find_transitions(from,to,transinf,prob_spec,prob_dep,democatpos,democat);
	
	for(auto c = 0u; c < comp.size(); c++){                            // Adds internal Erlang transitions
		const auto &co = comp[c];
		if(co.shape != UNSET && co.num < co.shape - 1){
			auto cc = 0u; while(cc < comp.size() && !(comp[cc].name == co.name && comp[cc].num == co.num+1)) cc++;
			if(cc == comp.size()) emsgEC("Input Comps",2);
			add_internal_transition(trans,comp,c,cc);
		}
	}
	
	for(auto tr = 0u; tr < from.size(); tr++){                         // Adds transitions between compartments
		add_transition(trans,comp,from[tr],to[tr],transinf[tr],prob_spec[tr],prob_dep[tr]);
	}
}


/// Adds an internal transition to the model (for the Erlang distribution)
void Inputs::add_internal_transition(vector <Transition> &trans, vector <Compartment> &comp, const unsigned int c_from, const unsigned int c_to)
{
	Transition tr;
	
	string st =  comp[c_from].name; if(comp[c_from].num > 0) st += to_string(comp[c_from].num);
	st += "->";
	st += comp[c_to].name; if(comp[c_to].num > 0) st += to_string(comp[c_to].num);
	tr.name = st;
	tr.name_file = "";
	tr.from = c_from;
	tr.to = c_to;
	tr.inf = TRANS_NOTINFECTION;	
	tr.comptrans_ref = comp[tr.from].trans.size(); comp[tr.from].trans.push_back(trans.size());
	tr.prob_dep = "";
	
	trans.push_back(tr);
}


/// Adds a transition to the model
void Inputs::add_transition(vector <Transition> &trans, vector <Compartment> &comp, const string &from, const string &to, const TransInf inf, const vector <ParamSpec> &prob_spec, const string prob_dep)
{
	Transition tr;
	
	//tr.name = from+"→"+to;
	tr.name = from+"->"+to;
	tr.name_file = from+"-"+to;
	
	tr.from = get_compartment(comp,from,LAST);
	if(tr.from == UNSET) emsgroot("In 'trans' cannot find the 'from' compartment '"+from+"' for the transition");

	tr.inf = inf;
	if(from == to) emsgroot("In 'trans' the 'from' and 'to' values cannot be the same");
		
	if(inf == TRANS_INFECTION) tr.comptrans_ref = UNSET;
	else{
		tr.comptrans_ref = comp[tr.from].trans.size(); comp[tr.from].trans.push_back(trans.size());
	}
	
	tr.to = get_compartment(comp,to,FIRST);
	if(tr.to == UNSET) emsgroot("In 'trans' cannot find the 'to' compartment '"+to+"' for the transition");	
	
	tr.prob_spec = prob_spec;
	tr.prob_dep = prob_dep;
	
	trans.push_back(tr);
}


/// For ML approaches it is neccesary that brahches do not merge into the same compartment
void Inputs::split_merged(vector <Compartment> &comp, vector <Transition> &trans)
{
	bool flag;
	do{
		flag = false;
		
		for(const auto &tra : trans){
			for(auto &tra2 : trans){
				if(tra.to == tra2.to && tra.from != tra2.from){
					flag = true;
					
					auto c = tra.to;
					
					tra2.to = comp.size();
					comp.push_back(comp[c]);
									
					vector <CompCopy> cc_list;
					
					CompCopy ccopy;
					ccopy.c = c;
					ccopy.c_copy = comp.size()-1;
					ccopy.index = 0;
					cc_list.push_back(ccopy);
					
					while(cc_list.size() > 0){
						auto &cc = cc_list[cc_list.size()-1];
						
						auto c = cc.c;
						auto c_copy = cc.c_copy;
						
						if(cc.index < comp[c].trans.size()){
							auto tr = comp[c].trans[cc.index];
							auto c_to = trans[tr].to;
					
							auto c_new = comp.size();
							auto tr_new = trans.size();
								
							trans.push_back(trans[tr]);
							
							comp.push_back(comp[c_to]);
				
							trans[tr_new].from = c_copy;
							trans[tr_new].to = c_new;
							comp[c_copy].trans[cc.index] = tr_new;

							cc.index++;
							
							CompCopy ccopy;
							ccopy.c = c_to;
							ccopy.c_copy = c_new;
							ccopy.index = 0;
							cc_list.push_back(ccopy);
						}
						else{
							cc_list.pop_back();
						}
					}
				}
				if(flag == true) break;
			}
			if(flag == true) break;
		}
		//break;
	}while(flag == true);
}


/// Gerenrates a list of
void Inputs::add_comp_name(vector <Compartment> &comp, vector <CompartmentName> &comp_name)
{
	for(auto co = 0u; co < comp.size(); co++){
		auto name = comp[co].name;
		
		auto j = 0u; while(j < comp_name.size() && comp_name[j].name != name) j++;
		if(j == comp_name.size()){
			CompartmentName cn; cn.name = name;
			comp_name.push_back(cn);
		}
		comp_name[j].comp.push_back(co);
	}
}


/// Gets a comparment from it's name 
unsigned int Inputs::get_compartment(const vector <Compartment> &comp, const string compname, const ErlangPos pos) const
{
	for(auto c = 0u; c < comp.size(); c++){
		if(compname == comp[c].name){
			if(comp[c].num == UNSET || (pos == FIRST && comp[c].num == 0) || (pos == LAST && comp[c].num == comp[c].shape-1)){
				return c;
			}
		}
	}
	return UNSET;
}

	
/// Returns a list of compartment names along with the susceptible compartment (used for reading 'init_pop' file)
vector <string> Inputs::find_compartment_names(unsigned int &co_sus, const vector <DemographicCategory> &democat, const vector < vector<unsigned int> > &democatpos)
{
	vector <Transition> trans;                          // Stores model transitions
	vector <Compartment> comp;	                        // Stores model compartments
	vector <CompartmentName> comp_name;                 // Groups compartments with the same name
		 
	load_compartmental_model(comp,trans,comp_name,democat,democatpos,mode());
		
	co_sus = UNSET;
	for(auto i = 0u; i < trans.size(); i++){
		if(trans[i].inf == TRANS_INFECTION){
			if(co_sus != UNSET) emsgroot("Cannot have multiple infection transitions");
			co_sus = trans[i].from;
		}
	}
	
	if(co_sus == UNSET){
		emsgroot("An infection transition must exist");
	}
	
	vector <string> name;
	for(const auto &co : comp) name.push_back(co.name);

	return name;
}


/// Finds 'trans'
void Inputs::find_transitions(vector <string> &from, vector <string> &to, vector <TransInf> &transinf, vector < vector <ParamSpec> > &prob_spec, vector <string> &prob_dep, const vector < vector<unsigned int> > &democatpos, const vector <DemographicCategory> &democat)
{
	if(basedata->contains("trans")){
		const auto transin = basedata->open("trans",used);
		for(auto j = 0u; j < transin.size(); j++){
			auto trans = transin[j]; trans.set_used();
			
			auto dep = get_dep(trans,"trans","prob_param","prob_value","prob_prior");
			                  
			auto fr_temp = trans.stringfield("from","In 'trans'");
			auto to_temp = trans.stringfield("to","In 'trans'");
			//auto name = fr_temp+"→"+to_temp;
			auto name = fr_temp+"->"+to_temp;
			
			vector <ParamSpec> probsp;
			if(trans.contains("prob_param") || trans.contains("prob_value") || trans.contains("prob_prior")){
				probsp = find_ct_param_spec(trans,"trans","prob_param","prob_value","prob_prior",dep,democatpos,democat,name+" branch prob.");
			}
				
			TransInf tinf = TRANS_NOTINFECTION; 
			auto inf = trans.stringfield("infection","");
			if(inf != ""){
				if(inf=="yes") tinf = TRANS_INFECTION;
				else{
					if(inf != "no") emsgroot("In 'comps' the value of 'infection' can either be 'yes' or 'no' or omitted.");
				}
			}
			
			from.push_back(fr_temp); 
			to.push_back(to_temp);
			transinf.push_back(tinf);
			prob_spec.push_back(probsp);
			prob_dep.push_back(dep);
			trans.check_used("trans");
		}
	}
	else emsgroot("The input file must contain transition definitions through 'trans'.");
}

/// Finds 'comps'
void Inputs::find_compartments(vector <string> &name, vector < vector <ParamSpec> > &mean_spec, vector <string> &mean_dep, vector <unsigned int> &k, vector <ParamSpec> &infectivity, const vector < vector<unsigned int> > &democatpos, const vector <DemographicCategory> &democat)
{
	if(basedata->contains("comps")) {
		auto compsin = basedata->open("comps",used);
		for(auto j = 0u; j < compsin.size(); j++){
			auto comps = compsin[j]; comps.set_used();

			auto nam = comps.stringfield("name","In 'comps'");

			auto root = "comps' name='"+nam;
			
			auto par = comps.stringfield("mean_param","");
			
			auto dep = get_dep(comps,root,"mean_param","mean_value","mean_prior");

			auto k_temp = UNSET;
			
			auto dist = comps.stringfield("dist","");
			if(dist != ""){
				if(dist == "Exp") k_temp = 1;
				else{
					if(dist == "Erlang"){
						auto k_str = comps.stringfield("k","");
						
						if(k_str == "") emsgroot("In 'comps' compartment '"+nam+"' an integer shape parameter 'k' must be set.");
						k_temp = get_int(k_str,"In 'comps' compartment '"+nam+"' the integer 'k'");
					}
					else{
						emsgroot("In 'comps' for the '"+nam+"' compartment the value 'dist="+dist+"' is not recognised.");
					}
				}
			}
			
			auto meansp = find_ct_param_spec(comps,root,"mean_param","mean_value","mean_prior",dep,democatpos,democat,nam+" mean occupancy");
			
			auto inf = ps_zero();
			string inf_tr = "";
			
			if(comps.contains("inf_param") || comps.contains("inf_value") || comps.contains("inf_prior")){		
				auto infvec = get_paramspec(comps,"comps","inf_param","inf_value","inf_prior","","",1,false);
				
				if(param_not_set(infvec)) infvec[0].name = nam+" Infectivity";    // Sets parameter name (if uninitialised)
				inf = infvec[0];
			}
	
			name.push_back(nam);
			mean_spec.push_back(meansp);
			mean_dep.push_back(dep);
			k.push_back(k_temp);
			infectivity.push_back(inf);
			
			comps.check_used("comps");
		}
	}
	else{ emsgroot("The input file must contain compartment definitions through 'comps'");}
}

