// This file defines the compartmental model and is used to simulate compartmental dynamics after an individual becomes infected

#include <string>
#include <sstream>
#include <iostream>
#include <fstream>
#include <cmath>
 
#include "math.h"

using namespace std;

#include "model.hh"
#include "matrix.hh"

/// Initialises the model 
Model::Model(Inputs &inputs, const Details &details, Data &data, Mpi &mpi) : details(details), data(data), inputs(inputs), mpi(mpi)
{
	load_model();

	complete_datatables();
	
	setup_modification();
	
	prior_order();
	
	set_parameter_type();
	
	if(false) print_parameter_types();
	
	if(mpi.core == 0){
		cout << "Number of  model parameters: " << param.size()-2 << endl; 
		if(details.siminf == INFERENCE) cout << "Number of parameters to be inferred: " << param_not_fixed.size() << endl; 
	}
}


/// Loads the model from the input file
void Model::load_model()
{
	inputs.load_compartmental_model(comp,trans,comp_name,data.democat,data.democatpos,details.mode);
	
	add_compartmental_model_parameters();                            // Adds compartmental model parameters
	
	compartmental_model_check();                                     // Checks model correctly specified
	
	setup_infectivity();                                             // Sets up infection transition
	
	add_splines();                                                   // Adds splines to the model

	add_region_effect();                                             // Adds regional effects

	add_probreach();                                                 // Adds output giving e.g. ifr
	
	set_dirichlet();                                                 // Sets up the dirichlet distributions

	for(unsigned int c = 0; c < data.ncovar; c++){                   // Adds parameters for area covariates
		covariate_param.push_back(add_parameter(data.covar[c].ps,COVAR_PARAM));
	}

	for(auto &str : data.strain){                                    // Adds parameters giving Rfactors for strains
		str.Rfactor_param = add_parameter(str.Rfactor_ps,RFACTOR_PARAM);
	}

	set_infected_uninfected();                                       // Divides model into infected and uninfected compartments
	
	for(auto th = 0u; th < param.size(); th++){                      // Lists all parameters which are not fixed
		if(param[th].priortype != FIXED_PRIOR) param_not_fixed.push_back(th);
	}
	
	if(param_not_fixed.size() == 0 && details.siminf == INFERENCE){
		emsg("For inference the model must contain unspecified parameters");
	}
	
	if(false){                                                       // Print model parameters
		for(auto th = 0u; th < param.size(); th++){      
			cout <<th << " " << param[th].name << endl;
		}
		emsg("Model Parameters");
	}
}


/// Adds parameters associated with the compartmental model
void Model::add_compartmental_model_parameters()
{
	for(auto &co : comp){
		if(co.mean_spec.size() > 0){
			if(co.mean_spec.size() != data.ndemocatpos) emsgEC("Model",1);
			
			for(auto dp = 0u; dp < data.ndemocatpos; dp++){
				co.param_mean.push_back(add_parameter(co.mean_spec[dp],DISTVAL_PARAM));
			}
		}
		co.infectivity_param = add_parameter(co.inf_spec,INF_PARAM);
	}
	
	for(auto &tr : trans){
		if(tr.prob_spec.size() > 0){                                       // Adds branching probabilties
			if(tr.prob_spec.size() != data.ndemocatpos) emsgEC("Input Comps",3);
			
			for(auto dp = 0u; dp < data.ndemocatpos; dp++){
				tr.param_prob.push_back(add_parameter(tr.prob_spec[dp],BRANCHPROB_PARAM));
			}
		}
	}
}


/// Sets the uninfected and infect compartments in the model.
void Model::set_infected_uninfected()
{ 
	vector <bool> flag_forward(comp.size());
	for(auto co = 0u; co < comp.size(); co++) flag_forward[co] = false;
		
	vector <unsigned int> list;                                      // List all states which link to infection
	auto c = trans[infection_trans].to;	
	list.push_back(c); flag_forward[c] = true; 
		
	for(auto i = 0u; i < list.size(); i++){
		c = list[i];	
		for(auto tr : comp[c].trans){
			if(trans[tr].inf == TRANS_NOTINFECTION){
				auto cc = trans[tr].to;
				if(flag_forward[cc] == false && cc != trans[infection_trans].from){
					list.push_back(cc); flag_forward[cc] = true; 
				}
			}
		}
	}

	vector <bool> flag_back(comp.size());                             // Lists all states which go on become infectious
	for(auto co = 0u; co < comp.size(); co++) flag_back[co] = false;

	for(auto co = 0u; co < comp.size(); co++){
		if(param[comp[co].infectivity_param].name != "zero"){
			if(flag_forward[co] != true){
				emsgroot("The infectious compartment '"+comp[co].name+"' derives from the wrong infection transition."); 
			}

			vector <unsigned int> list;
			list.push_back(co); flag_back[co] = true;
			
			for(auto i = 0u; i < list.size(); i++){
				auto c = list[i]; 
				for(auto tr = 0u; tr < trans.size(); tr++){
					auto tra = trans[tr];
					if(tra.inf == TRANS_NOTINFECTION && tra.to == c){
						auto cc = tra.from;
						if(flag_back[cc] == false){
							list.push_back(cc); flag_back[cc] = true;
						}
					}
				}
			}
		}
	}

	inf_state_ref.resize(comp.size());
	for(auto co = 0u; co < comp.size(); co++){
		if(flag_forward[co] == true && flag_back[co] == true){
			inf_state_ref[co] = inf_state.size();
			inf_state.push_back(co);
		}
		else inf_state_ref[co] = UNSET;
	}

	if(false){
		cout << endl << "INFECT " << ":" << endl;
		for(auto co : inf_state){
			cout << comp[co].name_num  << " infectious state" << endl;
		}	

		for(auto co = 0u; co < comp.size(); co++){
			cout << comp[co].name_num  << " " << inf_state_ref[co];
			cout << "  Forward:"; if(flag_forward[co] == true) cout << "True"; else cout << "False";
			cout << "  Back:"; if(flag_back[co] == true) cout << "True"; else cout << "False";
			cout << "    ref" << endl;
		}
		emsg("Done");
	}
}


/// Sets up infection transitions and which compartments cause infectivity
void Model::setup_infectivity()
{
	infection_trans = UNSET;
	for(auto tr = 0u; tr < trans.size(); tr++){                      // Calculates all the infection transitions
		if(trans[tr].inf == TRANS_INFECTION){
			if(infection_trans == UNSET) infection_trans = tr;
			else emsgroot("In 'trans' only one transition can have 'infection=\"yes\"'");
		}
	}
	if(infection_trans == UNSET) emsgroot("In 'trans' one transition must be set to 'infection=\"yes\"'");
	
	start_compartment = trans[infection_trans].from;                 // Determines the start compartment for individuals
}


/// Add parameters for regional effects
void Model::add_region_effect()
{	
	vector <ParamSpec> ps_vec; 
	ParamSpec sigma;
	
	inputs.find_region_effect(ps_vec,sigma,data.area);
	
	if(ps_vec.size() == 0){ region_effect.on = false; return;}

	region_effect.on = true; 

	if(ps_vec.size() != data.narea) emsgEC("Model",4);
	
	if(ps_vec[0].value == "sample"){                                // Samples regional effects
		string file = "simulated_region_effect.csv";
		auto file_full = details.output_directory+"/"+file;
		vector <double> RE(data.narea);
				
		if(mpi.core == 0){
			if(details.mode == SIM || details.mode == MULTISIM){
				if(sigma.value == "") emsgroot("In 'region_effect' when using 'sample' the value of 'sigma_value' must be set.");
				auto sd = get_double(sigma.value,"In 'region_effect' for 'sigma_value'");
			
				cout << "The values for sampled regional effects are shown in '" << file_full << "'" << endl;
				ofstream reout(file_full.c_str());
				if(!reout) emsg("The file '"+file_full+"' could not be openned");
				
				reout << "Area,Parameter,Region effect" << endl;
				for(auto c = 0u; c < data.narea; c++){
					RE[c] = normal_sample(0,sd);
					reout << data.area[c].code << "," << ps_vec[c].name << "," << RE[c] << endl;
				}
			}
			else{
				ifstream file_check(file_full);
				if(file_check){
					file_check.close();
		
					auto tab = data.get_table(file,details.output_directory);
					auto column = data.get_table_column("Region effect",tab);
					
					if(column.size() != data.narea) emsgEC("Model",5);
					for(auto c = 0u; c < data.narea; c++){
						RE[c] = get_double(column[c],"In 'region_effect' from file '"+file+"'");
					}
				}
			}
		}
	
		mpi.bcast(RE);
		for(auto c = 0u; c < data.narea; c++) ps_vec[c].value = to_string(RE[c]); 
	}
	
	if(sigma.value == "") region_effect.sigma_param = UNSET;
	else region_effect.sigma_param = add_parameter(sigma,SIGMA_PARAM);
	
	auto param_beg = param.size();
	for(auto c = 0u; c < data.narea; c++){
		region_effect.area_param.push_back(add_parameter(ps_vec[c],RE_PARAM));
	}
	
	for(auto th = param_beg; th < param.size(); th++) region_effect.param_list.push_back(th);
}


/// Add outputs which give the probability of reaching a certain compartment
void Model::add_probreach()
{
	vector <string> name, comp;
	
	inputs.find_probreach(name,comp);

	for(auto i = 0u; i < name.size(); i++){
		ProbReach pr;
		pr.name = name[i];
		emsg("TO DO");
		//pr.comp = get_compartment(comp[i],FIRST);		
		prob_reach.push_back(pr);
	}
}


/// Finds a parameter from a string
unsigned int Model::find_param(string name, string em) const 
{
	auto th = 0u; while(th < param.size() && param[th].name != name) th++;
	if(th == param.size()) emsg(em+" cannot find the parameter '"+name+"'");
	return th;
}


/// Sets any prior order in parameters e.g. so one parameter must be bigger than another
void Model::prior_order() 
{
	auto prior_order = inputs.find_string("prior_order","");
	if(prior_order == "") return;
	
	auto spl = split(prior_order,'|');
	
	for(auto st : spl){
		if(st == "") emsgroot("In 'prior_order' the expression '"+prior_order+"' is not recognised");
		
		auto spl2 = split(st,'>');
		if(spl2.size() > 2) emsgroot("In 'prior_order' the expression '"+st+"' is not recognised");
		
		if(spl2.size() == 2){
			auto th1 = find_param(spl2[0],"In 'prior_order'");
			auto th2 = find_param(spl2[1],"In 'prior_order'");
			if(th1 == th2) emsgroot("In 'prior_order' cannot have the expression '"+st+"'.");
			
			param[th1].greaterthan.push_back(th2);
			param[th2].lessthan.push_back(th1);
		}		
		else{
			auto spl2 = split(st,'<');
			if(spl2.size() > 2) emsgroot("In 'prior_order' the expression '"+st+"' is not recognised");
		
			if(spl2.size() == 2){
				auto th1 = find_param(spl2[0],"In 'prior_order'");
				auto th2 = find_param(spl2[1],"In 'prior_order'");
				if(th1 == th2) emsgroot("In 'prior_order' cannot have the expression '"+st+"'.");
			
				param[th1].lessthan.push_back(th2);
				param[th2].greaterthan.push_back(th1);
			}
			else emsgroot("In 'prior_order' the expression '"+st+"' is not recognised");
		}
	}
}


// Checks ordering specified in 'prior_order' is correct
bool Model::prior_order_correct(const vector <double> &paramv) const
{
	for(auto th = 0u; th < param.size(); th++){
		for(auto th2 : param[th].greaterthan){
			if(paramv[th] < paramv[th2]) return false; 
		}
							
		for(auto th2 : param[th].lessthan){
			if(paramv[th] > paramv[th2]) return false; 
		}
	}
	return true;
}					
							
	
/// Sets the types for different types of parameters
void Model::set_parameter_type() 
{
	parameter_type.resize(PARAMTYPEMAX);
	for(auto th = 0u; th < param.size(); th++){
		if(param[th].priortype == FIXED_PRIOR) parameter_type[param[th].type].push_back(th);
	}
	
	for(const auto &Rpl : Rspline_info){
		const auto &spl = spline[Rpl.spline_ref];
		for(const auto &po :spl.p){
			auto th = po.param;
			if(param[th].priortype != FIXED_PRIOR) add_vec(R_param,th);
		}
	}
}


/// Updates the references to the spline
void Model::update_R_ref(string &name, string &desc, string area_str, vector <unsigned int>& spl_ref, vector <SplineInfo> &spline_info) const
{	
	vector <unsigned int> area;
	if(area_str != ""){
		auto c = 0u; while(c < data.narea && area_str != data.area[c].code) c++;
		if(c == data.narea) emsgroot("In '"+name+"' cannot find 'area="+area_str+"'");
		
		area.push_back(c);
		
		if(data.narea > 1){
			name += " "+area_str;
			desc += " for area "+area_str;
		}
	}
	else{
		for(auto c = 0u; c < data.narea; c++) area.push_back(c);
	}
	
	for(auto c : area){
		if(spl_ref[c] != UNSET){
			emsgroot("For '"+name+"' the area '"+data.area[c].code+"' already has a spline specified");			
		}
		spl_ref[c] = spline_info.size();
	}	
	
	SplineInfo si;
	si.name = name;
	si.spline_ref = spline.size();
	si.area = area;
	si.democatpos_dist.resize(data.ndemocatpos_per_strain);
	auto total = 0.0;
	for(auto c : area){
		for(auto co = 0u; co < comp.size(); co++){
			for(auto dp = 0u; dp < data.ndemocatpos; dp++){
				si.democatpos_dist[dp%data.ndemocatpos_per_strain] += data.area[c].pop_init[co][dp];
				total += data.area[c].pop_init[co][dp];
			}
		}
	}

	for(auto dp = 0u; dp < data.ndemocatpos_per_strain; dp++) si.democatpos_dist[dp] /= total;
	si.total_pop = total;	

	spline_info.push_back(si);
}


/// Updates the references to the spline
void Model::update_efoi_ref(string &name, string &desc, string area_str, string strain_str, const vector <double> efoi_ad, vector < vector <unsigned int> > &spl_ref, vector <SplineInfo> &spline_info) const
{	
	vector <unsigned int> area;
	if(area_str != ""){
		auto c = 0u; while(c < data.narea && area_str != data.area[c].code) c++;
		if(c == data.narea) emsgroot("In '"+name+"' cannot find 'area="+area_str+"'");
		
		area.push_back(c);
		
		if(data.narea > 1){
			name += " "+area_str;
			desc += " for area "+area_str;
		}
	}
	else{
		for(auto c = 0u; c < data.narea; c++) area.push_back(c);
	}
	
	vector <unsigned int> sta;
	if(strain_str != ""){
		auto st = 0u; while(st < data.nstrain && strain_str != data.strain[st].name) st++;
		if(st == data.nstrain) emsgroot("In '"+name+"' cannot find 'strain="+strain_str+"'");
		
		sta.push_back(st);
		
		if(data.nstrain > 1){
			name += " "+strain_str;
			desc += " for strain "+strain_str;
		}
	}
	else{
		for(auto st = 0u; st < data.nstrain; st++) sta.push_back(st);
	}
	
	for(auto st : sta){
		for(auto c : area){
			if(spl_ref[st][c] != UNSET){
				emsgroot("For '"+name+"' the area '"+data.area[c].code+"' already has a spline specified");			
			}
			spl_ref[st][c] = spline_info.size();
		}
	}		
	
	SplineInfo si;
	si.name = name;
	si.spline_ref = spline.size(); si.efoi_agedist = efoi_ad;
	si.area = area;
	si.democatpos_dist.resize(data.ndemocatpos_per_strain);
	auto total = 0.0;
	for(auto c : area){
		for(auto co = 0u; co < comp.size(); co++){
			for(auto dp = 0u; dp < data.ndemocatpos; dp++){
				si.democatpos_dist[dp%data.ndemocatpos_per_strain] += data.area[c].pop_init[co][dp];
				total += data.area[c].pop_init[co][dp];
			}
		}
	}

	for(auto dp = 0u; dp < data.ndemocatpos_per_strain; dp++) si.democatpos_dist[dp] /= total;
	si.total_pop = total;	

	spline_info.push_back(si);
}
	
	
/// Adds splines to the model
void Model::add_splines()
{	
	vector <string> name_vec;
	vector < vector <ParamSpec> > ps_vec; 
	vector <ParamSpec> factor_param_vec; 
	vector < vector <unsigned int> > bp_vec;  
	vector <SmoothType> smooth_type_vec;
	vector <string> strain;
	vector <string> area;
	vector < vector <double> > efoi_agedist_vec;
	vector <double> empty;
	
	R_spl_ref.resize(data.narea);                                      // Sets up splines for R(t)
	for(auto c = 0u; c < data.narea; c++) R_spl_ref[c] = UNSET;
	
	inputs.find_spline("R_spline",name_vec,ps_vec,factor_param_vec,bp_vec,smooth_type_vec,strain,area,efoi_agedist_vec,empty,details);  
	
	if(ps_vec.size() == 0) emsgroot("A value for 'R_spline' must be set.");
	
	for(auto num = 0u; num < ps_vec.size(); num++){
		string desc = "The reproduction number";
		if(name_vec[num] == "") name_vec[num] = "R_t"; 
		if(efoi_agedist_vec[num].size() != 0) emsgroot("In 'R_spline' a value for 'age_dist' cannot be set");
		if(strain[num].size() != 0) emsgroot("In 'R_spline' a value for 'strain' cannot be set");

		update_R_ref(name_vec[num],desc,area[num],R_spl_ref,Rspline_info);
		
		create_spline(name_vec[num],desc,ps_vec[num],bp_vec[num],factor_param_vec[num],smooth_type_vec[num],R_EFOI_PARAM);
	}
	
	for(auto c = 0u; c < data.narea; c++){
		if(R_spl_ref[c] == UNSET) emsgroot("'R_spline' must be set for all areas and infection transitions.");
	}
	
	efoi_spl_ref.resize(data.nstrain);                                // Sets up splines for efoi
	for(auto st = 0u; st < data.nstrain; st++){
		efoi_spl_ref[st].resize(data.narea);
		for(auto c = 0u; c < data.narea; c++) efoi_spl_ref[st][c] = UNSET; 
	}

	auto efoi_agedist = data.agedist;	
	inputs.find_spline("efoi_spline",name_vec,ps_vec,factor_param_vec,bp_vec,smooth_type_vec,strain,area,efoi_agedist_vec,data.agedist,details);  // Sets up efoi(t)
	
	if(ps_vec.size() == 0) emsgroot("A value for 'efoi_spline' must be set.");

	for(auto num = 0u; num < ps_vec.size(); num++){
		string desc = "External force of infection";
		auto name = name_vec[num];
		if(name == "") name = "efoi(t)";
		
		auto str_st = strain[num];
		if(str_st != ""){
			if(data.nstrain > 1){ name += " "+str_st; desc += " for strain "+str_st;}
		}
		
		/*
		vector <unsigned int> str;
		auto str_st = strain[num];
		if(str_st != ""){
			auto i = 0u; while(i < data.nstrain && str_st != data.strain[i].name) i++;
			if(i == data.nstrain) emsgroot("The strain '"+str_st+"' is not recognised");
			str.push_back(i);
			if(data.nstrain > 1){ name += " "+str_st; desc += " for strain "+str_st;}
		}
		else{
			for(auto i = 0u; i < data.nstrain; i++) str.push_back(i);
		}
	
		for(auto st : str){
			*/
			
		update_efoi_ref(name,desc,area[num],str_st,efoi_agedist_vec[num],efoi_spl_ref,efoispline_info);
		
		create_spline(name,desc,ps_vec[num],bp_vec[num],factor_param_vec[num],smooth_type_vec[num],R_EFOI_PARAM);
	}
	
	for(auto st = 0u; st < data.nstrain; st++){
		for(auto c = 0u; c < data.narea; c++){
			if(efoi_spl_ref[st][c] == UNSET) emsgroot("'efoi_spline' must be set for all areas and strains.");
		}
	}
	
	geo_spline_ref = spline.size();                                     // Sets up spline for geograpic mixing
	inputs.find_spline("geo_mixing_modify",name_vec,ps_vec,factor_param_vec,bp_vec,smooth_type_vec,strain,area,efoi_agedist_vec,empty,details);
	
	if(ps_vec.size() == 0){                                             // If not set then sets up a spline with the value unity
		vector <unsigned int> bp;
		bp.push_back(0); bp.push_back(details.end-details.start);
		bp_vec.push_back(bp);
		
		vector <ParamSpec> ps_list; ps_list.push_back(ps_one()); ps_list.push_back(ps_one()); 
		ps_vec.push_back(ps_list);
			
		factor_param_vec.push_back(ps_one());
		
		smooth_type_vec.push_back(LOGSMOOTH);
		name_vec.push_back("UNSET");
	}
	
	if(ps_vec.size() != 1) emsgroot("Only one spline can represent 'geo_mixing'.");
	if(name_vec[0] == "") name_vec[0] = "m(t)";
	create_spline(name_vec[0],"Relative within / between region mixing.",ps_vec[0],bp_vec[0],factor_param_vec[0],smooth_type_vec[0],R_EFOI_PARAM);
		
	for(auto &mm : data.genQ.matmod){
		switch(mm.type){
		case ROW_COLUMN: case PERTURB:
			mm.spline_ref = spline.size();	
			
			create_spline(mm.name,mm.desc,mm.ps,mm.bp,ps_one(),mm.smoothtype,MODIFY_PARAM);
			break;
		}
	}
				
	// Spline mapping observation to data
	inputs.find_spline("obs_spline",name_vec,ps_vec,factor_param_vec,bp_vec,smooth_type_vec,strain,area,efoi_agedist_vec,empty,details);  
	
	for(auto num = 0u; num < ps_vec.size(); num++){
		if(name_vec[num] == "") emsgroot("Each spline in 'obs_spline' must have 'name' set");
		
		auto flag = false;
		for(const auto &dt : data.datatable){
			if(dt.factor_spline == name_vec[num]) flag = true;
		}
		
		if(flag == true){ 		
			string desc = "*"+name_vec[num]+"* is a factor mapping observation to data";
			if(efoi_agedist_vec[num].size() != 0) emsgroot("In 'obs_spline' a value for 'age_dist' cannot be set");
			if(strain[num].size() != 0) emsgroot("In 'obs_spline' a value for 'strain' cannot be set");

			create_spline(name_vec[num],desc,ps_vec[num],bp_vec[num],factor_param_vec[num],smooth_type_vec[num],OBS_PARAM);
		}
	}
	
	set_disc_spline_grad();
}


/// Creates a spline based on a given specification
void Model::create_spline(const string name, const string desc, const vector <ParamSpec> &ps, const vector<unsigned int> &bp, const ParamSpec param_fac, const SmoothType smooth_type, const ParamType type)
{
	if(false){
		cout << name << " name" << endl;
		for(auto i = 0u; i < bp.size(); i++){	
			cout <<  bp[i] << " " << ps[i].name << " " << ps[i].value << "  pro:" << ps[i].prior << " " << ps[i].smooth << " " << ps[i].factor << endl;
		}
	}
	
	Spline spli;	
	for(auto i = 0u; i < bp.size(); i++){	
		SplinePoint spl;
		spl.t = bp[i];
		spl.param = add_parameter(ps[i],type); 
		spl.multfac = ps[i].factor;
		spl.smooth_value = ps[i].smooth;
		spli.p.push_back(spl);
	}
	
	spli.param_factor = add_parameter(param_fac,type); 
	spli.name = name;
	spli.desc = desc;
	spli.smoothtype = smooth_type;
	
	spline.push_back(spli);
}

		
/// Samples from the prior for a given parameter
double Model::param_prior_sample(const unsigned int th, const vector <double> &paramv) const
{	
	switch(param[th].priortype){
			case FIXED_PRIOR:
				return get_val1(th,paramv);
				
			case UNIFORM_PRIOR:	
				{
					auto min = get_val1(th,paramv), max = get_val2(th,paramv);
					return min + ran()*(max - min);
				}
				
			case EXP_PRIOR:
				{
					auto rate = get_val1(th,paramv);
					return -log(ran())/rate;
				}
				
			case NORMAL_PRIOR:
				return normal_sample(get_val1(th,paramv),get_val2(th,paramv));
			
			case DIRICHLET_PRIOR: case DIRICHLET_ALPHA_PRIOR: case MDIRICHLET_PRIOR:
				return gamma_sample(get_val1(th,paramv),get_val2(th,paramv));
			
			case DIRICHLET_FLAT_PRIOR:
				return exp_sample(1);
						
			default:			
				emsgEC("model",6);
				break;
		}
}


/// Randomly samples the parameter values from the prior
vector <double> Model::sample_from_prior() const
{
	vector <double> paramv(param.size());
	
	auto loop = 0u;
	do{
		for(auto th = 0u; th < param.size(); th++) paramv[th] = param_prior_sample(th,paramv); 
	
		if(false){ for(auto th = 0u; th < param.size(); th++) cout << param[th].name << " " << paramv[th] << "sample" << endl;}
		
		spline_sample(paramv);                                        // Smoothing on the spline is accounted for 

		loop++;
	}while(loop < sample_try && prior_order_correct(paramv) == false);
	
	if(loop == sample_try) emsgroot("Could not get the correct ordering of parameters (please check 'prior_order')");
	
	return paramv;
}


/// Samples from the smoothing spline
void Model::spline_sample(vector <double> &paramv) const 
{
	for(const auto &spl : spline){  
		if(spl.smoothtype != NOSMOOTH){
			auto spli = spl.p;
			bool flag;
			auto loop = 0u;
			do{
				vector <bool> param_set(param.size());
				for(auto th = 0u; th < param.size(); th++) param_set[th] = false;
			
				flag = true;
				auto th1 = spli[0].param;

				if(param[th1].priortype != FIXED_PRIOR && param[th1].priortype != UNIFORM_PRIOR){
					if(details.mode != ML_INF){ 
						emsgroot("For smoothed splines only 'Fixed' and 'Uniform' priors can be used, otherwise the prior can't be sampled from.");
					}
				}
				
				paramv[th1] = param_prior_sample(th1,paramv);
				
				param_set[th1] = true;
				
				for(auto i = 0u; i < spli.size()-1; i++){
					auto th1 = spli[i].param;
					auto th2 = spli[i+1].param;
					
					if(th1 != th2 && param[th2].priortype != FIXED_PRIOR){
						if(spli[i].t == spli[i+1].t){
							paramv[th2] = param_prior_sample(th2,paramv);
						}
						else{
							switch(spl.smoothtype){
							case LOGSMOOTH:
								{
									auto sm = spli[i].smooth_value;
									if(sm == UNSET) emsgroot("In spline '"+spl.name+"' a value for 'smooth' must be set");
									paramv[th2] = (spli[i+1].multfac/spli[i].multfac)*paramv[th1]*exp(normal_sample(0,sm)); 
								}
								break;
								
							case SMOOTH:
								{
									auto sm = spli[i].smooth_value;
									if(sm == UNSET) emsgroot("In spline '"+spl.name+"' a value for 'smooth' must be set");
							
									if(spli[i+1].multfac != spli[i].multfac) emsgEC("Model",7);
									paramv[th2] = paramv[th1] + normal_sample(0,sm); 
								}
								break;
								
							default: emsgEC("Model",8);
							}
						
							switch(param[th1].priortype){
								case UNIFORM_PRIOR:
									{
										auto min = param[th2].val1, max = param[th2].val2; 
							
										if(paramv[th2] <= min || paramv[th2] >= max) flag = false;
									}
									break;

								case MDIRICHLET_PRIOR:
									if(paramv[th2] <= 0) flag = false;
									break;
							
								default:
									emsgroot("Prior not supported");
									break;
							}
							
							for(auto th : param[th2].greaterthan){
								if(param_set[th] == true && paramv[th2] < paramv[th]) flag = false; 
							}
							
							for(auto th : param[th2].lessthan){
								if(param_set[th] == true && paramv[th2] > paramv[th]) flag = false; 
							}
							
							param_set[th2] = true;
						}
					}
					
					if(flag == false) break;
				}
			
				loop++;
			}while(flag == false && loop < spline_sample_try);	
	
			if(loop == spline_sample_try) emsg("Spline smoothing is found to be not possible");
		}
	}
}


/// Gets the infectivity difference when going down a transition
double Model::get_infectivity_dif(const unsigned int tr, const vector <double> &paramv) const
{
	return paramv[comp[trans[tr].to].infectivity_param] - paramv[comp[trans[tr].from].infectivity_param];
}


/// Gets either a number or a reference to another parameter
void Model::get_prior_val(const string name, const string st, double &val, unsigned int &valparam)
{
	val = get_double_with_unset(st,"For parameter '"+name+"'"); valparam = UNSET;
	if(val == UNSET){ 
		auto th = 0u; while(th < param.size() && param[th].name != st) th++;
		if(th < param.size()){ val = UNSET; valparam = th;}
		else emsgroot("The expression '"+st+"' is not a number or the name of a parameter");
	}
}
	
	
/// Adds a parameter to the model
unsigned int Model::add_parameter(const ParamSpec ps, const ParamType type)
{
	PriorType pt = FIXED_PRIOR;
	double val1=UNSET, val2=UNSET, value=UNSET;
	unsigned int val1_param=UNSET, val2_param=UNSET;
	double dir_mean=UNSET, dir_sd=UNSET;
	
	if(ps.name == "") emsgroot("The parameter must have a name");
	
	switch(details.siminf){
	case SIMULATE:
		if(ps.value == "") emsgroot("The value for parameter '"+ps.name+"' must be set.");
		pt = FIXED_PRIOR;
		val1 = get_double_with_tobeset(ps.value,"For parameter '"+ps.name+"'");
		value = val1;
		break;
		
	case INFERENCE: case DATAVIEW:
		if(ps.value != "") value = get_double_with_tobeset(ps.value,"For parameter '"+ps.name+"'");
		
		if(ps.prior == "") emsgroot("The prior for parameter '"+ps.name+"' must be set.");
		
		auto spl = split(ps.prior,'(');
		if(spl.size() != 2){
			emsgroot("For parameter '"+ps.name+"' there was a problem with the prior expression '"+ps.prior+"'");
		}
		auto spl2 = split(spl[1],')');
		
		if(spl2.size() > 2){
			emsgroot("For parameter '"+ps.name+"' there was a problem with the prior expression '"+ps.prior+"'");
		}
		auto vals = split(spl2[0],',');
		
		if(toLower(spl[0]) == "fixed"){
			pt = FIXED_PRIOR;
			if(vals.size() != 1){
				emsgroot("For parameter '"+ps.name+"' there was a problem with the prior expression '"+ps.prior+"'"); 
			}
		
			if(vals[0] == "*") val1 = TOBESET;
			else get_prior_val(ps.name,vals[0],val1,val1_param);
		}
		else{
			if(toLower(spl[0]) == "uniform"){
				pt = UNIFORM_PRIOR;
				if(vals.size() != 2) emsgroot("For parameter '"+ps.name+"' there was a problem with the prior expression '"+ps.prior+"'"); 
				get_prior_val(ps.name,vals[0],val1,val1_param);	
				get_prior_val(ps.name,vals[1],val2,val2_param);
			}
			else{
				if(toLower(spl[0]) == "exp"){
					pt = EXP_PRIOR;
					if(vals.size() != 1) emsgroot("For parameter '"+ps.name+"' there was a problem with the prior expression '"+ps.prior+"'"); 
					get_prior_val(ps.name,vals[0],val1,val1_param);
				}
				else{
					if(toLower(spl[0]) == "normal"){
						pt = NORMAL_PRIOR;
						if(vals.size() != 2) emsgroot("For parameter '"+ps.name+"' there was a problem with the prior expression '"+ps.prior+"'"); 
						get_prior_val(ps.name,vals[0],val1,val1_param);	
						get_prior_val(ps.name,vals[1],val2,val2_param);
					}
					else{
						if(toLower(spl[0]) == "dir"){
							if(vals.size() == 1){
								if(vals[0] == "*"){ 
									pt = DIRICHLET_FLAT_PRIOR;
									dir_mean = TOBESET;
									dir_sd = TOBESET;
								}
								else{
									pt = DIRICHLET_ALPHA_PRIOR;
									val1 = get_double(vals[0],"For parameter '"+ps.name+"' the prior '"+ps.prior+"'");
									val2 = 1;
									dir_mean = UNSET;
									dir_sd = UNSET;
								}
							}
							else{
								pt = DIRICHLET_PRIOR;
								if(vals.size() != 2) emsgroot("For parameter '"+ps.name+"' the prior '"+ps.prior+"' is not recognised.");
								dir_mean = get_double_with_tobeset(vals[0],"For parameter '"+ps.name+"' the prior '"+ps.prior+"'");
								dir_sd = get_double_with_tobeset(vals[1],"For parameter '"+ps.name+"' the prior '"+ps.prior+"'");
							}
						}
						else{
							if(toLower(spl[0]) == "mdir"){
								if(vals.size() != 1){
									emsgroot("The prior '"+ps.prior+"' should only have one argument.");
								}			
								pt = MDIRICHLET_PRIOR;
								dir_mean = get_double(vals[0],"For parameter '"+ps.name+"' the prior '"+ps.prior+"'");
								dir_sd = TOBESET;
							}
							else{
								emsgroot("For parameter '"+ps.name+"' for the prior '"+spl[0]+"' is not recognised");
							}
						}
					}
				}
			}
		}
		break;
	}
	
	Param par;
	par.name = ps.name; par.priortype = pt;
	par.val1 = val1; par.val2 = val2; par.val1_param = val1_param; par.val2_param = val2_param; 
	par.fixed = val1; 
	
	par.dir_mean = dir_mean; par.dir_sd = dir_sd;
	par.value = value;
	par.type = type;
	par.ps = ps;

	auto th = 0u; while(th < param.size() && param[th].name != par.name) th++;
	
	if(th == param.size()){
		if(par.val1_param != UNSET) param[par.val1_param].param_dep.push_back(param.size());
		if(par.val2_param != UNSET) param[par.val2_param].param_dep.push_back(param.size());
		
		param.push_back(par);
	}
	else{
		if(param[th].priortype != par.priortype || param[th].val1 != par.val1 || param[th].val2 != par.val2 
     || param[th].value != par.value){
			 if(param[th].value != par.value){
				 if(mpi.core == 0){
					cout << "'" << ps.name << "' is defined to have a value of both ";
					cout << param[th].value << " and " <<  par.value << endl;
				 }
			 }
			 else{
				 if(param[th].priortype != par.priortype){
					  if(mpi.core == 0) cout << "'" << ps.name << "' is defined to have two different prior types." << endl;
				 }
			 }
				 
			 emsgroot("Parameter '"+ps.name+"' has multiple definitions.");	
		 }
	}
	
	return th;
}


/// Returns the parameters used to simulate the data
vector <double> Model::simulated_values() const
{
	vector <double> paramv(param.size());
	for(auto th = 0u; th < param.size(); th++){
		paramv[th] = param[th].value;
	}	

	for(const auto &dir : dirichlet){
		if(dir.dirtype == DIR_MODIFIED){
			for(auto par : dir.param) paramv[par.th] *= par.frac;
		}
	}		
	
	if(false){
		for(auto th = 0u; th < param.size(); th++){
			cout << param[th].name << " " << paramv[th] << " Parameter Value" << endl;
		}
	}
	
	return paramv;
}


/// Creates all the splines
vector < vector <double> > Model::create_disc_spline(const vector<double> &paramv_dir) const
{
	vector < vector <double> > disc_spline; disc_spline.resize(spline.size());
	for(auto sp = 0u; sp < spline.size(); sp++) disc_spline[sp] = create_disc_spline(sp,paramv_dir);
	
	return disc_spline;
}


/// Creates a discrete spline for a given reference ref
vector <double> Model::create_disc_spline(const unsigned int ref, const vector<double> &paramv_dir) const
{
	vector <double> disc_spline(details.ndivision);
	
	auto factor = paramv_dir[spline[ref].param_factor];
	auto spl = spline[ref].p;
	
	auto p = 0u;
	for(auto s = 0u; s < details.ndivision; s++){	
		auto t = double((s+0.5)*details.period)/details.ndivision;
		
		while(p+1 < spl.size() && t > spl[p+1].t) p++;
		
		if(p+1 == spl.size()){ 
			if(details.siminf != SIMULATE) emsgEC("Model",9);
			if(s == 0) emsgEC("Model",10);
			disc_spline[s] = disc_spline[s-1];
		}
		else{
			auto fac = (t-spl[p].t)/(spl[p+1].t-spl[p].t);
			disc_spline[s] = factor*(paramv_dir[spl[p].param]*spl[p].multfac*(1-fac) + paramv_dir[spl[p+1].param]*spl[p+1].multfac*fac);
		}
	}
	
	if(details.mode == PREDICTION){
		for(auto s = 0u; s < details.ndivision; s++){	
			disc_spline[s] *= modelmod.spline_mult[ref][s];
			auto set = modelmod.spline_set[ref][s];
			if(set != UNSET) disc_spline[s] = set;
		}
	}			
	
	return disc_spline;
}


/// Sets gradient of discrete spline w.r.t. a variable
void Model::set_disc_spline_grad()
{
	if(details.siminf == SIMULATE) return;
	
	disc_spline_grad.resize(spline.size());                       // The gradient in a spline w.r.t. a variable
	
	for(auto sp = 0u; sp < spline.size(); sp++){
		auto spl = spline[sp].p;
		
		disc_spline_grad[sp].resize(param.size());
		for(auto th = 0u; th < param.size(); th++){
			disc_spline_grad[sp][th].resize(details.ndivision);
			for(auto s = 0u; s < details.ndivision; s++) disc_spline_grad[sp][th][s] = 0;
		}		
	
		auto p = 0;
		for(auto s = 0u; s < details.ndivision; s++){	
			auto t = double((s+0.5)*details.period)/details.ndivision;
			
			while(p < int(spl.size())-1 && t > spl[p+1].t) p++;
			
			auto fac = (t-spl[p].t)/(spl[p+1].t-spl[p].t);
			
			disc_spline_grad[sp][spl[p].param][s] += (1-fac)*spl[p].multfac;
			disc_spline_grad[sp][spl[p+1].param][s] += fac*spl[p].multfac;
		}
	}
	
	for(auto sp = 0u; sp < spline.size(); sp++){ 
		auto spl = spline[sp];
		if(spl.smoothtype != NOSMOOTH){			
			auto spli = spl.p;
			for(auto i = 0u; i < spli.size()-1; i++){
				if(spli[i].t != spli[i+1].t){
					auto par1 = spli[i].param;
					auto par2 = spli[i+1].param;
					
					if(par1 != par2){
						ParamSplineRef psref;	psref.spline = sp; psref.i = i; 
						param[par1].paramsplineref.push_back(psref);
						param[par2].paramsplineref.push_back(psref);
					}
				}
			}
		}
	}
	
	if(false){
		for(auto par: param){
			cout << par.name << " ";
			for(auto psref : par.paramsplineref) cout << spline[psref.spline].name << " " << psref.i << ", ";
			cout << endl;
		}
		emsg("Done");
	}
}


/// In the case of a popultion model, this calculates the transition rates
vector < vector <double> > Model::create_transrate(const vector<double> &paramv_dir) const
{
	vector< vector <double> > transrate;
	transrate.resize(trans.size());
	for(auto c = 0u; c < comp.size(); c++){
		auto shape = comp[c].shape;
		
		auto kmax = comp[c].trans.size();
		for(auto k = 0u; k < kmax; k++){
			auto tr = comp[c].trans[k];
		
			transrate[tr].resize(data.ndemocatpos);
			if(trans[tr].inf == TRANS_NOTINFECTION){	
				for(auto dp = 0u; dp < data.ndemocatpos; dp++){
					auto mean = paramv_dir[comp[c].param_mean[dp]]/shape; if(mean == 0) emsgEC("Model",11);
					
					transrate[tr][dp] = 1.0/mean;
					if(kmax > 1) transrate[tr][dp] *= paramv_dir[trans[tr].param_prob[dp]];
				}
			}
		}
	}
	
	return transrate;
}


/// Defines the relative susceptibility of individuals
vector <double> Model::create_susceptibility(const vector<double> &paramv_dir) const 
{
	vector <double> susceptibility(data.ndemocatpos);

	vector < vector <double> > sus;
	sus.resize(data.ndemocat);
	for(auto c = 0u; c < data.ndemocat; c++){                        // Checks susceptibility is correctly averaged
		auto nval = data.democat[c].value.size(); 
		sus[c].resize(nval);
		if(data.democat[c].sus_vari == true){		
			for(auto j = 0u; j < nval; j++) sus[c][j] = paramv_dir[data.democat[c].sus_param[j]];
		}
		else{
			for(auto j = 0u; j < nval; j++) sus[c][j] = 1;
		}
	}
	
	for(auto dp = 0u; dp < data.ndemocatpos; dp++){
		auto val = 1.0; for(auto c = 0u; c < data.ndemocat; c++) val *= sus[c][data.democatpos[dp][c]];
		susceptibility[dp] = val;
	}
	
	if(false){
		for(auto c = 0u; c < data.ndemocat; c++){                       // Checks the susceptibility is correctly averaged
			auto sum = 0.0;
		
			auto nval = data.democat[c].value.size(); 
	
			for(auto j = 0u; j < nval; j++) sum += sus[c][j]*data.democat_dist[c][j];
			cout << sum << "sum check" << endl;
		}
	}
	
	return susceptibility;
}


/// Defines the relative transmission rate for different areas
vector < vector <double> > Model::create_areafactor(const vector<double> &paramv_dir) const
{
	vector <double> af(data.narea);
	for(auto c = 0u; c < data.narea; c++){
		if(region_effect.on == true) af[c] = exp(paramv_dir[region_effect.area_param[c]]);
		else af[c] = 1;
		
		if(data.area_effect.on == true) af[c] *= paramv_dir[area_effect_param_list[c]];
			
		if(std::isnan(af[c])) emsgEC("Model",12);
	}
	
	vector < vector <double> > areafactor;
	areafactor.resize(details.period);
	for(auto t = 0u; t < details.period; t++) areafactor[t] = af;
	
	for(auto j = 0u; j < data.ncovar; j++){
		for(auto t = 0u; t < details.period; t++){
			for(auto c = 0u; c < data.narea; c++){
				areafactor[t][c] *= exp(paramv_dir[covariate_param[j]]*data.covar[j].value[c][t]);
			}
		}
	}
	
	if(data.level_effect.on == true){
		auto nlev_ef = level_effect_param_list.size();
		vector <double> lev_ef(nlev_ef);
		for(auto i = 0u; i < nlev_ef; i++) lev_ef[i] = paramv_dir[level_effect_param_list[i]];
		
		for(auto t = 0u; t < details.period; t++){
			for(auto c = 0u; c < data.narea; c++){
				areafactor[t][c] *= lev_ef[data.level_effect.param_map[t][c]];
			}				
		}
	}
	
	auto area_av = calculate_area_av(areafactor);
	for(auto t = 0u; t < details.period; t++){            // Normalises the result
		for(auto c = 0u; c < data.narea; c++){
			areafactor[t][c] /= area_av;
		}
	}
	
	if(false){
		for(auto t = 0u; t < details.period; t++){
			cout << t << " ";
			for(auto c = 0u; c < data.narea; c++) cout << areafactor[t][c] << " ";
			cout << "areafactor" << endl;
		}
		emsg("Done");
	}
	
	return areafactor;
}


/// Calcualtes the average of the area factor
double Model::calculate_area_av(const vector < vector <double> > &areafactor) const
{
	auto area_av = 0.0; 
	for(auto c = 0u; c < data.narea; c++){
		auto av = 0.0; for(auto t = 0u; t < details.period; t++) av += areafactor[t][c];
		av /= details.period;		
		area_av += av*data.area[c].total_pop;
	}
	area_av /= data.popsize;
	
	return area_av;
}


/// Defines the time variation in the age mixing matrix
vector < vector < vector <double> > > Model::create_Ntime(const vector < vector<double> > &disc_spline) const
{	
	timer[TIME_CREATEN].start();
	
	vector < vector < vector <double> > > Ntime;

	Ntime.resize(details.ndivision);
	for(auto sett = 0u; sett < details.ndivision; sett++){
		auto &N = Ntime[sett];
		N = data.genQ.N[0].ele;
		
		for(auto &mm : data.genQ.matmod){                             // These modify the basic contact matrix
			switch(mm.type){
				case ROW_COLUMN: case PERTURB:
					auto si = N.size();
					vector <bool> flag(si);
					//for(auto i = 0u; i < si; i++) flag[i] = false;
					
					//for(auto a : mm.ages) flag[a] = true;
				
					auto fac = disc_spline[mm.spline_ref][sett];
				
					for(auto a : mm.ages){
						for(auto i = 0u; i < si; i++){
							N[a][i] *= fac;
							N[i][a] *= fac; 
							//if(flag[i] == false) N[i][a] *= fac; 
						}
					}
					break;
			}
		}
		
		if(false){
			for(auto j = 0u; j < N.size(); j++){
				for(auto i = 0u; i < N.size(); i++){
					cout << N[j][i] << " ";
					
				}
				cout << "   " << sett << " Mixing matrix" << endl;
			}				
		}
	}	
	
	timer[TIME_CREATEN].stop();
		
	return Ntime;
}


/// Determines if two matrices are equal or not
bool Model::equal(const vector < vector <double> > &Ntime_1, const vector < vector <double> > &Ntime_2) const  
{
	for(auto i = 0u; i < Ntime_1.size(); i++){
		for(auto j = 0u; j < Ntime_1[i].size(); j++){
			auto val1 = Ntime_1[i][j], val2 = Ntime_2[i][j];
			if(val1 < val2-TINY || val1 > val2+TINY) return false;
		}
	}
	return true;
}
	

/// Calculates the transmission rate beta from R
vector < vector < vector <double> > > Model::calculate_beta_from_R(const vector <double> &susceptibility, const vector <double> &paramv_dir, const vector < vector < vector <double> > > &Ntime, const vector < vector <double> > &transrate,	const vector < vector <double> > &disc_spline) const 
{
	timer[TIME_BETA_FROM_R].start();
		
	vector < vector < vector <double> > > beta;

	beta.resize(data.nstrain);
	for(auto st = 0u; st < data.nstrain; st++){
		beta[st].resize(data.narea);
		
		auto Vinv = calculate_Vinv(transrate,st);

		for(const auto &info : Rspline_info){			
			auto spl = info.spline_ref;
			auto c = info.area[0];
			beta[st][c].resize(details.ndivision);
			
			auto ratio = 0.0;
			for(auto sett = 0u; sett < details.ndivision; sett++){
				if(sett == 0 || equal(Ntime[sett-1],Ntime[sett]) == false){
					ratio = calculate_R_beta_ratio_using_NGM(paramv_dir,susceptibility,Ntime[sett],Vinv,st,info.democatpos_dist);
				}		
				beta[st][c][sett] = paramv_dir[data.strain[st].Rfactor_param]*disc_spline[spl][sett]/ratio;
			}

			for(auto i = 1u; i < info.area.size(); i++){                   // Copies to other areas
				beta[st][info.area[i]] = beta[st][c];
			}
		}
	}
	
	if(details.mode == PREDICTION){                                    // Potential model changes
		for(auto st = 0u; st < data.nstrain; st++){
			for(auto c = 0u; c < data.narea; c++){
				for(auto sett = 0u; sett < details.ndivision; sett++){	
					beta[st][c][sett] *= modelmod.beta_mult[sett][c][st];              
				}
			}
		}
	}
	
	timer[TIME_BETA_FROM_R].stop();
		
	return beta;
}

 
/// From a string gets the transition number
vector <unsigned int> Model::get_transition(const string transname)
{
	vector <unsigned int> tr_list;
	for(auto tr = 0u; tr < trans.size(); tr++){
		if(transname == comp[trans[tr].from].name + "->" + comp[trans[tr].to].name){
			tr_list.push_back(tr);
		}
	}
	if(tr_list.size() == 0)	emsgroot("The transition '"+transname+"' is not recognised");
	
	return tr_list;
}


/// Completes information for datatables based on the model
void Model::complete_datatables()
{
	for(auto d = 0u; d < data.datatable.size(); d++){
		auto &dt = data.datatable[d];
		
		auto error_desc = "for file '"+dt.file+"'"; if(dt.optype == OUTPUT) error_desc = "in 'state_outputs'";
		
		auto obs = dt.observation;
		auto sp = split(obs,',');
		
		for(auto &st : sp){
			auto fl = false;
			for(auto co = 0u; co < comp.size(); co++){
				if(st == comp[co].name){
					if(find_in(dt.complist,co) != UNSET){
						emsg("The observation '"+obs+"' "+error_desc+" should not contain multiple instances of the same compartment.");
					}
						
					dt.complist.push_back(co);
					if(dt.type != POP && dt.type != POPFRAC){
						emsg("The observation '"+obs+"' "+error_desc+" should not contain compartmental population data");
					}
					fl = true;
				}
			}
			
			if(fl == false){
				for(auto tr = 0u; tr < trans.size(); tr++){
					if(st == comp[trans[tr].from].name + "->" + comp[trans[tr].to].name){
						if(find_in(dt.translist,tr) != UNSET){
							emsg("The observation '"+obs+"' "+error_desc+" should not contain multiple instances of the same transition.");
						}
						
						dt.translist.push_back(tr);
						
						if(dt.type != TRANS && dt.type != MARGINAL){
							emsg("The observation '"+obs+"' "+error_desc+" should not contain compartmental transition data in 'observation'");
						}		
						fl = true;
					}
				}
			}
			
			if(fl == false){
				emsgroot("The observation '"+st+"' "+error_desc+" is not a compartment or a transition");
			}
		}
		
		if(dt.factor_spline != ""){
			auto sp = 0u; while(sp < spline.size() && spline[sp].name != dt.factor_spline) sp++;
			if(sp == spline.size()) emsgroot("In 'data_tables' cannot find '"+dt.factor_spline+"'");
			
			for(auto &ob : data.obs){
				if(ob.datatable == d) ob.factor_spline = sp;
			}
			
			for(auto &gr : data.graph){
				if(gr.datatable == d) gr.factor_spline = sp;
			}
		}
	}
}


/// Implements model modifications
void Model::setup_modification()
{
	modelmod.transmean_mult.resize(details.ndivision);
	for(auto sett = 0u; sett < details.ndivision; sett++){					
		modelmod.transmean_mult[sett].resize(data.narea);
		for(auto c = 0u; c < data.narea; c++){
			modelmod.transmean_mult[sett][c].resize(trans.size());
			for(auto tr = 0u; tr < trans.size(); tr++){
				modelmod.transmean_mult[sett][c][tr].resize(data.ndemocatpos);
				for(auto dp = 0u; dp < data.ndemocatpos; dp++){
					modelmod.transmean_mult[sett][c][tr][dp] = 1;
				}
			}
		}
	}
	
	modelmod.beta_mult.resize(details.ndivision);
	for(auto sett = 0u; sett < details.ndivision; sett++){					
		modelmod.beta_mult[sett].resize(data.narea);
		for(auto c = 0u; c < data.narea; c++){
			modelmod.beta_mult[sett][c].resize(data.nstrain);
			for(auto i = 0u; i < data.nstrain; i++){
				modelmod.beta_mult[sett][c][i] = 1;
			}
		}
	}
	
	modelmod.efoi_mult.resize(details.ndivision);
	for(auto sett = 0u; sett < details.ndivision; sett++){					
		modelmod.efoi_mult[sett].resize(data.narea);
		for(auto c = 0u; c < data.narea; c++){
			modelmod.efoi_mult[sett][c].resize(data.nstrain);
			for(auto i = 0u; i < data.nstrain; i++){
				modelmod.efoi_mult[sett][c][i] = 1;
			}
		}
	}

	modelmod.spline_mult.resize(spline.size()); modelmod.spline_set.resize(spline.size());
	for(auto sp = 0u; sp < spline.size(); sp++){
		modelmod.spline_mult[sp].resize(details.ndivision); modelmod.spline_set[sp].resize(details.ndivision);
		for(auto sett = 0u; sett < details.ndivision; sett++){
			modelmod.spline_mult[sp][sett] = 1;
			modelmod.spline_set[sp][sett] = UNSET;
		}
	}
	
	modelmod.pred_start = UNSET;
	for(const auto &cf : data.modification){
		if(modelmod.pred_start == UNSET || cf.start < modelmod.pred_start) modelmod.pred_start = cf.start;
		
		auto start = cf.start*details.division_per_time;
		auto end = cf.end*details.division_per_time;
		auto fac = cf.factor;
		
		switch(cf.type){
			case CF_TRANS_RATE:
				{
					auto tr_list = get_transition(cf.trans_str);
				
					for(auto tr : tr_list){
						for(auto sett = start; sett < end; sett++){					
							for(auto c : cf.area){
								for(auto dp : cf.dp_sel){
									modelmod.transmean_mult[sett][c][tr][dp] *= fac;
								}
							}
						}
					}
				}
				break;
						
			case CF_BETA:
				{
					auto st=0u;
					if(cf.strain_str == ""){
						if(data.nstrain != 1) emsgroot("In 'modification' a value for 'strain' must be set");
					}
					else{
						st = 0; while(st < data.nstrain && cf.strain_str != data.strain[st].name) st++;
						if(st == data.nstrain) emsgroot("In 'modification' the value 'strain=\""+cf.strain_str+"\"' is not recognised");
					}
					
					for(auto sett = start; sett < end; sett++){					
						for(auto c : cf.area){
							modelmod.beta_mult[sett][c][st] *= fac;
						}
					}
				}
				break;
				
			case CF_EFOI:
				{
					auto st=0u;
					if(cf.strain_str == ""){
						if(data.nstrain != 1) emsgroot("In 'modification' a value for 'strain' must be set");
					}
					else{
						st = 0; while(st < data.nstrain && cf.strain_str != data.strain[st].name) st++;
						if(st == data.nstrain) emsgroot("In 'modification' the value 'strain=\""+cf.strain_str+"\"' is not recognised");
					}
					
					for(auto sett = start; sett < end; sett++){					
						for(auto c : cf.area){
							modelmod.efoi_mult[sett][c][st] *= fac;
						}
					}
				}
				break;
				
			case CF_SPLINEFAC: case CF_SPLINESET:
				auto sp = 0u; while(sp < spline.size() && spline[sp].name != cf.spline_name_str) sp++;
				if(sp == spline.size()) emsgroot("In 'modification' cannot find spline with the name '"+cf.spline_name_str+"'");
				
				for(auto sett = start; sett < end; sett++){	
					if(cf.type == CF_SPLINEFAC){ 
						if(modelmod.spline_set[sp][sett] != UNSET){
							emsgroot("In 'modification' a spline cannot both be directly set and have a factor change");
						}
						modelmod.spline_mult[sp][sett] *= fac;
					}
					else{
						if(modelmod.spline_mult[sp][sett] != 1){
							emsgroot("In 'modification' a spline cannot both be directly set and have a factor change");
						}
						modelmod.spline_set[sp][sett] = fac;
					}
				}		
				break;
		}
	}
	
	if(details.pred_start == UNSET){
		if(modelmod.pred_start == UNSET) modelmod.pred_start = details.end - details.start;
	}
	else{
		if(modelmod.pred_start < details.pred_start) emsgroot("Dates in 'modification' cannot be before 'pred_start'"); 
		modelmod.pred_start = details.pred_start;
	}
}


/// Calculates the prior probability
double Model::prior(const vector<double> &paramv) const 
{
	double Pr = 0;
	
	for(auto th = 0u; th < param.size(); th++){
		switch(param[th].priortype){
			case FIXED_PRIOR: case UNIFORM_PRIOR:
				break;
			
			case EXP_PRIOR:
				{
					auto val = paramv[th];
					auto rate = get_val1(th,paramv);
					Pr += log(rate) - rate*val;
				}
				break;
				
			case NORMAL_PRIOR:
				{
					auto val = paramv[th];
					auto mu = get_val1(th,paramv);
					auto sd = get_val2(th,paramv);
					Pr += normal_probability(val,mu,sd*sd);
				}
				break;
				
			case DIRICHLET_PRIOR: case DIRICHLET_ALPHA_PRIOR: case MDIRICHLET_PRIOR:
				{
					auto val = paramv[th];
					auto alpha = get_val1(th,paramv);
					auto beta = get_val2(th,paramv);
					if(details.mode == ML_INF){
						if(beta != 1) emsgEC("Model",54);
						if(alpha != 1){
							Pr += (alpha-1)*log(val);
							//cout << param[th].name << " " << alpha << " "<< beta << " " << val << " hh\n";
						}
					}
					else{
						Pr += gamma_probability(val,alpha,beta);
					}
				}
				break;
				
			case DIRICHLET_FLAT_PRIOR:
				{
					if(details.mode != ML_INF){
						auto val = paramv[th];
						Pr += -val;
					}
				}
				break;
				
			default:
				emsgEC("model",13);
				break;
		}
	}
	
	Pr += spline_prior(paramv);

	return Pr;
}
	

/// Calculates the smoothing spline prior
double Model::spline_prior(const vector<double> &paramv) const
{
	auto Pr = 0.0;
	if(details.siminf == SIMULATE) return Pr;
	
	for(const auto &spl : spline){
		if(spl.smoothtype != NOSMOOTH){ 
			auto spli = spl.p;
			for(auto i = 0u; i < spli.size()-1; i++){
				auto par1 = spli[i].param;
				auto par2 = spli[i+1].param;
				
				if(par1 != par2 && spli[i].t != spli[i+1].t){
					if(param[par1].priortype != FIXED_PRIOR || param[par2].priortype != FIXED_PRIOR){
						auto value = spli[i].smooth_value;
						if(value == UNSET) emsgroot("For the spline '"+spl.name+"' a value for 'smooth' must be set");
							
						switch(spl.smoothtype){
							case LOGSMOOTH:                                    // This is the pdf for the log-normal distribution
								{
									double mu = log(paramv[par1]*spli[i].multfac);
									double logx = log(paramv[par2]*spli[i+1].multfac);							
									Pr += -(logx-mu)*(logx-mu)/(2*value*value) - logx; 
								}
								break;
								
							case SMOOTH:
								{
									double d = paramv[par1] - paramv[par2];
									Pr += -d*d/(2*value*value); 
								}
								break;
							default: emsgEC("Model",14);
						}
					}
				}
			}
		}
	}
	
	return Pr;
}


/// Gets val1 from the parameter
double Model::get_val1(const unsigned int th, const vector <double> &paramv) const 
{
	auto val1 = param[th].val1; if(val1 == UNSET) val1 = paramv[param[th].val1_param]; if(val1 == UNSET) emsgEC("Model",15);
	return val1;
}


/// Gets val2 from the parameter
double Model::get_val2(const unsigned int th, const vector <double> &paramv) const 
{
	auto val2 = param[th].val2; 
	if(val2 == UNSET) val2 = paramv[param[th].val2_param]; 
	if(val2 == UNSET) emsgEC("Model",16);
	return val2;
}


/// Returns the gradient in the prior
double Model::dPr_dth(const unsigned int th, const vector<double> &paramv) const 
{
	auto dPr_dth = 0.0;
	
	auto val = paramv[th];

	switch(param[th].priortype){                                      // The direct contribution due to the prior
		case FIXED_PRIOR: case UNIFORM_PRIOR: break;
		
		case EXP_PRIOR:
			{
				auto rate = get_val1(th,paramv);
				dPr_dth = -rate;
			}
			break;
			
		case NORMAL_PRIOR:
			{
				auto mu = get_val1(th,paramv), sd = get_val2(th,paramv);
				dPr_dth += -(val-mu)/(sd*sd);
			}
			break;
			
		case DIRICHLET_PRIOR: case DIRICHLET_ALPHA_PRIOR: case MDIRICHLET_PRIOR:
			{
				auto a = get_val1(th,paramv), b = get_val2(th,paramv);
				dPr_dth += (a-1)/val - b;
			}
			break;
			
		case DIRICHLET_FLAT_PRIOR:
			{
			}
			break;
			
		default:
			emsgEC("model",17);
			break;
	}
	
	for(auto th_dep : param[th].param_dep){                    // Includes other parameters which depend on this one
		switch(param[th_dep].priortype){
			case NORMAL_PRIOR:
				if(param[th_dep].val2_param == th){
					dPr_dth += -1.0/val + paramv[th_dep]*paramv[th_dep]/(val*val*val);
				}
				else emsgEC("Model",18);
				break;
			
			default: emsgEC("Model",19); break;
		}
	}
	
	for(auto psref : param[th].paramsplineref){                // Includes any splines which depend on this one
		auto i = psref.i;		 
		auto spl = spline[psref.spline];	
		auto spli = spl.p;
			
		if(spli[i].t != spli[i+1].t){
			auto par1 = spli[i].param;
			auto par2 = spli[i+1].param;
		
			if((par2 == th || par1 == th) && par2 != par1){
				auto value = spli[i].smooth_value; 
				if(value == UNSET) emsgroot("For the spline '"+spl.name+"' a value for 'smooth' must be set");

				switch(spl.smoothtype){
					case LOGSMOOTH:
						{
							double dval = log(paramv[par2]*spli[i+1].multfac/(paramv[par1]*spli[i].multfac))/value;
						
							if(par2 == th){
								auto dval_dth = 1.0/(value*paramv[th]);
								dPr_dth += -dval*dval_dth - 1.0/paramv[th];
							}
							
							if(par1 == th){
								auto dval_dth = -1.0/(value*paramv[th]);
								dPr_dth += -dval*dval_dth;
							}
						}
						break;
						
					case SMOOTH:
						{
							double d = paramv[par1] - paramv[par2];
							if(par2 == th) dPr_dth += d/(value*value);
							if(par1 == th) dPr_dth -= d/(value*value);
						}
						break;
					
					default: emsgEC("Model",20); break;
				}
			}
		}
	}
		
	if(false){
		auto par = paramv;
		auto val = prior(par);
		auto d = 0.00001;
		par[th] += d;
		auto up = prior(par);
		par[th] -= d;
		cout << param[th].name << " " << (up-val)/d << " " << dPr_dth << " grad" << endl;
	}
	
	return dPr_dth;
}


/// Based on transition probabilites, this calculates the probability of an passing through a given compartment
vector <CompProb> Model::create_compprob(const vector<double> &paramv_dir, const unsigned int st) const
{
	vector <CompProb> compprob(comp.size());
	
	auto dpmax = data.ndemocatpos_per_strain;

	for(auto &co : compprob){
		co.value.resize(dpmax);
		for(auto &value : co.value) value = 0;
	}

	for(auto dp = 0u; dp < dpmax; dp++){
		vector <unsigned int> cst, kst;
		vector <double> probst;
	
		auto prob = 1.0;
		auto c = trans[infection_trans].to;  
		do{
			compprob[c].value[dp] += prob;
			
			unsigned int k;
			if(comp[c].trans.size() == 0){
				if(cst.size() == 0) break;
				
				do{
					auto j = cst.size()-1;
					c = cst[j];
					prob = probst[j]; 
					kst[j]++;
					k = kst[j];
					if(k < comp[c].trans.size()) break;
					
					cst.pop_back();
					probst.pop_back();
					kst.pop_back();
				}while(cst.size() > 0);
				if(k == comp[c].trans.size()) break;
			}
			else{
				k = 0;
				if(comp[c].trans.size() > 1){
					cst.push_back(c);
					probst.push_back(prob);
					kst.push_back(0);
				}
			}
			if(comp[c].trans.size() > 1){
				prob *= paramv_dir[trans[comp[c].trans[k]].param_prob[st*dpmax+dp]];
			}

			c = trans[comp[c].trans[k]].to;		
		}while(1 == 1);
		
		if(cst.size() > 0 || kst.size() > 0 || probst.size() > 0) emsgEC("Model",21);
	}
	
	return compprob;
}


/// Calculates the number of infections caused externally per strain
vector <double> Model::calculate_external_ninf(const vector<double> &paramv_dir) const
{	
	vector < vector<double> > spline_val;
	spline_val.resize(spline.size());
	for(auto sp = 0u; sp < spline.size(); sp++) spline_val[sp] = create_disc_spline(sp,paramv_dir);
	
	auto susceptibility = create_susceptibility(paramv_dir);   
	auto dt = double(details.period)/details.ndivision;
	
	auto dpmax = data.ndemocatpos_per_strain;
	
	vector <double> ninf(data.nstrain);
	for(auto st = 0u; st < data.nstrain; st++){
		ninf[st] = 0;
		for(auto sett = 0u; sett < details.ndivision; sett++){
			for(auto c = 0u; c < data.narea; c++){
				auto sp_info = efoispline_info[efoi_spl_ref[st][c]];
				auto val = spline_val[sp_info.spline_ref][sett];
		
				for(auto dp = 0u; dp < dpmax; dp++){
					auto popsum = 0u; for(auto co = 0u; co < comp.size(); co++) popsum += data.area[c].pop_init[co][dp];
					auto a = data.democatpos[dp][0];
					ninf[st] += susceptibility[st*dpmax + dp]*val*popsum*sp_info.efoi_agedist[a]*dt;
				}
			}
		}
	}
	
	return ninf;
}


/// Calculates the probability of reaching a specified compartment
vector <double> Model::calculate_probreach(const vector<double> &paramv_dir, const unsigned int st) const
{
	vector <double> vec;
	for(auto pr = 0u; pr < prob_reach.size(); pr++){
		auto probr = prob_reach[pr];
		auto compprob = create_compprob(paramv_dir,st);
		
		auto sum = 0.0, wsum = 0.0;
		for(auto dp = 0u; dp < data.ndemocatpos_per_strain; dp++){
			auto w = data.democatpos_dist[dp];
			sum += compprob[probr.comp].value[dp]*w;
			wsum += w;
		}
		sum /= wsum;
		
		vec.push_back(sum);
	}
	
	return vec; 
}

	
/// Returns a sample of all the splines that will be output
vector <SplineOutput> Model::get_spline_output(const vector <double> &paramv_dir, const vector < vector < vector < vector <double> > > > &pop) const
{
	vector <SplineOutput> splineout;
	
	vector <bool> flag(spline.size()); for(auto sp = 0u; sp < spline.size(); sp++) flag[sp] = false;
	
	for(auto st = 0u; st < data.nstrain; st++){                       // Adds splines for R
		auto R_eff = calculate_R_eff(paramv_dir,pop,st);
		
		string name_add; if(data.nstrain > 1) name_add = " for "+data.strain[st].name;
		
		for(auto i = 0u; i < Rspline_info.size(); i++){
			auto sp = Rspline_info[i].spline_ref; flag[sp] = true;
			
			if(pop.size() > 0){                                    // This outputs effective Rt
				SplineOutput so;
				so.name = "Effective "+spline[sp].name+name_add;
				so.desc = "The effective reproduction number"+name_add;
				so.fulldesc = "Effective reproduction number &R^{eff}&: This graph shows the temporal variation in the effective reproduction number &R^{eff}&.  &R^{eff}& is defined to be the expected number of cases directly caused by an infected individual. This takes into account the fact that as the epidemic progresses the fraction of susceptible individuals reduces (causing effective transmission upon contact of individuals to become less and less common). &R^{eff}& is always less than &R& and if it reduces below 1 herd immunity in reached (i.e. the disease naturally dies out over time).";
			
				so.tab = details.analysis_type; 
				so.tab2 = "Transmission";
				so.tab3 = "&R^{eff}&"; 
				if(data.nstrain > 1) so.tab4 += data.strain[st].name+" Time Variation"; 
				else{
					if(data.narea > 1) so.tab4 = "Time Variation"; else so.tab4 = "";
				}
			
				so.splineval = R_eff[i];
				so.spline_param_JSON = "";
				
				splineout.push_back(so);
			}
			else{                                                       // For simulated parameters there is no output
				SplineOutput so;
				splineout.push_back(so);
			}
			
			SplineOutput so;                                            // This outputs Rt
			so.name = spline[sp].name+name_add;
			so.desc = spline[sp].desc+name_add; 
			so.fulldesc = "Reproduction number &R_t&: This graph shows the temporal variation in the reproduction number &R&.  &R_t& is defined as the expected number of cases directly caused by an infected individual, assuming an otherwise suceptible population (note, this does not account for the fact that in reality the susceptible population reduces with time). As such, &R_t& is interpreted as a quantity proportional to the rate at which individuals come into contact with each other (with each contact allowing for the possibility of disease transmission).";
			
			if(Rspline_info.size() > 1 || data.nstrain > 1){
				so.fulldesc += " Based on the mathematical description given in the model section, this graph is calculated using &";
				if(Rspline_info.size() > 1) so.fulldesc += "R_{a,t}"; else so.fulldesc += "R_{t}";
				if(data.nstrain > 1) so.fulldesc += "Ψ_{"+data.strain[st].name+"}";
				so.fulldesc += "&."; 
			}
			
			so.tab = details.analysis_type;
			so.tab2 = "Transmission";
			so.tab3 = "&R&"; 
			if(data.nstrain > 1) so.tab4 += data.strain[st].name+" Time Variation";
			else{
				if(data.narea > 1) so.tab4 = "Time Variation"; else so.tab4 = "";
			}
				
			so.splineval = create_disc_spline(sp,paramv_dir);
			auto fac = paramv_dir[data.strain[st].Rfactor_param];
			for(auto &val : so.splineval) val *= fac;
			
			so.spline_param_JSON = spline_param_JSON(sp);
			
			splineout.push_back(so);
		}
	}
	
	if(false){   // This has temporarily been turned off, but may be put on again when we look at time-varying N
		if(data.nage > 1){ 
			auto R_age = calculate_R_age(paramv_dir);
			for(auto a = 0u; a < data.nage; a++){
				SplineOutput so;
				so.name = "&R_t& for age:"+data.democat[0].value[a];
				so.desc = "The reproduction number for age: "+data.democat[0].value[a];
				so.splineval = R_age[a];
				splineout.push_back(so);
			}
		}
	}
	
	for(auto i = 0u; i < efoispline_info.size(); i++){                      // Adds splines for efoi
		auto sp = efoispline_info[i].spline_ref; flag[sp] = true;
			
		SplineOutput so;
		so.name = spline[sp].name;
		so.desc = spline[sp].desc;
		so.fulldesc = "External force of infection: The external force of infection sets the average number of infections per unit time (for every "+to_string(details.efoi_factor)+" individuals) caused by contacts from outside the population under investigation.";
		so.tab = details.analysis_type; 
		so.tab2 = "Transmission";
		so.tab3 = "External FOI"; 
		if(efoispline_info.size() > 1) so.tab4 = efoispline_info[i].name;
				
		so.splineval = create_disc_spline(sp,paramv_dir);
		for(auto &val : so.splineval) val *= details.efoi_factor;
		so.desc += " (per "+to_string(details.efoi_factor)+")";
		
		so.spline_param_JSON = spline_param_JSON(sp);	
		
		splineout.push_back(so);		
	}
	
	for(auto sp = 0u; sp < spline.size(); sp++){                            // Any other splines          
		if(flag[sp] == false){
			SplineOutput so;
			so.name = spline[sp].name;
			so.desc = spline[sp].desc;
			if(sp == geo_spline_ref){
				so.fulldesc = "Time variation in spatial mixing: This plot shows the time variation in spatial mixing. A value of &m_{t}&=1 implies that the geographical mixing matrix &M_{a,a′}& is as specifed by the file given in the *geo-mixing-matrix* command. A value of &m_{t}&=0, however, implies no contacts between areas.";
				so.tab = details.analysis_type; so.tab2 = "Spatial Mixing"; so.tab3 = "Time series"; so.tab4 = "&m_{t}&";
			}
			else{
				auto mm = 0u; while(mm < data.genQ.matmod.size() && data.genQ.matmod[mm].spline_ref != sp) mm++;
				if(mm < data.genQ.matmod.size()){
					const auto &mmod = data.genQ.matmod[mm];
					
					string st;
					
					so.fulldesc = "Age mixing modification: ";
					switch(mmod.type){
						case ROW_COLUMN: case PERTURB:
							so.fulldesc += "This spline represents a time-varying factor that multiplies the ";
							if(mmod.ages.size() == 1){
								st = "*"+data.democat[0].value[mmod.ages[0]]+"*";
								so.fulldesc += st+" row and column";
							}
							else{
								for(auto k = 0u; k < mmod.ages.size(); k++){
									if(k > 0){
										if(k < mmod.ages.size()-1) st += ", ";
										else st += " and ";
									}
									st += "*"+data.democat[0].value[mmod.ages[k]]+"*";
								}
								so.fulldesc += st+" rows and columns";
							}
					
							so.fulldesc += " of the age mixing matrix.";
							break;
					}
					so.tab = details.analysis_type; so.tab2 = "Age Mixing"; so.tab3 = "Time series"; so.tab4 = "Modification "+to_string(mm+1); 
				}
				else{
					so.fulldesc = "Spline: "+spline[sp].desc;
					so.tab = details.analysis_type; so.tab2 = spline[sp].name; 
				}
			}
	
			so.splineval = create_disc_spline(sp,paramv_dir);
			
			so.spline_param_JSON = spline_param_JSON(sp);
		
			splineout.push_back(so);
		}
	}
	
	return splineout;
}	
	

/// Returns a sample of derived parameters
vector <DerivedParam> Model::get_derived_param(const vector<double> &paramv_dir, const vector <double> &susceptibility, const vector <vector <double> > &A, const vector < vector <double> > &transrate) const
{
	vector <DerivedParam> derpar;
	
	for(auto st = 0u; st < data.nstrain; st++){
		DerivedParam dp;
		dp.name = "Generation time";
		dp.desc = "Generation time";
		dp.file = "generation_time";
		if(data.nstrain > 1){
			auto name = data.strain[st].name;
			dp.name += " "+name;
			dp.desc += " for "+name;
			dp.file += "_"+name;
		}
		dp.file += ".csv"; 
			
		dp.value = calculate_generation_time(paramv_dir,susceptibility,A,transrate,data.democatpos_dist,st);
		derpar.push_back(dp);
	}
	
	for(auto st = 0u; st < data.nstrain; st++){
		auto prvec = calculate_probreach(paramv_dir,st);                // The probability of reaching a given compartment
		for(auto pr = 0u; pr < prob_reach.size(); pr++){
			auto &probr = prob_reach[pr];
			DerivedParam dp;
			dp.name = probr.name;
			dp.file = probr.name;
			dp.desc = "The probability of reaching "+comp[probr.comp].name;
			if(data.nstrain > 0){
				auto name = data.strain[st].name;
				dp.name += " for strain "+name;
				dp.desc += " for strain "+name;
				dp.file += " for strain "+name;
			}
			dp.value = prvec[pr];
			derpar.push_back(dp);
		}
	}
	
	auto exf_ninf = calculate_external_ninf(paramv_dir);
	for(auto st = 0u; st < data.nstrain; st++){                       // The generation time
		DerivedParam dp;
		dp.name = "External Infections";
		dp.desc = "Number of external infections";
		dp.file = dp.name+".csv";
		if(data.nstrain > 1){
			auto name = data.strain[st].name;
			dp.name += " "+name;
			dp.desc += " for "+name;
		}
		dp.value = exf_ninf[st];
		derpar.push_back(dp);
	}
	
	return derpar;
}


/// Returns the distribution in the susceptible individuals for a set of areas
vector <double> Model::get_sus_dist(const unsigned int sett, const vector <unsigned int> &area, const vector < vector < vector < vector <double> > > > &pop) const
{
	vector <double> dist(data.ndemocatpos_per_strain);
	for(auto dp = 0u; dp < data.ndemocatpos_per_strain; dp++) dist[dp] = 0;
	
	auto total_pop = 0.0;
	
	for(auto c : area){
		for(auto dp = 0u; dp < data.ndemocatpos; dp++){
			auto p = pop[sett][c][start_compartment][dp];
			if(p < 1) p = 1;
			dist[dp%data.ndemocatpos_per_strain] += p;
			
			for(auto co = 0u; co < comp.size(); co++) total_pop += data.area[c].pop_init[co][dp];
		}
	}
	for(auto dp = 0u; dp < data.ndemocatpos_per_strain; dp++) dist[dp] /= total_pop;
	
	return dist;
}


/// Gets the initial susceptible distribution
vector <double> Model::get_sus_dist_init(const vector <unsigned int> &area) const
{
	vector <double> dist(data.ndemocatpos_per_strain);
	auto total_pop = 0.0;
	for(auto c : area){
		for(auto co = 0u; co < comp.size(); co++){
			for(auto dp = 0u; dp < data.ndemocatpos; dp++){
				auto pop = data.area[c].pop_init[co][dp];
				dist[dp%data.ndemocatpos_per_strain] += pop;
				total_pop += pop;
			}
		}
	}
	
	for(auto dp = 0u; dp < data.ndemocatpos_per_strain; dp++) dist[dp] /= total_pop;
	
	return dist;
}
			

/// Stores maps for the reproduction number
vector <RMap> Model::get_Rmap(const vector<double> &paramv_dir, const vector < vector <double> > &disc_spline, const vector < vector <double> > &areafactor, const vector <double> &susceptibility, const vector < vector < vector <double> > > &Ntime, const vector < vector <double> > &transrate, const vector < vector < vector < vector <double> > > > &pop) const
{
	vector <RMap> rmap_list;
	
	for(auto st = 0u; st < data.nstrain; st++){
		RMap rmap;
			
		auto &map = rmap.map;
		
		map.resize(details.period);
		for(auto t = 0u; t < details.period; t++) map[t].resize(data.narea);
		
		for(const auto &info : Rspline_info){
			const auto &Rt = disc_spline[info.spline_ref];
			auto Rfac = paramv_dir[data.strain[st].Rfactor_param];
		
			for(auto t = 0u; t < details.period; t++){
				auto sett = t*details.division_per_time;
				for(auto c : info.area) map[t][c] = Rfac*Rt[sett]*areafactor[t][c];
			}
		}
	
		rmap.file = "R_map"; if(data.nstrain > 1) rmap.file += "_"+data.strain[st].name; rmap.file += ".csv";			
		
		rmap.fulldesc = "Reproduction number &R&: This shows a map giving the spatial and temporal variation in the reproduction number &R& (the logarithmic colour scale is set such that green / red indicates that &R& is less than / greater than one, respectivly).  &R& is defined as the expected number of cases directly caused by an infected individual, assuming an otherwise suceptible population (note, this does not account for the fact that in reality the susceptible population reduces with time). Mathematically this map is calculated using &";
		if(Rspline_info.size() > 1) rmap.fulldesc += "R_{a,t}"; else rmap.fulldesc += "R_{t}";
		if(areaeffect_timevary() == false) rmap.fulldesc += "α_{a}"; else rmap.fulldesc += "α_{a,t}";
		if(data.nstrain > 1) rmap.fulldesc += "Ψ_{"+data.strain[st].name+"}";
		rmap.fulldesc += "&."; 
		rmap.fulldesc += "  Clicking on the play button animates &R& over time (note, the left and right arrow keys can be used to step single time units).";
		
		rmap.tab = details.analysis_type;
		rmap.tab2 = "Transmission";
		rmap.tab3 = "&R&"; 
		if(data.nstrain > 1) rmap.tab4 = data.strain[st].name+" Map"; 
		else rmap.tab4 = "Map";
		
		rmap_list.push_back(rmap);
		
		RMap effrmap;
		effrmap.map = rmap.map;
		auto &effmap = effrmap.map;
		
		auto Vinv = calculate_Vinv(transrate,st);
		
		for(auto c = 0u; c < data.narea; c++){  
			vector <unsigned int> ar; ar.push_back(c);
			auto dist = get_sus_dist_init(ar);
		
			for(auto t = 0u; t < details.period; t++){ 
				auto sett = t*details.division_per_time;
				effmap[t][c] /= calculate_R_beta_ratio_using_NGM(paramv_dir,susceptibility,Ntime[sett],Vinv,st,dist);
			}
			
			for(auto t = 0u; t < details.period; t++){ 
				auto sett = t*details.division_per_time;
				
				auto dist = get_sus_dist(sett,ar,pop);
			
				effmap[t][c] *= calculate_R_beta_ratio_using_NGM(paramv_dir,susceptibility,Ntime[sett],Vinv,st,dist);
			}
		}
		
		effrmap.file = "Eff_R_map"; if(data.nstrain > 1) effrmap.file += "_"+data.strain[st].name; effrmap.file += ".csv";			
		
		effrmap.fulldesc = "Effective reproduction number &R^{eff}&: This shows a map giving the spatial and temporal variation in the effective reproduction number &R^{eff}& (the logarithmic colour scale is set such that green / red indicates that &R^{eff}& is less than / greater than one, respectivly).  &R^{eff}& is defined to be the expected number of cases directly caused by an infected individual. This takes into account the fact that as the epidemic progresses the fraction of susceptible individuals reduces (causing effective transmission upon contact of individuals to become less and less common). &R^{eff}& is always less than &R& and if it reduces below 1 herd immunity in reached (i.e. the disease naturally dies out over time). ";
		
		effrmap.fulldesc += "  Clicking on the play button animates &R^{eff}& over time (note, the left and right arrow keys can be used to step single time units).";
		
		effrmap.tab = details.analysis_type;
		effrmap.tab2 = "Transmission";
		effrmap.tab3 = "&R^{eff}&"; 
		if(data.nstrain > 1) effrmap.tab4 = data.strain[st].name+" Map";
		else effrmap.tab4 = "Map";
		
		rmap_list.push_back(effrmap);
	}
	
	
	return rmap_list;
}


/// Calculates the generation time
double Model::calculate_generation_time(const vector<double> &paramv_dir, const vector <double> &susceptibility, const vector <vector <double> > &A, const vector < vector <double> > &transrate, const vector <double> &democatpos_dist, const unsigned int st) const
{
	vector <double> vec(data.nstrain);
	
	auto Q = data.ndemocatpos_per_strain;
	auto S = inf_state.size();
	auto N = Q*S;
	
	vector <double> comp_tsi(N);                                  // Calculate the compartment time since infection time
	for(auto i = 0u; i < N; i++) comp_tsi[i] = UNSET;
	
	bool flag; 
	do{
		flag = false;
		for(auto tr : trans){
			auto cfr = tr.from;
			auto from = inf_state_ref[cfr];   
			auto to = inf_state_ref[tr.to];   
			if(to != UNSET){
				for(auto dp = 0u; dp < Q; dp++){
					auto ito = Q*to+dp, ifrom = Q*from+dp;
						
					if(tr.inf == TRANS_INFECTION){
						if(comp_tsi[ito] == UNSET){ comp_tsi[ito] = 0; flag = true;}
					}
					else{
						if(from != UNSET){
							if(comp_tsi[ifrom] != UNSET && comp_tsi[ito] == UNSET){
								auto shape = comp[cfr].shape; if(shape == UNSET) emsgEC("Model",22);
								auto mean = paramv_dir[comp[cfr].param_mean[st*Q+dp]]/shape;
								comp_tsi[ito] = comp_tsi[ifrom] + mean;
								flag = true;
							}
						}
					}
				}
			}
		}
	}while(flag == true);
	
	auto Vinv = calculate_Vinv(transrate,st);

	auto F = calculate_F(paramv_dir,susceptibility,A,st,democatpos_dist);

	auto NGM = calculate_NGM(F,Vinv);
	
	for(auto j = 0u; j < NGM.size(); j++){
		for(auto i = 0u; i < NGM.size(); i++) if(NGM[i][j] < 0) emsg("NGM prob 2");
	}
	
	vector <double> eigenvector;
	largest_eigenvalue(NGM,eigenvector);
	
	if(false){
		for(auto c= 0u; c < comp.size(); c++){
			cout << comp[c].name_num << " " << inf_state_ref[c] << " Infectious states" << endl; 
		}
		
		for(auto th = 0u; th < param.size(); th++){
			cout << param[th].name << " " << paramv_dir[th] << "param" << endl;
		}
			
		for(auto j = 0u; j < N; j++){ 
			auto s = j/Q, dp = j%Q;
			cout << comp[inf_state[s]].name << " " << dp << " " << comp_tsi[j] << " " << Vinv[j][j] << endl;
		}
		emsg("Done");
	}
	
	auto sum = 0.0, sum2 = 0.0;
	for(auto i = 0u; i < Q; i++){
		auto num = eigenvector[i];
		for(auto j = 0u; j < N; j++){ 
			for(auto q = 0u; q < Q; q++){
				sum += F[q][j]*(comp_tsi[j]+0.5*Vinv[j][j])*Vinv[j][i]*num;
				sum2 += F[q][j]*Vinv[j][i]*num;
			}
		}
	}
	
	return sum/sum2;
}


/// Works out effective R_t
vector < vector <double> > Model::calculate_R_eff(const vector <double> &paramv_dir, const vector < vector < vector < vector <double> > > > &pop, const unsigned int st) const
{
	vector < vector <double> > R_eff; 
	if(pop.size() == 0) return R_eff;
	
	auto transrate = create_transrate(paramv_dir);
	auto susceptibility = create_susceptibility(paramv_dir);   
	auto disc_spline = create_disc_spline(paramv_dir);
	auto Ntime = create_Ntime(disc_spline);	
	auto Vinv = calculate_Vinv(transrate,st);
	auto fac = paramv_dir[data.strain[st].Rfactor_param];

	R_eff.resize(Rspline_info.size());
	for(auto i = 0u; i < Rspline_info.size(); i++){			
		const auto &info = Rspline_info[i];
	
		R_eff[i] = disc_spline[info.spline_ref];
		for(auto &val : R_eff[i]) val *= fac;
			
		auto dist = get_sus_dist_init(info.area);
		for(auto sett = 0u; sett < details.ndivision; sett++){
			R_eff[i][sett] /= calculate_R_beta_ratio_using_NGM(paramv_dir,susceptibility,Ntime[sett],Vinv,st,info.democatpos_dist);
		}

		for(auto sett = 0u; sett < details.ndivision; sett++){
			auto dist = get_sus_dist(sett,info.area,pop);
			R_eff[i][sett] *= calculate_R_beta_ratio_using_NGM(paramv_dir,susceptibility,Ntime[sett],Vinv,st,dist);
		}
	}
	
	return R_eff;
}


/// Works out R as a function of age for each of the infection transitions
vector < vector <double> > Model::calculate_R_age(const vector <double> &paramv_dir) const
{	
	auto st = 0u;
	auto transrate = create_transrate(paramv_dir);
	auto susceptibility = create_susceptibility(paramv_dir);   
	auto disc_spline = create_disc_spline(paramv_dir);
	auto Ntime = create_Ntime(disc_spline);
	auto Vinv = calculate_Vinv(transrate,st);
	auto dpmax = data.ndemocatpos_per_strain;
	
	vector < vector <double> > R_age;
	R_age.resize(data.nage);
	for(auto a = 0u; a < data.nage; a++){
		R_age[a].resize(details.ndivision);
		for(auto sett = 0u; sett < details.ndivision; sett++) R_age[a][sett] = 0;
	}	
		
	auto Rfac = paramv_dir[data.strain[st].Rfactor_param];
	for(const auto &info : Rspline_info){                          // Goes through all Rsplines
		auto frac_pop = 0.0; for(auto c : info.area) frac_pop += data.area[c].total_pop/data.popsize;
		
		double ratio;
		vector <double> vec(dpmax);
		for(auto sett = 0u; sett < details.ndivision; sett++){
			if(sett == 0 || equal(Ntime[sett-1],Ntime[sett]) == false){
				auto F = calculate_F(paramv_dir,susceptibility,Ntime[sett],st,info.democatpos_dist);
				auto NGM = calculate_NGM(F,Vinv);

				for(auto j = 0u; j < NGM.size(); j++){
					for(auto i = 0u; i < NGM.size(); i++) if(NGM[i][j] < 0) emsg("NGM prob 3");
				}
	
				vector <double> eigenvector;
				ratio = largest_eigenvalue(NGM,eigenvector);
	
				if(NGM.size() != dpmax){
					emsg("The NGM is not the right size. It is size '"+to_string(NGM.size())+"' instead of '"+to_string(data.ndemocatpos_per_strain)+"'");
				}
				
				for(auto dp = 0u; dp < dpmax; dp++){
					vec[dp] = 0; for(auto dpp = 0u; dpp < dpmax; dpp++) vec[dp] += NGM[dpp][dp];
				}
			}
			
			
			for(auto dp = 0u; dp < dpmax; dp++){
				auto a = data.democatpos[dp][0];
				R_age[a][sett] += frac_pop*Rfac*disc_spline[info.spline_ref][sett]*vec[dp]/ratio;
			}
		}
	}
	
	return R_age;
}


/// Calculate Vinv (used to calculate next generation matrix)
vector < vector <double> > Model::calculate_Vinv(const vector < vector <double> > &transrate, const unsigned int st) const
{
	auto Q = data.ndemocatpos_per_strain;
	auto N = Q*inf_state.size();
	
	vector < vector <double> > V;
	V.resize(N);
	for(auto j = 0u; j < N; j++){
		V[j].resize(N);
		for(auto i = 0u; i < N; i++) V[j][i] = 0;
	}
	
	for(auto tr = 0u; tr < trans.size(); tr++){
		auto from = inf_state_ref[trans[tr].from], to = inf_state_ref[trans[tr].to];
		if(from != UNSET){
			for(auto dp = 0u; dp < Q; dp++){
				auto rate = transrate[tr][st*Q+dp];
				auto j = from*Q + dp;
				V[j][j] += rate;
				if(to != UNSET){
					auto k = to*Q + dp;
					V[k][j] -= rate;
				}
			}
		}
	}
	
	if(false){
		cout << st << " strain" << endl;
		for(auto j = 0u; j < N; j++){
			for(auto i = 0u; i < N; i++) cout << V[j][i] << ",";
			cout << "V" << endl;
		}
		emsg("V");
	}
	
	return invert_matrix(V);
}


/// Calculates the matrix giving the source of new infections based on individuals in different states
vector < vector <double> > Model::calculate_F(const vector <double> &paramv_dir, const vector <double> &susceptibility, const vector <vector <double> > &A, const unsigned int st, const vector <double> &democatpos_dist) const
{
	vector < vector <double> > F;
	auto Q = data.ndemocatpos_per_strain;
	auto S = inf_state.size();
	auto N = Q*S;

	F.resize(Q);
	for(auto q = 0u; q < Q; q++){
		F[q].resize(N);
		for(auto i = 0u; i < N; i++){
			F[q][i] = 0;
		}
	}
	
	for(auto s = 0u; s < S; s++){
		auto inf = paramv_dir[comp[inf_state[s]].infectivity_param];
		if(inf != 0){
			for(auto dp = 0u; dp < Q; dp++){
				auto i = s*Q + dp;
				auto a = data.democatpos[dp][0];
				for(auto q = 0u; q < Q; q++){
					auto aa = data.democatpos[q][0];
					F[q][i] += democatpos_dist[q]*susceptibility[Q*st + q]*inf*A[aa][a];
				}
			}
		}		
	}
	
	return F;
}


/// Calculates the next generation matrix
vector < vector <double> > Model::calculate_NGM(const vector < vector<double> > &F, const vector < vector<double> > &Vinv) const
{
	vector < vector <double> > NGM;	
	
	auto Q = F.size();
	auto N = Vinv.size();
	
	NGM.resize(Q);
	for(auto q = 0u; q < Q; q++){
		NGM[q].resize(Q);
		for(auto i = 0u; i < Q; i++){
			auto sum = 0.0;	for(auto j = 0u; j < N; j++) sum += F[q][j]*Vinv[j][i];
			NGM[q][i] = sum;
		}
	}
	return NGM;
}


/// Estimates the ratio between R and beta based on largest eigenvector of next generation matrix
double Model::calculate_R_beta_ratio_using_NGM(const vector<double> &paramv_dir, const vector <double> &susceptibility, const vector <vector <double> > &A, const vector < vector <double> > &Vinv, const unsigned int st, const vector <double> &democatpos_dist) const
{
	auto F = calculate_F(paramv_dir,susceptibility,A,st,democatpos_dist);
	
	auto NGM = calculate_NGM(F,Vinv);
	
	if(checkon == true){
		for(auto j = 0u; j < NGM.size(); j++){
			for(auto i = 0u; i < NGM.size(); i++){
				if(NGM[i][j] < 0) emsg("NGM prob 1");
			}
		}
	}
	
	vector <double> eigenvector;
	auto ratio = largest_eigenvalue(NGM,eigenvector);
	
	if(false){ for(auto val: eigenvector) cout << val << ","; cout << "eigen" << endl;}
	
	if(false){
		auto N = Vinv.size();
		for(auto j = 0u; j < N; j++){
			for(auto i = 0u; i < N; i++) cout << Vinv[j][i] << ",";
			cout << "Vinv" << endl;
		}
		
		for(auto dp = 0u; dp < data.ndemocatpos_per_strain; dp++){
			for(auto i = 0u; i < N; i++) cout << F[dp][i] << ",";
			cout << " FF" << endl;
		}
		
		for(auto q = 0u; q < NGM.size(); q++){
			for(auto qq = 0u; qq < NGM.size(); qq++) cout << NGM[q][qq] << ",";
			cout << " NGM" << endl;
		}
		
		cout << ratio << "Ratio" << endl;
		emsg("done");		
	}
	
	return ratio;
}


/// Returns true if all the parameters are within the ranges permited by the prior
bool Model::inbounds(const vector <double> &paramv) const
{
	for(auto th = 0u; th < paramv.size(); th++){
		auto val = paramv[th];
		
		switch(param[th].priortype){
			case FIXED_PRIOR:
				if(val != param[th].val1) return false;
				break;
				
			case UNIFORM_PRIOR:
				if(val < param[th].val1 || val > param[th].val2) return false;
				break;
				
			case EXP_PRIOR:
				if(val < 0) return false;
				break;
				
			case NORMAL_PRIOR:
				break;
				
			case DIRICHLET_PRIOR: case DIRICHLET_ALPHA_PRIOR: case MDIRICHLET_PRIOR:
				if(val < 0) return false;
				break;
				
			case DIRICHLET_FLAT_PRIOR:
				if(val < 0) return false;
				break;
				
			default:
				emsgEC("model",23);
				break;
		}
		
		for(auto th2 : param[th].greaterthan){
			if(val < paramv[th2]) return false;
		}
		
		for(auto th2 : param[th].lessthan){
			if(val > paramv[th2]) return false;
		}
	}	
	
	return true;
}


/// Determines if it is necessary to do MBPs for compartmental transitions
bool Model::do_mbp_events(const vector <double> &parami, const vector <double> &paramp) const
{
	for(auto th = 0u; th < param.size(); th++){
		if(parami[th] != paramp[th] && (param[th].type == DISTVAL_PARAM || param[th].type == BRANCHPROB_PARAM)){
			return true;
		}
	}
	return false;
}


/// Takes the raw MCMC parameters and converts them into Dirichlet distributed variables 
vector <double> Model::dirichlet_correct(const vector <double> &paramv) const
{
	auto paramv_dir = paramv;
	
	for(const auto &dir : dirichlet){
		auto sum = 0.0; for(auto par : dir.param) sum += paramv[par.th];
		for(auto par : dir.param) paramv_dir[par.th] = paramv[par.th]/(par.frac*sum);
	}	
	
	return paramv_dir;
}

	
/// Adds a (potential) Dirichlet distribution by adding a list of parameters and a list of fractions for each
void Model::add_dirichlet(unsigned int th_min, const vector <unsigned int> &th_list, const vector <double> &frac, const DirType dirtype) 
{	 
	Dirichlet dir; dir.th_min = th_min; dir.dirtype = dirtype;
	for(auto i = 0u; i < th_list.size(); i++){
		auto th = th_list[i];
		
		auto j = 0u; while(j < dir.param.size() && th != dir.param[j].th) j++;
		if(j == dir.param.size()){
			DirichletParam par; par.th = th; par.frac = frac[i];
			dir.param.push_back(par);
		}
		else dir.param[j].frac += frac[i];
	}
			
	bool flag = true;
	auto imax = dir.param.size();
	for(const auto &dir_prev : dirichlet){                       // Checks if Dirichlet already exists
		if(dir_prev.param.size() == imax){
			auto i = 0u; while(i < imax && dir_prev.param[i].th == dir.param[i].th 
			                            && dir_prev.param[i].frac == dir.param[i].frac) i++;
			if(i == imax){ flag = false; break;} 
		}
	}
	
	if(flag == true){
		for(const auto &par : dir.param){                          // Checks parameters are unique to dirichlet distribution
			if(par.th < dir.th_min) emsgroot("The parameter '"+param[par.th].name+"' is duplicated and must be unique.");
		}
	
		dirichlet.push_back(dir);
	}
}


/// Generate a string listing a a series of Dirichlet parameters
string Model::list_dir_param(const vector <DirichletParam> &list) const 
{
	string st;
	bool flag = false;
	for(const auto &par : list){
		if(flag == true) st += ","; 
		flag = true;
		st += "'"+param[par.th].name+"' ";
	}
	return st;
}


/// Sets up the dirichlet distributions
void Model::set_dirichlet()
{
	const auto &le = data.level_effect;                               // Level effects
	if(le.on == true){
		auto th_min = param.size();
		for(auto i = 0u; i < le.ps.size(); i++){
			level_effect_param_list.push_back(add_parameter(le.ps[i],LEVEL_EFFECT_PARAM));
		}
		add_dirichlet(th_min,level_effect_param_list,le.frac,DIR_MODIFIED);
	}
	
	const auto &ae = data.area_effect;                                // Area effects
	if(ae.on == true){
		auto th_min = param.size();
		for(auto i = 0u; i < ae.ps.size(); i++){
			area_effect_param_list.push_back(add_parameter(ae.ps[i],AREA_EFFECT_PARAM));
		}
		add_dirichlet(th_min,area_effect_param_list,ae.frac,DIR_MODIFIED);
	}
	
	
	for(auto c = 0u; c < data.ndemocat; c++){                         // Susceptibility
		auto &dc = data.democat[c];
		if(dc.sus_vari == true){	
			auto si = dc.value.size();
			
			auto th_min = param.size(); vector <double> frac;
			dc.sus_param.resize(si);
			for(auto fi = 0u; fi < si; fi++){	
				dc.sus_param[fi] = add_parameter(dc.ps[fi],SUSCEPTIBILITY_PARAM);
				if(c == data.ndemocat-1) frac.push_back(1);
				else frac.push_back(data.democat_dist[c][fi]);
			}
		
			add_dirichlet(th_min,dc.sus_param,frac,DIR_MODIFIED);
		}
	}
	
	auto i=0u;                                                        // Matrix modify perturb
	const auto &mm = data.genQ.matmod;
	while(i < mm.size() && mm[i].type != PERTURB) i++;
	if(i < mm.size()){
		auto np = spline[mm[i].spline_ref].p.size(); 
	
		auto ii = i; while(ii < mm.size() && mm[ii].type == PERTURB) ii++;
		auto num = ii-i;
		
		vector <bool> ageuse(data.nage); 
		for(auto a = 0u; a < data.nage; a++) ageuse[i] = false;
	
		vector <double> frac(num);
		for(auto j = 0u; j < num; j++){
			frac[j] = 0.0;
			for(auto a : mm[i+j].ages){
				frac[j] += data.democat_dist[0][a];
				if(ageuse[a] != false){
					emsg("'age_matrix_perturb' should not use age category '"+data.democat[0].value[a]+"' more than once.");
				}
				ageuse[a] = true;
			}
		}
		
		for(auto a = 0u; a < data.nage; a++){
			if(ageuse[a] == false){
				emsg("'age_matrix_perturb' should use age category '"+data.democat[0].value[a]+"'.");
			}
		}		
		
		for(auto p = 0u; p < np; p++){
			auto th_min = 0;
			vector <unsigned int> param_list;
			for(auto j = 0u; j < num; j++){
				auto th = spline[mm[i+j].spline_ref].p[p].param;
				param_list.push_back(th);
			}
			add_dirichlet(th_min,param_list,frac,DIR_MODIFIED);
		}
	}
	
	for(const auto &co : comp){                                       // Branching probabilities
		auto kmax = co.trans.size();
		if(kmax > 1){
			for(auto dp = 0u; dp < data.ndemocatpos; dp++){
				auto th_min = 0; vector <unsigned int> th_list; vector <double> frac;
				for(auto k = 0u; k < kmax; k++){ th_list.push_back(trans[co.trans[k]].param_prob[dp]); frac.push_back(1);}
				add_dirichlet(th_min,th_list,frac,DIR_NORM);
			}			
		}
	}	
	
	for(const auto &dir : dirichlet){                                 // This adjusts the value used in simulation
		if(param[dir.param[0].th].value != UNSET){
			vector <DirichletParam> tobeset;
			auto sum = 0.0;
			for(auto i = 0u; i < dir.param.size(); i++){
				const auto &par = dir.param[i];
				auto th = par.th;
				if(param[th].value == UNSET) emsgroot("The parameter '"+param[th].name+"' must be set");
				if(param[th].value == TOBESET) tobeset.push_back(par);
				else{
					if(param[th].priortype == FIXED_PRIOR) param[th].val1 *= par.frac;
					sum += param[th].value*par.frac; 
				}
			}
			
			if(tobeset.size() == 0){
				emsgroot("One of the parameters: "+list_dir_param(dir.param)+"must have a value set to '*'");
			}

			if(tobeset.size() > 1){
				emsgroot("All but one of the parameters: "+list_dir_param(tobeset)+"must have their value set");
			}

			if(tobeset.size() == 1){
				if(sum > 1) emsgroot("The value of the parameter '"+param[tobeset[0].th].name+"' cannot be set (because the average is over one)");
			
				auto th = tobeset[0].th;
				auto val = (1-sum)/tobeset[0].frac;
				
				param[th].value = val; 
				if(param[th].priortype == FIXED_PRIOR){
					param[th].val1 = 1-sum;
					param[th].fixed = val;
				}
		
				if(details.siminf == SIMULATE && mpi.core == 0){
					cout << "Parameter '"+param[tobeset[0].th].name+"' has been set to "+to_string(val) << endl;
				}
			}
		}
	}
	
	if(details.siminf == INFERENCE || details.siminf == DATAVIEW){
		for(const auto &dir : dirichlet){                         // Sets the prior distributions used by the Dirichtlet
			switch(dir.dirtype){
			case DIR_NORM:
				{
					auto nfixed = 0u, ndir = 0u, ndiral = 0u, ndirflat = 0u;
					for(auto par : dir.param){
						auto th = par.th;
						if(param[th].priortype == FIXED_PRIOR) nfixed++;
						if(param[th].priortype == DIRICHLET_PRIOR) ndir++;
						if(param[th].priortype == DIRICHLET_ALPHA_PRIOR) ndiral++;
						if(param[th].priortype == DIRICHLET_FLAT_PRIOR) ndirflat++;
					}			
					
					if(nfixed == dir.param.size()){             // The case of fixed priors
						auto sum = 0.0;
						vector <unsigned int> unset_list;
						for(auto par : dir.param){
							auto th = par.th; 
							if(param[th].val1 == TOBESET) unset_list.push_back(th);
							else{
								param[th].value = param[th].val1;
								sum += param[th].value;
							}
						}
						
						auto flag = false;
						switch(unset_list.size()){
							case 0: if(sum < 1-TINY || sum > 1+TINY) flag = true; break;
							case 1:
								if(sum > 1) flag = true;
								param[unset_list[0]].value = 1-sum;
								param[unset_list[0]].val1 = 1-sum;
								break;
							default: flag = true; break;
						}
						
						if(flag == true){
							stringstream ss; ss << "The prior for parameters: ";
							for(auto par : dir.param) ss << "'" << param[par.th].name << "' ";
							ss << "is not valid";
							emsgroot(ss.str());
						}
					}	
					else{
						if(ndir != dir.param.size() && ndirflat != dir.param.size() && ndiral != dir.param.size()){
							stringstream ss; ss << "The parameters: ";
							for(auto par : dir.param) ss << param[par.th].name << " ";
							ss << "must either all be fixed (by setting a value and leaving the prior unspecified)";
							ss << ", all set to a prior value of 'Dir(*)', which implies a flat Dirichlet distribution, or set to 'Dir(α) or 'Dir(mean,sd)' for a constraining prior (see manual).";
							emsgroot(ss.str());
						}
				
						if(ndiral == dir.param.size()){    // This is the case when the values of alpha are specified							
						}
						else{                              // This is the case when mu and sd are specified
							if(ndir == dir.param.size()){
								vector <DirichletParam> tobeset;
								auto sum = 0.0;
								for(auto i = 0u; i < dir.param.size(); i++){
									const auto &par = dir.param[i];
									auto th = par.th;
									if(param[th].dir_mean == UNSET){
										emsgroot("The Dirichlet mean for the parameter '"+param[th].name+"' must be set");
									}

									if(param[th].dir_mean == TOBESET) tobeset.push_back(par);
									else sum += param[th].dir_mean*par.frac; 
								}
							
								if(tobeset.size() == 0){
									emsgroot("The mean of the Dirichlet prior for one of the parameters: "+list_dir_param(dir.param)+"must have a value set to '*'");
								}
							
								if(tobeset.size() > 1){
									emsgroot("For the means of the Dirichlet prior, all but one of the parameters: "+list_dir_param(tobeset)+"must have their value set");
								}
							
								if(tobeset.size() == 1){
									if(sum > 1){
										emsgroot("The value for the Dirichlet prior mean of the parameter '"+param[tobeset[0].th].name+"' cannot be set");
									}
									
									auto th = tobeset[0].th;
									auto val = (1-sum)/tobeset[0].frac;
								
									param[th].dir_mean = val; if(param[th].priortype == FIXED_PRIOR) param[th].val1 = 1-sum;
						
									if(details.siminf == SIMULATE && mpi.core == 0){
										cout << "The Dirichlet mean for parameter '"+param[tobeset[0].th].name+"' has been set to "+to_string(val) << endl;
									}
								}
					
								vector <DirichletParam> sdset;
								for(auto i = 0u; i < dir.param.size(); i++){
									const auto &par = dir.param[i];
									auto th = par.th;
									if(param[th].dir_sd == UNSET) emsgroot("The Dirichlet standard deviation '"+param[th].name+"' must be set");
									if(param[th].dir_sd != TOBESET) sdset.push_back(par);
								}
								
								if(sdset.size() == 0){
									emsgroot("The standard deviation of the Dirichlet prior, for one of the parameters: "+list_dir_param(dir.param)+"must be set");
								}
						
								if(sdset.size() > 1){
									emsgroot("The standard deviation of the Dirichlet prior, for only one of the parameters: "+list_dir_param(dir.param)+"must be set");
								}
								
								auto th = sdset[0].th;
								auto mu = param[th].dir_mean*sdset[0].frac;
								auto var = param[th].dir_sd*param[th].dir_sd*sdset[0].frac*sdset[0].frac;
								auto alpha0 = mu*(1-mu)/var - 1;
								
								sum = 0.0;
								for(auto i = 0u; i < dir.param.size(); i++){
									auto th = dir.param[i].th;
									param[th].val1 = param[th].dir_mean*dir.param[i].frac;
									param[th].val2 = 1; 
									sum += param[th].val1;
								}
								
								for(auto i = 0u; i < dir.param.size(); i++){
									auto th = dir.param[i].th;
									param[th].val1 *= alpha0/sum;
								}
								
								double a0=0;
								for(auto i = 0u; i < dir.param.size(); i++){
									auto th = dir.param[i].th;
									a0 += param[th].val1;
								}
								
								for(auto i = 0u; i < dir.param.size(); i++){
									auto th = dir.param[i].th;
									if(param[th].dir_sd == TOBESET){
										param[th].dir_sd = sqrt((param[th].val1/a0)*(1-(param[th].val1/a0))/(a0+1))/dir.param[i].frac;
									}
								}
							}
						}
					}
				}
				break;
					
			case DIR_MODIFIED:
				{
					auto nfixed = 0u, nmdir = 0u;
					for(auto par : dir.param){
						auto th = par.th;
						if(param[th].priortype == FIXED_PRIOR) nfixed++;
						if(param[th].priortype == MDIRICHLET_PRIOR) nmdir++;
					}
					
					if(nfixed != dir.param.size()){
						if(nmdir != dir.param.size()){
							stringstream ss; ss << "The parameters: ";
							for(auto par : dir.param) ss << param[par.th].name << " ";
							ss << "must either all be fixed (by setting a value and leaving the prior unspecified)";
							ss << " or set to 'MDir(k)' where hyperparameter set how constraining the prior is (see manual).";
							emsgroot(ss.str());
						}
					}
							
					if(nmdir == dir.param.size()){ // This applies a modified Dirichlet prior
						auto kappa = param[dir.param[0].th].dir_mean;
						// cout << kappa << "kappa\n";
						auto beta = ((dir.param.size()-1)/(kappa*kappa))-1;
						
						for(auto i = 0u; i < dir.param.size(); i++){
							auto th = dir.param[i].th;
							
							param[th].val1 = beta*dir.param[i].frac;
							if(param[th].val1 < 0.1){
								emsgroot("The prior '"+param[th].ps.prior+"' is too diffuse. Please reduce this number.");
							}
							param[th].val2 = 1; 
						}
					}
				}
				break;
			}
		}	
	}
	
	if(false){
		for(const auto &dir : dirichlet){  
			for(auto dp : dir.param){
				cout << param[dp.th].name << ", ";
			}
			cout << " Dirichlet params\n";
		}
	}
}


/// If N-1 direchlet varaibles are set then this calculates the last missing value
void Model::calculate_dirichlet_missing(vector <double> &param) const
{
	for(const auto &dir : dirichlet){                       // Fills in dependent diriclet values
		auto sum = 0.0;				
		for(auto j = 1u; j < dir.param.size(); j++){
			sum += param[dir.param[j].th];
		}
		
		if(sum > 1) emsg("Dirichlet varaibles out of range");
		param[dir.param[0].th] = 1-sum;
	}
}
	

/// This checks the prior by sampling from it and 
void Model::check_prior() const 
{
	const auto imax = 1000000u;
	
	vector <double> av(param.size()), av2(param.size());
	
	for(auto th = 0u; th < param.size(); th++){ av[th] = 0; av2[th] = 0;}
		
	for(auto i = 0u; i < imax; i++){
		if(i%1000 == 0) cout << i << "i" << endl;
		auto paramv = sample_from_prior();
		auto paramv_dir = dirichlet_correct(paramv);
	
		for(auto th = 0u; th < param.size(); th++){ 
			av[th] += paramv_dir[th];
			av2[th] += paramv_dir[th]*paramv_dir[th];
		}
	}

	for(const auto &dir : dirichlet){
		cout << "Dir" << endl;
		for(auto par : dir.param){
			auto th = par.th;
			cout << th << " " << param[th].name << " " << av[th]/imax << " " << (av2[th]/imax) - (av[th]/imax)*(av[th]/imax) << " : " << param[th].dir_mean*param[th].dir_mean*(1-par.frac)/((dir.param.size()-1)*par.frac) << endl;
		}
	}	

	emsg("Done");	
}


/// Provides a normal approximation to the prior
NormalApprox Model::prior_normal_approx(unsigned int th) const 
{
	NormalApprox na;
	const auto &pa = param[th];
	
	switch(pa.priortype){
		case FIXED_PRIOR:
			if(pa.val1 != UNSET){
				if(pa.fixed == TOBESET) na.mean = pa.val1;
				else na.mean = pa.fixed;
			}
			else{
				na = prior_normal_approx(pa.val1_param);
			}
			na.var = 0;
			break;			
		
		case UNIFORM_PRIOR: 
			{
				auto val1 = pa.val1;
				if(val1 == UNSET){ auto naa = prior_normal_approx(pa.val1_param); val1 = naa.mean;}
				
				auto val2 = pa.val2;
				if(val2 == UNSET){ auto naa = prior_normal_approx(pa.val2_param); val2 = naa.mean;}
				
				na.mean = (val1 + val2)/2;
				na.var = (val2 - val1)*(val2 - val1)/4;
			}
			break;
	
		case DIRICHLET_PRIOR: case DIRICHLET_ALPHA_PRIOR: case MDIRICHLET_PRIOR:
			{
				auto val1 = pa.val1;
				if(val1 == UNSET){ auto naa = prior_normal_approx(pa.val1_param); val1 = naa.mean;}
				
				auto val2 = pa.val2;
				if(val2 == UNSET){ auto naa = prior_normal_approx(pa.val2_param); val2 = naa.mean;}
					
				auto alpha = val1;
				auto beta = val2;
				na.mean = alpha/beta;
				na.var = alpha/(beta*beta);
						
				if(details.mode == ML_INF){                 // This account for the fact that diricehlets are defined differently
					auto inside = false;
					for(const auto &dir : dirichlet){	
						for(auto j = 0u; j < dir.param.size(); j++){
							if(dir.param[j].th == th) inside = true;
						}
						if(inside == true){
							auto mean_sum = 0.0;
							for(auto j = 0u; j < dir.param.size(); j++){
								const auto &pa2 = param[dir.param[j].th];
								
								auto val1 = pa2.val1;
								if(val1 == UNSET){ auto naa = prior_normal_approx(pa2.val1_param); val1 = naa.mean;}
				
								auto val2 = pa2.val2;
								if(val2 == UNSET){ auto naa = prior_normal_approx(pa2.val2_param); val2 = naa.mean;}
		
								mean_sum += val1/val2;
							}
							
							na.mean /= mean_sum;
							na.var /= (mean_sum*mean_sum);
							break;
						}
					}
					if(inside == false) emsgEC("Model",45);
				}
			}
			break;
			
		case EXP_PRIOR:
			{
				auto val1 = pa.val1;
				if(val1 == UNSET){ auto naa = prior_normal_approx(pa.val1_param); val1 = naa.mean;}
				
				na.mean = 1.0/val1; na.var = 1.0/(val1*val1);
			}
			break;
			
		case NORMAL_PRIOR:
			{
				auto val1 = pa.val1;
				if(val1 == UNSET){ auto naa = prior_normal_approx(pa.val1_param); val1 = naa.mean;}
				
				auto val2 = pa.val2;
				if(val2 == UNSET){ auto naa = prior_normal_approx(pa.val2_param); val2 = naa.mean;}
				
				na.mean = val1; na.var = val2*val2;
			}
			break;
		
		case DIRICHLET_FLAT_PRIOR:
			na.mean = 1; na.var = 1;
			break;
	
		default: emsg("Not done yet"); break;
	}
	
	return na;
}


/// Prints a set of parameters to the terminal
void Model::print_param(const string name, const vector <double> &paramv) const 
{
	cout << endl << name << ":" << endl;
	for(auto th = 0u; th < param.size(); th++){
		cout << param[th].name << " " << paramv[th] << endl;
	}
}


/// Prints the prior for a parameter
string Model::print_prior(const unsigned int p) const
{
	stringstream ss;
	switch(param[p].priortype){
		case FIXED_PRIOR:
			ss << "Fixed(";
			if(param[p].val1 != UNSET){
				if(param[p].fixed == TOBESET) ss << param[p].val1;
				else ss << param[p].fixed;
			}
			else{
				ss << param[param[p].val1_param].name;
			}
			ss << ")";
			break;
		
		case UNIFORM_PRIOR:
			ss << "Uniform(";
			if(param[p].val1 != UNSET) ss << param[p].val1; else ss << param[param[p].val1_param].name;
			ss << ","; 
			if(param[p].val2 != UNSET) ss << param[p].val2; else ss << param[param[p].val2_param].name;
			ss << ")";
			break;
			
		case EXP_PRIOR:
			ss << "Exp(";
			if(param[p].val1 != UNSET) ss << param[p].val1; else ss << param[param[p].val1_param].name;
			ss << ")";
			break;
			
		case NORMAL_PRIOR:
			ss << "Normal(";
			if(param[p].val1 != UNSET) ss << param[p].val1; else ss << param[param[p].val1_param].name;
			ss << ","; 
			if(param[p].val2 != UNSET) ss << param[p].val2; else ss << param[param[p].val2_param].name;
			ss << ")";
			break;
			
		case DIRICHLET_PRIOR:
			ss << "Dir(";
			ss << param[p].dir_mean;
			ss << ","; 
			ss << param[p].dir_sd;
			ss << ")";
			break;
			
		case DIRICHLET_ALPHA_PRIOR:
			ss << "Dir(";
			ss << param[p].val1;
			ss << ")";
			break;
			
		case MDIRICHLET_PRIOR:
			ss << "MDir(";
			ss << param[p].dir_mean;
			ss << ")";
			break;
			
		case DIRICHLET_FLAT_PRIOR:
			ss << "Dir(*)";
			break;
			
		default:
			emsgEC("model",24);
			break;
	}
						
	return ss.str();
}


/// Outputs a summary of the model
string Model::print() const
{                      
	stringstream ss;
                
	switch(details.siminf){
		case SIMULATE:
			ss << "Parameters:" << endl;
			for(auto p = 0u; p < param.size(); p++){
				if(param[p].name != "zero" && param[p].name != "one"){
					ss << "  " << param[p].name << " = " << param[p].val1 << endl;
				}
			}
			break;
		
		case INFERENCE: case DATAVIEW:
			ss << "Priors:" << endl;
			for(auto p = 0u; p < param.size(); p++){
				ss << "  " << param[p].name << " = ";
				ss << print_prior(p);
			
				if(false && param[p].param_dep.size() > 0){
					ss << " Dependencies: "; 	for(auto th : param[p].param_dep) ss << param[th].name << " ";
				}
					
				for(auto th : param[p].greaterthan) ss << " Greater than " << param[th].name << " ";
				for(auto th : param[p].lessthan) ss << " Less than " << param[th].name << " ";
				
				ss << endl;
			}
			
			if(dirichlet.size() > 0){
				ss << endl << "Dirichlet distributions:" << endl;
				for(const auto &dir : dirichlet){
					ss << "  ";
					for(auto par : dir.param) ss << param[par.th].name << " " <<  par.frac << "   ";
					ss << endl;
				}
			}
			break;
		
		default: 
			emsgEC("Model",25);
			break;
	}
	ss << endl;
		
	ss << "Compartments:" << endl; 
	auto num = 0u;
	for(const auto &co : comp){
		ss << num << "  "; num++;
		for(auto tr : co.trans) ss << tr << ","; ss << "   ";
			
		ss << "  " << co.name;
		if(co.num != UNSET && co.num > 0) ss << co.num;
		if(co.shape != UNSET){
			if(co.shape == 1) ss << "   Exponential";
			else ss << "   Erlang";
			
			ss <<	"    mean = ";
			for(auto i = 0u; i < co.param_mean.size(); i++){
				auto th = co.param_mean[i];
				if(i != 0) ss << " | ";
				if(co.shape > 1) ss << 1.0/co.shape << "*";
				ss << "'" << param[th].name << "'";
			}
		}
		
		auto th = co.infectivity_param;
		if(param[th].name != "zero") ss << "  Infectivity: " << param[th].name;
		if(co.mean_dep != "") ss << "  Dependency: " << co.mean_dep;
		ss << endl; 		
	}
	ss << endl;
	
	ss << "Transitions:" << endl; 
	auto num2 = 0u;
	for(const auto &tr : trans){
		ss << num2 << "   " << tr.from << "->" << tr.to << "   "; num2++;
		const auto &coi = comp[tr.from];
		const auto &cof = comp[tr.to];
		
		ss << "  " << coi.name; if(coi.num != UNSET && coi.num != 0) ss << coi.num;
		ss << "->" << cof.name; if(cof.num != UNSET && cof.num != 0) ss << cof.num;
		
		if(tr.inf == TRANS_INFECTION) ss << " Infection";
		
		if(tr.param_prob.size() > 0){
			ss << " with probability ";
			for(auto j = 0u; j < tr.param_prob.size(); j++){
				if(j > 0) ss << " | ";
				ss << "'" << param[tr.param_prob[j]].name << "'";
			}
			ss << endl;
		}
		
		if(tr.prob_dep != "") ss << "  Dependency: " << tr.prob_dep;
		ss << endl;
	}
	ss << endl;
	
	if(prob_reach.size() > 0){
		ss << "Probability of reaching:" << endl; 
		for(auto pr = 0u; pr < prob_reach.size(); pr++){
			auto probr = prob_reach[pr];
			ss << "  Name = " << probr.name << "   Compartment = " << comp[probr.comp].name << endl;
		}
	}
	
	return ss.str();
}


/// Prints all the parmaeter types
void Model::print_parameter_types()
{
	for(auto &pa: param){
		cout << pa.name << ": ";
		switch(pa.type){
			case R_EFOI_PARAM: cout << "R / efoi" << endl; break;
			case GEOMIX_PARAM: cout << "Geographic mixing" << endl; break;
			case MODIFY_PARAM: cout << "Modify contact matrix" << endl; break;
			case SIGMA_PARAM: cout << "Sigma" << endl; break;
			case COVAR_PARAM: cout << "Covariate" << endl; break;
			case DISTVAL_PARAM: cout << "Distribution value" << endl; break;
			case BRANCHPROB_PARAM: cout << "Branching probability" << endl; break;
			case RE_PARAM: cout << "Regional effect" << endl; break;
			case OBS_PARAM: cout << "Observation parameter" << endl; break;
			case INF_PARAM: cout << "Infectivity parameter" << endl; break;
			case SUSCEPTIBILITY_PARAM: cout << "Susceptibility parameter" << endl; break; 
			case LEVEL_EFFECT_PARAM: cout << "Level effect parameter" << endl; break; 
			default: emsgEC("Model",26); break;
		}
	}
}


/// Looks at different models to estimate eigenvectors (this is used in the raw data analysis)
void Model::eignevector_compare_models(const vector <double> &susceptibility, const vector <double> &paramv_dir, const vector < vector < vector <double> > > &Ntime, const vector < vector <double> > &transrate) const
{
	auto Vinv = calculate_Vinv(transrate,0);
		
	auto F = calculate_F(paramv_dir,susceptibility,Ntime[0],0,data.democat_dist[0]);
	 
	auto NGM = calculate_NGM(F,Vinv);
	
	vector <double> eigenvector;
	largest_eigenvalue(NGM,eigenvector);
	
	auto nage = susceptibility.size();
	vector <double> sus(nage);
	for(auto a = 0u; a < nage; a++) sus[a] = 1;
	sus[nage-1] = 1.8;
	
	// This is the vector we are aiming for 
	vector <double> aim = {0.041777881,0.083194264,0.091183693,0.089258079,0.086831252,0.077545859,0.064406066,0.056545308,0.056864718,0.068969452,0.078228807,0.075130564,0.042264548,0.026687841,0.022323116,0.019930333,0.018858221};
	
	// This is for getting the relative susceptibility
	for(auto loop = 0u; loop < 10000; loop++){
		auto F = calculate_F(paramv_dir,sus,Ntime[0],0,data.democat_dist[0]);
		auto NGM = calculate_NGM(F,Vinv);
		
		largest_eigenvalue(NGM,eigenvector);
			
		for(auto a = 0u; a < nage-1; a++){
			sus[a] += 0.1*(aim[a]-eigenvector[a]);
		}
		
		auto sum = 0.0;
		for(auto a = 0u; a < nage; a++){
			sum += data.democat_dist[0][a]*sus[a]; 
		}
		
		for(auto a = 0u; a < nage; a++) sus[a] /= sum;
		sus[nage-1] = 1.8;
	}
	for(auto val : sus) cout << val << endl; 
	 cout << " relative susceptibility" << endl;	
	
	vector <double> v(nage);
	for(auto a = 0u; a < nage; a++) v[a] = 1;
	v[nage-1] = 1.6;
	
	for(auto a = 0u; a < nage; a++) sus[a] = 1;
	
	// This is for getting the adjusted age mixing matrix
	for(auto loop = 0u; loop < 10000; loop++){
		auto A = Ntime[0];
		
		for(auto j = 0u; j < nage; j++){
			for(auto i = 0u; i < nage; i++){
				A[j][i] *= v[j]*v[i];
			}
		}
		
		auto F = calculate_F(paramv_dir,sus,A,0,data.democat_dist[0]);
		auto NGM = calculate_NGM(F,Vinv);
		
		largest_eigenvalue(NGM,eigenvector);
			
		for(auto a = 0u; a < nage-1; a++){
			v[a] += 0.1*(aim[a]-eigenvector[a]);
		}
		
		auto sum = 0.0;
		for(auto a = 0u; a < nage; a++){
			sum += data.democat_dist[0][a]*v[a]; 
		}
		
		for(auto a = 0u; a < nage; a++) v[a] /= sum;
		v[nage-1] = 1.6;
	}
	for(auto val : v) cout << val << " "; 
	cout << " age contact factors" << endl;	

	emsg("done");
}


/// Divides compartments in the model based on a transition and returns the comparmtents in a given direction 
vector <unsigned int> Model::trans_comp_divide(unsigned int tr, unsigned int dir) const 
{
	vector <unsigned int> map(comp.size());
	for(auto c = 0u; c < comp.size(); c++) map[c] = UNSET;
	
	map[trans[tr].from] = 0; map[trans[tr].to] = 1;
	
	bool flag;
	do{
		flag = false;
		for(const auto &tra : trans){
			auto from = map[tra.from], to = map[tra.to];
			auto num = 0u;
			if(from == UNSET) num++; if(to == UNSET) num++;
			switch(num){
			case 0:
				if(from != to && !(tra.from == trans[tr].from && tra.to == trans[tr].to)){
					emsgroot("Under this inference scheme transitions cannot merge into the same comparment");
				}
				break;
				
			case 1:
				if(from == UNSET) map[tra.from] = to;
				else map[tra.to] = from;
				flag = true;
				break;
			}
		}
	}while(flag == true);
	
	if(false){
		cout << trans[tr].name << ":   ";
		for(auto c = 0u; c < comp.size(); c++) cout << comp[c].name << ":" << map[c] << ", ";
		cout << "\n"; 
	}
	
	vector <unsigned int> clist;
	for(auto c = 0u; c < comp.size(); c++){
		if(map[c] == dir) clist.push_back(c);
	}
	
	return clist;
}

/// Checks that transition have been correctly specified;
void Model::compartmental_model_check() const
{
	for(const auto &co : comp){                                      // Checks if prob_dep is specified for all branches
		for(auto tr : co.trans){
			for(auto tr2 : co.trans){
				auto tr_dep = trans[tr].prob_dep;
				auto tr_dep2 = trans[tr2].prob_dep;
				if(tr_dep != tr_dep2){
					string em = "In 'trans' the transitions '"+trans[tr].name+"' and '"+trans[tr2].name+"' must have a consistent dependency. ";
					em += "Here '"+trans[tr].name+"' ";
					if(tr_dep == "") em += "has no dependency ";
					else em += "is dependent on '"+tr_dep+"' ";
					
					em += "and '"+trans[tr2].name+"' ";
					if(tr_dep2 == "") em += "has no dependency.";
					else em += "is dependent on '"+tr_dep2+"'.";
		
					emsgroot(em);
				}
			}
		}			
	}	
	
	for(auto tr = 0u; tr < trans.size(); tr++){                      // Checks branching probabilities are correctly specified  
		if(comp[trans[tr].from].trans.size() > 1){
			if(trans[tr].param_prob.size() == 0){
				emsgroot("In 'trans' transition '"+trans[tr].name+"' a value for 'prob' must be set.");
			}
		}
		else{
			if(trans[tr].param_prob.size() != 0){
				emsgroot("In 'trans' transition '"+trans[tr].name+"' a value for 'prob' should not be set.");
			}
		}
	}
	
	for(const auto &co : comp){                                      // Checks distributions are correctly identified
		if(co.trans.size() > 0){
			if(co.param_mean.size() == 0){
				if(details.siminf == SIMULATE){
					emsgroot("In 'comps' for compartment '"+co.name+"' a value for 'mean_value' must be set.");
				}
				else{
					emsgroot("In 'comps' for compartment '"+co.name+"' a value for 'mean_value' or 'mean_prior' must be set.");
				}
			}
		}
		else{
			if(co.param_mean.size() != 0){
				emsgroot("In 'comps' for compartment '"+co.name+"' a value for 'mean_value' or 'mean_prior' should not be set.");
			}
			if(co.num != UNSET){
				emsgroot("In 'comps' for compartment '"+co.name+"' a value for 'dist' should not be set.");
			}
		}
	}
	
	auto c = 0u;                                                     // Checks infectious states exist
	while(c < comp.size() && param[comp[c].infectivity_param].name == "zero") c++;
	if(c == comp.size()){
		emsgroot("At least one compartment must have non-zero infectivity");
	}
	
	vector <bool> flag(comp.size());                                 // Checks compartments are connected
	for(auto co = 0u; co < comp.size(); co++) flag[co] = false;
	for(auto tr : trans){ flag[tr.from] = true; flag[tr.to] = true;}
	for(auto co = 0u; co < comp.size(); co++){
		if(flag[co] == false){
			emsgroot("The compartment '"+comp[co].name+"' is not connected to the other compartments by a transition.");
		}	
	}
}
