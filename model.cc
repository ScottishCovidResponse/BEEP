// This file defines the compartmental model and is used to simulate compartmental dynamics after an individual becomes infected

#include <string>
#include <sstream>
#include <iostream>
#include <fstream>
#include <cmath>
 
#include "math.h"

using namespace std;

#include "model.hh"

/// Initialises the model 
Model::Model(const Inputs &inputs, const Details &details, Data &data) : details(details), data(data), inputs(inputs)
{
	maximum_infected = inputs.find_integer("infected_max",UNSET);

	if(modeltype == IND_MODEL){
		if((details.mode == ABC_SMC || details.mode == ABC_MBP || details.mode == ABC_MBP_GR) && maximum_infected == LARGE){
			emsgroot("Input file must contain a limit on the maximum number of individuals through 'maximum_infected'.");
		}
	}		

	load_model();

	complete_datatables();
	
	counterfactual_modification();
	
	if(false) print_parameter_types();
}


/// Loads the model from the input file
void Model::load_model()
{
	add_trans_comps();                                               // Adds compartments and transitions to the model

	add_splines();                                                   // Adds splines to the model

	add_region_effect();                                             // Adds regional effects

	add_probreach();                                                 // Adds output giving e.g. ifr
	
	set_susceptibility_param();

	for(unsigned int c = 0; c < data.ncovar; c++){                   // Adds parameters for area covariates
		covariate_param.push_back(add_parameter(data.covar[c].ps,COVAR_PARAM));
	}
	
	set_infected_uninfected();                                       // Divides model into infected and uninfected compartments
	
	for(auto th = 0u; th < param.size(); th++){                      // Lists all parameters which are not fixed
		if(param[th].priortype != FIXED_PRIOR) param_not_fixed.push_back(th);
	}
}


/// Sets the uninfected and infect compartments in the model.
void Model::set_infected_uninfected()
{
	if(ninfection_trans == 0) emsg("There must be at least one infection transition.");
	 
	inf_state.resize(ninfection_trans);       
	inf_state_ref.resize(ninfection_trans);       
	
	for(auto inft = 0u; inft < ninfection_trans; inft++){
		vector <bool> flag_forward(comp.size());
		for(auto co = 0u; co < comp.size(); co++) flag_forward[co] = false;
		
		vector <unsigned int> list;                           // List all states which link to infection
		auto c = trans[infection_trans[inft]].to;	
		list.push_back(c); flag_forward[c] = true; 
		
		for(auto i = 0u; i < list.size(); i++){
			c = list[i];
			for(auto tr : comp[c].trans){
				if(trans[tr].type != INFECTION_DIST){
					auto cc = trans[tr].to;
					if(flag_forward[cc] == false){
						list.push_back(cc); flag_forward[cc] = true; 
					}
				}
			}
		}
	
		vector <bool> flag_back(comp.size());                   // Lists all states which go to infectious individuals
		for(auto co = 0u; co < comp.size(); co++) flag_back[co] = false;
	
		for(auto co = 0u; co < comp.size(); co++){
			if(param[comp[co].infectivity_param].name != "zero"){
				if(comp[co].infection_transition != inft && flag_forward[co] == true){
					emsg("The infectious state '"+comp[co].name+"' derives from the wrong infection transition."); 
				}
		
				if(comp[co].infection_transition == inft){
					if(flag_forward[co] != true) emsg("The  infectious state '"+comp[co].name+"' derives from the wrong infection transition."); 
					
					vector <unsigned int> list;
					list.push_back(co); flag_back[co] = true;
					
					for(auto i = 0u; i < list.size(); i++){
						auto c = list[i]; 
						for(auto tr = 0u; tr < trans.size(); tr++){
							auto tra = trans[tr];
							if(tra.type != INFECTION_DIST && tra.to == c){
								auto cc = tra.from;
								if(flag_back[cc] == false){
									list.push_back(cc); flag_back[cc] = true;
								}
							}
						}
					}
				}
			}
		}
	
		inf_state_ref[inft].resize(comp.size());
		for(auto co = 0u; co < comp.size(); co++){
			if(flag_forward[co] == true && flag_back[co] == true){
				inf_state_ref[inft][co] = inf_state[inft].size();
				inf_state[inft].push_back(co);
			}
			else inf_state_ref[inft][co] = UNSET;
		}
	
		if(false){
			cout << "\nINFECT " << inft << ": \n";
			for(auto co : inf_state[inft]){
				cout << comp[co].name  << " infectious state\n";
			}	
	
			for(auto co = 0u; co < comp.size(); co++){
				cout << comp[co].name  << " " << inf_state_ref[inft][co] << "  Forward:" << flag_forward[co] << "   Back:" << flag_back[co] << " ref\n";
			}
		}
	}
}

	
/// Sets the susceptibility parameters
void Model::set_susceptibility_param()
{
	sus_param_list.resize(data.ndemocat);                       // Adds parameters for susceptibility difference
	sus_param_popfrac.resize(data.ndemocat); 
	sus_ref.resize(data.ndemocat); 	
	for(auto c = 0u; c < data.ndemocat; c++){
		if(data.democat[c].sus_vari == true){		
			for(auto fi = 0u; fi < data.democat[c].value.size(); fi++){	
				auto th = add_parameter(data.democat[c].ps[fi],SUSCEPTIBILITY_PARAM);
			
				auto j = 0u; while(j < sus_param_list[c].size() && th != sus_param_list[c][j]) j++;
				if(j == sus_param_list[c].size()){
					sus_param_list[c].push_back(th);
					sus_param_popfrac[c].push_back(data.democat_dist[c][fi]);
				}
				else sus_param_popfrac[c][j] += data.democat_dist[c][fi];
				
				sus_ref[c].push_back(j);
			}
			
			switch(details.siminf){                                           // Accounts for fractions of population in different groups
				case SIMULATE:
					for(auto j = 0u; j < sus_param_list[c].size(); j++){
						param[sus_param_list[c][j]].val1 *= sus_param_popfrac[c][j];
					}
					break;
					
				default: break;
			}
		}
	}
}

	
/// Adds transitions and comaprtments to the model
void Model::add_trans_comps()
{
	vector <string> name;
	vector <ParamSpec> inf;
	vector <string> inf_trans;
	inputs.find_compartments(name,inf,inf_trans);        // Adds compartments
	for(auto c = 0u; c < name.size(); c++) add_compartment(name[c],inf[c],inf_trans[c]);	

	vector <string> from, to;                            // Adds transitions
	vector < vector <ParamSpec> > distspec, probspec;

	vector <unsigned int>	shape;
	vector <int> type;
	
	inputs.find_transitions(from,to,type,distspec,probspec,shape,data.ndemocatpos);
	for(auto tr = 0u; tr < from.size(); tr++){
		if(type[tr] == ERLANG_DIST){                       // Adds in the intermediary steps for the Erlang distribution
			unsigned int k = shape[tr];
			if(k == 1){
				add_transition(from[tr],from[tr],to[tr],distspec[tr],probspec[tr],type[tr],1);
			}
			else{
				auto c = 0u; auto cmax = comp.size(); while(c < cmax && from[tr] != comp[c].name) c++;
				if(c == cmax) emsg("Cannot find the 'from' compartment '"+from[tr]+"' for the transition");
				ParamSpec inf = param[comp[c].infectivity_param].ps;
				string inf_trans = comp[c].inf_trans; 
				
				vector <string> inter;
				for(auto i = 1u; i < k; i++){                  // First generates intermediary compartments				
					stringstream ss; ss << from[tr] << "_" << i;
					inter.push_back(ss.str());
					
					add_compartment(ss.str(),inf,inf_trans);
				}
				
				add_transition(from[tr],from[tr],inter[0],distspec[tr],probspec[tr],EXP_DIST,1.0/k);
				vector <ParamSpec> empty;
				for(auto i = 0u; i < k-2; k++) add_transition(from[tr],inter[i],inter[i+1],empty,empty,EXP_DIST,1.0/k);
				add_transition(from[tr],inter[k-2],to[tr],distspec[tr],probspec[tr],EXP_DIST,1.0/k);			
			}
		}
		else{
			add_transition(from[tr],from[tr],to[tr],distspec[tr],probspec[tr],type[tr],1);
		}
	}
	
	for(auto tr = 0u; tr < trans.size(); tr++){          // Checks all branckhing probabilities are correctly specified  
		if(comp[trans[tr].from].trans.size() > 1){
			if(trans[tr].probparam.size() == 0) emsg("Transition "+trans_str(tr)+" must have 'prob' specified.");
		}
		else{
			if(trans[tr].probparam.size() != 0) emsg("Transition "+trans_str(tr)+" should not have 'prob' specified.");
		}
	}
	
	for(auto tr = 0u; tr < trans.size(); tr++){                                // Calculates all the infection transitions
		if(trans[tr].type == INFECTION_DIST) infection_trans.push_back(tr);
	}
	ninfection_trans = infection_trans.size();
	
	start_compartment = UNSET;
	for(auto tr : infection_trans){
		if(start_compartment == UNSET){
			start_compartment = trans[tr].from;
		}
		else{
			if(start_compartment != trans[tr].from) emsg("More than one potential start compartments"); 
		}
	}
	
	for(auto &co : comp){                                             // Works out where the infectivity is applied
		auto th = co.infectivity_param;
		if(param[th].name != "zero"){
			if(co.inf_trans == ""){
				if(infection_trans.size() == 1){
					co.infection_transition = 0;
				}
				else{
					emsg("For the compartment '"+co.name+"' the user must use 'inf_trans' to disambiguate with transition the infectivty applies.");
				}
			}
			else{
				auto tr = get_transition(co.inf_trans);
				auto i = find_in(infection_trans,tr);
				if(i == UNSET) emsg("The transition '"+co.inf_trans+"' could not be found");
				
				co.infection_transition = i;
			}
		}
		else{
			if(co.inf_trans != ""){
				emsg("inf_trans='"+co.inf_trans+"' should not be set."); 
			}				
		}
	}
	
	vector <bool> flag(comp.size());                                  // Checks compartments are connected
	for(auto co = 0u; co < comp.size(); co++) flag[co] = false;
	for(auto tr : trans){ flag[tr.from] = true; flag[tr.to] = true;}
	for(auto co = 0u; co < comp.size(); co++){
		if(flag[co] == false) emsg("Compartment '"+comp[co].name+"' not connected by a transition.");
	}
}
		

/// Add parameters for regional effects
void Model::add_region_effect()
{		
	region_effect = 0;
	sigma_param = UNSET;
	if(data.region_effect.size() > 0){
		region_effect = 1;
		add_parameter(data.sigma,SIGMA_PARAM);
		
		region_effect_param.resize(data.region_effect.size());
		for(auto r = 0u; r < data.region_effect.size(); r++){
			region_effect_param[r] = param.size();
			stringstream ss; ss << "reff_" << data.region_effect[r].region;
			
			ParamSpec ps;
			ps.name = ss.str(); ps.value = ""; ps.prior="Uniform(-1,1)";
			
			add_parameter(data.sigma,SIGMA_PARAM);
		}
	}
	
	if(region_effect != 0){                                              // Calculates the shift in region effect 
		region_effect_shift.resize(data.ncovar);
		for(auto j = 0u; j < data.ncovar; j++){    
			region_effect_shift[j].resize(data.region_effect.size());
			for(auto r = 0u; r < data.region_effect.size(); r++){
				auto sum = 0.0, sum2 = 0.0; 
				for(auto c : data.region_effect[r].area){ 
					sum += data.area[c].covar[j]; sum2++;
				}
				region_effect_shift[j][r] = -sum/sum2;
			}
		}
	}	
}


/// Add outputs which give the probability of reaching a certain compartment
void Model::add_probreach()
{
	vector <string> name, comp, inf_trans;
	
	inputs.find_probreach(name,comp,inf_trans);

	for(auto i = 0u; i < name.size(); i++){
		ProbReach pr;
		pr.name = name[i];
		pr.comp = get_compartment(comp[i]);

		auto tr = get_transition(inf_trans[i]);
		pr.inft = find_in(infection_trans,tr);
		if(pr.inft == UNSET) emsg("The transition '"+inf_trans[i]+"' could not be found");
		
		prob_reach.push_back(pr);
	}
}


/// Updates the references to the spline
void Model::update_spline_ref(string &name, string &desc, string inft_str, string area_str, const vector <double> efoi_ad, vector < vector <unsigned int> >& spline_ref)
{
	auto inftmin = 0u, inftmax = ninfection_trans, areamin = 0u, areamax = data.narea;
		
	if(inft_str != ""){
		auto tr = get_transition(inft_str);
		inftmin = find_in(infection_trans,tr);
		if(inftmin == UNSET) emsg("The transition '"+inft_str+"' is not an infection transition.");
		inftmax = inftmin+1;

		if(ninfection_trans > 1){
			name += " "+inft_str;
			desc += " for transition "+inft_str;
		}
	}
		
	if(area_str != ""){
		auto c = 0u; while(c < data.narea && area_str != data.area[c].code) c++;
		if(c == data.narea) emsg("Cannot find area '"+area_str+"'");
		
		areamin = c; areamax = areamin+1;

		if(data.narea > 1){
			name += " "+area_str;
			desc += " for area "+area_str;
		}
	}

	for(auto inft = inftmin; inft < inftmax; inft++){
		for(auto c = areamin; c < areamax; c++){
			spline_ref[inft][c] = spline.size();
			if(efoi_ad.size() != 0) efoi_agedist[inft][c] = efoi_ad;
		}
	}	
}
	
	
/// Adds splines to the model
void Model::add_splines()
{	
	vector < vector <ParamSpec> > ps_vec; 
	vector <ParamSpec> factor_param_vec; 
	vector < vector <unsigned int> > bp_vec;  
	vector <Smooth> smooth_type_vec;
	vector <string> inft;
	vector <string> area;
	vector < vector <double> > efoi_agedist_vec;
	vector <double> empty;
	
	R_spline_ref.resize(ninfection_trans); efoi_spline_ref.resize(ninfection_trans); efoi_agedist.resize(ninfection_trans);
	for(auto inft = 0u; inft < ninfection_trans; inft++){
		R_spline_ref[inft].resize(data.narea); efoi_spline_ref[inft].resize(data.narea); efoi_agedist[inft].resize(data.narea);
		for(auto c = 0u; c < data.narea; c++){
			R_spline_ref[inft][c] = UNSET; efoi_spline_ref[inft][c] = UNSET; 
		}
	}

	inputs.find_spline("R_spline",ps_vec,factor_param_vec,bp_vec,smooth_type_vec,inft,area,efoi_agedist_vec,empty,details);  // Sets up splines for R(t)
	
	if(ps_vec.size() == 0) emsg("The quantity 'R_spline' must be specified.");
	
	for(auto num = 0u; num < ps_vec.size(); num++){
		string name = "R(t)", desc = "The reproduction number";
		if(efoi_agedist_vec[num].size() != 0) emsg("'age_dist' cannot be applied to R_spline");

		update_spline_ref(name,desc,inft[num],area[num],empty,R_spline_ref);
		create_spline(name,desc,ps_vec[num],bp_vec[num],factor_param_vec[num],smooth_type_vec[num],R_EFOI_PARAM);
	}
	
	auto efoi_agedist = data.agedist;	
	inputs.find_spline("efoi_spline",ps_vec,factor_param_vec,bp_vec,smooth_type_vec,inft,area,efoi_agedist_vec,data.agedist,details);  // Sets up efoi(t)
	
	if(ps_vec.size() == 0) emsg("The quantity 'efoi_spline' must be specified.");
	
	for(auto num = 0u; num < ps_vec.size(); num++){
		string name = "efoi(t)", desc = "External force of infection";
		update_spline_ref(name,desc,inft[num],area[num],efoi_agedist_vec[num],efoi_spline_ref);
		create_spline(name,desc,ps_vec[num],bp_vec[num],factor_param_vec[num],smooth_type_vec[num],R_EFOI_PARAM);
	}
	
	for(auto inft = 0u; inft < ninfection_trans; inft++){
		for(auto c = 0u; c < data.narea; c++){
			if(R_spline_ref[inft][c] == UNSET) emsg("'R_spline' must be set for all areas and infection transitions.");
			if(efoi_spline_ref[inft][c] == UNSET) emsg("'efoi_spline' must be set for all areas and infection transitions.");
		}
	}
	
	geo_spline_ref = spline.size();                                       // Sets up spline for geograpic mixing
	inputs.find_spline("geo_mixing_modify",ps_vec,factor_param_vec,bp_vec,smooth_type_vec,inft,area,efoi_agedist_vec,empty,details);
	
	if(ps_vec.size() == 0){                                               // If not set then sets up a spline with the value unity
		ParamSpec one; one.name = "one"; one.value = "1"; one.prior = "Fixed(1)";
		
		vector <unsigned int> bp;
		bp.push_back(0); bp.push_back(details.end-details.start);
		bp_vec.push_back(bp);
		
		vector <ParamSpec> ps_list; ps_list.push_back(one); ps_list.push_back(one); 
		ps_vec.push_back(ps_list);
			
		factor_param_vec.push_back(one);
		
		smooth_type_vec.push_back(LOGSMOOTH);
	}
	
	if(ps_vec.size() != 1) emsg("Only one spline can represent 'geo_mixing'.");
	create_spline("d(t)","Relative within / between region mixing.",ps_vec[0],bp_vec[0],factor_param_vec[0],smooth_type_vec[0],R_EFOI_PARAM);
		
	for(auto& mm : data.genQ.matmod){
		switch(mm.type){
		case ROW_COLUMN:
			mm.spline_ref = spline.size();
			ParamSpec one; one.name = "one"; one.value = "1"; one.prior = "Fixed(1)";
			create_spline(mm.name,mm.desc,mm.ps,mm.bp,one,LOGSMOOTH,MODIFY_PARAM);
			break;
		}
	}
						
	set_disc_spline_grad();
}


/// Creates a spline based on a given specification
void Model::create_spline(const string name, const string desc, const vector <ParamSpec> &ps, const vector<unsigned int> &bp, const ParamSpec param_fac, Smooth smooth_type, ParamType type)
{
	if(false){
		cout << name << " name\n";
		for(auto i = 0u; i < bp.size(); i++){	
			cout <<  bp[i] << " " << ps[i].name << " " << ps[i].value << "  pro:" << ps[i].prior << " " << ps[i].smooth << " " << ps[i].factor << " \n";
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
	spli.smooth = smooth_type;
	
	spline.push_back(spli);
}

		
/// Samples from the prior for a given parameter
double Model::param_prior_sample(unsigned int th) const
{
	switch(param[th].priortype){
			case FIXED_PRIOR: 
				return param[th].val1;
				
			case UNIFORM_PRIOR:	
				return param[th].val1 + ran()*(param[th].val2 - param[th].val1);
			
			case EXP_PRIOR:
				return -log(ran())/param[th].val1;
				
			default:			
				emsgEC("model",58);
				break;
		}
}


/// Randomly samples the parameter values from the prior
vector <double> Model::sample_from_prior() const
{
	vector <double> paramv(param.size());
	for(auto th = 0u; th < param.size(); th++){	     // Most parameters have a uniform distribution.
		paramv[th] = param_prior_sample(th); 
	}
	
	if(region_effect == 1){                          // Regional effects are random effects
		for(auto th : region_effect_param){
			if(param[th].priortype != FIXED_PRIOR){
				bool flag;
				do{
					paramv[th] = normal_sample(0,paramv[sigma_param]);
					
					flag = false;
					switch(param[th].priortype){
						case UNIFORM_PRIOR: 
							if(paramv[th] < param[th].val1 || paramv[th] > param[th].val2) flag = true;
							break;
				
						default:
							emsgEC("model",57);
							break;
					}
				}while(flag == true);
			}
		}
	}
	
	spline_sample(paramv);                          // Smoothing on the spline is accounted for 

	return paramv;
}


/// Samples from the smoothing spline
void Model::spline_sample(vector <double> &paramv) const 
{
	for(const auto &spl : spline){  
		if(spl.smooth != NOSMOOTH){
			auto spli = spl.p;
			bool flag;
			auto loop = 0u;
			do{
				vector <bool> param_set(param.size());
				for(auto th = 0u; th < param.size(); th++) param_set[th] = false;
			
				flag = true;
				auto th1 = spli[0].param;

				if(param[th1].priortype != FIXED_PRIOR && param[th1].priortype != UNIFORM_PRIOR) emsgEC("model",60);
					
				paramv[th1] = param_prior_sample(th1);
				
				param_set[th1] = true;
				
				for(auto i = 0u; i < spli.size()-1; i++){
					auto th1 = spli[i].param;
					auto th2 = spli[i+1].param;
					
					if(th1 != th2 && param[th2].priortype != FIXED_PRIOR){
						auto min = param[th2].val1, max = param[th2].val2; 
						
						if(spli[i].t == spli[i+1].t){
							paramv[th2] = param_prior_sample(th2);
						}
						else{
							switch(spl.smooth){
							case LOGSMOOTH:
								{
									auto sm = spli[i].smooth_value; if(sm == UNSET) emsg("The smoothed value in the spline shoud be set");
									paramv[th2] = (spli[i+1].multfac/spli[i].multfac)*paramv[th1]*exp(normal_sample(0,sm)); 
								}
								break;
								
							case SMOOTH:
								{
									auto sm = spli[i].smooth_value; if(sm == UNSET) emsg("The smoothed value in the spline shoud be set");
							
									if(spli[i+1].multfac != spli[i].multfac) emsgEC("Model",56);
									paramv[th2] = paramv[th1] + normal_sample(0,sm); 
								}
								break;
								
							default: emsgEC("Model",43);
							}
						
							if(paramv[th2] <= min || paramv[th2] >= max){ flag = false; break;}
							
							for(auto th : param[th2].greaterthan){
								if(param_set[th] == true && paramv[th2] < paramv[th]) flag = false; 
							}
							
							for(auto th : param[th2].lessthan){
								if(param_set[th] == true && paramv[th2] > paramv[th]) flag = false; 
							}
							
							param_set[th2] = true;
						}
					}
				}
			
				loop++;
			}while(flag == false && loop < spline_sample_try);	
	
			if(loop == spline_sample_try) emsg("Spline smoothing not possible");
		}
	}
}


/// Gets the infectivity difference when going down a transition
double Model::get_infectivity_dif(unsigned int inft, unsigned int tr, const vector <double> paramv) const
{
	double di = 0;
	auto &compto = comp[trans[tr].to];
	auto &compfrom = comp[trans[tr].from];
	if(compto.infection_transition == inft) di += paramv[compto.infectivity_param];
	if(compfrom.infection_transition == inft) di -= paramv[compfrom.infectivity_param];
	return di;
	//return paramv[comp[trans[tr].to].infectivity_param] - paramv[comp[trans[tr].from].infectivity_param];
}


/// Adds a compartment to the model
void Model::add_compartment(const string& name, const ParamSpec &inf, const string &inf_trans)
{
	Compartment co;
	co.name = name;
	co.infectivity_param = add_parameter(inf,INF_PARAM); 
	co.inf_trans = inf_trans;
	co.infection_transition = UNSET;
	
	comp.push_back(co);	
}


/// Adds a parameter to the model
unsigned int Model::add_parameter(const ParamSpec ps, ParamType type)
{
	PriorType pt = FIXED_PRIOR;
	double val1=UNSET, val2=UNSET, value=UNSET;
	
	if(ps.name =="") emsg("The parameter must have a name");
	
	switch(details.siminf){
	case SIMULATE:
		if(ps.value == "") emsg("The simulation value should be set for parameter '"+ps.name+"'.");
		pt = FIXED_PRIOR;
		val1 = atof(ps.value.c_str());
		value = val1;
		break;
		
	case INFERENCE:
		if(ps.value != "") value = atof(ps.value.c_str());
		
		if(ps.prior == "") emsg("The prior should be set for parameter '"+ps.name+"'.");
		
		auto spl = split(ps.prior,'(');
		if(spl.size() != 2) emsg("There was a problem with the expression '"+ps.prior+"'");
		auto spl2 = split(spl[1],')');
		
		if(spl2.size() > 2) emsg("There was a problem with the expression '"+ps.prior+"'");
		auto vals = split(spl2[0],',');
		
		if(spl[0] == "Fixed"){
			pt = FIXED_PRIOR;
			if(vals.size() != 1) emsg("There was a problem with the expression '"+ps.prior+"'"); 
			val1 = atof(vals[0].c_str());
		}
		else{
			if(spl[0] == "Uniform"){
				pt = UNIFORM_PRIOR;
				if(vals.size() != 2) emsg("There was a problem with the expression '"+ps.prior+"'"); 
				val1 = atof(vals[0].c_str());
				val2 = atof(vals[1].c_str());
			}
			else{
				if(spl[0] == "Exp"){
					pt = EXP_PRIOR;
					if(vals.size() != 1) emsg("There was a problem with the expression '"+ps.prior+"'"); 
					val1 = atof(vals[0].c_str());
				}
			}
		}
		break;
	}
	
	Param par;
	par.name = ps.name; par.priortype = pt, par.val1 = val1; par.val2 = val2; par.value = value;
	par.type = type;
	par.ps = ps;

	auto th = 0u; while(th < param.size() && param[th].name != par.name) th++;
	
	if(th == param.size()) param.push_back(par);
	else{
		if(param[th].priortype != par.priortype || param[th].val1 != par.val1 || param[th].val2 != par.val2 
     || param[th].value != par.value  || param[th].type != par.type) emsg("Parameter '"+ps.name+"' has multiple definitions.");	
	}
	
	return th;
}



/// Adds a transition to the model

void Model::add_transition(const string& from_origin, const string& from, const string& to, const vector <ParamSpec>& distspec, const vector <ParamSpec>& probspec, unsigned int type, double mean_factor)
{
	Transition tr;
	
	tr.from_origin = get_compartment(from_origin);
	if(tr.from_origin == UNSET) emsg("Cannot find the 'from' compartment '"+from+"' for the transition");

	tr.from = get_compartment(from);
	if(tr.from == UNSET) emsg("Cannot find the 'from' compartment '"+from+"' for the transition");

	if(from == to) emsg("'from' cannot be the same as 'to'");
		
	if(type == INFECTION_DIST) tr.comptrans_ref = UNSET;
	else{
		tr.comptrans_ref = comp[tr.from].trans.size(); comp[tr.from].trans.push_back(trans.size());
	}
	
	tr.to = get_compartment(to);
	if(tr.to == UNSET) emsg("Cannot find the 'to' compartment '"+to+"' for the transition");	
	
	if(probspec.size() > 0){       // Adds pranching probabilties
		if(probspec.size() != data.ndemocatpos) emsgEC("Model",48);
		
		for(auto dp = 0u; dp < data.ndemocatpos; dp++){
			tr.probparam.push_back(add_parameter(probspec[dp],BRANCHPROB_PARAM));
		}
	}
	
	tr.type = type;
	tr.mean_factor = mean_factor;
	
	if(distspec.size() > 0){
		if(distspec.size() != data.ndemocatpos) emsgEC("Model",48);
		
		for(auto dp = 0u; dp < data.ndemocatpos; dp++){
			tr.param_mean.push_back(add_parameter(distspec[dp],DISTVAL_PARAM));
		}
	}	
	
	
	if(modeltype == POP_MODEL && (type == GAMMA_DIST || type == LOGNORM_DIST)){
		emsg("Only exponential transitions can be used with a population model");
	}
	
	trans.push_back(tr);
}


/// Creates a discrete spline for a given reference ref
vector <double> Model::create_disc_spline(unsigned int ref, const vector<double> &paramv) const
{
	vector <double> disc_spline(details.ndivision);
	
	auto factor = paramv[spline[ref].param_factor];
	auto spl = spline[ref].p;
	
	auto p = 0;
	for(auto s = 0u; s < details.ndivision; s++){	
		auto t = double((s+0.5)*details.period)/details.ndivision;
		
		while(p < int(spl.size())-1 && t > spl[p+1].t) p++;
		
		auto fac = (t-spl[p].t)/(spl[p+1].t-spl[p].t);
		disc_spline[s] = factor*(paramv[spl[p].param]*spl[p].multfac*(1-fac) + paramv[spl[p+1].param]*spl[p+1].multfac*fac);
	}
	
	return disc_spline;
}


/// Sets gradient of discrete spline w.r.t. a variable
void Model::set_disc_spline_grad()
{
	disc_spline_grad.resize(spline.size()); // The gradient in a spline w.r.t. a variable
	
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
	
	for(auto sp = 0u; sp < spline.size(); sp++){   // Works out which splines 
		auto spl = spline[sp];
		if(spl.smooth != NOSMOOTH){			
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
			cout << "\n";
		}
		emsg("P");
	}
}


/// Sets branching transition probabilies based on the parameters
vector <CompTransProb> Model::create_comptransprob(const vector<double> &paramv) const
{
	vector <CompTransProb> comptransprob;
	
	comptransprob.resize(comp.size());
	for(auto c = 0u; c < comp.size(); c++){
		auto kmax = comp[c].trans.size();
		if(kmax > 1){
			comptransprob[c].prob.resize(data.ndemocatpos);
			comptransprob[c].probsum.resize(data.ndemocatpos);
			for(auto dp = 0u; dp < data.ndemocatpos; dp++){
				comptransprob[c].prob[dp].resize(kmax);
				comptransprob[c].probsum[dp].resize(kmax);
				
				auto sum = 0.0;
				for(auto k = 0u; k < kmax; k++){
					auto p = trans[comp[c].trans[k]].probparam[dp]; if(p == UNSET) emsgEC("model",1);	
					if(paramv[p] < 0)  emsgEC("model",2);	
					sum += paramv[p];
				}
				
				for(auto k = 0u; k < kmax; k++){
					auto p = trans[comp[c].trans[k]].probparam[dp]; 
					comptransprob[c].prob[dp][k] = paramv[p]/sum;
				}
				
				sum = 0.0;
				for(auto k = 0u; k < kmax; k++){
					sum += comptransprob[c].prob[dp][k];
					comptransprob[c].probsum[dp][k] = sum;
				}
			}
		}
	}
	
	return comptransprob; 
}


/// In the case of a popultion model, this calculates the transition rates
vector < vector <double> > Model::create_transrate(vector <CompTransProb> &comptransprob, const vector<double> &paramv) const
{
	vector< vector <double> > transrate;
	transrate.resize(trans.size());
	for(auto c = 0u; c < comp.size(); c++){
		auto kmax = comp[c].trans.size();
		for(auto k = 0u; k < kmax; k++){
			auto tr = comp[c].trans[k];
		
			transrate[tr].resize(data.ndemocatpos);
			if(trans[tr].type != INFECTION_DIST){	
				if(trans[tr].type != EXP_DIST) emsgEC("Model",64);
				
				for(auto dp = 0u; dp < data.ndemocatpos; dp++){
					auto mean = paramv[trans[tr].param_mean[dp]]*trans[tr].mean_factor; if(mean == 0) emsgEC("Model",65);
					auto rate = 1.0/mean;
			
					if(kmax == 1) transrate[tr][dp] = rate;
					else transrate[tr][dp] = rate*comptransprob[c].prob[dp][k];
				}
			}
		}
	}
	
	for(auto tr : infection_trans){   // These are setup to be filled in later
		transrate[tr].resize(data.ndemocatpos);
	}

	return transrate;
}


/// Defines the relative susceptibility of individuals
vector <double> Model::create_susceptibility(const vector<double> &paramv) const 
{
	vector <double> susceptibility(data.ndemocatpos);

	vector < vector <double> > sus;
	sus.resize(data.ndemocat);
	for(auto c = 0u; c < data.ndemocat; c++){                          // Checks the susceptibility is correctly averaged
		auto nval = data.democat[c].value.size(); 
		sus[c].resize(nval);
			
		if(data.democat[c].sus_vari == true){		
			auto npa = sus_param_list[c].size();
			
			auto sum = 0.0; for(auto j = 0u; j < npa; j++) sum += paramv[sus_param_list[c][j]];
			double store[npa];
			for(auto j = 0u; j < npa; j++) store[j] = paramv[sus_param_list[c][j]]/(sum*sus_param_popfrac[c][j]);
			
			for(auto j = 0u; j < nval; j++) sus[c][j] = store[sus_ref[c][j]];
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
		for(auto c = 0u; c < data.ndemocat; c++){                          // Checks the susceptibility is correctly averaged
			auto sum = 0.0;
		
			auto nval = data.democat[c].value.size(); 
	
			for(auto j = 0u; j < nval; j++) sum += sus[c][j]*data.democat_dist[c][j];
			cout << sum << "sum check\n";
		}
	}
	
	return susceptibility;
}


/// Defines the relative transmission rate for different areas
vector <double> Model::create_areafactor(const vector<double> &paramv) const
{
	vector <double> areafactor(data.narea);
	for(auto c = 0u; c < data.narea; c++){
		auto sum = 0.0;
		for(auto j = 0u; j < data.ncovar; j++) sum += paramv[covariate_param[j]]*data.area[c].covar[j];
		
		if(region_effect != 0) sum += paramv[region_effect_param[data.area[c].region_effect]];
		
		areafactor[c] = exp(sum);
		if(std::isnan(areafactor[c])) emsgEC("Model",90);
	}
	
	return areafactor;
}


/// Defines the time variation in the age mixing matrix
void Model::create_Ntime(vector < vector < vector <double> > > &Ntime, const vector < vector<double> > &disc_spline) const
{	
	for(auto sett = 0u; sett < details.ndivision; sett++){
		auto& N = Ntime[sett];
			if(data.genQ.N.size() == 0) emsg("l");
		N = data.genQ.N[0].ele;
		
		for(auto& mm : data.genQ.matmod){                         // These modify the basic contact matrix
			switch(mm.type){
				case ROW_COLUMN:	
					auto si = N.size();
					vector <bool> flag(si);
					for(auto i = 0u; i < si; i++) flag[i] = false;
					
					for(auto a : mm.ages) flag[a] = true;
				
					auto fac = disc_spline[mm.spline_ref][sett];
				
					for(auto a : mm.ages){
						for(auto i = 0u; i < si; i++){
							N[a][i] *= fac;
							if(flag[i] == false)	N[i][a] *= fac; 
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
				cout << "before\n";
			}				
		}
	}	
}


/// Determines if two matrices are equal or not
bool Model::equal(const vector < vector <double> > &Ntime_1, const vector < vector <double> > &Ntime_2) const  
{
	for(auto i = 0u; i < Ntime_1.size(); i++){
		for(auto j = 0u; j < Ntime_1[i].size(); j++){
			if(Ntime_1[i][j] != Ntime_2[i][j]) return false;
		}
	}
	return true;
}
	
	
/// Creates the time variation in beta from the time variation in R
vector < vector <double> > Model::create_beta_R_ratio(const vector <double> &susceptibility, const vector <double> &areafactor, const vector <double> &paramv, vector < vector < vector <double> > > &Ntime, const vector < vector <double> > &transrate) const 
{
	vector < vector <double> > beta_R_ratio;  
	
	auto area_av = calculate_area_av(areafactor);
		
	beta_R_ratio.resize(ninfection_trans);
	 
	for(auto inft = 0u; inft < ninfection_trans; inft++){
		beta_R_ratio[inft].resize(details.ndivision);
		
		auto Vinv = calculate_Vinv(transrate, inft);
	
		auto ratio = 0.0;
		for(auto sett = 0u; sett < details.ndivision; sett++){
			if(sett == 0 || equal(Ntime[sett-1],Ntime[sett]) == false){
				ratio = calculate_R_beta_ratio_using_NGM(paramv,susceptibility,Ntime[sett],Vinv,inft);
			}		
			beta_R_ratio[inft][sett] = 1.0/(area_av*ratio);
		}
	}
	
	return beta_R_ratio;
}


/// Gets the comparment from it's name 
unsigned int Model::get_compartment(string compname)
{
	for(auto c = 0u; c < comp.size(); c++){ if(compname == comp[c].name) return c;}
	return UNSET;
}

 
/// From a string gets the transition number
unsigned int Model::get_transition(string transname)
{
	for(auto tr = 0u; tr < trans.size(); tr++){
		if(transname == comp[trans[tr].from_origin].name + "->" + comp[trans[tr].to].name) return tr;
	}
	emsg("Transition '"+transname+"' not recognised");
}


/// Completes information for datatables based on the model
void Model::complete_datatables()
{
	for(auto d = 0u; d < data.datatable.size(); d++){
		auto& dt = data.datatable[d];
		
		auto sp = split(dt.observation,',');
		
		for(auto& st : sp){
			auto co = 0u; while(co < comp.size() && st != comp[co].name) co++;
			if(co < comp.size()){
				dt.complist.push_back(co);
				if(dt.type != POP) emsg("'"+dt.file+"' should not contain compartmental population data");
			}
			else{
				dt.translist.push_back(get_transition(st));
				if(dt.type != TRANS && dt.type != MARGINAL) emsg("'"+dt.file+"' should not contain compartmental transition data");			
			}
		}
		
		if(dt.fraction_spline != ""){
			emsg("P");
			auto sp = 0u; while(sp < spline.size() && spline[sp].name != dt.fraction_spline) sp++;
			if(sp == spline.size()) emsgroot("Cannot find '"+dt.fraction_spline+"'");
			
			for(auto& ob : data.obs){
				if(ob.datatable == d) ob.fraction_spline_ref = sp;
			}
			
			for(auto& gr : data.graph){
				if(gr.datatable == d) gr.fraction_spline_ref = sp;
			}
		}
	}
}


/// Completes information for counterfactual based on the model
void Model::counterfactual_modification()
{
	countermod.transmean_mult.resize(details.ndivision);
	for(auto sett = 0u; sett < details.ndivision; sett++){					
		countermod.transmean_mult[sett].resize(data.narea);
		for(auto c = 0u; c < data.narea; c++){
			countermod.transmean_mult[sett][c].resize(trans.size());
			for(auto tr = 0u; tr < trans.size(); tr++){
				countermod.transmean_mult[sett][c][tr].resize(data.ndemocatpos);
				for(auto dp = 0u; dp < data.ndemocatpos; dp++){
					countermod.transmean_mult[sett][c][tr][dp] = 1;
				}
			}
		}
	}
	
	countermod.efoi_mult.resize(details.ndivision);
	for(auto sett = 0u; sett < details.ndivision; sett++){					
		countermod.efoi_mult[sett].resize(data.narea);
		for(auto c = 0u; c < data.narea; c++){
			countermod.efoi_mult[sett][c].resize(ninfection_trans);
			for(auto i = 0u; i < ninfection_trans; i++){
				countermod.efoi_mult[sett][c][i] = 1;
			}
		}
	}

	countermod.sett_start = details.ndivision;
	for(const auto &cf : data.counterfact){
		auto start = cf.start*details.division_per_time;
		auto end = cf.end*details.division_per_time;
		auto fac = cf.factor;
		
		if(start < countermod.sett_start) countermod.sett_start = start;
		
		switch(cf.type){
			case CF_PPC: break;
			
			case CF_TRANS_RATE:
				{
					auto tr = get_transition(cf.trans_str);
				
					for(auto sett = start; sett < end; sett++){					
						for(auto c : cf.area){
							for(auto dp : cf.dp_sel){
								countermod.transmean_mult[sett][c][tr][dp] *= fac;
							}
						}
					}
				}
				break;
				
			case CF_EFOI:
				{
					auto i=0u;
					if(cf.trans_str == ""){
						if(ninfection_trans != 1) emsg("In 'counterfactual' the infection transition must be specified");
					}
					else{
						auto tr = get_transition(cf.trans_str);
						i = find_in(infection_trans,tr);
						if(i == UNSET) emsg("The transition '"+cf.trans_str+"' could not be found");
					}
					
					for(auto sett = start; sett < end; sett++){					
						for(auto c : cf.area){
							countermod.efoi_mult[sett][c][i] *= fac;
						}
					}
				}
				break;
		}
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
					auto rate = param[th].val1;
					Pr += log(rate) - rate*val;
				}
				break;
				
			default:
				emsgEC("model",67);
				break;
		}
	}
	
	if(region_effect == 1){
		auto sd = paramv[sigma_param];
		for(auto th : region_effect_param) Pr += normal_probability(paramv[th],0,sd*sd);
	}
	
	Pr += spline_prior(paramv);

	return Pr;
}
	

/// Calculates the smoothing spline prior
double Model::spline_prior(const vector<double> &paramv) const
{
	auto Pr = 0.0;

	for(const auto& spl : spline){
		if(spl.smooth != NOSMOOTH){ 
			auto spli = spl.p;
			for(auto i = 0u; i < spli.size()-1; i++){
				auto par1 = spli[i].param;
				auto par2 = spli[i+1].param;
				
				if(par1 != par2 && spli[i].t != spli[i+1].t){
					auto value = spli[i].smooth_value; if(value == UNSET) emsg("The spline smooth value must be set");
						
					switch(spl.smooth){
						case LOGSMOOTH:       // This is the pdf for the log-normal distribution
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
						default: emsgEC("Model",43);
					}
				}
			}
		}
	}
	
	return Pr;
}


/// Returns the gradient in the prior
double Model::dPr_dth(unsigned int th, const vector<double> &paramv) const 
{
	auto dPr_dth = 0.0;
	
	switch(param[th].priortype){
		case FIXED_PRIOR: case UNIFORM_PRIOR: break;
		
		case EXP_PRIOR:
			dPr_dth = -param[th].val1;
			break;
			
		default:
			emsgEC("model",77);
			break;
	}
	
	switch(param[th].type){
		case SIGMA_PARAM:
			{
				auto sd = paramv[sigma_param];
				for(auto th : region_effect_param){
					dPr_dth += -1.0/sd + paramv[th]*paramv[th]/(sd*sd*sd);
				}
			}
			break;
	
		case RE_PARAM:
			if(region_effect == 1){
				for(auto r = 0u; r < data.region_effect.size(); r++){ 
					if(region_effect_param[r] == th){ 
						auto sd = paramv[sigma_param];
						dPr_dth += -paramv[th]/(sd*sd);
					}
				}
			}
			break;
			
		case R_EFOI_PARAM: case GEOMIX_PARAM: case MODIFY_PARAM:
			for(auto psref : param[th].paramsplineref){
				auto i = psref.i;		 
				auto spl = spline[psref.spline];	
				auto spli = spl.p;
					
				if(spli[i].t != spli[i+1].t){
					auto par1 = spli[i].param;
					auto par2 = spli[i+1].param;
				
					if((par2 == th || par1 == th) && par2 != par1){
						auto value = spli[i].smooth_value; if(value == UNSET) emsg("The spline smooth value must be set");
		
						switch(spl.smooth){
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
							
							default: emsgEC("Model",33); break;
						}
					}
				}
			}
			break;
			
		default: break;
	}
	
	if(false){
		auto par = paramv;
		auto val = prior(par);
		auto d = 0.001;
		par[th] += d;
		auto up = prior(par);
		par[th] -= d;
		cout << param[th].name << " " << (up-val)/d << " " << dPr_dth << " grad\n";
	}
	
	return dPr_dth;
}


/// Based in transition probabilites, this calculates the probability of an passing through a given compartment
vector <CompProb> Model::create_compprob(const vector <CompTransProb> &comptransprob, unsigned int inft) const
{
	vector <CompProb> compprob(comp.size());
	
	for(auto& co : compprob){
		co.value.resize(data.ndemocatpos);
		for(auto& value : co.value) value = 0;
	}

	for(auto dp = 0u; dp < data.ndemocatpos; dp++){
		vector <unsigned int> cst, kst;
		vector <double> probst;
	
		auto prob = 1.0;
		auto c = trans[infection_trans[inft]].to;  
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
			if(comp[c].trans.size() > 1) prob *= comptransprob[c].prob[dp][k];
			c = trans[comp[c].trans[k]].to;		
		}while(1 == 1);
		
		if(cst.size() > 0 || kst.size() > 0 || probst.size() > 0) emsgEC("Model",8);
	}
	
	return compprob;
}


/// Calculates the number of infections caused externally [inft][area]
vector <double> Model::calculate_external_ninf(const vector<double> &paramv) const
{	
	vector < vector<double> > spline_val;
	spline_val.resize(spline.size());
	for(auto sp = 0u; sp < spline.size(); sp++) spline_val[sp] = create_disc_spline(sp,paramv);
	
	auto susceptibility = create_susceptibility(paramv);   
	auto dt = double(details.period)/details.ndivision;
	
	vector <double> ninf;
	ninf.resize(ninfection_trans);	
	for(auto inft = 0u; inft < ninfection_trans; inft++){
		ninf[inft] = 0;
		for(auto sett = 0u; sett < details.ndivision; sett++){
			for(auto c = 0u; c < data.narea; c++){
				auto val = spline_val[efoi_spline_ref[inft][c]][sett];
		
				for(auto dp = 0u; dp < data.ndemocatpos; dp++){	
					auto a = data.democatpos[dp][0];
					ninf[inft] += susceptibility[dp]*val*data.area[c].pop[dp]*efoi_agedist[inft][c][a]*dt;
				}
			}
		}
	}
	
	return ninf;
}


/// Calculates the probability of reaching a specified compartment
vector <double> Model::calculate_probreach(const vector<double> &paramv) const
{
	auto comptransprob = create_comptransprob(paramv);
	
	vector <double> vec;
	for(auto pr = 0u; pr < prob_reach.size(); pr++){
		auto probr = prob_reach[pr];
		auto compprob = create_compprob(comptransprob,probr.inft);
		
		auto sum = 0.0, wsum = 0.0;
		for(auto dp = 0u; dp < data.ndemocatpos; dp++){
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
vector <SplineOutput> Model::get_spline_output(const vector <double> &paramv) const
{
	vector <SplineOutput> splineout;
	
	for(auto sp = 0u; sp < spline.size(); sp++){
		SplineOutput so;
		so.name = spline[sp].name;
		so.desc = spline[sp].desc;
		so.splineval = create_disc_spline(sp,paramv);;
		if(so.name.length() >= 7 && so.name.substr(0,7) == "efoi(t)"){
			for(auto& val : so.splineval) val *= details.efoi_factor;
			so.desc += " (per "+to_string(details.efoi_factor)+")";
		}
		splineout.push_back(so);
	}
	
	if(data.nage > 1){
		auto R_age = calculate_R_age(paramv);
		
		for(auto inft = 0u; inft < ninfection_trans; inft++){ 
			for(auto a = 0u; a < data.nage; a++){
				SplineOutput so;
				so.name = "R(t) age:"+data.democat[0].value[a];
				so.desc = "The reproduction number for age:"+data.democat[0].value[a];
				if(ninfection_trans > 1){
					auto trname = trans_str(infection_trans[inft]);
					so.name += " "+trname;
					so.desc += " for "+trname;
				}
				so.splineval = R_age[inft][a];
				splineout.push_back(so);
			}
		}
	}
	
	return splineout;
}	
	

/// Returns a sample of derived parameters
vector <DerivedParam> Model::get_derived_param(const vector<double> &paramv, const vector <double> &susceptibility, const vector <vector <double> > &A, const vector < vector <double> > &transrate) const
{
	vector <DerivedParam> derpar;
	for(auto inft = 0u; inft < ninfection_trans; inft++){   // Generation times
		DerivedParam dp;
		dp.name = "GT";
		dp.desc = "Generation time";
		dp.filename = "GT"; 
		if(ninfection_trans > 1){
			dp.name += " "+trans_str(infection_trans[inft]);
			dp.desc += " for "+trans_str(infection_trans[inft]);
		}
		dp.value = calculate_generation_time(paramv,susceptibility,A,transrate,inft);
		derpar.push_back(dp);
	}
	
	auto prvec = calculate_probreach(paramv);                // The probability of reaching a given compartment
	for(auto pr = 0u; pr < prob_reach.size(); pr++){
		auto &probr = prob_reach[pr];
		DerivedParam dp;
		dp.name = probr.name;
		dp.desc = "The probability of reaching "+comp[probr.comp].name+" from "+trans_str(infection_trans[probr.inft]);
		dp.filename = probr.name;
		dp.value = prvec[pr];
		derpar.push_back(dp);
	}
	
	auto exf_ninf = calculate_external_ninf(paramv);
	for(auto inft = 0u; inft < ninfection_trans; inft++){   // Generation times
		DerivedParam dp;
		dp.name = "External Infections";
		dp.desc = "Number of external infections";
		if(ninfection_trans > 1){
			dp.name += " "+trans_str(infection_trans[inft]);
			dp.desc += " for "+trans_str(infection_trans[inft]);
		}
		dp.value = exf_ninf[inft];
		derpar.push_back(dp);
	}
	
	return derpar;
}


/// Calculates the generation time
double Model::calculate_generation_time(const vector<double> &paramv, const vector <double> &susceptibility, const vector <vector <double> > &A, const vector < vector <double> > &transrate, const unsigned int inft) const
{
	vector <double> vec(ninfection_trans);
	
	auto Q = data.ndemocatpos;
	auto S = inf_state[inft].size();
	auto N = Q*S;
	
	vector <double> comp_tsi(N);                 // Calculate the compartment time since infection
	for(auto i = 0u; i < N; i++) comp_tsi[i] = UNSET;
	
	bool flag; 
	do{
		flag = false;
		for(auto tr : trans){
			auto from = inf_state_ref[inft][tr.from];   
			auto to = inf_state_ref[inft][tr.to];   
			if(to != UNSET){
				for(auto dp = 0u; dp < data.ndemocatpos; dp++){
					if(tr.type == INFECTION_DIST) comp_tsi[Q*to+dp] = 0;
					else{
						if(from != UNSET){
							auto ito = Q*to+dp, ifrom = Q*from+dp;
							if(comp_tsi[ifrom] != UNSET && comp_tsi[ito] == UNSET){
								comp_tsi[ito] = comp_tsi[ifrom] + paramv[tr.param_mean[dp]]*tr.mean_factor;
								flag = true;
							}
						}
					}
				}
			}
		}
	}while(flag == true);
	
	auto Vinv = calculate_Vinv(transrate,inft);

	auto F = calculate_F(paramv,susceptibility,A,inft);

	auto NGM = calculate_NGM(F,Vinv);
	
	vector <double> eigenvector;
	largest_eigenvalue(NGM,eigenvector);
	
	if(false){
		for(auto th = 0u; th < param.size(); th++){
			cout << param[th].name << " " << paramv[th] << "param\n";
		}
			
		for(auto j = 0u; j < N; j++){ 
			auto s = j/Q, dp = j%Q;
			cout << comp[inf_state[inft][s]].name << " " << dp << " " << comp_tsi[j] << " " << Vinv[j][j] << " \n";
		}
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


/// Works out R as a function of age for each of the infection transitions
vector < vector <vector <double> > > Model::calculate_R_age(const vector <double> &paramv) const
{
	vector < vector < vector <double> > > R_age;
	
	auto comptransprob = create_comptransprob(paramv);
	auto transrate = create_transrate(comptransprob,paramv);
	
	vector < vector <double> > disc_spline;
	disc_spline.resize(spline.size());
	for(auto sp = 0u; sp < spline.size(); sp++) disc_spline[sp] = create_disc_spline(sp,paramv);

	auto susceptibility = create_susceptibility(paramv);    
	auto areafactor = create_areafactor(paramv);  

	vector < vector < vector <double> > > Ntime;  
	Ntime.resize(details.ndivision);
	for(auto sett = 0u; sett < details.ndivision; sett++){
		Ntime[sett].resize(data.nage);
		for(auto a = 0u; a < data.nage; a++) Ntime[sett][a].resize(data.nage);
	}	
	create_Ntime(Ntime,disc_spline);
 
	R_age.resize(ninfection_trans);
	for(auto inft = 0u; inft < ninfection_trans; inft++){
		auto Vinv = calculate_Vinv(transrate,inft);
	
		R_age[inft].resize(data.nage);
		for(auto a = 0u; a < data.nage; a++) R_age[inft][a].resize(details.ndivision);
		
		double ratio;
		vector <double> vec(data.nage);
		for(auto sett = 0u; sett < details.ndivision; sett++){
			if(sett == 0 || equal(Ntime[sett-1],Ntime[sett]) == false){
				auto F = calculate_F(paramv,susceptibility,Ntime[sett],inft);
				auto NGM = calculate_NGM(F,Vinv);

				vector <double> eigenvector;
				ratio = largest_eigenvalue(NGM,eigenvector);
	
				if(NGM.size() != data.nage) emsg("NGM not the right size");
				
				for(auto i = 0u; i < data.nage; i++){
					vec[i] = 0; for(auto j = 0u; j < data.nage; j++) vec[i] += NGM[j][i];
				}
			}
			
			auto Rav = 0.0;  // Gets the values of R by averaging over areas
			for(auto c = 0u; c < data.narea; c++) Rav += disc_spline[R_spline_ref[inft][c]][sett]; 
			Rav /= data.narea;
			
			for(auto a = 0u; a < data.nage; a++) R_age[inft][a][sett] = Rav*vec[a]/ratio;
		}
	}
	
	return R_age;
}


/// Calculate Vinv (used to calculate next generation matrix)
vector < vector <double> > Model::calculate_Vinv(const vector < vector <double> > &transrate, const unsigned int inft) const
{
	auto N = data.ndemocatpos*inf_state[inft].size();
	
	vector < vector <double> > V;
	V.resize(N);
	for(auto j = 0u; j < N; j++){
		V[j].resize(N);
		for(auto i = 0u; i < N; i++) V[j][i] = 0;
	}
	
	for(auto tr = 0u; tr < trans.size(); tr++){
		auto from = inf_state_ref[inft][trans[tr].from], to = inf_state_ref[inft][trans[tr].to];
		if(from != UNSET){
			for(auto dp = 0u; dp < data.ndemocatpos; dp++){
				auto rate = transrate[tr][dp];
				auto j = from*data.ndemocatpos + dp;
				V[j][j] += rate;
				if(to != UNSET){
					auto k = to*data.ndemocatpos + dp;
					V[k][j] -= rate;
				}
			}
		}
	}
	
	if(false){
		cout << inft << " inft\n";
		for(auto j = 0u; j < N; j++){
			for(auto i = 0u; i < N; i++) cout << V[j][i] << ",";
			cout << "V\n";
		}
		emsg("V");
	}
	
	return invert_matrix(V);
}


/// Calcualtes the average of the area factor
double Model::calculate_area_av(const vector <double> &areafactor) const
{
	auto area_av = 0.0; 
	for(auto c = 0u; c < data.narea; c++){
		auto popsum = 0.0; for(auto pop : data.area[c].pop) popsum += pop;
		area_av += areafactor[c]*popsum;
	}
	area_av /= data.popsize;
	
	return area_av;
}


/// Calculates the matrix giving the source of new infections based on individuals in different states
vector < vector <double> > Model::calculate_F(const vector<double> &paramv, const vector <double> &susceptibility, const vector <vector <double> > &A, const unsigned int inft) const
{
	vector < vector <double> > F;
	auto Q = data.ndemocatpos;
	auto S = inf_state[inft].size();
	auto N = Q*S;
	
	F.resize(Q);
	for(auto q = 0u; q < Q; q++){
		F[q].resize(N);
		for(auto i = 0u; i < N; i++){
			F[q][i] = 0;
		}
	}
	
	for(auto s = 0u; s < S; s++){
		auto inf = paramv[comp[inf_state[inft][s]].infectivity_param];
		if(inf != 0){
			for(auto dp = 0u; dp < data.ndemocatpos; dp++){
				auto i = s*data.ndemocatpos + dp;
				auto a = data.democatpos[dp][0];
				for(auto q = 0u; q < Q; q++){
					auto aa = data.democatpos[q][0];
					F[q][i] += data.democatpos_dist[q]*susceptibility[q]*inf*A[aa][a];
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
	
	
/// Estimates the ration between R and beta based on largest eigenvector of next generation matrix
double Model::calculate_R_beta_ratio_using_NGM(const vector<double> &paramv, const vector <double> &susceptibility, const vector <vector <double> > &A, const vector < vector <double> > &Vinv, const unsigned int inft) const
{
	auto F = calculate_F(paramv,susceptibility,A,inft);

	auto NGM = calculate_NGM(F,Vinv);
	
	vector <double> eigenvector;
	auto ratio = largest_eigenvalue(NGM,eigenvector);
	
	if(false){ for(auto val: eigenvector) cout << val << ","; cout << "eig\n";}
	
	if(false){
		auto N = Vinv.size();
		for(auto j = 0u; j < N; j++){
			for(auto i = 0u; i < N; i++) cout << Vinv[j][i] << ",";
			cout << "Vinv\n";
		}
		
		for(auto dp = 0u; dp < data.ndemocatpos; dp++){
			for(auto i = 0u; i < N; i++) cout << F[dp][i] << ",";
			cout << " FF\n";
		}
		
		for(auto q = 0u; q < NGM.size(); q++){
			for(auto qq = 0u; qq < NGM.size(); qq++) cout << NGM[q][qq] << ",";
			cout << " NGM\n";
		}
		
		cout << ratio << "Ratio\n"; 
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
				
			default:
				emsgEC("model",62);
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
		if(parami[th] != paramp[th] && (param[th].type == DISTVAL_PARAM || param[th].type == BRANCHPROB_PARAM)) return true;
	}
	return false;
}


/// prints a transition to the terminal
void Model::print_transition(unsigned int tr) const 
{
	cout << comp[trans[tr].from].name << " -> " <<  comp[trans[tr].to].name;
}


/// Generates a string from a transition number 
string Model::trans_str(unsigned int tr) const
{ 
	return comp[trans[tr].from].name+"->"+comp[trans[tr].to].name;
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
		
		case INFERENCE:
			ss << "Priors:" << endl;
			for(auto p = 0u; p < param.size(); p++){
				ss << "  " << param[p].name << " = ";
				switch(param[p].priortype){
					case FIXED_PRIOR:
						ss << "Fixed(" << param[p].val1 << ")";
						break;
					
					case UNIFORM_PRIOR:
						ss << "Uniform(" << param[p].val1 << " - " << param[p].val2 << ")";
						break;
						
					case EXP_PRIOR:
						ss << "Exp(" << param[p].val1 << ")";
						break;
						
					default:
						emsgEC("model",65);
						break;
				}
				
				for(auto th : param[p].greaterthan) ss << " Greater than " << param[th].name << " ";
				for(auto th : param[p].lessthan) ss << " Less than " << param[th].name << " ";
				
				ss << endl;
			}
			break;
		
		default: 
			emsg("Mode not recognised");
			break;
	}
	ss << endl;
		
	ss << "Compartments:" << endl; 
	for(const auto& co : comp){
		ss << "  " << co.name;
		auto th = co.infectivity_param;
		if(param[th].name != "zero"){
			ss << "  Infectivity: " << param[th].name << " acting on ";
			auto tr = co.infection_transition;
			ss << comp[trans[tr].from].name << "->" << comp[trans[tr].to].name;
		}
		ss << endl; 		
	}
	ss << endl;
	
	ss << "Transitions:" << endl; 
	for(const auto& tr : trans){
		if(tr.from != tr.to){
			ss << "  " << comp[tr.from].name << " → " << comp[tr.to].name;
			
			switch(tr.type){
				case INFECTION_DIST: 
					ss << " Infection";
					break;
				case EXP_DIST:
					ss << " Exponential" << endl;
					ss <<	"    mean = ";
					for(auto i = 0u; i < tr.param_mean.size(); i++){
						auto th = tr.param_mean[i];
						if(i != 0) ss << " | ";
						if(tr.mean_factor != 1) ss << tr.mean_factor << "*";
						ss << param[th].name;
					}
					ss << endl;
					break;
				case GAMMA_DIST:
					ss << " Gamma mean="; for(auto th : tr.param_mean) ss << param[th].name << " "; 
					ss << " cv="; for(auto th : tr.param_cv) ss << param[th].name; 
					break;
				case LOGNORM_DIST:
					ss << " Lognormal mean="; for(auto th : tr.param_mean) ss << param[th].name << " " ; 
					ss << " cv="; for(auto th : tr.param_cv) ss << param[th].name << " "; 
					break;
				default:
					break;
			}
			
			if(tr.probparam.size() > 0){
				ss << "    probability ";
				for(auto j = 0u; j < tr.probparam.size(); j++){
					if(j > 0) ss << " | ";
					ss << param[tr.probparam[j]].name;
				}
				ss << endl;
			}
		}
	}
	ss << endl;
	
	if(prob_reach.size() > 0){
		ss << "Probability of reaching:" << endl; 
		for(auto pr = 0u; pr < prob_reach.size(); pr++){
			auto probr = prob_reach[pr];
			auto tr = infection_trans[probr.inft];
			ss << "  Name = " << probr.name << "   Compartment = " << comp[probr.comp].name << "    Infection = ";
			ss << comp[trans[tr].from].name << "->" <<  comp[trans[tr].to].name << "\n";
		}
	}
	
	return ss.str();
}


/// Prints all the parmaeter types
void Model::print_parameter_types()
{
	for(auto& pa: param){
		cout << pa.name << ": ";
		switch(pa.type){
			case R_EFOI_PARAM: cout << "R / efoi\n"; break;
			case GEOMIX_PARAM: cout << "Geographic mixing\n"; break;
			case MODIFY_PARAM: cout << "Modify contact matrix\n"; break;
			case SIGMA_PARAM: cout << "Sigma\n"; break;
			case COVAR_PARAM: cout << "Covariate\n"; break;
			case DISTVAL_PARAM: cout << "Distribution value\n"; break;
			case BRANCHPROB_PARAM: cout << "Branching probability\n"; break;
			case RE_PARAM: cout << "Regional effect\n"; break;
			case OBS_PARAM: cout << "Observation parameter\n"; break;
			case INF_PARAM: cout << "Infectivity parameter\n"; break;
			case SUSCEPTIBILITY_PARAM: cout << "Susceptibility parameter\n"; break; 
			default:emsg("error"); break;
		}
	}
}

// Gets a vector of ParamSample from the generation 
vector <ParamSample> Generation::get_psamp() const
{
	vector <ParamSample> psamp;
	for(const auto& samp : param_samp){
		ParamSample psa; 
		psa.paramval = samp;
		psamp.push_back(psa);
	}
	return psamp;
}
