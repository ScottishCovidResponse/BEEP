// This includes JSON outputs used in the graphical interface

#include <string>
#include <sstream>
#include <iostream>
#include <fstream>
#include <cmath>
 
#include "math.h"

using namespace std;

#include "model.hh"

/// Outputs a summary of the compartmental model in JSON format
string Model::comparmtental_model_JSON() const 
{
	stringstream vis;
	
	vis << "{\"comp\":[";
	
	auto fl = false;
	for(auto i = 0u; i < comp.size(); i++){
		const auto &co = comp[i];
	
		if(co.num == UNSET || co.num == 0){
			if(fl == true) vis << ",";
			fl = true;
			
			vis << "{\"name\":\"" << co.name;
			vis << "\"";
				
			if(co.shape != UNSET){
				vis << ",\"mean\":[";
				
				for(auto i = 0u; i < co.param_mean.size(); i++){
					auto th = co.param_mean[i];
					if(i != 0) vis << ",";
					vis << "{\"param\":\"&m&^{" << co.name << "}_{" << data.democatpos_name[i] << "}\"";
					vis << ",\"name\":\"";
					vis  << param[th].name;
					vis << "\"}";
				}
				vis << "]";
				
				vis << ",\"shape\":" << co.shape;
			}
			if(co.mean_dep != "") vis << ",\"mean_dep\":\"" << co.mean_dep << "\"";
			auto th = co.infectivity_param;
			if(param[th].name != "zero"){
				vis << ",\"inf\":\"" << param[th].name << "\"";
			}
			vis << "}";
		}
	}
	vis << "]";
	
	auto flag = false;
	
	vis << ",\"trans\":["; 
	for(auto i = 0u; i < trans.size(); i++){
		const auto &tr = trans[i];
		if(comp[tr.from].name != comp[tr.to].name){
			if(flag == true) vis << ",";
			flag = true;
			
			vis << "{\"name\":\"" << tr.name << "\"";
			if(tr.inf == TRANS_INFECTION) vis << ",\"infection\":\"yes\"";
			
			vis << ",\"from\":\"" << comp[tr.from].name << "\"";
			vis << ",\"to\":\"" << comp[tr.to].name << "\"";
			
			if(tr.param_prob.size() > 0){
				vis << ",\"prob\":[";
				for(auto j = 0u; j < tr.param_prob.size(); j++){
					if(j != 0) vis << ",";
					vis << "{\"param\":\"&b&^{" << comp[tr.from].name << "→" << comp[tr.to].name << "}_{" << data.democatpos_name[j] << "}\"";
					vis << ",\"name\":\"";
					vis  << param[tr.param_prob[j]].name;
					vis << "\"}";
				}
				vis << "]";
			}
			
			if(tr.prob_dep != "") vis << ",\"prob_dep\":\"" << tr.prob_dep << "\"";
			vis << "}";
		}
	}
	vis << "]";
		
	vis << "}";

	return vis.str();
}


/// Works out if the area effect has time variation
bool Model::areaeffect_timevary() const
{
	for(const auto &cov : data.covar){ if(cov.timevary == true) return true;}
	if(level_effect_param_list.size() > 0) return true;
	return false;
}


/// This convert a vector of parameters into JSON format
string Model::param_JSON(const vector <unsigned int> &vec) const
{
	stringstream ss;
	bool fl = false;
	ss << ",\"param\":[";
	for(auto th : vec){
		auto name = param[th].name;
		if(name != "zero" && name != "one"){
			if(fl == true) ss << ","; fl = true;
			ss << "\"" << name << "\"";
		}
	}
	if(fl == false) ss << "\"No parameters\"";
	ss << "]";
	
	return ss.str();
}


/// This outputs an expression for the force of infection
string Model::foi_model_JSON() const 
{
	stringstream vis;
		
	vis << "{";
	vis << "\"name\":\"&λ_{";
	
	auto fl = false;
	bool descfl= false;
	
	stringstream ss;
	if(data.nstrain > 1){ 
		if(fl == true) vis << ","; fl = true;
		vis << "s";
		
		if(descfl == true) ss << ","; descfl = true;
		ss << "\"&s& indexes pathogen strain.\"";
	}

	if(data.narea > 1){ 
		if(fl == true) vis << ","; fl = true;
		vis << "a";
		
		if(descfl == true) ss << ","; descfl = true;
		ss << "\"&a& indexes area.\"";
	}
	
	if(data.ndemocatpos_per_strain > 1){ 
		if(fl == true) vis << ","; fl = true;
		vis << "d";
		
		if(descfl == true) ss << ","; descfl = true;
		ss << "\"&d& indexes demographic group.\"";
	}
	
	if(fl == true) vis << ","; fl = true;
	vis << "t";
	
	if(descfl == true) ss << ","; descfl = true;
	ss << "\"&t& indexes time.\"";

	vis << "}&\"";
	
	vis << ",\"eq\":[";

	bool eqfl = false;
	
	bool susvar = false;
	for(auto c = 0u; c < data.ndemocat; c++){
		if(data.democat[c].sus_vari == true) susvar = true;
	}
	
	if(susvar == true){                                   // sigma
		if(descfl == true) ss << ","; descfl = true;
		ss << "\"&σ_d& sets the relative susceptiblity of demographic group &d&.\"";
		
		if(eqfl == true) vis << ","; eqfl = true;
		vis << "{\"text\":\"&σ_d&\"";
		
		vector <unsigned int> sus_param_list;
		for(auto c = 0u; c < data.ndemocat; c++){
			if(data.democat[c].sus_vari == true){		
				for(auto th : data.democat[c].sus_param) add_vec(sus_param_list,th);
			}
		}
		vis << param_JSON(sus_param_list) << "}";
	}
	
	if(susvar == true){
		if(eqfl == true) vis << ","; eqfl = true;
		vis << "{\"text\":\"[\"}";
	}
	
	string r;                                           // Adds r
	if(data.nstrain > 1) r = "&r_{s,t}&";
	else r = "&r_{t}&";
		
	if(descfl == true) ss << ","; descfl = true;
	ss << "\"" << r << " is an automatically calculated propotionality contant.\"";
	
	if(eqfl == true) vis << ","; eqfl = true;
	vis << "{\"text\":\"" << r << "\"}";
	
	string R;                                           // Adds R
	if(Rspline_info.size() > 1) R = "&R_{a,t}&";
	else R = "&R_{t}&";
	
	if(descfl == true) ss << ","; descfl = true;
	ss << "\"" << R << " represents the ";
	if(data.narea > 1) ss << "spatially averaged ";
	ss << "reproduction number.\"";
	
	if(eqfl == true) vis << ","; eqfl = true;
	vis << "{\"text\":\"" << R << "\"";
	
	vector <unsigned int> param_list;
	for(const auto &info : Rspline_info){
		auto sp = info.spline_ref;
		for(const auto &p : spline[sp].p)	add_vec(param_list,p.param);
	}
	
	vis << param_JSON(param_list) << "}";
	
	vector <unsigned int> area_param_list;                           // Area effect
	
	for(auto th : area_effect_param_list) add_vec(area_param_list,th);
	for(auto th : covariate_param) add_vec(area_param_list,th);
	for(auto th : level_effect_param_list) add_vec(area_param_list,th);
	
	if(area_param_list.size() > 0){
		if(descfl == true) ss << ","; descfl = true;
		
		string alpha;
		if(areaeffect_timevary() == false){
			alpha = "&α_a&"; ss << "\"" << alpha << " gives a factor change in &R& for different areas.\"";
		}
		else{
			alpha = "&α_{a,t}&"; ss << "\"" << alpha << " gives a time-varying factor change in &R& for different areas.\"";
		}			
		
		if(eqfl == true) vis << ","; eqfl = true;
		vis << "{\"text\":\"" << alpha << "\"";
	
		vis << param_JSON(area_param_list) << "}";
	}
	
	if(data.nstrain > 1){                              // psi		
		if(descfl == true) ss << ","; descfl = true;
		ss << "\"&Ψ_s& gives a factor change in &R& for different strains.\"";
	
		if(eqfl == true) vis << ","; eqfl = true;
		vis << "{\"text\":\"&Ψ_s&\"";
	
		vector <unsigned int> psi_param_list;
		for(const auto &str : data.strain) add_vec(psi_param_list,str.Rfactor_param);
	
		vis << param_JSON(psi_param_list) << "}";
	}
	
	
	if(eqfl == true) vis << ","; eqfl = true;
	vis << "{\"text\":\"(\"}";
	
	bool sumfl = false;                                        // Sum
	string sum_desc = "sums over ";
	string sum = "Σ&_{";
	if(data.narea > 1){
		if(sumfl == true){ sum += ","; sum_desc += ", ";}; sumfl = true;
		sum += "a′";
		sum_desc += "area &a′&";
	}
	
	if(data.ndemocatpos_per_strain > 1){
		if(sumfl == true){ sum += ","; sum_desc += ", ";}; sumfl = true;
		sum += "d′";
		sum_desc += "demographic group &d′&";
	}
	
	if(sumfl == true){ sum += ","; sum_desc += ", ";}; sumfl = true;
	sum += "c";
	sum_desc += "compartment &c&.";

	sum += "}&";
	
	if(descfl == true) ss << ","; descfl = true;
	ss << "\"" << sum << sum_desc << "\"";
	
	
	if(eqfl == true) vis << ","; eqfl = true;
	vis << "{\"text\":\"" << sum << "\"}";
	
	if(data.narea > 1){                                       // Adds M
		vector <unsigned int> M_param_list;      
	
		auto sp = geo_spline_ref;
		for(const auto &p : spline[sp].p)	add_vec(M_param_list,p.param);

		string M;                         
		if(spline[sp].name != "UNSET") M = "&M_{a,a′,t}&";
		else M = "&M_{a,a′}&";
		
		if(descfl == true) ss << ","; descfl = true;
		ss << "\"" << M << " represents the geographical mixing of individuals.\"";
	
		if(eqfl == true) vis << ","; eqfl = true;
		vis << "{\"text\":\"" << M << "\"" << param_JSON(M_param_list) << "}";
	}
	
	if(data.nage > 1){                                       // Adds A
		vector <unsigned int> A_param_list;              
		for(const auto &mm : data.genQ.matmod){
			auto sp = mm.spline_ref;
			for(const auto &p : spline[sp].p)	add_vec(A_param_list,p.param);
		}
	
		string A;                         
		if(data.genQ.matmod.size() > 0) A = "&A_{d,d′,t}&";
		else A = "&A_{d,d′}&";
		
		if(descfl == true) ss << ","; descfl = true;
		ss << "\"" << A << " represents the demographic mixing of individuals.\"";
	
		if(eqfl == true) vis << ","; eqfl = true;
		vis << "{\"text\":\"" << A << "\"" << param_JSON(A_param_list) << "}";
	}
	
	                        
	vector <unsigned int> i_param_list;                    // Adds i   
	for(const auto &co : comp) add_vec(i_param_list,co.infectivity_param);
	
	if(descfl == true) ss << ","; descfl = true;
	ss << "\"&i^c& gives the relative infectivity of compartments.\"";
	
	if(eqfl == true) vis << ","; eqfl = true;
	vis << "{\"text\":\"&i^c&\"" << param_JSON(i_param_list) << "}";
	
	bool Nfl = false;                                      // Adds N
	string N = "&N^c_{";
	string Ndesc = " gives the give the number of individuals in &c&";
	if(data.narea > 1){
		if(Nfl == true) N += ","; Nfl = true;
		Ndesc += ", ";
		N += "a′"; Ndesc += "&a′&";
	}
	if(data.ndemocatpos_per_strain > 1){
		if(Nfl == true) N += ","; Nfl = true;
		Ndesc += ", ";
		N += "d′"; Ndesc += "&d′&";
	}
	N += "}&"; Ndesc += ".";
	
	if(descfl == true) ss << ","; descfl = true;
	ss << "\"" << N << Ndesc << "\"";
	
	if(eqfl == true) vis << ","; eqfl = true;
	vis << "{\"text\":\"" << N << "\"" << "}";
	
	if(eqfl == true) vis << ","; eqfl = true;
	vis << "{\"text\":\")\"}";
		
	if(eqfl == true) vis << ","; eqfl = true;
	vis << "{\"text\":\"+\"}";
	
	if(eqfl == true) vis << ","; eqfl = true;       // Adds EFOI
	vis << "{\"text\":\"(\"}";
	
	string efoi;                                         
	if(efoispline_info.size() > 1){
		if(data.nstrain > 1) efoi = "&η_{s,a,d,t}&";
		else efoi = "&η_{a,d,t}&";
	}
	else{
		if(data.nstrain > 1) efoi = "&η_{s,t}&";
		else efoi = "&η_{t}&";
	}
	
	if(descfl == true) ss << ","; descfl = true;
	ss << "\"" << efoi << " represents the external force of infection (FOI).\"";
	
	if(eqfl == true) vis << ","; eqfl = true;
	vis << "{\"text\":\"" << efoi << "\"";
	
	vector <unsigned int> efoi_param_list;
	for(const auto &info : efoispline_info){
		auto sp = info.spline_ref;
		for(const auto &p : spline[sp].p){
			add_vec(efoi_param_list,p.param);
		}
	}
	
	vis << param_JSON(efoi_param_list) << "}";
	
	if(eqfl == true) vis << ","; eqfl = true;               // f
	vis << "{\"text\":\"/\"}";
	
	if(descfl == true) ss << ","; descfl = true;
	ss << "\"&f& converts the external FOI such that it is per " << details.efoi_factor << " individuals.\"";
	
	if(eqfl == true) vis << ","; eqfl = true;
	vis << "{\"text\":\"&f&\"}";
	
	if(eqfl == true) vis << ","; eqfl = true;              // Adds EFOI
	vis << "{\"text\":\")\"}";
	
	if(susvar == true){
		if(eqfl == true) vis << ","; eqfl = true;
		vis << "{\"text\":\"]\"}";
	}
	
	vis << "]";
	
	vis << ",\"desc\":[" << ss.str() << "]";

	vis << "}";

	return vis.str();
}


/// Returns information about parameters used on splines
string Model::spline_param_JSON(const unsigned int sp) const
{
	stringstream ss;
	const auto &spl = spline[sp];
	
	ss << "[";
	for(auto i = 0u; i < spl.p.size(); i++){
		if(i > 0) ss << ",";
		ss << "{\"param\":\"" << param[spl.p[i].param].name << "\"";
		ss << ",\"time\":" << spl.p[i].t << "}";
	}	
	ss << "]";
	
	return ss.str();
}

