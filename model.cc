// This file defines the compartmental model and is used to simulate compartmental dynamics after an individual becomes infected

#include <string>
#include <sstream>
#include <iostream>
#include <cmath>
 
#include "timers.hh"
#include "math.h"

using namespace std;

#include "consts.hh"
#include "utils.hh"
#include "model.hh"
#include "utils.hh"
#include "inputs.hh"

/// Initialises the model 
Model::Model(const Inputs &inputs, const Details &details, Data &data) : details(details), data(data)
{
	maximum_infected = inputs.find_integer("infmax",LARGE);
	if((details.mode == MCMCMC || details.mode == ABC_SMC || details.mode == ABC_MBP) && maximum_infected == LARGE){
		emsgroot("Input file must contain a limit on the maximum number of individuals through 'maximum_infected'.");
	}	
	
	time_period = data.time_period; ntime_period = time_period.size();

	switch(details.mode){
		case SIM: case MULTISIM:
			{
				vector <string> name;
				vector <double> val;
				inputs.find_parameter(name,val);
				for(auto th = 0u; th < name.size(); th++){
					add_parameter(name[th],val[th],val[th]);
				}
			}
			break;
		
		case MCMCMC: case ABC_SMC: case ABC_MBP:
			{
				vector <string> name;
				vector <double> min,max;
				inputs.find_prior(name,min,max);
				for(auto th = 0u; th < name.size(); th++){
					add_parameter(name[th],min[th],max[th]);
				}
			}
			break;
	
		default:
			emsgEC("Model",56);
			break;
	}
	
	region_effect = 0;
	sigma_param = UNSET;
	for(auto th = 0u; th < param.size(); th++){
		if(param[th].name == "regeff_sigma"){
			param[th].used = true;
			
			region_effect = 1;
			sigma_param = th;	
			
			regioneffect_param.resize(data.nregion);
			for(auto r = 0u; r < data.nregion; r++){
				regioneffect_param[r] = param.size();
				stringstream ss; ss << "reff_" << data.region[r].code;
				add_parameter(ss.str(),-LARGE,LARGE);
				param[param.size()-1].used = true;
			}
		}
	}
	
	if(region_effect == 0){
		region_effect = 2;
		
		regioneffect_param.resize(data.nregion);
		for(auto r = 0u; r < data.nregion; r++){
			regioneffect_param[r] = param.size();
			stringstream ss; ss << "reff_" << data.region[r].code;
			add_parameter(ss.str(),-0.2,0.2);
			param[param.size()-1].used = true;
		}
	}
	
	add_parameter("zero",TINY,TINY);

	vector <string> name;
	vector <double> infectivity;
	inputs.find_compartments(name,infectivity);
	for(auto c = 0u; c < name.size(); c++) add_compartment(name[c],infectivity[c]);
	
	vector <string> from, to, prpar, mean, cv;
	vector <int> type;
	inputs.find_transitions(from,to,prpar,type,mean,cv);
	for(auto tr = 0u; tr < from.size(); tr++) add_transition(from[tr],to[tr],prpar[tr],type[tr],mean[tr],cv[tr]);
	
	if(details.mode == MCMCMC || details.mode == ABC_SMC || details.mode == ABC_MBP){
		priorcomps = inputs.find_priorcomps(comp);
	}
			
	vector <int> time;
	vector <string> pname; 
	string splinetype;

	for(auto loop = 0u; loop < 2; loop++){        // Loads up the different splines
		double fac;
		switch(loop){
		case 0: fac = 1; splinetype = "betaspline"; break;
		case 1: fac = 1.0/data.popsize; splinetype = "phispline"; break;
		}
		
		vector <SplinePoint> spli;
		inputs.find_spline(details,splinetype,time,pname);
		for(auto i = 0u; i < time.size(); i++){
			SplinePoint spl;
			spl.t = time[i];
			spl.param = find_parameter(pname[i]);
			spl.multfac = fac;
			
			spli.push_back(spl);
		}	
		spline.push_back(spli);
	}
	betaspline_ref = 0;
	phispline_ref = 1;
	
	suscetibility_param.resize(data.ndemocat);
	for(auto c = 0u; c < data.ndemocat; c++){
		for(auto fi = 0u; fi < data.democat[c].value.size(); fi++){
			suscetibility_param[c].push_back(find_parameter(data.democat[c].param[fi]));
		}
	}

	for(unsigned int c = 0; c < data.ncovar; c++){
		covariate_param.push_back(find_parameter(data.covar[c].param));
	}

	for(auto p = 0u; p < param.size(); p++){
		if(param[p].used == false) emsg("The [arameter '"+param[p].name+"' is not used in the model.");
	}

	for(const auto& tr : trans){
		switch(tr.type){
		case EXP_DIST: param[tr.param_mean].type = DISTVAL_PARAM; break;
		case GAMMA_DIST: case LOGNORM_DIST: param[tr.param_mean].type = DISTVAL_PARAM; param[tr.param_cv].type = DISTVAL_PARAM; break;
		}
		
		if(tr.istimep == false){
			for(const auto th : tr.probparam) param[th].type = BRANCHPROB_PARAM; 
		}
	}
	
	add_Q();
	check_data_files();
}

/// Outputs a summary of the model
void Model::print_to_terminal() const
{
	cout << endl;                                           
	switch(details.mode){
		case SIM: case MULTISIM:
			cout << "Parameters:" << endl;
			for(auto p = 0u; p < param.size()-1; p++){
				cout << "  " << param[p].name << " = " << param[p].min << endl;
			}
			break;
		
		case MCMCMC: case ABC_SMC: case ABC_MBP:
			cout << "Priors:" << endl;
			for(auto p = 0u; p < param.size()-1; p++){
				cout << "  " << param[p].name << " = ";
				if(param[p].min ==  param[p].max) cout << param[p].min << endl;
				else cout << "Uniform(" << param[p].min << " - " << param[p].max << ")" << endl;
			}
			break;
	}
	cout << endl;
		
	cout << "Compartments:" << endl; 
	for(const auto& co : comp){
		cout << "  " << co.name << "  Infectivity: " << co.infectivity << endl; 			
	}
	cout << endl;
	
	cout << "Transitions:" << endl; 
	for(const auto& tr : trans){
		if(tr.from != tr.to){
			cout << "  " << comp[tr.from].name << " → " << comp[tr.to].name;
			if(tr.probparam.size() > 0){
				cout << " with probability ";
				for(auto j = 0u; j < tr.probparam.size(); j++){
					if(j > 0) cout << ", ";
					cout << param[tr.probparam[j]].name;
				}
				cout << "  ";
			}
			
			switch(tr.type){
				case INFECTION_DIST: 
					cout << " Infection";
					break;
				case EXP_DIST:
					cout << " Exponential  mean=" << param[tr.param_mean].name;
					break;
				case GAMMA_DIST:
					cout << " Gamma mean=" << param[tr.param_mean].name 
							 << " cv=" << param[tr.param_cv].name; 
					break;
				case LOGNORM_DIST:
					cout << " Lognormal mean=" << param[tr.param_mean].name 
							 << " cv=" << param[tr.param_cv].name; 
					break;
				default:
					break;
			}
			cout << endl;
		}
	}
	cout << endl;
}

/// Adds in the tensor Q to the model
void Model::add_Q()
{
	for(const auto& co : comp) add_transition(co.name,co.name,"",TIMEP_DIST,"","");  
	
	for(const auto& QQ : data.Q){
		unsigned int c;
		for(c = 0; c < comp.size(); c++) if(QQ.comp == comp[c].name) break;
		if(c == comp.size()) emsg("The compartment '"+QQ.comp+"' in '"+QQ.name+"' is not recognised.");
	}
	
	DQinfo dq;
	dq.q.resize(2); dq.fac.resize(2);
	for(auto& tr : trans){
		auto ci = tr.from;
		auto cf = tr.to;
		auto compi = comp[ci].name;
		auto compf = comp[cf].name;
	
		for(auto timep = 0u; timep < ntime_period; timep++){
			auto timepi = timep, timepf = timep;
			if(compi == compf) timepf++;
			
			if(timepf < ntime_period){
				auto qi = 0u; while(qi < data.Q.size() && !(data.Q[qi].comp == compi && data.Q[qi].timep == timepi)) qi++;
				if(qi == data.Q.size()) qi = UNSET;
				
				auto qf = 0u; while(qf < data.Q.size() && !(data.Q[qf].comp == compf && data.Q[qf].timep == timepf)) qf++;
				if(qf == data.Q.size()) qf = UNSET;
				
				if(qi == UNSET && qf == UNSET){
					tr.DQ.push_back(UNSET);
				}
				else{
					tr.DQ.push_back(DQ.size());
					dq.q[0] = qi; dq.fac[0] = -comp[ci].infectivity;
					dq.q[1] = qf; dq.fac[1] = comp[cf].infectivity; 
					DQ.push_back(dq);
				}
			}
		}
	}
}

/// Randomly samples the initial parameter values from the prior (which are uniform distributions
vector <double> Model::sample_from_prior() const
{
	vector <double> paramv(param.size());
	for(auto th = 0u; th < param.size(); th++){	
		paramv[th] = param[th].min + ran()*(param[th].max - param[th].min);
	}
	
	if(region_effect == 1){
		for(auto th : regioneffect_param) paramv[th] = normal(0,paramv[sigma_param]);
	}
	
	if(smooth_spline == 1){
		
		for(const auto &spli : spline){
			for(auto i = 0u; i < spli.size()-1; i++){
				auto th1 = spli[i].param;
				auto th2 = spli[i+1].param;
				
				if(spli[i].t != spli[i+1].t && th1 != th2 && param[th2].min != param[th2].max){
					double val;
					do{
						double fac = exp(normal(0,smooth));
						val = paramv[th1]*fac;
					}while(val <= param[th2].min || val >= param[th2].max);
					paramv[th2] = val;
				}
			}
		}
	}
	//for(auto th = 0u; th < param.size(); th++) cout << "paramv[" << th <<"] = " << paramv[th] << ";" << endl;

	return paramv;
}

/// Gets the infectivity of a compartment
double Model::get_infectivity(const string& name) const
{
	auto c = 0u; while(c < comp.size() && comp[c].name != name) c++;
	if(c == comp.size()) emsg("Cannot find the compartment '"+name+"'");
	return comp[c].infectivity;
}

/// Adds a compartment to the model
void Model::add_compartment(const string& name, double infectivity)
{
	Compartment co;
	co.name = name;
	co.infectivity = infectivity;
	co.transtimep = UNSET;

	comp.push_back(co);	
}

/// Finds a parameter from a string
unsigned int Model::find_parameter(const string& name)
{
	auto p = 0u; auto pmax = param.size(); while(p < pmax && name != param[p].name) p++;
	if(p == pmax) emsg("Cannot find the parameter '"+name+"'");	
	param[p].used = true;
	
	return p;
}

/// Adds a parameter to the model
void Model::add_parameter(const string& name, double min, double max)
{
	Param par;
	par.name = name; par.min = min; par.max = max; par.ntr = 0; par.nac = 0; par.jump = 0.5*(min+max)/10; if(par.jump == 0) par.jump = 0.1;
	par.used = false; par.type = OTHER_PARAM;

	param.push_back(par);
}

/// Adds a transition to the model
void Model::add_transition(const string& from, const string& to, const string& prpar, unsigned int type, const string& mean, const string& cv)
{
	Transition tr;
	auto c = 0u; auto cmax = comp.size(); while(c < cmax && from != comp[c].name) c++;
	if(c == cmax) emsg("Cannot find the 'from' compartment '"+from+"' for the transition");
	tr.from = c;	
	
	if(from != to){ tr.istimep = false; comp[c].trans.push_back(trans.size());}
	else{ tr.istimep = true; comp[c].transtimep = trans.size();}
		
	c = 0; cmax = comp.size(); while(c < cmax && to != comp[c].name) c++;
	if(c == cmax) emsg("Cannot find the 'to' compartment '"+to+"' for the transition");	
	tr.to = c;
	
	if(prpar != ""){
		vector <string> probparam = split(prpar,',');
	
		if(probparam.size() != data.nage){
			emsg("For the transition '"+from+"→"+to+"' the number of parameters in expression '"+prpar+"' should equal the number of age groups.");
		}
		
		for(const auto& parname : probparam) tr.probparam.push_back(find_parameter(parname));
	}
	
	tr.type = type;
	if(mean != "") tr.param_mean = find_parameter(mean);	else tr.param_mean = UNSET;
	if(cv != "") tr.param_cv = find_parameter(cv); else tr.param_cv = UNSET;
	
	trans.push_back(tr);
}

vector <double> Model::create_disc_spline(unsigned int ref, const vector<double> &paramv) const
{
	vector <double> disc_spline(details.ndivision);
	
	auto p = 0;
	for(auto s = 0u; s < details.ndivision; s++){	
		auto t = double((s+0.5)*details.period)/details.ndivision;
		
		while(p < int(spline[ref].size())-1 && t > spline[ref][p+1].t) p++;
		
		auto fac = (t-spline[ref][p].t)/(spline[ref][p+1].t-spline[ref][p].t);
		disc_spline[s] = (paramv[spline[ref][p].param]*(1-fac) + paramv[spline[ref][p+1].param]*fac)*spline[ref][p].multfac;
	}
	
	return disc_spline;
}

/// Sets the transition probabilies based on the parameters
unsigned int Model::create_comptransprob(vector <CompTransProb> &comptransprob, const vector<double> &paramv) const
{
	comptransprob.resize(comp.size());
	for(auto c = 0u; c < comp.size(); c++){
		auto kmax = comp[c].trans.size();
		if(kmax > 1){
			comptransprob[c].prob.resize(data.nage);
			comptransprob[c].probsum.resize(data.nage);
			for(auto a = 0u; a < data.nage; a++){
				comptransprob[c].prob[a].resize(kmax);
				comptransprob[c].probsum[a].resize(kmax);
				
				auto sum = 0.0;
				for(auto k = 0u; k < kmax-1; k++){
					auto p = trans[comp[c].trans[k]].probparam[a]; if(p == UNSET) emsgEC("model",1);		
					auto prob = paramv[p]; comptransprob[c].prob[a][k] = prob; sum += prob; if(prob < 0) return 1;
				} 
				auto prob = 1-sum; comptransprob[c].prob[a][kmax-1] = prob; if(prob < 0) return 1;
				
				sum = 0.0;
				for(auto k = 0u; k < kmax; k++){
					sum += comptransprob[c].prob[a][k];
					comptransprob[c].probsum[a][k] = sum;
				}
			}
		}
	}
	
	return 0;
}

 
/// Defines the relative susceptibility of individuals
vector <double> Model::create_susceptibility(const vector<double> &paramv) const 
{
	vector <double> susceptibility(data.ndemocatpos);
	for(auto dp = 0u; dp < data.ndemocatpos; dp++){
		auto val = 1.0;
		for(auto c = 0u; c < data.ndemocat; c++){
			auto j = data.democatpos[dp][c];
			val *= exp(paramv[suscetibility_param[c][j]]);
		}
		susceptibility[dp] = val;
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
		
		if(region_effect != 0) sum += paramv[regioneffect_param[data.area[c].region]];
		
		areafactor[c] = exp(sum);
		if(std::isnan(areafactor[c])) emsgEC("Model",90);
	}
	
	return areafactor;
}

/// Checks that the transition and population data is correct
void Model::check_data_files() const
{
	for(auto& td : data.transdata){
		auto from = td.fromstr, to = td.tostr; 
		
		unsigned int tra;
		for(tra = 0u; tra < trans.size(); tra++){ if(comp[trans[tra].from].name == from && comp[trans[tra].to].name == to) break;}
		if(tra == trans.size()) emsg("Cannot find the transition "+from+"→"+to+" in file '"+td.file+"'.");
		
		td.trans = tra;
	}
	
	for(auto& pd : data.popdata){
		auto compstr = pd.compstr;
	
		unsigned int c;
		for(c = 0u; c < comp.size(); c++){ if(comp[c].name == compstr) break;}
		if(c == comp.size()) emsg("Cannot find the compartment '"+compstr+"' in file '"+pd.file+"'.");
		
		pd.comp = c;
	}
	
	for(auto& md : data.margdata){
		auto from = md.fromstr, to = md.tostr; 
		
		unsigned int tra;
		for(tra = 0; tra < trans.size(); tra++){ if(comp[trans[tra].from].name == from && comp[trans[tra].to].name == to) break;}
		if(tra == trans.size()) emsg("Cannot find the transition "+from+"→"+to+" in file '"+md.file+"'.");
	
		md.trans = tra;
	}
}

/// Calculates the prior probability
double Model::prior(const vector<double> &paramv) const 
{
	double Pr = 0;
	
	/*
	unsigned int pc, a, c;
	double prob;
	if(priorcomps.size() > 0){
		calcprobin();
		for(pc = 0; pc < priorcomps.size(); pc++){
			c = priorcomps[pc].comp;
			prob = 0; for(a = 0; a < data.nage; a++) prob += data.agedist[a]*comp[c].probin[a];
			Pr += normalprob(prob,priorcomps[pc].value,priorcomps[pc].sd*priorcomps[pc].sd);
		}
	}
	*/
	
	if(region_effect == 1){
		auto sd = paramv[sigma_param];
		for(auto th : regioneffect_param) Pr += normalprob(paramv[th],0,sd*sd);
	}
	
	if(smooth_spline == 1){
		for(const auto& spli : spline){
			for(auto i = 0u; i < spli.size()-1; i++){
				if(spli[i].t != spli[i+1].t){
					double dval = log(paramv[spli[i+1].param]/paramv[spli[i].param])/smooth;
				
					Pr += -dval*dval/2;
				}
			}
		}
	}
	
	return Pr;
}
	
/// Calculate compartmental probabilities
vector <CompProb> Model::create_compprob(const vector <CompTransProb> &comptransprob) const
{
	vector <CompProb> compprob(comp.size());
	
	for(auto& co : compprob){
		co.value.resize(data.nage);
		for(auto& value : co.value) value = 0;
	}

	for(auto a = 0u; a < data.nage; a++){
		vector <unsigned int> cst, kst;
		vector <double> probst;
	
		auto prob = 1.0;
		auto c = 0u;
		do{
			compprob[c].value[a] += prob;
			
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
			if(comp[c].trans.size() > 1) prob *= comptransprob[c].prob[a][k];
			c = trans[comp[c].trans[k]].to;		
		}while(1 == 1);
		
		if(cst.size() > 0 || kst.size() > 0 || probst.size() > 0) emsgEC("Model",8);
	}
	
	return compprob;
}

/// Calculates R0
vector <double> Model::calculate_R_vs_time(const vector<double> &paramv) const
{
	vector <CompTransProb> comptransprob;
	if(create_comptransprob(comptransprob,paramv) == 1) emsgEC("Model",81);
	
	vector <CompProb> compprob = create_compprob(comptransprob);
	
	vector < vector<double> > infint;   // Stores the time integrated infectivity
	
	infint.resize(comp.size());
	for(auto c = 0u; c < comp.size(); c++){
		infint[c].resize(data.nage);
		for(auto a = 0u; a < data.nage; a++){
			infint[c][a] = 0;
		
			auto kmax = comp[c].trans.size();
			for(auto k = 0u; k < kmax; k++){
				auto tra = comp[c].trans[k];
				
				double dt;
				switch(trans[tra].type){
				case INFECTION_DIST: dt = 0; break;
				case EXP_DIST: dt = paramv[trans[tra].param_mean]; break;
				case GAMMA_DIST: dt = paramv[trans[tra].param_mean]; break;
				case LOGNORM_DIST: dt = paramv[trans[tra].param_mean]; break;
				default: emsgEC("Model",9); break;
				}	
				if(kmax == 1) infint[c][a] += compprob[c].value[a]*comp[c].infectivity*dt;
				else infint[c][a] += compprob[c].value[a]*comp[c].infectivity*dt*comptransprob[c].prob[a][k];
			}
		}
	}
	
	vector <double> R0fac(ntime_period);
	for(auto& R0fa : R0fac) R0fa = 0;
		
	vector <double> sus = create_susceptibility(paramv);
	vector <double> areafactor = create_areafactor(paramv);
	
	for(const auto& QQ : data.Q){
		auto timep = QQ.timep;
		auto qt = QQ.Qtenref;
		
		unsigned int co;
		for(co = 0u; co < comp.size(); co++) if(QQ.comp == comp[co].name) break;
		if(co == comp.size()) emsg("Compartment '"+QQ.comp+"' in '"+QQ.name+"' not recognised.");
		
		for(auto c = 0u; c < data.narea; c++){
			for(auto a = 0u; a < data.nage; a++){
				auto fac = infint[co][a]*double(data.area[c].agepop[a])/data.popsize; 
				
				if(fac != 0){
					auto vi = c*data.nage + a;
					auto kmax = data.genQ.Qten[qt].ntof[vi];
					for(auto k = 0u; k < kmax; k++){
						auto cc = data.genQ.Qten[qt].tof[vi][k];
						for(auto dp = 0u; dp < data.ndemocatpos; dp++){
							auto aa = data.democatpos[dp][0];
							R0fac[timep] += fac*sus[dp]*areafactor[cc]*data.genQ.Qten[qt].valf[vi][k][aa]*data.area[cc].pop[dp];
						}
					}
				}
			}
		}
	}
	
	vector <double> R0(details.ndivision);

	vector <double> beta = create_disc_spline(betaspline_ref,paramv);
	auto timep = 0u; 
	for(auto st = 0u; st < details.ndivision; st++){
		auto t = details.division_time[st];
		while(timep < ntime_period && t > time_period[timep].tend) timep++;
		
		R0[st] = beta[st]*R0fac[timep];
	}
		
	return R0;
}

bool Model::inbounds(const vector <double> &paramv) const
{
	for(auto th = 0u; th < paramv.size(); th++){
		if(inbounds(paramv[th],th) == false) return false;
	}	
	return true;
}

bool Model::inbounds(double val, unsigned int th) const
{
	if(val >= param[th].min && val <= param[th].max) return true;
	return false;
}

/// Outputs an event sequence (used for debugging)
void Model::print_events(const string& name, const vector <Event> &ev) const 
{
	cout << name << ":" << endl;
	for(auto e = 0u; e < ev.size(); e++){
		auto tra = ev[e].trans;
		cout << comp[trans[tra].from].name << "->" << comp[trans[tra].to].name << "  " << ev[e].t << endl;
	}
}

/// Determines if it is necessary to do mbp for the exisiting event sequence
bool Model::do_mbp_events(const vector <double> &parami, const vector <double> &paramp) const
{
	for(auto th = 0u; th < param.size(); th++){
		if(parami[th] != paramp[th] && (param[th].type == DISTVAL_PARAM || param[th].type == BRANCHPROB_PARAM)) return true;
	}
	return false;
}
