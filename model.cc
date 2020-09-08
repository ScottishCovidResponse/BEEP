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
MODEL::MODEL(Inputs &inputs, const Details &details, DATA &data) : details(details), data(data)
{
	infmax = inputs.find_int("infmax",large);
	if((details.mode == inf || details.mode == abcsmc || details.mode == abcmbp) && infmax == large){
		emsgroot("Input file must contain a limit on the maximum number of individuals through 'infmax'.");
	}	
	
	timeperiod = data.timeperiod; ntimeperiod = timeperiod.size();

	if(details.mode == sim || details.mode == multisim){
		vector <string> name;
		vector <double> val;
		inputs.find_param(name,val);
		for(auto th = 0u; th < name.size(); th++){
			addparam(name[th],val[th],val[th]);
		}
	}

	if(details.mode == inf || details.mode == abcsmc || details.mode == abcmbp){
		vector <string> name;
		vector <double> min,max;
		inputs.find_prior(name,min,max);
		for(auto th = 0u; th < name.size(); th++){
			addparam(name[th],min[th],max[th]);
		}
	}
	
	regioneffect = 0;
	sigma_param = UNSET;
	for(auto th = 0u; th < param.size(); th++){
		if(param[th].name == "regeff_sigma"){
			param[th].used = 1;
			
			regioneffect = 1;
			sigma_param = th;	
			
			regioneff_param.resize(data.nregion);
			for(auto r = 0u; r < data.nregion; r++){
				regioneff_param[r] = param.size();
				stringstream ss; ss << "reff_" << data.region[r].code;
				addparam(ss.str(),-large,large);
				param[param.size()-1].used = 1;
			}
		}
	}
	
	if(regioneffect == 0){
		regioneffect = 2;
		
		regioneff_param.resize(data.nregion);
		for(auto r = 0u; r < data.nregion; r++){
			regioneff_param[r] = param.size();
			stringstream ss; ss << "reff_" << data.region[r].code;
			addparam(ss.str(),-0.2,0.2);
			param[param.size()-1].used = 1;
		}
	}
	
	addparam("zero",tiny,tiny);

	vector <string> name;
	vector <double> infectivity;
	inputs.find_comps(name,infectivity);
	for(auto c = 0u; c < name.size(); c++) addcomp(name[c],infectivity[c]);
	
	vector <string> from, to, prpar, mean, cv;
	vector <int> type;
	inputs.find_trans(from,to,prpar,type,mean,cv);
	for(auto tr = 0u; tr < from.size(); tr++) addtrans(from[tr],to[tr],prpar[tr],type[tr],mean[tr],cv[tr]);
	
	if(details.mode == inf || details.mode == abcsmc || details.mode == abcmbp){
		priorcomps = inputs.find_priorcomps(comp);
	}
			
	vector <int> time;
	vector <string> pname; 
	string splinetype;

	splinetype = "betaspline";
	inputs.find_spline(details,splinetype,time,pname);
	for(auto i = 0u; i < time.size(); i++){
		SPLINEP spl;
		spl.t = time[i];
		spl.param = findparam(pname[i]);
		betaspline.push_back(spl);
	}	
	
	splinetype = "phispline";
	inputs.find_spline(details,splinetype,time,pname);
	for(auto i = 0u; i < time.size(); i++){
		SPLINEP spl;
		spl.t = time[i];
		spl.param = findparam(pname[i]);
		phispline.push_back(spl);
	}	
		
	sus_param.resize(data.ndemocat);
	for(auto c = 0u; c < data.ndemocat; c++){
		for(auto fi = 0u; fi < data.democat[c].value.size(); fi++){
			sus_param[c].push_back(findparam(data.democat[c].param[fi]));
		}
	}

	for(unsigned int c = 0; c < data.ncovar; c++){
		covar_param.push_back(findparam(data.covar[c].param));
	}

	for(auto p = 0u; p < param.size(); p++){
		if(param[p].used == 0) emsg("The [arameter '"+param[p].name+"' is not used in the model.");
	}

	ntr = 0; nac = 0;
	parami.resize(param.size()); paramp.resize(param.size());

	beta.resize(details.nsettime); phi.resize(details.nsettime);
	
	for(const auto& tr : trans){
		switch(tr.type){
		case exp_dist: param[tr.param_mean].type = distval_paramtype; break;
		case gamma_dist: case lognorm_dist: param[tr.param_mean].type = distval_paramtype; param[tr.param_cv].type = distval_paramtype; break;
		}
		
		if(tr.istimep == 0){
			for(const auto th : tr.probparam) param[th].type = branchprob_paramtype; 
		}
	}
	
	addQ();
	checkdata();
}

/// Outputs a summary of the model
void MODEL::print_to_terminal() const
{
	cout << endl;                                           
	if(details.mode == sim || details.mode == multisim){
		cout << "Parameters:" << endl;
		for(auto p = 0u; p < param.size()-1; p++){
			cout << "  " << param[p].name << " = " << param[p].min << endl;
		}
	}
	else{
		cout << "Priors:" << endl;
		for(auto p = 0u; p < param.size()-1; p++){
			cout << "  " << param[p].name << " = ";
			if(param[p].min ==  param[p].max) cout << param[p].min << endl;
			else cout << "Uniform(" << param[p].min << " - " << param[p].max << ")" << endl;
		}
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
				case infection_dist: 
					cout << " Infection";
					break;
				case exp_dist:
					cout << " Exponential  mean=" << param[tr.param_mean].name;
					break;
				case gamma_dist:
					cout << " Gamma mean=" << param[tr.param_mean].name 
							 << " cv=" << param[tr.param_cv].name; 
					break;
				case lognorm_dist:
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

/// Sets up the model with a given set of parameters
unsigned int MODEL::setup(const vector <double> &paramv)
{
	if(settransprob(paramv) == 1) return 1;
	
	timevariation(paramv);
	setsus(paramv);
	setarea(paramv);
	return 0;
}

/// Copies values used for the initial state (used in MBPs)
void MODEL::copyi(const vector<double> &paramv)
{
	parami = paramv; betai = beta; phii = phi; susi = sus; areafaci = areafac;
	for(auto& co : comp) co.probi = co.prob;
}
	
/// Copies values used for the proposed state (used in MBPs)
void MODEL::copyp(const vector<double> &paramv)
{
	paramp = paramv; betap = beta; phip = phi; susp = sus; areafacp = areafac;
	for(auto& co : comp) co.probp = co.prob;
}

/// Adds in the tensor Q to the model
void MODEL::addQ()
{
	for(const auto& co : comp) addtrans(co.name,co.name,"",timep_dist,"","");  
	
	for(const auto& QQ : data.Q){
		unsigned int c;
		for(c = 0; c < comp.size(); c++) if(QQ.comp == comp[c].name) break;
		if(c == comp.size()) emsg("The compartment '"+QQ.comp+"' in '"+QQ.name+"' is not recognised.");
	}
	
	DQINFO dq;
	dq.q.resize(2); dq.fac.resize(2);
	for(auto& tr : trans){
		auto ci = tr.from;
		auto cf = tr.to;
		auto compi = comp[ci].name;
		auto compf = comp[cf].name;
	
		for(auto timep = 0u; timep < ntimeperiod; timep++){
			auto timepi = timep, timepf = timep;
			if(compi == compf) timepf++;
			
			if(timepf < ntimeperiod){
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
vector <double> MODEL::priorsamp()
{
	vector <double> paramv(param.size());
	for(auto th = 0u; th < param.size(); th++){	
		paramv[th] = param[th].min + ran()*(param[th].max - param[th].min);
	}
	
	if(regioneffect == 1){
		for(auto th : regioneff_param) paramv[th] = normal(0,paramv[sigma_param]);
	}
	
	if(smooth_spline == 1){
		for(auto i = 0u; i < betaspline.size()-1; i++){
			auto th1 = betaspline[i].param;
			auto th2 = betaspline[i+1].param;
			
			if(betaspline[i].t != betaspline[i+1].t && th1 != th2 && param[th2].min != param[th2].max){
				double val;
				do{
					double fac = exp(normal(0,smooth));
					val = paramv[th1]*fac;
				}while(val <= param[th2].min || val >= param[th2].max);
			  paramv[th2] = val;
			}
		}
		
		for(auto i = 0u; i < phispline.size()-1; i++){
			auto th1 = phispline[i].param;
			auto th2 = phispline[i+1].param;
			
			if(phispline[i].t != phispline[i+1].t && th1 != th2 && param[th2].min != param[th2].max){
				double val;
				do{
					double fac = exp(normal(0,smooth));
					val = paramv[th1]*fac;
				}while(val <= param[th2].min || val >= param[th2].max);
			  paramv[th2] = val;
			}
		}
	}
	//for(auto th = 0u; th < param.size(); th++) cout << "paramv[" << th <<"] = " << paramv[th] << ";" << endl;

	return paramv;
}

/// Gets the infectivity of a compartment
double MODEL::getinfectivity(const string& name)
{
	auto c = 0u; while(c < comp.size() && comp[c].name != name) c++;
	if(c == comp.size()) emsg("Cannot find the compartment '"+name+"'");
	return comp[c].infectivity;
}

/// Adds a compartment to the model
void MODEL::addcomp(const string& name, double infectivity)
{
	COMP co;
	co.name = name;
	co.infectivity = infectivity;
	co.transtimep = UNSET;

	comp.push_back(co);	
}

/// Finds a parameter from a string
unsigned int MODEL::findparam(const string& name)
{
	auto p = 0u; auto pmax = param.size(); while(p < pmax && name != param[p].name) p++;
	if(p == pmax) emsg("Cannot find the parameter '"+name+"'");	
	param[p].used = 1;
	
	return p;
}

/// Adds a parameter to the model
void MODEL::addparam(const string& name, double min, double max)
{
	PARAM par;
	par.name = name; par.min = min; par.max = max; par.ntr = 0; par.nac = 0; par.jump = 0.5*(min+max)/10; if(par.jump == 0) par.jump = 0.1;
	par.used = 0; par.type = other_paramtype;

	param.push_back(par);
}

/// Adds a transition to the model
void MODEL::addtrans(const string& from, const string& to, const string& prpar, unsigned int type, const string& mean, const string& cv)
{
	TRANS tr;
	auto c = 0u; auto cmax = comp.size(); while(c < cmax && from != comp[c].name) c++;
	if(c == cmax) emsg("Cannot find the 'from' compartment '"+from+"' for the transition");
	tr.from = c;	
	
	if(from != to){ tr.istimep = 0; comp[c].trans.push_back(trans.size());}
	else{ tr.istimep = 1; comp[c].transtimep = trans.size();}
		
	c = 0; cmax = comp.size(); while(c < cmax && to != comp[c].name) c++;
	if(c == cmax) emsg("Cannot find the 'to' compartment '"+to+"' for the transition");	
	tr.to = c;
	
	if(prpar != ""){
		vector <string> probparam = split(prpar,',');
	
		if(probparam.size() != data.nage){
			emsg("For the transition '"+from+"→"+to+"' the number of parameters in expression '"+prpar+"' should equal the number of age groups.");
		}
		
		for(const auto& parname : probparam) tr.probparam.push_back(findparam(parname));
	}
	
	tr.type = type;
	if(mean != "") tr.param_mean = findparam(mean);	else tr.param_mean = UNSET;
	if(cv != "") tr.param_cv = findparam(cv); else tr.param_cv = UNSET;
	tr.numvisittot = UNSET;
	tr.dtsum = 0.;
	
	trans.push_back(tr);
}
	
/// Generates the time variation in beta and phi from the parameters
void MODEL::timevariation(const vector<double> &paramv)
{
  // Uses a linear spline for beta
	auto p = 0;
	for(auto s = 0u; s < details.nsettime; s++){	
		auto t = double((s+0.5)*details.period)/details.nsettime;
		
		while(p < int(betaspline.size())-1 && t > betaspline[p+1].t) p++;
		
		auto fac = (t-betaspline[p].t)/(betaspline[p+1].t-betaspline[p].t);
		beta[s] = (paramv[betaspline[p].param]*(1-fac) + paramv[betaspline[p+1].param]*fac);
	}
	
	// Uses a linear spline for phi
	p = 0;
	for(auto s = 0u; s < details.nsettime; s++){	
		auto t = double((s+0.5)*details.period)/details.nsettime;
		
		while(p < int(phispline.size())-1 && t > phispline[p+1].t) p++;
		
		auto fac = (t-phispline[p].t)/(phispline[p+1].t-phispline[p].t);
		phi[s] = (paramv[phispline[p].param]*(1-fac) + paramv[phispline[p+1].param]*fac)/data.popsize;
	}
}

/// Sets the transition probabilies based on the parameters
unsigned int MODEL::settransprob(const vector<double> &paramv)
{
	for(auto& co : comp){
		auto kmax = co.trans.size();
		if(kmax > 1){
			co.prob.resize(data.nage);
			co.probsum.resize(data.nage);
			for(auto a = 0u; a < data.nage; a++){
				co.prob[a].resize(kmax);
				co.probsum[a].resize(kmax);
				
				auto sum = 0.0;
				for(auto k = 0u; k < kmax-1; k++){
					auto p = trans[co.trans[k]].probparam[a]; if(p == UNSET) emsgEC("model",1);		
					auto prob = paramv[p]; co.prob[a][k] = prob; sum += prob; if(prob < 0) return 1;
				} 
				auto prob = 1-sum; co.prob[a][kmax-1] = prob; if(prob < 0) return 1;
				
				sum = 0.0;
				for(auto k = 0u; k < kmax; k++){
					sum += co.prob[a][k];
					co.probsum[a][k] = sum;
				}
			}
		}
	}
	
	return 0;
}

/// This simulates from the model and generates an event list
void MODEL::simmodel(const vector<double> &paramv, vector <FEV> &evlist, unsigned int i, unsigned int c, double t)
{
	auto timep = 0u; while(timep < ntimeperiod && t > timeperiod[timep].tend) timep++;

	FEV ev;
	ev.ind = i; ev.timep = timep; 
		
	auto a = data.democatpos[data.ind[i].dp][0];
		
	unsigned int tra;
	if(c == 0){
		evlist.clear();	
		tra = 0;
		ev.trans = tra; ev.ind = i; ev.t = t; 
		evlist.push_back(ev);
	
		c = trans[tra].to;
	}
	
	double dt, z;
	do{
		auto kmax = comp[c].trans.size();
		if(kmax == 0) break;
		
		if(kmax == 1) tra = comp[c].trans[0];
		else{
			z = ran(); auto k = 0u; while(k < kmax && z > comp[c].probsum[a][k]) k++;
			if(k == kmax) emsgEC("Model",2);
			tra = comp[c].trans[k];
		}
		
		switch(trans[tra].type){
		case exp_dist:
			dt = -log(ran())*paramv[trans[tra].param_mean];
			break;
		
		case gamma_dist:
			{
				auto mean = paramv[trans[tra].param_mean]; auto sd = paramv[trans[tra].param_cv]*mean;
				dt = gammasamp(mean*mean/(sd*sd),mean/(sd*sd));
			}
			break;
			
		case lognorm_dist:
			{
				auto mean_ns = paramv[trans[tra].param_mean], cv_ns = paramv[trans[tra].param_cv];
				auto sd = sqrt(log((1+cv_ns*cv_ns))), mean = log(mean_ns) - sd*sd/2;
				dt = exp(normal(mean,sd));
			}
			break;
			
		default: emsgEC("MODEL",3); break;
		}

		if(dt < tiny) dt = tiny;
		t += dt;
		
		while(timep < ntimeperiod-1 && timeperiod[timep].tend < t){    // Adds in changes in time period
			ev.trans = comp[c].transtimep; ev.t = timeperiod[timep].tend;
			evlist.push_back(ev);
			timep++;
			ev.timep = timep; 
		}
		
		ev.trans = tra; ev.t = t;
		evlist.push_back(ev);

		c = trans[tra].to; 
	}while(1 == 1);
}

/// This does an equivelent MBP for the compartmental model
void MODEL::mbpmodel(vector <FEV> &evlisti, vector <FEV> &evlistp)
{
	evlistp.clear();
	
	FEV ev = evlisti[0];
	auto tra = ev.trans;
	auto ti = ev.t; 
	auto i = ev.ind; 
	auto timep = ev.timep;
	evlistp.push_back(ev);
	
	auto a = data.democatpos[data.ind[i].dp][0];
		
	auto c = trans[tra].to;
	
	auto tp = ti;
	unsigned int emax = evlisti.size();
	for(auto e = 1u; e < emax; e++){
		tra = evlisti[e].trans;
		if(trans[tra].istimep == 0){		
			auto dt = evlisti[e].t - ti;
			ti = evlisti[e].t;
			
			unsigned int kmax = comp[c].trans.size();
			if(kmax == 0) break;
			
			if(kmax > 1){
				auto k = 0u; while(k < kmax && tra != comp[c].trans[k]) k++;
				if(k == kmax) emsgEC("model",4);
				
				if(comp[c].probp[a][k] < comp[c].probi[a][k]){  // Looks at switching to another branch
					if(ran() < 1 - comp[c].probp[a][k]/comp[c].probi[a][k]){
						auto sum = 0.0; 
						vector <double> sumst;
						for(auto k = 0u; k < kmax; k++){
							auto dif = comp[c].probp[a][k] - comp[c].probi[a][k];
							if(dif > 0) sum += dif;
							sumst.push_back(sum);
						}
						
						auto z = ran()*sum; auto k = 0u; while(k < kmax && z > sumst[k]) k++;
						if(k == kmax) emsgEC("Model",5);
						tra = comp[c].trans[k];
					}
				}	
			}
			
			double dtnew;
			switch(trans[tra].type){
			case exp_dist:
				{
					auto p = trans[tra].param_mean;
					dtnew = dt*paramp[p]/parami[p];
				}
				break;
			
			case gamma_dist:
				emsgEC("model",6);
				break;
				
			case lognorm_dist:
				{
					auto p = trans[tra].param_mean, p2 = trans[tra].param_cv;
					
					auto mean_nsi = parami[p], cv_nsi = parami[p2];
					auto mean_nsp = paramp[p], cv_nsp = paramp[p2];
					
					if(mean_nsi == mean_nsp && cv_nsi == cv_nsp) dtnew = dt;
					else{
						auto sdi = sqrt(log(1+cv_nsi*cv_nsi)), meani = log(mean_nsi) - sdi*sdi/2;
						auto sdp = sqrt(log(1+cv_nsp*cv_nsp)), meanp = log(mean_nsp) - sdp*sdp/2;
						dtnew = exp(meanp + (log(dt) - meani)*sdp/sdi); 
					}
				}
				break;
				
			default:
				emsgEC("Model",7);
				break;
			}
	
			if(dtnew < tiny) dtnew = tiny;
			tp += dtnew;
	
			while(timep < ntimeperiod-1 && timeperiod[timep].tend < tp){    // Adds in changes in time period
				ev.trans = comp[c].transtimep; ev.t = timeperiod[timep].tend;
				evlistp.push_back(ev);
				timep++;
				ev.timep = timep; 
			}
			
			ev.trans = tra; ev.t = tp;
			evlistp.push_back(ev);
			
			c = trans[tra].to; 
			if(evlisti[e].trans != tra) break;
		}
	}
	
	if(comp[c].trans.size() != 0)	simmodel(paramp,evlistp,i,c,tp);
}
  
/// Defines the relative susceptibility of individuals
void MODEL::setsus(const vector<double> &paramv)     
{
	sus.resize(data.ndemocatpos);
	for(auto dp = 0u; dp < data.ndemocatpos; dp++){
		auto val = 1.0;
		for(auto c = 0u; c < data.ndemocat; c++){
			auto j = data.democatpos[dp][c];
			val *= exp(paramv[sus_param[c][j]]);
		}
		sus[dp] = val;
	}
}

/// Defines the relative transmission rate for different areas
void MODEL::setarea(const vector<double> &paramv) 
{
	areafac.resize(data.narea);
	for(auto c = 0u; c < data.narea; c++){
		auto sum = 0.0;
		for(auto j = 0u; j < data.ncovar; j++) sum += paramv[covar_param[j]]*data.area[c].covar[j];
		
		if(regioneffect != 0) sum += paramv[regioneff_param[data.area[c].region]];
		
		areafac[c] = exp(sum);
		if(std::isnan(areafac[c])) emsgEC("Model",90);
	}
}

/// Checks that the transition and population data is correct
void MODEL::checkdata()
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
double MODEL::prior(const vector<double> &paramv)
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
	
	if(regioneffect == 1){
		auto sd = paramv[sigma_param];
		for(auto th : regioneff_param) Pr += normalprob(paramv[th],0,sd*sd);
	}
	
	if(smooth_spline == 1){
		for(auto i = 0u; i < betaspline.size()-1; i++){
			if(betaspline[i].t != betaspline[i+1].t){
				double dval = log(paramv[betaspline[i+1].param]/paramv[betaspline[i].param])/smooth;
			
				Pr += -dval*dval/2;
			}
		}

		for(auto i = 0u; i < phispline.size()-1; i++){
			if(phispline[i].t != phispline[i+1].t){
				double dval = log(paramv[phispline[i+1].param]/paramv[phispline[i].param])/smooth;
					
				Pr += -dval*dval/2;
			}
		}
	}
	
	return Pr;
}
	
/// Calculate compartmental probabilities
void MODEL::calcprobin()
{
	for(auto& co : comp){
		co.probin.resize(data.nage);
		for(auto& probin : co.probin) probin = 0;
	}

	for(auto a = 0u; a < data.nage; a++){
		vector <unsigned int> cst, kst;
		vector <double> probst;
	
		auto prob = 1.0;
		auto c = 0u;
		do{
			comp[c].probin[a] += prob;
			
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
			if(comp[c].trans.size() > 1) prob *= comp[c].prob[a][k];
			c = trans[comp[c].trans[k]].to;		
		}while(1 == 1);
		
		if(cst.size() > 0 || kst.size() > 0 || probst.size() > 0) emsgEC("Model",8);
	}
}

/// Calculates R0
vector <double> MODEL::R0calc(const vector<double> &paramv)
{
	setup(paramv);
	calcprobin();
	
	for(auto& co : comp){
		co.infint.resize(data.nage);
		for(auto a = 0u; a < data.nage; a++){
			co.infint[a] = 0;
		
			auto kmax = co.trans.size();
			for(auto k = 0u; k < kmax; k++){
				auto tra = co.trans[k];
				
				double dt;
				switch(trans[tra].type){
				case infection_dist: dt = 0; break;
				case exp_dist: dt = paramv[trans[tra].param_mean]; break;
				case gamma_dist: dt = paramv[trans[tra].param_mean]; break;
				case lognorm_dist: dt = paramv[trans[tra].param_mean]; break;
				default: emsgEC("MODEL",9); break;
				}	
				if(kmax == 1) co.infint[a] += co.probin[a]*co.infectivity*dt;
				else co.infint[a] += co.probin[a]*co.infectivity*dt*co.prob[a][k];
			}
		}
	}
	
	vector <double> R0fac(ntimeperiod);
	for(auto& R0fa : R0fac) R0fa = 0;
		
	for(const auto& QQ : data.Q){
		auto timep = QQ.timep;
		auto qt = QQ.Qtenref;
		
		unsigned int co;
		for(co = 0u; co < comp.size(); co++) if(QQ.comp == comp[co].name) break;
		if(co == comp.size()) emsg("Compartment '"+QQ.comp+"' in '"+QQ.name+"' not recognised.");
		
		for(auto c = 0u; c < data.narea; c++){
			for(auto a = 0u; a < data.nage; a++){
				auto fac = comp[co].infint[a]*double(data.area[c].agepop[a])/data.popsize; 
				if(fac != 0){
					auto vi = c*data.nage + a;
					auto kmax = data.genQ.Qten[qt].ntof[vi];
					for(auto k = 0u; k < kmax; k++){
						auto cc = data.genQ.Qten[qt].tof[vi][k];
						for(auto dp = 0u; dp < data.ndemocatpos; dp++){
							auto aa = data.democatpos[dp][0];
							R0fac[timep] += fac*sus[dp]*areafac[cc]*data.genQ.Qten[qt].valf[vi][k][aa]*data.area[cc].pop[dp];
						}
					}
				}
			}
		}
	}
	
	vector <double> R0(details.nsettime);

	auto timep = 0u; 
	for(auto st = 0u; st < details.nsettime; st++){
		auto t = details.settime[st];
		while(timep < ntimeperiod && t > timeperiod[timep].tend) timep++;
		
		R0[st] = beta[st]*R0fac[timep];
	}
		
	return R0;
}

/// Makes proposal to compartmental paramters
void MODEL::compparam_prop(unsigned int samp, unsigned int burnin, vector <EVREF> &x, vector <vector <FEV> > &indev, vector <double> &paramv,
												   vector <float> &paramjumpxi, vector <unsigned int> &ntrxi,  vector <unsigned int> &nacxi, double &Pri)
{	
	timers.timecompparam -= clock();

	for(auto& tr : trans){
		if(tr.istimep == 0){
			tr.num.resize(data.nage); for(auto& num : tr.num) num = 0;
		}
		tr.numvisittot = 0; tr.dtsum = 0; tr.dtlist.clear();
	}

	for(const auto& xx : x){                  // Extracts values based on the event sequence
		auto i = xx.ind;
		auto dp = data.ind[i].dp;
		auto a = data.democatpos[dp][0];
				
		auto t = 0.0;
		for(const auto& ev : indev[i]){
			auto tra = ev.trans;
			if(trans[tra].istimep == 0){
				trans[tra].num[a]++;
	
				auto dt = ev.t-t;
				if(dt == 0) emsgEC("Model",10);
				t += dt;
				
				trans[tra].numvisittot++;
				switch(trans[tra].type){
				case exp_dist: trans[tra].dtsum += dt; break;
				case gamma_dist: case lognorm_dist: trans[tra].dtlist.push_back(dt); break;
				case infection_dist: break;
				default: emsgEC("model",99); break;
				}
			}
		}
	}

	if(settransprob(paramv) == 1) emsgEC("Model",11);
	auto Li_dt = likelihood_dt(paramv);
	auto Li_prob = likelihood_prob();

	for(auto th = 0u; th < param.size(); th++){
		if((param[th].type == distval_paramtype || param[th].type == branchprob_paramtype) && param[th].min != param[th].max){
			vector <double> paramst = paramv;	
		
			paramv[th] += normal(0,paramjumpxi[th]);               // Makes a change to a parameter
			
			auto Lp_prob = Li_prob;
			auto dL=0.0;
			auto Prp = Pri;
			auto flag=0u;
	
			if(paramv[th] < param[th].min || paramv[th] > param[th].max) flag = 0;
			else{
				if(param[th].type == branchprob_paramtype){
					if(settransprob(paramv) == 1) Lp_prob = -large;
					else Lp_prob = likelihood_prob();
				}
				
				dL = dlikelihood_dt(paramst,paramv);
				Prp = prior(paramv);
				
				auto al = exp(dL+ Lp_prob - Li_prob + Prp-Pri);
				if(ran() < al) flag = 1; else flag = 0;
			}
			
			ntrxi[th]++;
			if(flag == 1){
				Li_dt += dL;
				Li_prob = Lp_prob;
				Pri = Prp;
				
				nacxi[th]++;
				if(samp < burnin) paramjumpxi[th] *= 1.01;
			}
			else{
				paramv = paramst;
				if(samp < burnin) paramjumpxi[th] *= 0.995;
			}	
		}
	}
	
	if(settransprob(paramv) == 1) emsgEC("Model",12);
	setup(paramv);
	
	if(checkon == 1){
		double dd;
		dd = likelihood_dt(paramv)-Li_dt; if(dd*dd > tiny) emsgEC("Model",13);
		dd = likelihood_prob()-Li_prob; if(dd*dd > tiny) emsgEC("Model",14);
	}
	
	timers.timecompparam += clock();
}

/// Calculates the likelihood relating to branching probabilities
double MODEL::likelihood_prob()
{
	auto L = 0.0;
	for(const auto& co : comp){	
		auto kmax = co.trans.size();
		if(kmax > 1){
			for(auto k = 0u; k < kmax; k++){
				for(auto a = 0u; a < data.nage; a++){
					auto num = trans[co.trans[k]].num[a];
					if(num > 0) L += num*log(co.prob[a][k]);
				}
			}
		}
	}
	return L;
}
			
/// Calculates the likelihood for the timings of the transitions
double MODEL::likelihood_dt(vector <double> &paramv)
{
	auto L = 0.0;
	for(const auto& tr : trans){
		switch(tr.type){
		case exp_dist:
			{
				auto r = 1.0/paramv[tr.param_mean];
				L += tr.numvisittot*log(r) - r*tr.dtsum;
			}
			break;
		
		case gamma_dist:
			{
				auto mean = paramv[tr.param_mean], sd = paramv[tr.param_cv]*mean;
				for(auto dt : tr.dtlist) L += gammaprob(dt,mean*mean/(sd*sd),mean/(sd*sd));
			}
			break;
			
		case lognorm_dist:
			{
				auto mean_ns = paramv[tr.param_mean], cv_ns = paramv[tr.param_cv];
				auto sd = sqrt(log(1+cv_ns*cv_ns)), mean = log(mean_ns) - sd*sd/2;
				for(auto dt : tr.dtlist) L += lognormprob(dt,mean,sd*sd);
			}
			break;
		}
	}
	
	return L;
}

/// Calculates the change in likelihood for a given change in parameters
double MODEL::dlikelihood_dt(vector <double> &paramvi, vector <double> &paramvf)
{
	auto L = 0.0;
	for(const auto& tr : trans){
		switch(tr.type){
		case exp_dist:
			{
				auto ri = 1.0/paramvi[tr.param_mean], rf = 1.0/paramvf[tr.param_mean];
				if(ri != rf){
					L += tr.numvisittot*log(rf) - rf*tr.dtsum;
					L -= tr.numvisittot*log(ri) - ri*tr.dtsum;
				}
			}
			break;
		
		case gamma_dist:
			{
				auto meani = paramvi[tr.param_mean], sdi = paramvi[tr.param_cv]*meani;
				auto meanf = paramvf[tr.param_mean], sdf = paramvf[tr.param_cv]*meanf;
		
				if(meani != meanf || sdi != sdf){
					for(auto dt : tr.dtlist){
						L += gammaprob(dt,meanf*meanf/(sdf*sdf),meanf/(sdf*sdf));
						L -= gammaprob(dt,meani*meani/(sdi*sdi),meani/(sdi*sdi));
					}
				}				
			}
			break;
			
		case lognorm_dist:
			{
				auto mean_nsi = paramvi[tr.param_mean], cv_nsi = paramvi[tr.param_cv];
				auto mean_nsf = paramvf[tr.param_mean], cv_nsf = paramvf[tr.param_cv];
				if(mean_nsi != mean_nsf || cv_nsi != cv_nsf){
					auto sdi = sqrt(log(1+cv_nsi*cv_nsi)), meani = log(mean_nsi) - sdi*sdi/2;
					auto sdf = sqrt(log(1+cv_nsf*cv_nsf)), meanf = log(mean_nsf) - sdf*sdf/2;
		
					for(auto dt : tr.dtlist){
						L += lognormprob(dt,meanf,sdf*sdf);
						L -= lognormprob(dt,meani,sdi*sdi);
					}
				}
			}
			break;
		}
	}
	
	return L;
}

/// Outputs an event sequence (used for debugging)
void MODEL::oe(const string& name, const vector <FEV> &ev)
{
	cout << name << ":" << endl;
	for(auto e = 0u; e < ev.size(); e++){
		auto tra = ev[e].trans;
		cout << comp[trans[tra].from].name << "->" << comp[trans[tra].to].name << "  " << ev[e].t << endl;
	}
}

/// Determines if it is necessary to do mbp for the exisiting event sequence
unsigned int MODEL::dombpevents()
{
	for(auto th = 0u; th < param.size(); th++){
		if(parami[th] != paramp[th] && (param[th].type == distval_paramtype || param[th].type == branchprob_paramtype)) return 1;
	}
	return 0;
}
