// This file defines the compartmental model and is used to simulate compartmental dynamics after an individual becomes infected

#include <string>
#include <sstream>
#include <iostream>

#include "timers.hh"
#include "math.h"
#include "consts.hh"
#include "utils.hh"
#include "model.hh"
#include "utils.hh"
#include "toml11/toml.hpp"

using namespace std;

MODEL::MODEL(DATA &data) : data(data)
{
}

/// Defines the compartmental model
void MODEL::definemodel(unsigned int core, double /* period */, unsigned int /* popsize*/, const toml::basic_value<::toml::discard_comments, std::unordered_map, std::vector> &tomldata)
{
	unsigned int p, c, j, fi, tra, a;
	SPLINEP spl;
	PRIORCOMP pricomp;

	timeperiod = data.timeperiod; ntimeperiod = timeperiod.size();

	if(data.mode == MODE_SIM || data.mode == MODE_MULTISIM){
		if(tomldata.contains("params")){
			string name;
			double value;
			
			const auto paramsin = toml::find(tomldata,"params");
			for(j = 0; j < paramsin.size(); j++){
				const auto params = toml::find(paramsin,j);
				if(!params.contains("name")) emsg("Parameter must contain a 'name' definition.");
				name = toml::find<std::string>(params,"name");
				
				if(!params.contains("value")) emsg("Parameter must contain a 'value' definition.");
				value = toml::find<double>(params,"value");
				
				addparam(name,value,value);
			}
		}
		else{ emsg("The input file must contain parameter values through 'params'.");}
	}

	if(data.mode == MODE_INF){
		if(tomldata.contains("priors")){
			string name;
			double value;
			
			const auto paramsin = toml::find(tomldata,"priors");
			for(j = 0; j < paramsin.size(); j++){
				const auto params = toml::find(paramsin,j);
				if(!params.contains("name")) emsg("Parameter must contain a 'name' definition.");
				name = toml::find<std::string>(params,"name");

				if(params.contains("value")){
					value = toml::find<double>(params,"value");
					addparam(name,value,value);
				}
				else{
					if(!params.contains("type")) emsg("The prior must have a 'value' or a 'type'");
					
					string type = toml::find<std::string>(params,"type");
					if(type == "uniform"){
						if(!params.contains("min")) emsg("A uniform prior must contain a 'min' definition.");
						double min = toml::find<double>(params,"min");
				
						if(!params.contains("max")) emsg("A uniform prior must contain a 'max' definition.");
						double max = toml::find<double>(params,"max");
				
						addparam(name,min,max);
					}
					else emsg("The prior type '"+type+"' is not recognised.");
				}
			}
		}
		else{ emsg("The input file must contain 'priors'.");}
	}
	addparam("zero",tiny,tiny);

	if(tomldata.contains("comps")) {
		string name;
		double inf;
		
		const auto compsin = toml::find(tomldata,"comps");
		for(j = 0; j < compsin.size(); j++){
			const auto comps = toml::find(compsin,j);
			if(!comps.contains("name")) emsg("Compartment must contain a 'name' definition.");

			name = toml::find<std::string>(comps,"name");
			if(!comps.contains("inf")) emsg("Compartment must contain an 'inf' definition.");
			inf = toml::find<double>(comps,"inf");

			addcomp(name,inf); 
		}
	}
	else{ emsg("The input file must contain compartment definition through 'comps'");}

	if(tomldata.contains("trans")){
		string dist, mean="", cv="";
		unsigned int distval;
		
		const auto transin = toml::find(tomldata,"trans");
		for(j = 0; j < transin.size(); j++){
			const auto trans = toml::find(transin,j);
			
			if(!trans.contains("from")) emsg("Transition must contain a 'from' definition.");
			const auto from = toml::find<std::string>(trans, "from");
			
			if(!trans.contains("to")) emsg("Transition must contain a 'to' definition.");
			const auto to = toml::find<std::string>(trans, "to");
		
			if(!trans.contains("dist")) emsg("'dist' must be set in the transitions in 'trans'");
			
			const auto dist = toml::find<std::string>(trans, "dist");
			
			distval = UNSET;
			
			if(dist == "infection"){
				distval = INFECTION;
			}
			
			if(dist == "exp"){
				distval = EXP_DIST;
				if(!trans.contains("mean")) emsg("Transition distribution must contain a 'mean' definition.");
				mean = toml::find<std::string>(trans, "mean");
			}
			
			if(dist == "lognorm"){
				distval = LOGNORM_DIST;
				if(!trans.contains("mean")) emsg("Transition distribution must contain a 'mean' definition.");
				mean = toml::find<std::string>(trans, "mean");
				
				if(!trans.contains("cv")) emsg("Transition distribution must contain an 'cv' definition.");
				cv = toml::find<std::string>(trans, "cv");
			}
			
			if(dist == "gamma"){
				distval = GAMMA_DIST;
				if(!trans.contains("mean")) emsg("Transition distribution must contain a 'mean' definition.");
				mean = toml::find<std::string>(trans, "mean");
				
				if(!trans.contains("cv")) emsg("Transition distribution must contain an 'cv' definition.");
				cv = toml::find<std::string>(trans, "cv");
			}
			
			if(distval == UNSET) emsg("Distribution '"+dist+"' not recognised.");
			
			if(trans.contains("prob")) {
				const auto prob = toml::find<std::string>(trans, "prob");
				addtrans(from,to,prob,distval,mean,cv);  
			}
			else{
				addtrans(from,to,"",distval,mean,cv);  
			}
		}
	}
	else{ emsg("The input file must contain transition definitions through 'trans'.");}

	if(data.mode == MODE_INF){
		if(tomldata.contains("priorcomps")){
			string co;
		
			const auto prcomps = toml::find(tomldata,"priorcomps");
			for(j = 0; j < prcomps.size(); j++){
				const auto prcomp = toml::find(prcomps,j);
				
				if(!prcomp.contains("comp")) emsg("'priorcomps' must contain a 'comp' definition.");
				co = toml::find<std::string>(prcomp,"comp");
				c = 0; while(c < comp.size() && comp[c].name != co) c++;
				if(c == comp.size()) emsg("Cannot find '"+co+"' in 'priorcomps'");
				pricomp.comp = c;
				
				if(!prcomp.contains("value")) emsg("'priorcomps' must contain a 'value' definition.");
				double val = toml::find<double>(prcomp,"value");
				pricomp.value = val;
				
				if(!prcomp.contains("sd")) emsg("'priorcomps' must contain a 'sd' definition.");
				double sd = toml::find<double>(prcomp,"sd");
				pricomp.sd = sd;
		
				priorcomps.push_back(pricomp);
			}
		}
	}

	if(tomldata.contains("betaspline")) {
		const auto bespin = toml::find(tomldata,"betaspline");
		for(j = 0; j < bespin.size(); j++){
			const auto besp = toml::find(bespin,j);
			
			if(!besp.contains("param")) emsg("Beta spline definition must contain a 'param' definition.");
			const auto name = toml::find<std::string>(besp,"param");
			
			if(!besp.contains("time")) emsg("Beta spline definition must contain a 'time' definition.");
			const auto timstr = toml::find<string>(besp,"time");
			int tim = data.gettime(timstr) - data.start;
			
			if(j == 0 && tim != 0) emsg("The first beta spline point must be at t=0.");
			if(j == bespin.size()-1 && tim != (int)data.period) emsg("The last beta spline point must be at t=period.");
			if(tim < 0 || tim > (int)data.period) emsg("The beta spline points must be within the time period.");
			
			spl.t = tim;
			spl.param = findparam(name);
			betaspline.push_back(spl);
		}
	}

	if(tomldata.contains("phispline")) {
		const auto bespin = toml::find(tomldata,"phispline");
		for(j = 0; j < bespin.size(); j++){
			const auto besp = toml::find(bespin,j);
			
			if(!besp.contains("param")) emsg("Phi spline definition must contain a 'param' definition.");
			const auto name = toml::find<std::string>(besp,"param");
			
			if(!besp.contains("time")) emsg("Phi spline definition must contain a 'time' definition.");
			const auto timstr = toml::find<string>(besp,"time");
			int tim = data.gettime(timstr) - data.start;
			
			if(j == 0 && tim != 0) emsg("The first phi spline point must be at t=0.");
			if(j == bespin.size()-1 && tim != (int)data.period) emsg("The last phi spline point must be at t=period.");
			if(tim < 0 || tim > (int)data.period) emsg("The phi spline points must be within the time period.");
			
			spl.t = tim;
			spl.param = findparam(name);
			phispline.push_back(spl);
		}
	}
	
	sus_param.resize(data.ndemocat);
	for(c = 0; c < data.ndemocat; c++){
		for(fi = 0; fi < data.democat[c].value.size(); fi++){
			sus_param[c].push_back(findparam(data.democat[c].param[fi]));
		}
	}

	for(c = 0; c < data.ncovar; c++){
		covar_param.push_back(findparam(data.covar[c].param));
	}

	for(p = 0; p < param.size(); p++){
		if(param[p].used == 0) emsg("Parameter '"+param[p].name+"' not in the model.");
	}

	ntr = 0; nac = 0;
	parami.resize(param.size()); paramp.resize(param.size());

	paramval.resize(param.size());
	priorsamp();
 
	beta.resize(data.nsettime);	phi.resize(data.nsettime);
	
	setup(paramval);
	
	for(tra = 0; tra < trans.size(); tra++){
		switch(trans[tra].type){
		case EXP_DIST: param[trans[tra].param1].type = 1; break;
		case GAMMA_DIST: case LOGNORM_DIST: param[trans[tra].param1].type = 1; param[trans[tra].param2].type = 1; break;
		}
		
		if(trans[tra].istimep == 0){
			if(trans[tra].probparam.size() > 0){
				for(a = 0; a < data.nage; a++) param[trans[tra].probparam[a]].type = 2; 
			}
		}
	}
	
	if(core == 0){
		cout << endl;                                               // Outputs a summary of the model
		if(data.mode == MODE_SIM || data.mode == MODE_MULTISIM){
			cout << "Parameters:" << endl;
			for(p = 0; p < param.size()-1; p++){
				cout << "  " << param[p].name << " = " << param[p].min << endl;
			}
		}
		else{
			cout << "Priors:" << endl;
			for(p = 0; p < param.size()-1; p++){
				cout << "  " << param[p].name << " = ";
				if(param[p].min ==  param[p].max) cout << param[p].min << endl;
				else cout << "Uniform(" << param[p].min << " - " << param[p].max << ")" << endl;
			}
		}
		cout << endl;
			
		cout << "Compartments:" << endl; 
		for(c = 0; c < comp.size(); c++){
			cout << "  " << comp[c].name << "  Infectivity: " << comp[c].infectivity << endl; 			
		}
		cout << endl;
		
		cout << "Transitions:" << endl; 
		for(tra = 0; tra < trans.size(); tra++){
			cout << "  " << comp[trans[tra].from].name << " → " << comp[trans[tra].to].name;
			if(trans[tra].probparam.size() > 0){
				cout << "  with probability ";
				for(j = 0; j < trans[tra].probparam.size(); j++){
					if(j > 0) cout << ", ";
					cout << param[trans[tra].probparam[j]].name;
				}
			}
			
			switch(trans[tra].type){
				case INFECTION: 
					cout << " Infection";
					break;
				case EXP_DIST:
					cout << " Exponential  mean=" << param[trans[tra].param1].name;
					break;
				case GAMMA_DIST:
					cout << " Gamma mean=" << param[trans[tra].param1].name 
							 << " cv=" << param[trans[tra].param2].name; 
					break;
				case LOGNORM_DIST:
					cout << " Lognormal mean=" << param[trans[tra].param1].name 
							 << " cv=" << param[trans[tra].param2].name; 
					break;
				default:
					break;
			}
			cout << endl;
		}
		cout << endl;
	}
}

/// Sets up the model with a set of parameters
unsigned int MODEL::setup(vector <double> &paramv)
{
	paramval = paramv;
	if(settransprob() == 1) return 1;
	
	timevariation();
	setsus();
	setarea();
	return 0;
}

/// Copies values used for the initial state (MBPs)
void MODEL::copyi()
{
	unsigned int c;
	
	parami = paramval; betai = beta; phii = phi; susi = sus; areafaci = areafac;
	for(c = 0; c < comp.size(); c++) comp[c].probi = comp[c].prob;
}
	
/// Copies values used for the proposed state (MBPs)
void MODEL::copyp()
{
	unsigned int c;
		
	paramp = paramval; betap = beta; phip = phi; susp = sus; areafacp = areafac;
	for(c = 0; c < comp.size(); c++) comp[c].probp = comp[c].prob;
}

/// Adds in the tensor Q to the model
void MODEL::addQ()
{
	unsigned int q, qi, qf, c, ci, cf, tra, timep, timepi, timepf;
	string compi, compf;
	DQINFO dq;
	
	for(c = 0; c < comp.size(); c++) addtrans(comp[c].name,comp[c].name,"",TIMEP_DIST,"","");  
	
	for(q = 0; q < data.Q.size(); q++){
		for(c = 0; c < comp.size(); c++) if(data.Q[q].comp == comp[c].name) break;
		if(c == comp.size()) emsg( "Compartment "+data.Q[q].comp+" not recognised.");
	}
	
	dq.q.resize(2); dq.fac.resize(2);
	for(tra = 0; tra < trans.size(); tra++){
		ci = trans[tra].from;
		cf = trans[tra].to;
		compi = comp[ci].name;
		compf = comp[cf].name;
	
		for(timep = 0; timep < ntimeperiod; timep++){
			timepi = timep; timepf = timep;
			if(compi == compf) timepf++;
			
			if(timepf < ntimeperiod){
				qi = 0; while(qi < data.Q.size() && !(data.Q[qi].comp == compi && data.Q[qi].timep == timepi)) qi++;
				if(qi == data.Q.size()) qi = UNSET;
				
				qf = 0; while(qf < data.Q.size() && !(data.Q[qf].comp == compf && data.Q[qf].timep == timepf)) qf++;
				if(qf == data.Q.size()) qf = UNSET;
				
				if(qi == UNSET && qf == UNSET){
					trans[tra].DQ.push_back(UNSET);
				}
				else{
					trans[tra].DQ.push_back(DQ.size());
					dq.q[0] = qi; dq.fac[0] = -comp[ci].infectivity;
					dq.q[1] = qf; dq.fac[1] = comp[cf].infectivity; 
					DQ.push_back(dq);
				}
			}
		}
	}
}

/// Randomly samples the initial parameter values from the prior (which are uniform distributions
void MODEL::priorsamp()
{
	unsigned int th;

	for(th = 0; th < param.size(); th++){	
		paramval[th] = param[th].min + ran()*(param[th].max - param[th].min);
	}
	
//	paramval[1] = 0.09;
	
//	for(th = 0; th < param.size(); th++) cout << "paramval[" << th <<"] = " << paramval[th] << ";" << endl;
	/*
	paramval[0] = 0.35;
paramval[1] = 0.06;
paramval[2] = 2.9;
paramval[3] = 4.5;
paramval[4] = 0.533;
paramval[5] = 8;
paramval[6] = 1.5;
paramval[7] = 5;
paramval[8] = 20;
paramval[9] = 0.5;
paramval[10] = 20;
paramval[11] = 0.5;
paramval[12] = 20;
paramval[13] = 0.5;
paramval[14] = 0.75;
paramval[15] = 0.5;
paramval[16] = 0.12;
paramval[17] = 0.33;
paramval[18] = 0.2;
paramval[19] = 1e-08;
*/
}

/// Gets a parameter value
double MODEL::getparam(string name)
{	
	return paramval[findparam(name)];
}

/// Gets the infectivity of a compartment
double MODEL::getinfectivity(string name)
{
	unsigned int c;

	c = 0; while(c < comp.size() && comp[c].name != name) c++;
	if(c == comp.size()) emsg("Cannot find compartment '"+name+"'");
	return comp[c].infectivity;
}

/// Adds a compartment to the model
void MODEL::addcomp(string name, double infectivity)
{
	COMP co;
	co.name = name;
	co.infectivity = infectivity;
	comp.push_back(co);	
}

/// Finds a parameter from a string
unsigned int MODEL::findparam(string name)
{
	unsigned int p, pmax;
	
	p = 0; pmax = param.size(); while(p < pmax && name != param[p].name) p++;
	if(p == pmax) emsg("Cannot find parameter '"+name+"'");	
	param[p].used = 1;
	
	return p;
}

/// Adds a parameter to the model
void MODEL::addparam(string name, double min, double max)
{
	PARAM par;
	par.name = name; par.min = min; par.max = max; par.ntr = 0; par.nac = 0; par.jump = 0.5*(min+max)/10; if(par.jump == 0) par.jump = 0.1;
	par.used = 0; par.type = 0;

	param.push_back(par);
}

/// Adds a transition to the model
void MODEL::addtrans(string from, string to, string prpar, unsigned int type, string param1, string param2)
{
	unsigned int c, cmax, j;
	
	TRANS tr;
	
	c = 0; cmax = comp.size(); while(c < cmax && from != comp[c].name) c++;
	if(c == cmax) emsg("Cannot find compartment");
	tr.from = c;	
	
	if(from != to){ tr.istimep = 0; comp[c].trans.push_back(trans.size());}
	else{ tr.istimep = 1; comp[c].transtimep = trans.size();}
		
	c = 0; cmax = comp.size(); while(c < cmax && to != comp[c].name) c++;
	if(c == cmax) emsg("Cannot find compartment");	
	tr.to = c;
	
	if(prpar != ""){
		vector <string> probparam = split(prpar,',');
	
		if(probparam.size() != data.nage) emsg("Wrong number of parameters in expression '"+prpar+"'.");
		for(j = 0; j < probparam.size(); j++){
			tr.probparam.push_back(findparam(probparam[j]));
		}
	}
	
	tr.type = type;
	if(param1 != "") tr.param1 = findparam(param1);	
	if(param2 != "") tr.param2 = findparam(param2);
	
	trans.push_back(tr);
}
	
/// Generates the time variation in beta and phi from the parameters
void MODEL::timevariation()
{
  unsigned int s;
	int p;
	double t, fac;
	
	/*  // This uses a cubic spline for beta
	double t, fac, dt, a[n+1], b[n], c[n+1], d[n], h[n], alpha[n], l[n+1], mu[n+1], z[n+1];
	
		int n = nspline-1;
		for(p = 0; p <= n; p++) a[p] = log(paramval[p]);
		for(p = 0; p < n; p++) h[p] = splinet[p+1]-splinet[p];
		for(p = 1; p < n; p++) alpha[p] = (3/h[p])*(a[p+1]-a[p]) - (3/h[p-1])*(a[p]-a[p-1]);

		l[0]=1; mu[0]=0; z[0]=0;
		for(p = 1; p < n; p++){
			l[p] = 2*(splinet[p+1]-splinet[p-1]) - h[p-1]*mu[p-1];
			mu[p] = h[p]/l[p];
			z[p] = (alpha[p]-h[p-1]*z[p-1])/l[p];
		}
		l[n] = 1; z[n] = 0; c[n] = 0;
		for(p = n-1; p >= 0; p--){
			c[p] = z[p]-mu[p]*c[p+1];
			b[p] = (a[p+1]-a[p])/h[p] - h[p]*(c[p+1]+2*c[p])/3;
			d[p] = (c[p+1]-c[p])/(3*h[p]);
		}
		
		p = 0;
		for(s = 0; s < data.nsettime; s++){		
			t = double((s+0.5)*data.period)/data.nsettime;
			while(p < int(nspline)-1 && t > splinet[p+1]) p++;
			
			dt = t-splinet[p];	
			beta[s] = exp(a[p]+ b[p]*dt + c[p]*dt*dt + d[p]*dt*dt*dt);
		}
		*/
	
  // This uses a linear spline for beta
	p = 0;
	for(s = 0; s < data.nsettime; s++){	
		t = double((s+0.5)*data.period)/data.nsettime;
		
		while(p < int(betaspline.size())-1 && t > betaspline[p+1].t) p++;
		
		fac = (t-betaspline[p].t)/(betaspline[p+1].t-betaspline[p].t);
		beta[s] = (paramval[betaspline[p].param]*(1-fac) + paramval[betaspline[p+1].param]*fac);
	}
	
	// This uses a linear spline for phi
	p = 0;
	for(s = 0; s < data.nsettime; s++){	
		t = double((s+0.5)*data.period)/data.nsettime;
		
		while(p < int(phispline.size())-1 && t > phispline[p+1].t) p++;
		
		fac = (t-phispline[p].t)/(phispline[p+1].t-phispline[p].t);
		phi[s] = (paramval[phispline[p].param]*(1-fac) + paramval[phispline[p+1].param]*fac)/data.popsize;
	}
}

/// Sets the transition probabilies based on the parameters
unsigned int MODEL::settransprob()
{
	unsigned int c, k, kmax, tra, a;
	int p;
	double sum, prob;

	if(checkon == 2){
		for(c = 0; c < comp.size(); c++){
			cout << "Comp " << comp[c].name << "  "; cout << comp[c].transtimep << ",  ";
			for(k = 0; k <  comp[c].trans.size(); k++) cout << comp[c].trans[k] << ", ";
			cout << endl;
		}
	
		for(tra = 0; tra < trans.size(); tra++){
			cout << tra << " " << comp[trans[tra].from].name << "->" << comp[trans[tra].to].name << " trans" << endl;
			for(k = 0; k < trans[tra].DQ.size(); k++) {
				cout <<  trans[tra].DQ[k] << ", ";
			}
			cout << "DQ" << endl;
		}
	}
	
	for(c = 0; c < comp.size(); c++){
		kmax = comp[c].trans.size();
		
		if(kmax > 1){
			comp[c].prob.resize(data.nage);
			comp[c].probsum.resize(data.nage);
			for(a = 0; a < data.nage; a++){
				comp[c].prob[a].resize(kmax);
				comp[c].probsum[a].resize(kmax);
				
				sum = 0;
				for(k = 0; k < kmax-1; k++){
					p = trans[comp[c].trans[k]].probparam[a]; if(p == UNSET) emsg("model: EC1a");		
					prob = paramval[p]; comp[c].prob[a][k] = prob; sum += prob; if(prob < 0) return 1;
				} 
				prob = 1-sum; comp[c].prob[a][kmax-1] = prob; if(prob < 0) return 1;
				
				sum = 0;
				for(k = 0; k < kmax; k++){
					sum += comp[c].prob[a][k];
					comp[c].probsum[a][k] = sum;
				}
			}
		}
	}
	
	return 0;
}

/// This simulates from the model and generates an event list
void MODEL::simmodel(vector <FEV> &evlist, unsigned int i, unsigned int c, double t)
{
	unsigned int k, kmax, tra, timep, a;
	double mean, sd, mean_ns, cv_ns, z, dt;
	FEV ev;
	vector <double> probsum;

	timep = 0; while(timep < ntimeperiod && t > timeperiod[timep].tend) timep++;

	ev.ind = i; ev.timep = timep; 
		
	a = data.democatpos[data.ind[i].dp][0];
		
	if(c == 0){
		evlist.clear();	
		tra = 0;
		ev.trans = tra; ev.ind = i; ev.t = t; 
		evlist.push_back(ev);
	
		c = trans[tra].to;
	}
	 
	do{
		kmax = comp[c].trans.size();
		if(kmax == 0) break;
		
		if(kmax == 1) tra = comp[c].trans[0];
		else{
			z = ran(); k = 0; while(k < kmax && z > comp[c].probsum[a][k]) k++;
			if(k == kmax) emsg("Model: EC3");
			tra = comp[c].trans[k];
		}
		
		switch(trans[tra].type){
		case EXP_DIST:
			dt = -log(ran())*paramval[trans[tra].param1];
			break;
		
		case GAMMA_DIST:
			mean = paramval[trans[tra].param1]; sd = paramval[trans[tra].param2]*mean;
			dt = gammasamp(mean*mean/(sd*sd),mean/sd*sd);
			break;
			
		case LOGNORM_DIST:
			mean_ns = paramval[trans[tra].param1]; cv_ns = paramval[trans[tra].param2];
			sd = sqrt(log((1+cv_ns*cv_ns))); mean = log(mean_ns) - sd*sd/2;
			dt = exp(normal(mean,sd));
			break;
			
		default: emsg("MODEL: EC2b"); break;
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
	unsigned int c, k, kmax, tra, e, emax, i, p, p2, timep, a;
	double meani, sdi, meanp, sdp, mean_nsi, cv_nsi, mean_nsp, cv_nsp, z, ti, tp, dt, sum, dif;
	FEV ev;
	vector <double> sumst;
	double dtnew=0;
	
	evlistp.clear();
	
	ev = evlisti[0];
	tra = ev.trans; ti = ev.t; i = ev.ind; timep = ev.timep;
	evlistp.push_back(ev);
	
	a = data.democatpos[data.ind[i].dp][0];
		
	c = trans[tra].to;
	
	tp = ti;
	emax = evlisti.size();
	for(e = 1; e < emax; e++){
		tra = evlisti[e].trans;
		if(trans[tra].istimep == 0){		
			dt = evlisti[e].t - ti;
			ti = evlisti[e].t;
			
			kmax = comp[c].trans.size();
			if(kmax == 0) break;
			
			if(kmax > 1){
				k = 0; while(k < kmax && tra != comp[c].trans[k]) k++;
				if(k == kmax) emsg("model: EC32");
				
				if(comp[c].probp[a][k] < comp[c].probi[a][k]){  // Looks at switching to another branch
					if(ran() < 1 - comp[c].probp[a][k]/comp[c].probi[a][k]){
						sum = 0; sumst.clear();
						for(k = 0; k < kmax; k++){
							dif = comp[c].probp[a][k] - comp[c].probi[a][k];
							if(dif > 0) sum += dif;
							sumst.push_back(sum);
						}
						
						z = ran()*sum; k = 0; while(k < kmax && z > sumst[k]) k++;
						if(k == kmax) emsg("Model: EC12");
						tra = comp[c].trans[k];
					}
				}	
			}
			
			switch(trans[tra].type){
			case EXP_DIST:
				p = trans[tra].param1;
				dtnew = dt*paramp[p]/parami[p];
				break;
			
			case GAMMA_DIST:
				emsg("model: EC9");
				break;
				
			case LOGNORM_DIST:
				p = trans[tra].param1; p2 = trans[tra].param2;
				
				mean_nsi = parami[p]; cv_nsi = parami[p2];
				mean_nsp = paramp[p]; cv_nsp = paramp[p2];
				
				if(mean_nsi == mean_nsp && cv_nsi == cv_nsp) dtnew = dt;
				else{
					sdi = sqrt(log(1+cv_nsi*cv_nsi)); meani = log(mean_nsi) - sdi*sdi/2;
					sdp = sqrt(log(1+cv_nsp*cv_nsp)); meanp = log(mean_nsp) - sdp*sdp/2;
					dtnew = exp(meanp + (log(dt) - meani)*sdp/sdi); 
				}
				break;
				
			default:
				emsg("MODEL: EC2b");
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
	
	if(comp[c].trans.size() != 0)	simmodel(evlistp,i,c,tp);
}
  
/// Defines the relative susceptibility of individuals
void MODEL::setsus()     
{
	unsigned int dp, c, j;
	double val;

	sus.resize(data.ndemocatpos);
	for(dp = 0; dp < data.ndemocatpos; dp++){
		val = 1;
		for(c = 0; c < data.ndemocat; c++){
			j = data.democatpos[dp][c];
			val *= exp(paramval[sus_param[c][j]]);
		}
		sus[dp] = val;
	}
}

/// Defines the relative transmission rate for different areas
void MODEL::setarea() 
{
	unsigned int c, j;
	double sum;
	
	areafac.resize(data.narea);
	for(c = 0; c < data.narea; c++){
		sum = 0; for(j = 0; j < data.ncovar; j++) sum += paramval[covar_param[j]]*data.area[c].covar[j];
		areafac[c] = exp(sum);
	}
}

/// Checks that the transition and population data is correct
void MODEL::checkdata()
{
	unsigned int td, pd, md, c, tra;
	string from, to, compstr;
	TRANS tr;
	
	for(td = 0; td < data.transdata.size(); td++){
		from = data.transdata[td].fromstr; to = data.transdata[td].tostr; 
		for(tra = 0; tra < trans.size(); tra++){ if(comp[trans[tra].from].name == from && comp[trans[tra].to].name == to) break;}
		
		if(tra == trans.size()) emsg("Cannot find the transition "+from+"→"+to+" in file '"+data.transdata[td].file+"'.");
		data.transdata[td].trans = tra;
	}
	
	for(pd = 0; pd < data.popdata.size(); pd++){
		compstr = data.popdata[pd].compstr;
		for(c = 0; c < comp.size(); c++){ if(comp[c].name == compstr) break;}
		if(c == comp.size()) emsg("Cannot find the compartment '"+compstr+"' in file '"+data.popdata[pd].file+"'.");
		data.popdata[pd].comp = c;
	}
	
	for(md = 0; md < data.margdata.size(); md++){
		from = data.margdata[md].fromstr; to = data.margdata[md].tostr; 
		for(tra = 0; tra < trans.size(); tra++){ if(comp[trans[tra].from].name == from && comp[trans[tra].to].name == to) break;}
		
		if(tra == trans.size()) emsg("Cannot find the transition "+from+"→"+to+" in file '"+data.margdata[md].file+"'.");
		data.margdata[md].trans = tra;
	}
}

/// Calculates the prior probability
double MODEL::prior()
{
	unsigned int pc, a, c;
	double Pr, prob;
	
	Pr = 0;
	if(priorcomps.size() > 0){
		calcprobin();
		for(pc = 0; pc < priorcomps.size(); pc++){
			c = priorcomps[pc].comp;
			prob = 0; for(a = 0; a < data.nage; a++) prob += data.agedist[a]*comp[c].probin[a];
			Pr += normalprob(prob,priorcomps[pc].value,priorcomps[pc].sd*priorcomps[pc].sd);
		}
	}
	
	return Pr;
}
	
/// Calculate compartmental probabilities
void MODEL::calcprobin()
{
	unsigned int c, a, k, j;
	double prob;
	vector <unsigned int> cst, kst;
	vector <double> probst;
	
	for(c = 0; c < comp.size(); c++){
		comp[c].probin.resize(data.nage);
		for(a = 0; a < data.nage; a++) comp[c].probin[a] = 0;
	}
	
	for(a = 0; a < data.nage; a++){
		prob = 1;
		c = 0;
		do{
			comp[c].probin[a] += prob;
			
			if(comp[c].trans.size() == 0){
				if(cst.size() == 0) break;
				
				do {
					j = cst.size()-1;
					c = cst[j];
					prob = probst[j]; 
					kst[j]++;
					k = kst[j];
					if(k < comp[c].trans.size()) break;
					
					cst.pop_back();
					probst.pop_back();
					kst.pop_back();
				} while(cst.size() > 0);
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
		
		if(cst.size() > 0 || kst.size() > 0 || probst.size() > 0) emsg("Model: EC54");
	}
}

/// Calculates R0
vector <double> MODEL::R0calc()
{
	unsigned int c, cc, a, aa, k, kmax, st, timep, q, qt, co, vi, dp, tra;
	double t, dt, fac;
	vector <double> R0;
	vector <double> R0fac;

	setup(paramval);
	calcprobin();
	
	for(c = 0; c < comp.size(); c++){
		comp[c].infint.resize(data.nage);
		for(a = 0; a < data.nage; a++){
			comp[c].infint[a] = 0;
		
			kmax = comp[c].trans.size();
			for(k = 0; k < kmax; k++){
				tra = comp[c].trans[k];
				
				switch(trans[tra].type){
				case INFECTION: dt = 0; break;
				case EXP_DIST: dt = paramval[trans[tra].param1]; break;
				case GAMMA_DIST: dt = paramval[trans[tra].param1]; break;
				case LOGNORM_DIST: dt = paramval[trans[tra].param1]; break;
				default:
					//dt=0;
					emsg("MODEL: EC56"); break;
				}	
				if(kmax == 1) comp[c].infint[a] += comp[c].probin[a]*comp[c].infectivity*dt;
				else comp[c].infint[a] += comp[c].probin[a]*comp[c].infectivity*dt*comp[c].prob[a][k];
			}
		}
	}
	
	R0fac.resize(ntimeperiod);
	for(timep = 0; timep < ntimeperiod; timep++) R0fac[timep] = 0;
		
	for(q = 0; q < data.Q.size(); q++){
		timep = data.Q[q].timep;
		qt = data.Q[q].Qtenref;
		
		for(co = 0; co < comp.size(); co++) if(data.Q[q].comp == comp[co].name) break;
		if(co == comp.size()) emsg( "Compartment "+data.Q[q].comp+" not recognised.");
		
		for(c = 0; c < data.narea; c++){
			for(a = 0; a < data.nage; a++){
				fac = comp[co].infint[a]*double(data.area[c].agepop[a])/data.popsize; 
				if(fac != 0){
					vi = c*data.nage + a;
					kmax = data.genQ.Qten[qt].ntof[vi];
					for(k = 0; k < kmax; k++){
						cc = data.genQ.Qten[qt].tof[vi][k];
						for(dp = 0; dp < data.ndemocatpos; dp++){
							aa = data.democatpos[dp][0];
							R0fac[timep] += fac*sus[dp]*areafac[cc]*data.genQ.Qten[qt].valf[vi][k][aa]*data.area[cc].pop[dp];
						}
					}
				}
			}
		}
	}
	
	R0.resize(data.nsettime);
	timep = 0; 
	for(st = 0; st < data.nsettime; st++){
		t = data.settime[st];
		while(timep < ntimeperiod && t > timeperiod[timep].tend) timep++;
		
		R0[st] = beta[st]*R0fac[timep];
	}
		
	return R0;
}

/// Makes proposal to compartmental paramters
void MODEL::compparam_prop(unsigned int samp, unsigned int burnin, vector <EVREF> &x, vector <vector <FEV> > &indev, vector <double> &paramv,
												   vector <float> &paramjumpxi, vector <unsigned int> &ntrxi,  vector <unsigned int> &nacxi, double &Pri)
{	
	unsigned int a, i, j, jmax, dp, e, emax, tra, th, loop, loopmax = 1, flag;
	double t, dt, Li_dt, Li_prob, Lp_prob, Prp, al, dL, dd;
	vector <double> paramst;
	
	timers.timecompparam -= clock();

	for(tra = 0; tra < trans.size(); tra++){
		if(trans[tra].istimep == 0){
			trans[tra].num.resize(data.nage); for(a = 0; a < data.nage; a++) trans[tra].num[a] = 0;
		}
	}
		
	for(tra = 0; tra < trans.size(); tra++){
		trans[tra].numvisittot = 0; trans[tra].dtsum = 0; trans[tra].dtlist.clear();
	}

	jmax = x.size();                                                        // Extracts values based on the event sequence
	for(j = 0; j < jmax; j++){
		i = x[j].ind;
		dp = data.ind[i].dp;
		a = data.democatpos[dp][0];
				
		t = 0;
		emax = indev[i].size();
		for(e = 0; e < emax; e++){
			tra = indev[i][e].trans;
			if(trans[tra].istimep == 0){
				trans[tra].num[a]++;
	
				dt = indev[i][e].t-t;
				if(dt == 0) emsg("Model: EC78");
				t += dt;
				
				trans[tra].numvisittot++;
				switch(trans[tra].type){
				case EXP_DIST: trans[tra].dtsum += dt; break;
				case GAMMA_DIST: case LOGNORM_DIST: trans[tra].dtlist.push_back(dt); break;
				}
			}
		}
	}
	
	Li_dt = likelihood_dt(paramv);
	Li_prob = likelihood_prob();
	paramval = paramv; if(settransprob() == 1) emsg("Model: EC32");
	
	for(loop = 0; loop < loopmax; loop++){
		for(th = 0; th < param.size(); th++){
			if(param[th].type > 0 && param[th].min != param[th].max){
				paramst = paramv;	
				paramv[th] += normal(0,paramjumpxi[th]);               // Makes a change to a parameter
				
        Lp_prob = Li_prob;
				if(paramv[th] < param[th].min || paramv[th] > param[th].max) flag = 0;
				else{
					if(param[th].type == 2){
						paramval = paramv;
						if(settransprob() == 1) Lp_prob = -large;
						else Lp_prob = likelihood_prob();
					}
					
					dL = dlikelihood_dt(paramst,paramv);
					
					Prp = prior();
					
					al = exp(dL+ Lp_prob - Li_prob + Prp-Pri);
					if(ran() < al) flag = 1; else flag = 0;
				}
				
				ntrxi[th]++;
				if(flag == 1){
					Li_dt += dL;
					Li_prob = Lp_prob;
					Pri = prior();
					
					nacxi[th]++;
					if(samp < burnin) paramjumpxi[th] *= 1.01;
				}
				else{
					paramv = paramst;
					if(samp < burnin) paramjumpxi[th] *= 0.995;
				}
			}
		}
	}
	
	paramval = paramv; if(settransprob() == 1) emsg("Model: EC32");
		
	if(checkon == 1){
		dd = likelihood_dt(paramv)-Li_dt; if(dd*dd > tiny) emsg("Model: EC45");
		dd = likelihood_prob()-Li_prob; if(dd*dd > tiny) emsg("Model: EC46");
	}
	
	timers.timecompparam += clock();
}

/// Calculates the likelihood relating to branching probabilities
double MODEL::likelihood_prob()
{
	unsigned int c, k, kmax, a, num;
	double L;
	
	L = 0;
	for(c = 0; c < comp.size(); c++){	
		kmax = comp[c].trans.size();
		if(kmax > 1){
			for(k = 0; k < kmax; k++){
				for(a = 0; a < data.nage; a++){
					num = trans[comp[c].trans[k]].num[a];
					if(num > 0) L += num*log(comp[c].prob[a][k]);
				}
			}
		}
	}
	return L;
}
			
/// Calculates the likelihood for the timings of the transitions
double MODEL::likelihood_dt(vector <double> &paramv)
{
	unsigned int tra, j, jmax;
	double L, r, mean, sd, mean_ns, cv_ns;
	
	L = 0;
	for(tra = 0; tra < trans.size(); tra++){
		switch(trans[tra].type){
		case EXP_DIST:
			r = 1.0/paramv[trans[tra].param1];
			L += trans[tra].numvisittot*log(r) - r*trans[tra].dtsum;
			break;
		
		case GAMMA_DIST:
			mean = paramv[trans[tra].param1]; sd = paramv[trans[tra].param2]*mean;
			jmax = trans[tra].dtlist.size();
			for(j = 0; j < jmax; j++) L += gammaprob(trans[tra].dtlist[j],mean*mean/(sd*sd),mean/(sd*sd));
			break;
			
		case LOGNORM_DIST:
			mean_ns = paramval[trans[tra].param1]; cv_ns = paramval[trans[tra].param2];
			sd = sqrt(log(1+cv_ns*cv_ns)); mean = log(mean_ns) - sd*sd/2;
			jmax = trans[tra].dtlist.size();
			
			double ss = 0; for(j = 0; j < jmax; j++) ss += trans[tra].dtlist[j];
			for(j = 0; j < jmax; j++){
				if(trans[tra].dtlist[j] == 0) emsg("Time zero");
				L += lognormprob(trans[tra].dtlist[j],mean,sd*sd);
			}
			break;
		}
	}
	
	return L;
}

/// Calculates the change in likelihood for a given change in parameters
double MODEL::dlikelihood_dt(vector <double> &paramvi, vector <double> &paramvf)
{
	unsigned int tra, j, jmax;
	double L, ri, meani, sdi, rf, meanf, sdf, mean_nsi, cv_nsi, mean_nsf, cv_nsf;
	
	L = 0;
	for(tra = 0; tra < trans.size(); tra++){
		switch(trans[tra].type){
		case EXP_DIST:
			ri = 1.0/paramvi[trans[tra].param1]; rf = 1.0/paramvf[trans[tra].param1];
			if(ri != rf){
				L += trans[tra].numvisittot*log(rf) - rf*trans[tra].dtsum;
				L -= trans[tra].numvisittot*log(ri) - ri*trans[tra].dtsum;
			}
			break;
		
		case GAMMA_DIST:
			meani = paramvi[trans[tra].param1]; sdi = paramvi[trans[tra].param2]*meani;
			meanf = paramvf[trans[tra].param1]; sdf = paramvf[trans[tra].param2]*meanf;
	
			if(meani != meanf || sdi != sdf){
				jmax = trans[tra].dtlist.size();
				for(j = 0; j < jmax; j++){
					L += gammaprob(trans[tra].dtlist[j],meanf*meanf/(sdf*sdf),meanf/(sdf*sdf));
					L -= gammaprob(trans[tra].dtlist[j],meani*meani/(sdi*sdi),meani/(sdi*sdi));
				}
			}				
			break;
			
		case LOGNORM_DIST:
			mean_nsi = paramvi[trans[tra].param1]; cv_nsi = paramvi[trans[tra].param2];
			mean_nsf = paramvf[trans[tra].param1]; cv_nsf = paramvf[trans[tra].param2];
			if(mean_nsi != mean_nsf || cv_nsi != cv_nsf){
				sdi = sqrt(log(1+cv_nsi*cv_nsi)); meani = log(mean_nsi) - sdi*sdi/2;
				sdf = sqrt(log(1+cv_nsf*cv_nsf)); meanf = log(mean_nsf) - sdf*sdf/2;
	
				jmax = trans[tra].dtlist.size();
				for(j = 0; j < jmax; j++){
					L += lognormprob(trans[tra].dtlist[j],meanf,sdf*sdf);
					L -= lognormprob(trans[tra].dtlist[j],meani,sdi*sdi);
				}
			}
			break;
		}
	}
	
	return L;
}

/// Outputs an event sequence (used for debugging)
void MODEL::oe(string name, vector <FEV> &ev)
{
	unsigned int e, tra;
	
	cout << name << ":" << endl;
	for(e = 0; e < ev.size(); e++){
		tra = ev[e].trans;
		cout << comp[trans[tra].from].name << "->" << comp[trans[tra].to].name << "  " << ev[e].t << endl;
	}
}

/// Determines if it is necessary to do mbp for the exisiting event sequence
unsigned int MODEL::dombpevents()
{
	unsigned int th;
	for(th = 0; th < param.size(); th++){
		if(parami[th] != paramp[th] && param[th].type > 0) return 1;
	}
	return 0;
}
