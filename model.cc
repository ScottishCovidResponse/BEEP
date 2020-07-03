// This file defines the compartmental model and is used to simulate compartmental dynamics after an individual becomes infected

#include <string>
#include <sstream>
#include <iostream>

#include "math.h"
#include "consts.hh"
#include "utils.hh"
#include "model.hh"
#include "utils.hh"

using namespace std;

MODEL::MODEL(DATA &data) : data(data)
{
}

/// Defines the compartmental model
void MODEL::definemodel(unsigned int core, double period, unsigned int popsize, const toml::basic_value<::toml::discard_comments, std::unordered_map, std::vector> &tomldata)
{
	unsigned int p, c, t, j, fi;
	SPLINEP spl;
	
	timeperiod = data.timeperiod; ntimeperiod = timeperiod.size();

	if(data.mode == MODE_SIM){
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
				
				addparam(name,value,value,value);
			}
		}
		else{ emsg("The input file must contain parameter values through 'params'.");}
	}
	
	if(data.mode != MODE_SIM){
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
					addparam(name,value,value,value);
				}
				else{
					if(!params.contains("type")) emsg("The prior must have a 'value' or a 'type'");
					
					string type = toml::find<std::string>(params,"type");
					if(type == "uniform"){
						if(!params.contains("min")) emsg("A uniform prior must contain a 'min' definition.");
						double min = toml::find<double>(params,"min");
				
						if(!params.contains("max")) emsg("A uniform prior must contain a 'max' definition.");
						double max = toml::find<double>(params,"max");
				
						addparam(name,0.5+(max-min),min,max);
					}
					else emsg("The prior type '"+type+"' is not recognised.");
				}
			}
		}
		else{ emsg("The input file must contain parameter values through 'params'.");}
	}
	
	addparam("zero",tiny,tiny,tiny);
	
	if(tomldata.contains("comps")) {
		string name, dist, mean="", sd="";
		double inf;
		unsigned int distval = NO_DIST;
		
		const auto compsin = toml::find(tomldata,"comps");
		for(j = 0; j < compsin.size(); j++){
			const auto comps = toml::find(compsin,j);
			if(!comps.contains("name")) emsg("Compartment must contain a 'name' definition.");

			name = toml::find<std::string>(comps,"name");
			if(!comps.contains("inf")) emsg("Compartment must contain an 'inf' definition.");
			inf = toml::find<double>(comps,"inf");

			if(comps.contains("dist")){ 
				dist = toml::find<std::string>(comps, "dist");
				
				int fl = 0;
				
				if(dist == "infection"){
					distval = INFECTION;
					fl = 1;
				}
				
				if(dist == "exp"){
					distval = EXP_DIST;
					if(!comps.contains("mean")) emsg("Compartment distribution must contain a 'mean' definition.");
					mean = toml::find<std::string>(comps, "mean");
					fl = 1;
				}
				
				if(dist == "lognorm"){
					distval = LOGNORM_DIST;
					if(!comps.contains("mean")) emsg("Compartment distribution must contain a 'mean' definition.");
					mean = toml::find<std::string>(comps, "mean");
					
					if(!comps.contains("sd")) emsg("Compartment distribution must contain an 'sd' definition.");
					sd = toml::find<std::string>(comps, "sd");
					fl = 1;
				}
				
				if(dist == "gamma"){
					distval = GAMMA_DIST;
					if(!comps.contains("mean")) emsg("Compartment distribution must contain a 'mean' definition.");
					mean = toml::find<std::string>(comps, "mean");
					
					if(!comps.contains("sd")) emsg("Compartment distribution must contain an 'sd' definition.");
					sd = toml::find<std::string>(comps, "sd");
					fl = 1;
				}
				
				if(fl == 0) emsg("Distribution '"+dist+"' not recognised.");
			}

			addcomp(name,inf,distval,mean,sd); 
		}
	}
	else{ emsg("The input file must contain compartment definition through 'comps'");}
	
	if(tomldata.contains("trans")){
		const auto transin = toml::find(tomldata,"trans");
		for(j = 0; j < transin.size(); j++){
			const auto trans = toml::find(transin,j);
			
			if(!trans.contains("from")) emsg("Transition must contain a 'from' definition.");
			const auto from = toml::find<std::string>(trans, "from");
			
			if(!trans.contains("to")) emsg("Transition must contain a 'to' definition.");
			const auto to = toml::find<std::string>(trans, "to");
			
			if(trans.contains("prob")) {
				const auto prob = toml::find<std::string>(trans, "prob");
				addtrans(from,to,prob);  
			}
			else{
				addtrans(from,to,"");  
			}
		}
	}
	else{ emsg("The input file must contain transition definitions through 'trans'.");}
	
	if(tomldata.contains("betaspline")) {
		const auto bespin = toml::find(tomldata,"betaspline");
		for(j = 0; j < bespin.size(); j++){
			const auto besp = toml::find(bespin,j);
			
			if(!besp.contains("param")) emsg("Beta spline definition must contain a 'param' definition.");
			const auto name = toml::find<std::string>(besp,"param");
			
			if(!besp.contains("time")) emsg("Beta spline definition must contain a 'time' definition.");
			const auto tim = toml::find<double>(besp,"time");
			
			if(j == 0 && tim != 0) emsg("The first beta spline point must be at t=0.");
			if(j == bespin.size()-1 && tim != data.period) emsg("The last beta spline point must be at t=period.");
			if(tim < 0 || tim > data.period) emsg("The beta spline points must be within the time period.");
			
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
			const auto tim = toml::find<double>(besp,"time");
			
			if(j == 0 && tim != 0) emsg("The first phi spline point must be at t=0.");
			if(j == bespin.size()-1 && tim != data.period) emsg("The last phi spline point must be at t=period.");
			if(tim < 0 || tim > data.period) emsg("The phi spline points must be within the time period.");
			
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
		
	for(c = 0; c < comp.size(); c++){                             // Check distributions are set
		if(comp[c].trans.size() > 0 && comp[c].type == NO_DIST){
			emsg("A distribution must be put on '"+comp[c].name+"' because transitions leave it.");
		}
	}

	for(p = 0; p < param.size(); p++){
		if(param[p].used == 0){
			if(data.mode == MODE_SIM) emsg("Parameter '"+param[p].name+"' not in the model.");
		}
	}

	ntr = 0; nac = 0;
	parami.resize(param.size()); paramp.resize(param.size());

	paramval.resize(param.size()); for(p = 0; p < param.size(); p++) paramval[p] = param[p].valinit;
 
	beta.resize(data.nsettime);	phi.resize(data.nsettime);
	
	setup(paramval);
	
	if(core == 0){
		cout << endl;                                               // Outputs a summary of the model
		if(data.mode == MODE_SIM){
			cout << "Parameters:" << endl;
			for(p = 0; p < param.size()-1; p++){
				cout << param[p].name << " = " << param[p].valinit << endl;
			}
		}
		else{
			cout << "Priors:" << endl;
			for(p = 0; p < param.size()-1; p++){
				cout << param[p].name << " = ";
				if(param[p].min ==  param[p].max) cout << param[p].min << endl;
				else cout << "Uniform(" << param[p].min << " - " << param[p].max << ")" << endl;
			}
		}
		cout << endl;
			
		cout << "Compartments:" << endl; 
		for(c = 0; c < comp.size(); c++){
			cout << comp[c].name << "  Infectivity: " << comp[c].infectivity << "  "; 
			
			switch(comp[c].type){
				case INFECTION: 
					cout << " Infection" << endl;
					break;
					
				case EXP_DIST:
					cout << " Exponential  mean=" << param[comp[c].param1].name << endl;
					break;
				case GAMMA_DIST:
					cout << " Gamma mean=" << param[comp[c].param1].name 
							 << " sd=" << param[comp[c].param2].name  << endl; 
					break;
				case LOGNORM_DIST:
					cout << " Lognormal mean=" << param[comp[c].param1].name 
							 << " sd=" << param[comp[c].param2].name  << endl; 
					break;
				default:
					cout << endl;
					break;
			}
			
		}
		cout << endl;
		
		cout << "Transitions:" << endl; 
		for(t = 0; t < trans.size(); t++){
			cout << comp[trans[t].from].name << " → " << comp[trans[t].to].name;
			if(trans[t].probparam.size() > 0){
				cout << "  with probability ";
				for(j = 0; j < trans[t].probparam.size(); j++){
					if(j > 0) cout << ", ";
					cout << param[trans[t].probparam[j]].name;
				}
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
	vector <unsigned int> map;
	vector <unsigned int> nDQadd;               // Stores the mixing matrix between areas and ages 
	vector< vector <unsigned int> > DQtoadd;
	vector <vector< vector <double> > > DQvaladd;
	DQINFO dq;
	
	for(c = 0; c < comp.size(); c++) addtrans(comp[c].name,comp[c].name,"");  
	
	for(q = 0; q < data.Q.size(); q++){
		for(c = 0; c < comp.size(); c++) if(data.Q[q].comp == comp[c].name) break;
		if(c == comp.size()) emsg( "Compartment "+data.Q[q].comp+" not recognised.");
		
		if(data.Q[q].timep < 0 || data.Q[q].timep >= ntimeperiod) emsg("The time period for Q is out of range.");
	}
	
	map.resize(data.narea); for(c = 0; c < data.narea; c++) map[c] = UNSET;
	
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
		//paramval[th] = param[th].sim;
		paramval[th] = param[th].min + ran()*(param[th].max - param[th].min);
	}
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
	if(c ==  comp.size()) emsg("Cannot find compartment '"+name+"'");
	return comp[c].infectivity;
}

/// Adds a compartment to the model
void MODEL::addcomp(string name, double infectivity, unsigned int type, string param1, string param2)
{
	COMP co;
	co.name = name;
	co.infectivity = infectivity;
	co.type = type;
	
	if(param1 != "") co.param1 = findparam(param1);	
	if(param2 != "") co.param2 = findparam(param2);
	
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
void MODEL::addparam(string name, double val, double min, double max)
{
	PARAM par;
	par.name = name; par.valinit = val; par.sim = val; par.min = min; par.max = max; par.ntr = 0; par.nac = 0; par.jump = val/10;
	par.used = 0;

	param.push_back(par);
}

/// Adds a transition to the model
void MODEL::addtrans(string from, string to, string prpar)
{
	unsigned int c, cmax, j;
	
	TRANS tr;
	
	c = 0; cmax = comp.size(); while(c < cmax && from != comp[c].name) c++;
	if(c == cmax) emsg("Cannot find compartment");
	tr.from = c;	
	
	if(from != to) comp[c].trans.push_back(trans.size());
	else comp[c].transtimep = trans.size();
		
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
		for(c = 0; c <  comp.size(); c++){
			cout << "Comp " << comp[c].name << "  "; cout << comp[c].transtimep << ",  ";
			for(k = 0; k <  comp[c].trans.size(); k++) cout << comp[c].trans[k] << ", ";
			cout << endl;
		}
	
		for(tra = 0; tra < trans.size(); tra++){
			cout << tra << " " << comp[trans[tra].from].name << "->" << comp[trans[tra].to].name << " trans" << endl;
			for(k = 0; k < trans[tra].DQ.size(); k++) cout <<  trans[tra].DQ[k] << ", "; cout << "DQ" << endl;
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
	double mean, sd, z, dt;
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
		
		switch(comp[c].type){
		case EXP_DIST:
			t += -log(ran())*paramval[comp[c].param1];
			break;
		
		case GAMMA_DIST:
			mean = paramval[comp[c].param1]; sd = paramval[comp[c].param2];
			t += gammasamp(mean*mean/(sd*sd),mean/(sd*sd));
			break;
			
		case LOGNORM_DIST:
			mean = paramval[comp[c].param1]; sd = paramval[comp[c].param2];
			dt = exp(normal(mean,sd));
			t += dt;
			break;
			
		default: emsg("MODEL: EC2b"); break;
		}

		while(timep < ntimeperiod && timeperiod[timep].tend < t){    // Adds in changes in time period
			ev.trans = comp[c].transtimep; ev.t = timeperiod[timep].tend;
			evlist.push_back(ev);
			timep++;
			ev.timep = timep; 
		}
	
		if(kmax == 1) tra = comp[c].trans[0];
		else{
			z = ran(); k = 0; while(k < kmax && z > comp[c].probsum[a][k]) k++;
			if(k == kmax) emsg("Model: EC3");
			tra = comp[c].trans[k];
		}
		
		ev.trans = tra; ev.t = t;
		evlist.push_back(ev);

		c = trans[tra].to; 
	}while(1 == 1);
}

/// This does an equivelent MBP for the compartmental model
void MODEL::mbpmodel(vector <FEV> &evlisti, vector <FEV> &evlistp)
{
	unsigned int c, k, kmax, kk, tra, e, emax, i, p, p2, timep, a;
	double meani, sdi, meanp, sdp, z, t, dt;
	FEV ev;
	vector <double> probsum;

	evlistp.clear();
	
	ev = evlisti[0];
	tra = ev.trans; t = ev.t; i = ev.ind; timep = ev.timep;
	evlistp.push_back(ev);
	
	a = data.democatpos[data.ind[i].dp][0];
		
	c = trans[tra].to;
	
	emax = evlisti.size();
	for(e = 1; e < emax; e++){
		dt = evlisti[e].t-evlisti[e-1].t;
		
		kmax = comp[c].trans.size();
		if(kmax == 0) break;
		
		switch(comp[c].type){
		case EXP_DIST:
			p = comp[c].param1;
			t += dt*paramp[p]/parami[p];
			break;
		
		case GAMMA_DIST:
			emsg("model: EC9");
			break;
			
		case LOGNORM_DIST:
			p = comp[c].param1; p2 = comp[c].param2;
			meani = parami[p]; sdi = parami[p2];
			meanp = paramp[p]; sdp = paramp[p2];
			
			if(meani == meanp && sdi == sdp) t += dt;
			else t += exp(meanp + (log(dt) - meani)*sdp/sdi);
			break;
			
		default: emsg("MODEL: EC2b"); break;
		}
	
		tra = evlisti[e].trans;
		
		if(kmax > 1){
			k = 0; while(k < kmax && tra != comp[c].trans[k]) k++;
			if(k == kmax) emsg("model: EC32");
			
			if(comp[c].probp[k] < comp[c].probi[k]){  // Looks at switching to another branch
				if(ran() < 1 - comp[c].probp[a][k]/comp[c].probi[a][k]){
					do{
						z = ran(); kk = 0; while(kk < kmax && z > comp[c].probsum[a][kk]) kk++;
					}while(kk == k);
					tra = comp[c].trans[kk];
				}
			}			
		}
		
		ev.trans = tra; ev.t = t;
		evlistp.push_back(ev);

		c = trans[tra].to; 
		if(evlisti[e].trans != tra) break;
	}

	if(comp[c].trans.size() != 0) simmodel(evlistp,i,c,t);
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

/// Check that the transition data is correct
void MODEL::checktransdata()
{
	unsigned int td, tra;
	string from, to;
	TRANS tr;
	
	for(td = 0; td < data.transdata.size(); td++){
		from = data.transdata[td].from; to = data.transdata[td].to; 
		for(tra = 0; tra < trans.size(); tra++){
			tr = trans[tra];
			if(comp[tr.from].name == from && comp[tr.to].name == to) break;
		}
		
		if(tra == trans.size()){
			stringstream ss; ss << "Cannot find the transition " << from << "→" << to << ".";
			emsg(ss.str());
		}
	}
}

/// Calculates R0
vector <double> MODEL::R0calc()
{
	unsigned int c, cc, a, aa, flag, k, kmax, st, timep, q, co, vi, dp;
	double t, dt, mean, sd, fac;
	vector <double> R0;
	vector <double> R0fac;

	setup(paramval);
	
	for(c = 0; c < comp.size(); c++){
		comp[c].probin.resize(data.nage);
		comp[c].infint.resize(data.nage);
		for(a = 0; a < data.nage; a++){
			if(comp[c].name == "S") comp[c].probin[a] = 1;
			else comp[c].probin[a] = -1;
			comp[c].infint[a] = -1;
		}
	}
	
	do{
		flag = 0;
		for(c = 0; c < comp.size(); c++){
			if(comp[c].probin[0] >= 0 && comp[c].infint[0] == -1){
				switch(comp	[c].type){
				case NO_DIST: case INFECTION: dt = 0; break;
				case EXP_DIST: dt = paramval[comp[c].param1]; break;
				case GAMMA_DIST: dt = paramval[comp[c].param1]; break;
				case LOGNORM_DIST:
					mean = paramval[comp[c].param1]; sd = paramval[comp[c].param2];
					dt = exp(mean + sd*sd/2);
					break;
				default: dt = 0; emsg("MODEL: EC56"); break;
				}
				
				for(a = 0; a < data.nage; a++){
					comp[c].infint[a] = comp[c].probin[a]*comp[c].infectivity*dt;
				}					
			
				kmax = comp[c].trans.size();
				if(kmax == 1){
					cc = trans[comp[c].trans[0]].to; 
					for(a = 0; a < data.nage; a++) comp[cc].probin[a] = comp[c].probin[a];
				}
				else{
					for(k = 0; k < kmax; k++){
						cc = trans[comp[c].trans[k]].to; 
						for(a = 0; a < data.nage; a++) comp[cc].probin[a] = comp[c].probin[a]*comp[c].prob[a][k];
					}
				}
				flag = 1;
			}
		}
	}while(flag == 1);
	
	R0fac.resize(ntimeperiod);
	for(timep = 0; timep < ntimeperiod; timep++) R0fac[timep] = 0;
		
	for(q = 0; q < data.Q.size(); q++){
		timep = data.Q[q].timep;
		for(co = 0; co < comp.size(); co++) if(data.Q[q].comp == comp[co].name) break;
		if(co == comp.size()) emsg( "Compartment "+data.Q[q].comp+" not recognised.");
		
		//R0fac[timep] += comp[co].infint[0];
		
		for(c = 0; c < data.narea; c++){
			for(a = 0; a < data.nage; a++){
				fac = comp[co].infint[a]*double(data.area[c].agepop[a])/data.popsize; 
				if(fac != 0){
					vi = c*data.nage + a;
					kmax = data.Q[q].to[vi].size();
					for(k = 0; k < kmax; k++){
						cc = data.Q[q].to[vi][k];
						for(dp = 0; dp < data.ndemocatpos; dp++){
							aa = data.democatpos[dp][0];
							R0fac[timep] += fac*sus[dp]*areafac[cc]*data.Q[q].val[vi][k][aa]*data.area[cc].pop[dp];
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
