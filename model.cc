// This file defines the compartmental model and is used to simulate compartmental dynamics after an individual becomes infected

#include <string>
#include <sstream>
#include <iostream>

#include "math.h"
#include "consts.hh"
#include "utils.hh"
#include "model.hh"


using namespace std;

/// Defines the compartmental model
void MODEL::definemodel(DATA &data, unsigned int core, double period, unsigned int popsize, const toml::basic_value<::toml::discard_comments, std::unordered_map, std::vector> &tomldata)
{
	unsigned int p, c, t, j;
	int fi;
	SPLINEP spl;
	
	ntimeperiod = data.ntimeperiod; // Copies from data
	timeperiod = data.timeperiod;

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
	
	setup(data,paramval);
	
	if(core == 0){
		cout << endl;                                               // Outputs a summary of the model
		cout << "Parameters:" << endl;
		for(p = 0; p < param.size()-1; p++){
			cout << param[p].name << " " << param[p].valinit << " (" << param[p].min << " - " << param[p].max << ")" << endl;
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
			if(trans[t].probparam != UNSET) cout << "  with probability " << param[trans[t].probparam].name;
			cout << endl;
		}
		cout << endl;
	}
}

/// Sets up the model with a set of parameters
void MODEL::setup(DATA &data, vector <double> &paramv)
{
	paramval = paramv;
	timevariation(data);
	setsus(data);
}

/// Copies values used for the initial state (MBPs)
void MODEL::copyi()
{
	parami = paramval; betai = beta; phii = phi; susi = sus;
}
	
/// Copies values used for the proposed state (MBPs)
void MODEL::copyp()
{
	paramp = paramval; betap = beta; phip = phi; susp = sus;
}

/// Adds in the tensor Q to the model
void MODEL::addQ(DATA &data)
{
	unsigned int q, qi, qf, c, cc, ci, cf, tra, timep, timepi, timepf, loop, j, jmax, a, v;
	string compi, compf;
	double fac;
	vector <unsigned int> map;
	vector <unsigned int> nDQadd;               // Stores the mixing matrix between areas and ages 
	vector< vector <unsigned int> > DQtoadd;
	vector <vector< vector <double> > > DQvaladd;
	
	for(c = 0; c < comp.size(); c++) addtrans(comp[c].name,comp[c].name,"");  
	
	for(q = 0; q < data.Qnum; q++){
		for(c = 0; c < comp.size(); c++) if(data.Qcomp[q] == comp[c].name) break;
		if(c == comp.size()){ stringstream ss; ss << "Compartment " << data.Qcomp[q] << " not recognised."; emsg(ss.str());} 
		
		if(data.Qtimeperiod[q] < 0 || data.Qtimeperiod[q] >= ntimeperiod) emsg("The time period for Q is out of range.");
	}
	
	map.resize(data.narea); for(c = 0; c < data.narea; c++) map[c] = UNSET;
	
	for(tra = 0; tra < trans.size(); tra++){
		ci = trans[tra].from;
		cf = trans[tra].to;
		compi = comp[ci].name;
		compf = comp[cf].name;
	
		for(timep = 0; timep < ntimeperiod; timep++){
			timepi = timep; timepf = timep;
			if(compi == compf) timepf++;
			
			if(timepf < ntimeperiod){
				qi = 0; while(qi < data.Qnum && !(data.Qcomp[qi] == compi && data.Qtimeperiod[qi] == timepi)) qi++;
				if(qi == data.Qnum) qi = UNSET;
				
				qf = 0; while(qf < data.Qnum && !(data.Qcomp[qf] == compf && data.Qtimeperiod[qf] == timepf)) qf++;
				if(qf == data.Qnum) qf = UNSET;
				
				if(qi == UNSET && qf == UNSET){
					trans[tra].DQ.push_back(UNSET);
				}
				else{
					DQtoadd.clear(); DQvaladd.clear();
					nDQadd.resize(data.narage); DQtoadd.resize(data.narage); DQvaladd.resize(data.narage); 
					for(v = 0; v < data.narage; v++){
						nDQadd[v] = 0;
						
						for(loop = 0; loop < 2; loop++){
							switch(loop){
							case 0: q = qi; fac = -comp[ci].infectivity; break;
							case 1: q = qf; fac = comp[cf].infectivity; break;
							}
							if(q != UNSET){
								jmax = data.nQ[q][v];
								for(j = 0; j < jmax; j++){
									cc = data.Qto[q][v][j];
									if(map[cc] == UNSET){
										map[cc] = nDQadd[v];
										DQtoadd[v].push_back(cc);
										DQvaladd[v].push_back(vector <double> ());
										DQvaladd[v][nDQadd[v]].resize(data.nage);
										for(a = 0; a < data.nage; a++) DQvaladd[v][nDQadd[v]][a] = 0;
										nDQadd[v]++;
									}
									
									for(a = 0; a < data.nage; a++) DQvaladd[v][map[cc]][a] += fac*data.Qval[q][v][j][a];
								}
							}
						}
					
					
						for(j = 0; j < nDQadd[v]; j++) map[DQtoadd[v][j]] = UNSET;
					}
					
					trans[tra].DQ.push_back(nDQ.size());
					
					nDQ.push_back(nDQadd);
					DQto.push_back(DQtoadd);
					DQval.push_back(DQvaladd);
				}
			}
		}
	}
	
	DQnum = nDQ.size();
}

/// Randomly samples the initial parameter values from the prior (which are uniform distributions
void MODEL::priorsamp()
{
	unsigned int th;
	
	for(th = 0; th < param.size(); th++){	
		paramval[th] = param[th].sim;
		//paramval[th] = param[th].min + ran()*(param[th].max - param[th].min);
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
	unsigned int p, pmax;
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
void MODEL::addtrans(string from, string to, string probparam)
{
	unsigned int c, cmax, p, pmax;
	TRANS tr;
	
	c = 0; cmax = comp.size(); while(c < cmax && from != comp[c].name) c++;
	if(c == cmax) emsg("Cannot find compartment");
	tr.from = c;	
	
	if(from != to) comp[c].trans.push_back(trans.size());
	else comp[c].transtimep = trans.size();
		
	c = 0; cmax = comp.size(); while(c < cmax && to != comp[c].name) c++;
	if(c == cmax) emsg("Cannot find compartment");	
	tr.to = c;
	
	if(probparam == "") tr.probparam = UNSET;
	else{
		tr.probparam = findparam(probparam);
	}
	
	trans.push_back(tr);
}
	
/// Generates the time variation in beta and phi from the parameters
void MODEL::timevariation(DATA &data)
{
  unsigned int s, j;
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
		
		while(p < betaspline.size()-1 && t > betaspline[p+1].t) p++;
		
		fac = (t-betaspline[p].t)/(betaspline[p+1].t-betaspline[p].t);
		beta[s] = (paramval[betaspline[p].param]*(1-fac) + paramval[betaspline[p+1].param]*fac)/units;
	}
	
	// This uses a linear spline for phi
	p = 0;
	for(s = 0; s < data.nsettime; s++){	
		t = double((s+0.5)*data.period)/data.nsettime;
		
		while(p < phispline.size()-1 && t > phispline[p+1].t) p++;
		
		fac = (t-phispline[p].t)/(phispline[p+1].t-phispline[p].t);
		phi[s] = (paramval[phispline[p].param]*(1-fac) + paramval[phispline[p+1].param]*fac)/(units*data.popsize);
	}
}

/// Sets the transition probabilies based on the parameters
unsigned int MODEL::settransprob()
{
	unsigned int c, k, kmax, tra;
	int p;
	double sum, sumi, sump, prob;
	
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
			comp[c].prob.resize(kmax);
			comp[c].probsum.resize(kmax);
			comp[c].probi.resize(kmax);
			comp[c].probp.resize(kmax);
			
			sum = 0; sumi = 0; sump = 0; 
			for(k = 0; k < kmax-1; k++){
				p = trans[comp[c].trans[k]].probparam; if(p == UNSET) emsg("model: EC1a");
				
				prob = paramval[p]; comp[c].prob[k] = prob; sum += prob; if(prob < 0) return 0;
				prob = parami[p]; comp[c].probi[k] = prob; sumi += prob; if(prob < 0) return 0;
				prob = paramp[p]; comp[c].probp[k] = prob; sump += prob; if(prob < 0) return 0;
			} 
			prob = 1-sum; comp[c].prob[kmax-1] = prob; if(prob < 0) return 0;
			prob = 1-sumi; comp[c].probi[kmax-1] = prob; if(prob < 0) return 0;
			prob = 1-sump; comp[c].probp[kmax-1] = prob; if(prob < 0) return 0;
			
			sum = 0;
			for(k = 0; k < kmax; k++){
				sum += comp[c].prob[k];
				comp[c].probsum[k] = sum;
			}
		}
	}
	
	return 1;
}

/// This simulates from the model and generates an event list
void MODEL::simmodel(vector <FEV> &evlist, unsigned int i, unsigned int c, double t)
{
	unsigned int k, kmax, tra, timep;
	double mean, sd, z;

	FEV ev;
	vector <double> probsum;
	
	timep = 0; while(timep < ntimeperiod && t > timeperiod[timep]) timep++;

	ev.ind = i; ev.timep = timep; 
		
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
			t += -units*log(ran())/paramval[comp[c].param1];
			break;
		
		case GAMMA_DIST:
			mean = paramval[comp[c].param1]; sd = paramval[comp[c].param2];
			t += units*gammasamp(mean*mean/(sd*sd),mean/(sd*sd));
			break;
			
		case LOGNORM_DIST:
			mean = paramval[comp[c].param1]; sd = paramval[comp[c].param2];
			t += units*exp(normal(mean,sd));
			break;
			
		default: emsg("MODEL: EC2b"); break;
		}

		while(timeperiod[timep] < t){    // Adds in changes in time period
			ev.trans = comp[c].transtimep; ev.t = timeperiod[timep];
			evlist.push_back(ev);
			timep++;
			ev.timep = timep; 
		}
	
		if(kmax == 1) tra = comp[c].trans[0];
		else{
			z = ran(); k = 0; while(k < kmax && z > comp[c].probsum[k]) k++;
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
	unsigned int c, k, kmax, kk, tra, e, emax, i, p, p2, timep;
	double meani, sdi, meanp, sdp, z, t, dt;
	FEV ev;
	vector <double> probsum;

	evlistp.clear();
	
	ev = evlisti[0];
	tra = ev.trans; t = ev.t; i = ev.ind; timep = ev.timep;
	evlistp.push_back(ev);
	
	c = trans[tra].to;
	
	emax = evlisti.size();
	for(e = 1; e < emax; e++){
		dt = evlisti[e].t-evlisti[e-1].t;
		
		kmax = comp[c].trans.size();
		if(kmax == 0) break;
		
		switch(comp[c].type){
		case EXP_DIST:
			p = comp[c].param1;
			t += units*dt*parami[p]/paramp[p];
			break;
		
		case GAMMA_DIST:
			emsg("model: EC9");
			break;
			
		case LOGNORM_DIST:
			p = comp[c].param1; p2 = comp[c].param2;
			meani = parami[p]; sdi = parami[p2];
			meanp = paramp[p]; sdp = paramp[p2];
			
			if(meani == meanp && sdi == sdp) t += dt;
			else t += units*exp(meanp + (log(dt) - meani)*sdp/sdi);
			break;
			
		default: emsg("MODEL: EC2b"); break;
		}
	
		tra = evlisti[e].trans;
		
		if(kmax > 1){
			k = 0; while(k < kmax && tra != comp[c].trans[k]) k++;
			if(k == kmax) emsg("model: EC32");
			
			if(comp[c].probp[k] < comp[c].probi[k]){  // Looks at switching to another branch
				if(ran() < 1 - comp[c].probp[k]/comp[c].probi[k]){
					do{
						z = ran(); kk = 0; while(kk < kmax && z > comp[c].probsum[kk]) kk++;
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
void MODEL::setsus(DATA &data)     
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

/// Check that the transition data is correct
void MODEL::checktransdata(DATA &data)
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
