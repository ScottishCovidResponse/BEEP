// This file defines the compartmental model and is used to simulate compartmental dynamics after an individual becomes infected

#include <string>
#include <sstream>
#include <iostream>

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
	if((details.mode == inf || details.mode == abcsmc) && infmax == large){
		emsgroot("Input file must contain a limit on the maximum number of individuals through 'infmax'.");
	}	
	
	timeperiod = data.timeperiod; ntimeperiod = timeperiod.size();

	if(details.mode == sim || details.mode == multisim){
		vector <string> name;
		vector <double> val;
		inputs.find_param(name,val);
		for(unsigned int th = 0; th < name.size(); th++) addparam(name[th],val[th],val[th]);
	}

	if(details.mode == inf || details.mode == abcsmc){
		vector <string> name;
		vector <double> min,max;
		inputs.find_prior(name,min,max);
		for(unsigned int th = 0; th < name.size(); th++) addparam(name[th],min[th],max[th]);
	}
	
	regioneffect = 0;
	sigma_param = UNSET;
	for(unsigned int th = 0; th < param.size(); th++){
		if(param[th].name == "regeff_sigma"){
			param[th].used = 1;
			
			regioneffect = 1;
			sigma_param = th;	
			
			regioneff_param.resize(data.nregion);
			for(unsigned int r = 0; r < data.nregion; r++){
				regioneff_param[r] = param.size();
				stringstream ss; ss << "reff_" << data.region[r].code;
				addparam(ss.str(),-large,large);
				param[param.size()-1].used = 1;
			}
		}
	}
		
	unsigned int th;
	for(th = 0; th < param.size(); th++){
		if(param[th].name == "logbeta"){
			param[th].used = 1;
			logbeta_param = th;
			break;
		}
	}
	if(th == param.size()) emsg("In the input TOML file a 'logbeta' parameter must be specified");
	
	addparam("zero",tiny,tiny);

	vector <string> name;
	vector <double> infectivity;
	inputs.find_comps(name,infectivity);
	for(unsigned int c = 0; c < name.size(); c++) addcomp(name[c],infectivity[c]);
	
	vector <string> from, to, prpar, mean, cv;
	vector <int> type;
	inputs.find_trans(from,to,prpar,type,mean,cv);
	for(unsigned int tr = 0; tr < from.size(); tr++) addtrans(from[tr],to[tr],prpar[tr],type[tr],mean[tr],cv[tr]);
	
	if(details.mode == inf || details.mode == abcsmc){
		priorcomps = inputs.find_priorcomps(comp);
	}
			
	vector <int> time;
	vector <string> pname; 
	string splinetype;

	splinetype = "betaspline";
	inputs.find_spline(details,splinetype,time,pname);
	for(unsigned int i = 0; i < time.size(); i++){
		SPLINEP spl;
		spl.t = time[i];
		spl.param = findparam(pname[i]);
		betaspline.push_back(spl);
	}	
	
	splinetype = "phispline";
	inputs.find_spline(details,splinetype,time,pname);
	for(unsigned int i = 0; i < time.size(); i++){
		SPLINEP spl;
		spl.t = time[i];
		spl.param = findparam(pname[i]);
		phispline.push_back(spl);
	}	
		
	sus_param.resize(data.ndemocat);
	for(unsigned int c = 0; c < data.ndemocat; c++){
		for(unsigned int fi = 0; fi < data.democat[c].value.size(); fi++){
			sus_param[c].push_back(findparam(data.democat[c].param[fi]));
		}
	}

	for(unsigned int c = 0; c < data.ncovar; c++){
		covar_param.push_back(findparam(data.covar[c].param));
	}

	for(unsigned int p = 0; p < param.size(); p++){
		if(param[p].used == 0) emsg("The [arameter '"+param[p].name+"' is not used in the model.");
	}

	ntr = 0; nac = 0;
	parami.resize(param.size()); paramp.resize(param.size());

	paramval.resize(param.size());
	priorsamp();
 
	beta.resize(details.nsettime); phi.resize(details.nsettime);
	
	for(unsigned int tra = 0; tra < trans.size(); tra++){
		switch(trans[tra].type){
		case exp_dist: param[trans[tra].param_mean].type = 1; break;
		case gamma_dist: case lognorm_dist: param[trans[tra].param_mean].type = 1; param[trans[tra].param_cv].type = 1; break;
		}
		
		if(trans[tra].istimep == 0){
			if(trans[tra].probparam.size() > 0){
				for(unsigned int a = 0; a < data.nage; a++) param[trans[tra].probparam[a]].type = 2; 
			}
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
		for(unsigned int p = 0; p < param.size()-1; p++){
			cout << "  " << param[p].name << " = " << param[p].min << endl;
		}
	}
	else{
		cout << "Priors:" << endl;
		for(unsigned int p = 0; p < param.size()-1; p++){
			cout << "  " << param[p].name << " = ";
			if(param[p].min ==  param[p].max) cout << param[p].min << endl;
			else cout << "Uniform(" << param[p].min << " - " << param[p].max << ")" << endl;
		}
	}
	cout << endl;
		
	cout << "Compartments:" << endl; 
	for(unsigned int c = 0; c < comp.size(); c++){
		cout << "  " << comp[c].name << "  Infectivity: " << comp[c].infectivity << endl; 			
	}
	cout << endl;
	
	cout << "Transitions:" << endl; 
	for(unsigned int tra = 0; tra < trans.size(); tra++){
		if(trans[tra].from != trans[tra].to){
			cout << "  " << comp[trans[tra].from].name << " → " << comp[trans[tra].to].name;
			if(trans[tra].probparam.size() > 0){
				cout << " with probability ";
				for(unsigned int j = 0; j < trans[tra].probparam.size(); j++){
					if(j > 0) cout << ", ";
					cout << param[trans[tra].probparam[j]].name;
				}
				cout << "  ";
			}
			
			switch(trans[tra].type){
				case infection_dist: 
					cout << " Infection";
					break;
				case exp_dist:
					cout << " Exponential  mean=" << param[trans[tra].param_mean].name;
					break;
				case gamma_dist:
					cout << " Gamma mean=" << param[trans[tra].param_mean].name 
							 << " cv=" << param[trans[tra].param_cv].name; 
					break;
				case lognorm_dist:
					cout << " Lognormal mean=" << param[trans[tra].param_mean].name 
							 << " cv=" << param[trans[tra].param_cv].name; 
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
	paramval = paramv;
	if(settransprob() == 1) return 1;
	
	timevariation();
	setsus();
	setarea();
	return 0;
}

/// Copies values used for the initial state (used in MBPs)
void MODEL::copyi()
{
	unsigned int c;
	
	parami = paramval; betai = beta; phii = phi; susi = sus; areafaci = areafac;
	for(c = 0; c < comp.size(); c++) comp[c].probi = comp[c].prob;
}
	
/// Copies values used for the proposed state (used in MBPs)
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
	
	for(c = 0; c < comp.size(); c++) addtrans(comp[c].name,comp[c].name,"",timep_dist,"","");  
	
	for(q = 0; q < data.Q.size(); q++){
		for(c = 0; c < comp.size(); c++) if(data.Q[q].comp == comp[c].name) break;
		if(c == comp.size()) emsg("The compartment '"+data.Q[q].comp+"' in '"+data.Q[q].name+"' is not recognised.");
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
	unsigned int th, r;

	for(th = 0; th < param.size(); th++){	
		paramval[th] = param[th].min + ran()*(param[th].max - param[th].min);
	}
	
	if(regioneffect == 1){
		for(r = 0; r < data.nregion; r++) paramval[regioneff_param[r]] = normal(0,paramval[sigma_param]);
	}
	
	//for(th = 0; th < param.size(); th++) cout << "paramval[" << th <<"] = " << paramval[th] << ";" << endl;
}

/// Gets a parameter value
double MODEL::getparam(const string& name)
{	
	return paramval[findparam(name)];
}

/// Gets the infectivity of a compartment
double MODEL::getinfectivity(const string& name)
{
	unsigned int c;

	c = 0; while(c < comp.size() && comp[c].name != name) c++;
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
	unsigned int p, pmax;
	
	p = 0; pmax = param.size(); while(p < pmax && name != param[p].name) p++;
	if(p == pmax) emsg("Cannot find the parameter '"+name+"'");	
	param[p].used = 1;
	
	return p;
}

/// Adds a parameter to the model
void MODEL::addparam(const string& name, double min, double max)
{
	PARAM par;
	par.name = name; par.min = min; par.max = max; par.ntr = 0; par.nac = 0; par.jump = 0.5*(min+max)/10; if(par.jump == 0) par.jump = 0.1;
	par.used = 0; par.type = 0;

	param.push_back(par);
}

/// Adds a transition to the model
void MODEL::addtrans(const string& from, const string& to, const string& prpar, unsigned int type, const string& mean, const string& cv)
{
	unsigned int c, cmax, j;
	
	TRANS tr;
	
	c = 0; cmax = comp.size(); while(c < cmax && from != comp[c].name) c++;
	if(c == cmax) emsg("Cannot find the 'from' compartment '"+from+"' for the transition");
	tr.from = c;	
	
	if(from != to){ tr.istimep = 0; comp[c].trans.push_back(trans.size());}
	else{ tr.istimep = 1; comp[c].transtimep = trans.size();}
		
	c = 0; cmax = comp.size(); while(c < cmax && to != comp[c].name) c++;
	if(c == cmax) emsg("Cannot find the 'to' compartment '"+to+"' for the transition");	
	tr.to = c;
	
	if(prpar != ""){
		vector <string> probparam = split(prpar,',');
	
		if(probparam.size() != data.nage) emsg("For the transition '"+from+"→"+to+"' the number of parameters in expression '"+prpar+"' should equal the number of age groups.");
		for(j = 0; j < probparam.size(); j++){
			tr.probparam.push_back(findparam(probparam[j]));
		}
	}
	
	tr.type = type;
	if(mean != "") tr.param_mean = findparam(mean);	else tr.param_mean = UNSET;
	if(cv != "") tr.param_cv = findparam(cv); else tr.param_cv = UNSET;
	
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
		for(s = 0; s < details.nsettime; s++){		
			t = double((s+0.5)*details.period)/details.nsettime;
			while(p < int(nspline)-1 && t > splinet[p+1]) p++;
			
			dt = t-splinet[p];	
			beta[s] = exp(a[p]+ b[p]*dt + c[p]*dt*dt + d[p]*dt*dt*dt);
		}
	*/
	
  // Uses a linear spline for beta
	p = 0;
	for(s = 0; s < details.nsettime; s++){	
		t = double((s+0.5)*details.period)/details.nsettime;
		
		while(p < int(betaspline.size())-1 && t > betaspline[p+1].t) p++;
		
		fac = (t-betaspline[p].t)/(betaspline[p+1].t-betaspline[p].t);
		beta[s] = (paramval[betaspline[p].param]*(1-fac) + paramval[betaspline[p+1].param]*fac);
	}
	
	// Uses a linear spline for phi
	p = 0;
	for(s = 0; s < details.nsettime; s++){	
		t = double((s+0.5)*details.period)/details.nsettime;
		
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
					p = trans[comp[c].trans[k]].probparam[a]; if(p == UNSET) emsgEC("model",1);		
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
			if(k == kmax) emsgEC("Model",2);
			tra = comp[c].trans[k];
		}
		
		switch(trans[tra].type){
		case exp_dist:
			dt = -log(ran())*paramval[trans[tra].param_mean];
			break;
		
		case gamma_dist:
			mean = paramval[trans[tra].param_mean]; sd = paramval[trans[tra].param_cv]*mean;
			dt = gammasamp(mean*mean/(sd*sd),mean/(sd*sd));
			break;
			
		case lognorm_dist:
			mean_ns = paramval[trans[tra].param_mean]; cv_ns = paramval[trans[tra].param_cv];
			sd = sqrt(log((1+cv_ns*cv_ns))); mean = log(mean_ns) - sd*sd/2;
			dt = exp(normal(mean,sd));
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
				if(k == kmax) emsgEC("model",4);
				
				if(comp[c].probp[a][k] < comp[c].probi[a][k]){  // Looks at switching to another branch
					if(ran() < 1 - comp[c].probp[a][k]/comp[c].probi[a][k]){
						sum = 0; sumst.clear();
						for(k = 0; k < kmax; k++){
							dif = comp[c].probp[a][k] - comp[c].probi[a][k];
							if(dif > 0) sum += dif;
							sumst.push_back(sum);
						}
						
						z = ran()*sum; k = 0; while(k < kmax && z > sumst[k]) k++;
						if(k == kmax) emsgEC("Model",5);
						tra = comp[c].trans[k];
					}
				}	
			}
			
			switch(trans[tra].type){
			case exp_dist:
				p = trans[tra].param_mean;
				dtnew = dt*paramp[p]/parami[p];
				break;
			
			case gamma_dist:
				emsgEC("model",6);
				break;
				
			case lognorm_dist:
				p = trans[tra].param_mean; p2 = trans[tra].param_cv;
				
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
		sum = paramval[logbeta_param];
		for(j = 0; j < data.ncovar; j++) sum += paramval[covar_param[j]]*data.area[c].covar[j];
		
		if(regioneffect == 1) sum += paramval[regioneff_param[data.area[c].region]];
		
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
	unsigned int r;
	double Pr, sd;
	
	Pr = 0;
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
		sd = paramval[sigma_param];
		for(r = 0; r < data.nregion; r++){
			Pr += normalprob(paramval[regioneff_param[r]],0,sd*sd);
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
		
		if(cst.size() > 0 || kst.size() > 0 || probst.size() > 0) emsgEC("Model",8);
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
				case infection_dist: dt = 0; break;
				case exp_dist: dt = paramval[trans[tra].param_mean]; break;
				case gamma_dist: dt = paramval[trans[tra].param_mean]; break;
				case lognorm_dist: dt = paramval[trans[tra].param_mean]; break;
				default: emsgEC("MODEL",9); break;
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
		if(co == comp.size()) emsg("Compartment '"+data.Q[q].comp+"' in '"+data.Q[q].name+"' not recognised.");
		
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
	
	R0.resize(details.nsettime);
	timep = 0; 
	for(st = 0; st < details.nsettime; st++){
		t = details.settime[st];
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
				if(dt == 0) emsgEC("Model",10);
				t += dt;
				
				trans[tra].numvisittot++;
				switch(trans[tra].type){
				case exp_dist: trans[tra].dtsum += dt; break;
				case gamma_dist: case lognorm_dist: trans[tra].dtlist.push_back(dt); break;
				}
			}
		}
	}
	
	Li_dt = likelihood_dt(paramv);
	Li_prob = likelihood_prob();
	paramval = paramv; if(settransprob() == 1) emsgEC("Model",11);
	
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
	
	paramval = paramv; if(settransprob() == 1) emsgEC("Model",12);
		
	if(checkon == 1){
		dd = likelihood_dt(paramv)-Li_dt; if(dd*dd > tiny) emsgEC("Model",13);
		dd = likelihood_prob()-Li_prob; if(dd*dd > tiny) emsgEC("Model",14);
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
		case exp_dist:
			r = 1.0/paramv[trans[tra].param_mean];
			L += trans[tra].numvisittot*log(r) - r*trans[tra].dtsum;
			break;
		
		case gamma_dist:
			mean = paramv[trans[tra].param_mean]; sd = paramv[trans[tra].param_cv]*mean;
			jmax = trans[tra].dtlist.size();
			for(j = 0; j < jmax; j++) L += gammaprob(trans[tra].dtlist[j],mean*mean/(sd*sd),mean/(sd*sd));
			break;
			
		case lognorm_dist:
			mean_ns = paramval[trans[tra].param_mean]; cv_ns = paramval[trans[tra].param_cv];
			sd = sqrt(log(1+cv_ns*cv_ns)); mean = log(mean_ns) - sd*sd/2;
			jmax = trans[tra].dtlist.size();
			
			double ss = 0; for(j = 0; j < jmax; j++) ss += trans[tra].dtlist[j];
			for(j = 0; j < jmax; j++){
				if(trans[tra].dtlist[j] == 0) emsgEC("Model",15);
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
		case exp_dist:
			ri = 1.0/paramvi[trans[tra].param_mean]; rf = 1.0/paramvf[trans[tra].param_mean];
			if(ri != rf){
				L += trans[tra].numvisittot*log(rf) - rf*trans[tra].dtsum;
				L -= trans[tra].numvisittot*log(ri) - ri*trans[tra].dtsum;
			}
			break;
		
		case gamma_dist:
			meani = paramvi[trans[tra].param_mean]; sdi = paramvi[trans[tra].param_cv]*meani;
			meanf = paramvf[trans[tra].param_mean]; sdf = paramvf[trans[tra].param_cv]*meanf;
	
			if(meani != meanf || sdi != sdf){
				jmax = trans[tra].dtlist.size();
				for(j = 0; j < jmax; j++){
					L += gammaprob(trans[tra].dtlist[j],meanf*meanf/(sdf*sdf),meanf/(sdf*sdf));
					L -= gammaprob(trans[tra].dtlist[j],meani*meani/(sdi*sdi),meani/(sdi*sdi));
				}
			}				
			break;
			
		case lognorm_dist:
			mean_nsi = paramvi[trans[tra].param_mean]; cv_nsi = paramvi[trans[tra].param_cv];
			mean_nsf = paramvf[trans[tra].param_mean]; cv_nsf = paramvf[trans[tra].param_cv];
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
void MODEL::oe(const string& name, vector <FEV> &ev)
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
