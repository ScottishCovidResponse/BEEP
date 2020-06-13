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
void MODEL::definemodel(int core, double period, int popsize, int mod)
{
	double r, tinfav, facmin = 0.8, facmax = 1.2;
	int p, c, t, fi;
		
	modelsel = mod;
	
	switch(modelsel){
  case MOD_IRISH:  	// Irish model
		// R0 determines how many individuals an infected individual on average infects
		// These 8 values represent how R0 changes over time in the simulation (this captures the effect of lockdown) 
		//double R0sim[8] = {2.9,2.8,2.6,1.5,1.1,0.7,0.8,0.9};        
		//double R0sim[8] = {3.1,3.0,2.8,1.9,1.1,0.7,0.8,0.9};      
		double R0sim[8] = {3.4,3.4,3,2.3,1.1,0.7,0.8,0.9};      
		const double dpw = 7;
		//double R0sim[8] = {2.3,2.3,2.3,2.3,0.7,0.7,0.7,0.7};  	
		
		double probA=0.25, probH=0.5, probD=0.3, Ainf = 0.5, Pinf = 0.5;
		double muE = 1.63/dpw, sdE = 0.5/dpw, tA = 8/dpw;    // Estimates of transition times from literature
		double tP = 2/dpw, tI = 2/dpw, tH = 20/dpw;

		tinfav = probA*Ainf*tA + (1-probA)*(Pinf*tP+tI);          // Calculates the average integral of infectivity after infected
			
		nspline = 8;                                              // 8 spline points represent time variarion in beta 
		for(p = 0; p < nspline; p++){
			splinet.push_back(double(p*period)/(nspline-1));
			stringstream ss; ss << "beta_" << p;
			r = R0sim[p]/tinfav; addparam(ss.str(),r,0,3*r);
			param[int(param.size())-1].betachange = 1;
		}		
		
		phiparam = param.size();
		r = 4*7.0/popsize; addparam("phi",r,0,3*r); // Adds a small external force of infection

		addparam("muE",muE,facmin*muE,facmax*muE);              // Log-normal latent period
		addparam("sdE",sdE,facmin*sdE,facmax*sdE);

		r = 1.0/tA; addparam("rA",r,facmin*r,facmax*r);         // Recovery for asymptomatic
		r = 1.0/tP; addparam("rP",r,facmin*r,facmax*r); 
		r = 1.0/tI; addparam("rI",r,facmin*r,facmax*r);
		r = 1.0/tH; addparam("rH",r,facmin*r,facmax*r);
		
		addparam("Ainf",Ainf,Ainf,Ainf);  	          // The infectivity of asymptomatic
		addparam("Pinf",Pinf,Pinf,Pinf);  	          // The infectivity of presymptomatic

		addparam("probA",probA,probA,probA);  		    // The probability of being asymptomatic
		addparam("probH",probH,probH,probH);  		    // The probability of hospitalisation given I
		addparam("probD",probD,0,1);  	            	// The probability of death given H
				
		addcomp("S",0,NO_DIST,"","");                 // Different compartment in the model
		addcomp("E",0,LOGNORM_DIST,"muE","sdE");           
		addcomp("A",Ainf,EXP_DIST,"rA",""); 
		addcomp("R",0,NO_DIST,"","");        
  	addcomp("P",Pinf,EXP_DIST,"rP",""); 
		addcomp("I",1,EXP_DIST,"rI",""); 
		addcomp("H",0,EXP_DIST,"rH",""); 
		addcomp("D",0,NO_DIST,"",""); 
		
		addtrans("S","E","");                         // All the transition in the model with branching probabilities
		addtrans("E","A","probA");  
		addtrans("E","P","");		
		addtrans("A","R",""); 
		addtrans("P","I","");    		
		addtrans("I","H","probH");
		addtrans("I","R","");
		addtrans("H","D","probD");
		addtrans("H","R","");

		fix_sus_param.resize(nfix); fix_inf_param.resize(nfix);
		for(fi = 0; fi < nfix; fi++){                             // Adds fixed effects for susceptibility
			fix_sus_param[fi] = param.size(); 
			stringstream sssus; sssus << "fixsus_" << fi;
			addparam(sssus.str(),0,0,0);
			param[int(param.size())-1].suschange = 1;
			
			fix_inf_param[fi] = param.size(); 
			stringstream ssinf; ssinf << "fixinf_" << fi;
			addparam(ssinf.str(),0,0,0);
			param[int(param.size())-1].infchange = 1;
		}
		break;
	}

	ntr = 0; nac = 0;
	parami.resize(param.size()); paramp.resize(param.size());

	paramval.resize(param.size()); for(p = 0; p < param.size(); p++) paramval[p] = param[p].valinit;
	
	betaspline(period);

	if(core == 0){
		cout << endl;                                               // Outputs a summary of the model
		cout << "Parameters:" << endl;
		for(p = 0; p < param.size(); p++){
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
					cout << " Exponentially distributed with rate " << param[comp[c].param1].name << endl;
					break;
				case GAMMA_DIST:
					cout << " Gamma distributed with mean " << param[comp[c].param1].name 
							 << " and standard deviation " << param[comp[c].param2].name  << endl; 
					break;
				case LOGNORM_DIST:
					cout << " Lognormally distributed with mean " << param[comp[c].param1].name 
							 << " and standard deviation " << param[comp[c].param2].name  << endl; 
					break;
				default:
					cout << endl;
					break;
			}
			
		}
		cout << endl;
		
		cout << "Transitions:" << endl; 
		for(t = 0; t < trans.size(); t++){
			cout << "  From: " << comp[trans[t].from].name << "  To: " << comp[trans[t].to].name;
			if(trans[t].probparam != -1) cout << "  with probability " << param[trans[t].probparam].name;
			cout << endl;
		}
		cout << endl;
	}
}

/// Randomly samples the initial parameter values from the prior (which are uniform distributions
void MODEL::priorsamp()
{
	int th;
	
	for(th = 0; th < param.size(); th++){	
		paramval[th] = param[th].min + ran()*(param[th].max - param[th].min);
	}
}

/// Gets a parameter value
double MODEL::getparam(string name)
{
	int th;
	
	for(th = 0; th < param.size(); th++){ if(param[th].name == name) break;}
	if(th == param.size()) emsg("Cannot get parameter");
	return paramval[th];	
}

/// Adds a compartment to the model
void MODEL::addcomp(string name, double infectivity, int type, string param1, string param2)
{
	int p, pmax;
	COMP co;
	co.name = name;
	co.infectivity = infectivity;
	co.type = type;
	
	if(param1 != ""){
		p = 0; pmax = param.size(); while(p < pmax && param1 != param[p].name) p++;
		if(p == pmax) emsg("Cannot find parameter");	
		co.param1 = p;
	}
	
	if(param2 != ""){
		p = 0; pmax = param.size(); while(p < pmax && param2 != param[p].name) p++;
		if(p == pmax) emsg("Cannot find parameter");	
		co.param2 = p;
	}
	
	comp.push_back(co);	
}

/// Adds a parameter to the model
void MODEL::addparam(string name, double val, double min, double max)
{
	PARAM par;
	par.name = name; par.valinit = val; par.sim = val; par.min = min; par.max = max; par.ntr = 0; par.nac = 0; par.jump = val/10;
	par.betachange = 0;	par.suschange = 0; par.infchange = 0;

	param.push_back(par);
}

/// Adds a transition to the model
void MODEL::addtrans(string from, string to, string probparam)
{
	int c, cmax, p, pmax;
	TRANS tr;
	
	c = 0; cmax = comp.size(); while(c < cmax && from != comp[c].name) c++;
	if(c == cmax) emsg("Cannot find compartment");
	tr.from = c;	
	comp[c].trans.push_back(trans.size());
	
	c = 0; cmax = comp.size(); while(c < cmax && to != comp[c].name) c++;
	if(c == cmax) emsg("Cannot find compartment");	
	tr.to = c;
	
	if(probparam == "") tr.probparam = -1;
	else{
		p = 0; pmax = param.size(); while(p < pmax && probparam != param[p].name) p++;
		if(p == pmax) emsg("Cannot find parameter");	
		tr.probparam = p;
	}
	
	trans.push_back(tr);
}
	
// Converts the spline points to a finer timestep for use in simulations.
void MODEL::betaspline(double period)
{
  int p, s, n = nspline-1;
	double t, fac, dt, a[n+1], b[n], c[n+1], d[n], h[n], alpha[n], l[n+1], mu[n+1], z[n+1];
	
	settime.resize(nsettime); beta.resize(nsettime);
	
	if(1 == 0){   // This uses a cubic spline
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
		for(s = 0; s < nsettime; s++){
			settime[s] = double((s+1)*period)/nsettime;;
			
			t = double((s+0.5)*period)/nsettime;
			while(p < nspline-1 && t > splinet[p+1]) p++;
			
			dt = t-splinet[p];	
			beta[s] = exp(a[p]+ b[p]*dt + c[p]*dt*dt + d[p]*dt*dt*dt);
		}
	}
	else{  // This uses a linear spline
		p = 0;
		for(s = 0; s < nsettime; s++){
			settime[s] = double((s+1)*period)/nsettime;;
			
			t = double((s+0.5)*period)/nsettime;
			
			while(p < nspline-1 && t > splinet[p+1]) p++;
			
			fac = (t-splinet[p])/(splinet[p+1]-splinet[p]);
			beta[s] = paramval[p]*(1-fac) + paramval[p+1]*fac;
		}
	}
}

/// Sets the transition probabilies based on the parameters
int MODEL::settransprob()
{
	int c, p, k, kmax;
	double sum, sumi, sump, prob;
	
	for(c = 0; c < comp.size(); c++){
		kmax = comp[c].trans.size();
		
		if(kmax > 1){
			comp[c].prob.resize(kmax);
			comp[c].probsum.resize(kmax);
			comp[c].probi.resize(kmax);
			comp[c].probp.resize(kmax);
			
			sum = 0; sumi = 0; sump = 0; 
			for(k = 0; k < kmax-1; k++){
				p = trans[comp[c].trans[k]].probparam; if(p == -1) emsg("model: EC1a");
				
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
void MODEL::simmodel(vector <FEV> &evlist, int i, int c, double t)
{
	int k, kmax, tra;
	double mean, sd, z, prob, sum;
	FEV ev;
	vector <double> probsum;
	
	if(c == 0){
		evlist.clear();	
		tra = 0;
		ev.trans = tra; ev.ind = i; ev.t = t; ev.done = 1;
		evlist.push_back(ev);
	
		c = trans[tra].to;
	}
	else{
		ev.ind = i;
	}
	
	do{
		kmax = comp[c].trans.size();
		if(kmax == 0) break;
		
		switch(comp[c].type){
		case EXP_DIST:
			t += -log(ran())/paramval[comp[c].param1];
			break;
		
		case GAMMA_DIST:
			mean = paramval[comp[c].param1]; sd = paramval[comp[c].param2];
			t += gammasamp(mean*mean/(sd*sd),mean/(sd*sd));
			break;
			
		case LOGNORM_DIST:
			mean = paramval[comp[c].param1]; sd = paramval[comp[c].param2];
			t += exp(normal(mean,sd));
			break;
			
		default: emsg("MODEL: EC2b"); break;
		}

		if(kmax == 1) tra = comp[c].trans[0];
		else{
			z = ran(); k = 0; while(k < kmax && z > comp[c].probsum[k]) k++;
			if(k == kmax) emsg("Model: EC3");
			tra = comp[c].trans[k];
		}
		
		ev.trans = tra; ev.t = t; ev.done = 0;
		evlist.push_back(ev);

		c = trans[tra].to; 
	}while(1 == 1);
}

/// This does an equivelent MBP for the compartmental model
void MODEL::mbpmodel(vector <FEV> &evlisti, vector <FEV> &evlistp)
{
	int c, k, kmax, kk, tra, e, emax, i, p, p2;
	double meani, sdi, meanp, sdp, z, prob, sum, t, dt;
	FEV ev;
	vector <double> probsum;

	evlistp.clear();
	
	ev = evlisti[0];
	tra = ev.trans; t = ev.t; i = ev.ind;
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
			t += dt*parami[p]/paramp[p];
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
				if(ran() < 1 - comp[c].probp[k]/comp[c].probi[k]){
					do{
						z = ran(); kk = 0; while(kk < kmax && z > comp[c].probsum[kk]) kk++;
					}while(kk == k);
					tra = comp[c].trans[kk];
				}
			}			
		}
		
		ev.trans = tra; ev.t = t; ev.done = 0;
		evlistp.push_back(ev);

		c = trans[tra].to; 
		if(evlisti[e].trans != tra) break;
	}

	if(comp[c].trans.size() != 0) simmodel(evlistp,i,c,t);
}
