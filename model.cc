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
	
	double r, tinfav;
	int p, c, t, fi;
		
	modelsel = mod;
	
	if(modelsel == MOD_IRISH){  	// Irish model
		// R0 determines how many individuals an infected individual on average infects
		// These 8 values represent how R0 changes over time in the simulation (this captures the effect of lockdown) 
		//double R0sim[8] = {2.9,2.8,2.6,1.5,1.1,0.7,0.8,0.9};        
		double R0sim[8] = {3.1,3.0,2.8,1.9,1.1,0.7,0.8,0.9};      
		//double R0sim[8] = {2.3,2.3,2.3,2.3,0.7,0.7,0.7,0.7};  	
		
		double afrac = 0.25, muEA = 1.63/7.0, sdEA = 0.5/7.0, aI = 0.5;
		double tIaR = 8/7.0, tIH = 2/7.0, tHD = 100/7.0, tHR = 20/7.0; // Estimates of transition times from literature

		addcomp("S",0);
		addcomp("Ea",0); addcomp("E",0);                          // Different compartment in the model
		addcomp("Ia",aI); addcomp("R",0);          
  	addcomp("I",1); addcomp("H",0); addcomp("D",0); 
		
		
		tinfav = afrac*aI*tIaR + (1-afrac)*tIH;                   // Calculates the average integral of infectivity after infected
			
		nspline = 8;                                              // 8 spline points represent time variarion in beta 
		for(p = 0; p < nspline; p++){
			splinet.push_back(double(p*period)/(nspline-1));
			stringstream ss; ss << "beta_" << p;
			r = R0sim[p]/tinfav; addparam(ss.str(),r,0,3*r);
			param[int(param.size())-1].betachange = 1;
		}		
		
		phiparam = param.size();
		r = 7.0/popsize; addparam("phi",r,0,3*r);                 // Adds a small external force of infection
		
		afracparam = param.size();
		addparam("afrac",afrac,afrac,afrac);  		
		//addparam("afrac",afrac,0.2,0.3);  

		aIparam = param.size();
		addparam("aI",aI,aI,aI);  		

		fix_sus_param.resize(nfix); fix_inf_param.resize(nfix);
		for(fi = 0; fi < nfix; fi++){                             // Adds fixed effects for susceptibility
			fix_sus_param[fi] = param.size(); 
			stringstream sssus; sssus << "fixsus_" << fi;
			//addparam(sssus.str(),0.1,0.0,0.2);
			addparam(sssus.str(),0,0,0);
			param[int(param.size())-1].suschange = 1;
			
			fix_inf_param[fi] = param.size(); 
			stringstream ssinf; ssinf << "fixinf_" << fi;
			addparam(ssinf.str(),0,0,0);
			param[int(param.size())-1].infchange = 1;
		}
			
		addparam("muEA",muEA,muEA,muEA);                          // All the parameters in the model (with uniform priors)
		addparam("sdEA",sdEA,sdEA,sdEA);
		//r = 1.0/tIaR; addparam("rIaR",r,0.9*r,1.1*r);
		//r = 1.0/tIH; addparam("rIH",r,0.9*r,1.1*r);
		//r = 1.0/tHD; addparam("rHD",r,0.9*r,1.1*r);
		//r = 1.0/tHR; addparam("rHR",r,0.9*r,1.1*r);
	
		r = 1.0/tIaR; addparam("rIaR",r,r,r);
		r = 1.0/tIH; addparam("rIH",r,r,r);
		r = 1.0/tHD; addparam("rHD",r,r,r);
		r = 1.0/tHR; addparam("rHR",r,r,r);
	
		addtrans("S","Ea",INFECTION,"",""); 
		addtrans("S","E",INFECTION,"",""); 
		addtrans("Ea","Ia",LOGNORM_DIST,"muEA","sdEA");            // We define all the transition in the model
		addtrans("Ia","R",EXP_DIST,"rIaR","");
		addtrans("E","I",LOGNORM_DIST,"muEA","sdEA");    		
		addtrans("I","H",EXP_DIST,"rIH","");
		addtrans("H","R",EXP_DIST,"rHR","");
		addtrans("H","D",EXP_DIST,"rHD","");
	}
	else{                                                      	 // The previous model 
		double R0sim[5] = {2.3,2.0,1.5,0.5,0.5};       
	
		double tEA = 2.5/7.0, tAI = 3/7.0;                         // Estimates of transition times from literature
		double tIH =  3.1/7.0, tAR = 3.5/7.0, tIR = 16.7/7.0, tID = 12.9/7.0, tHD = 8/7.0, tHR = 13.1/7.0;
	
		double rAI, rAR, rIR, rIH, rID;

		addcomp("E",0); addcomp("A",0.2); addcomp("I",1);          // Different compartment in the model
		addcomp("H",0); addcomp("R",0); addcomp("D",0); 
		
		// This calculates the average integral of infectivity after an indvidual becomes infected
		rAI = 1.0/tAI; rAR = 1.0/tAR; rIH = 1.0/tIH; rIR = 1.0/tIR; rID = 1.0/tID;
		tinfav = 0.2/(rAI + rAR) + 1*(rAI/(rAI+rAR))*(1.0/(rIH + rIR + rID));
			 
		nspline = 5;                                               // 5 spline points represent time variarion in beta 
		for(p = 0; p < nspline; p++){
			splinet.push_back(double(p*period)/(nspline-1));
			stringstream ss; ss << "beta_" << p;
			//r = R0sim[p]/tinfav; addparam(ss.str(),r,0,3*r);
			r = R0sim[p]/tinfav; addparam(ss.str(),r,0,1.4*r);
			param[int(param.size())-1].betachange = 1;
		}		
	
		phiparam = param.size();
		r = 7.0/popsize; addparam("phi",r,r,r);                      // Adds a small external force of infection
		
		fix_sus_param.resize(nfix); fix_inf_param.resize(nfix);
		for(fi = 0; fi < nfix; fi++){                                // Adds fixed effects for susceptibility
			fix_sus_param[fi] = param.size(); 
			stringstream sssus; sssus << "fixsus_" << fi;
			addparam(sssus.str(),0.1,0.1,0.1);
			param[int(param.size())-1].suschange = 1;
			
			fix_inf_param[fi] = param.size(); 
			stringstream ssinf; ssinf << "fixinf_" << fi;
			addparam(ssinf.str(),0.1,0.1,0.1);
			param[int(param.size())-1].infchange = 1;
		}
		
		addparam("tEA",tEA,tEA,tEA);                                 // We define parameters in the model (with uniform priors)
		addparam("sdEA",tEA/2,tEA/2,tEA/2);
		r = 1.0/tEA; addparam("rEA",r,r,r);
		r = 1.0/tAI; addparam("rAI",r,r,r);
		r = 1.0/tAR; addparam("rAR",r,r,r);
		r = 1.0/tIR; addparam("rIR",r,r,r);
		r = 1.0/tIH; addparam("rIH",r,r,r);
		r = 1.0/tID; addparam("rID",r,r,r);
		r = 1.0/tHD; addparam("rHD",r,r,r);
		r = 1.0/tHR; addparam("rHR",r,r,r);

		addtrans("E","A",GAMMA_DIST,"tEA","sdEA");                   // We define all the transition in the model
		addtrans("A","I",EXP_DIST,"rAI","");
		addtrans("A","R",EXP_DIST,"rAR","");
		addtrans("I","R",EXP_DIST,"rIR","");
		addtrans("I","H",EXP_DIST,"rIH","");
		addtrans("I","D",EXP_DIST,"rID","");
		addtrans("H","D",EXP_DIST,"rHD","");
		addtrans("H","R",EXP_DIST,"rHR","");
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
			cout << comp[c].name << "  Infectivity: " << comp[c].infectivity << endl; 
		}
		cout << endl;
		
		cout << "Transitions:" << endl; 
		for(t = 0; t < trans.size(); t++){
			cout << "  From: " << comp[trans[t].from].name << "  To: " << comp[trans[t].to].name << "  ";
			switch(trans[t].type){
				case INFECTION: 
					cout << " Infection" << endl;
					break;
					
				case EXP_DIST:
					cout << " Exponentially distributed with rate " << param[trans[t].param1].name << endl;
					break;
				case GAMMA_DIST:
					cout << " Gamma distributed with mean " << param[trans[t].param1].name 
							 << " and standard deviation " << param[trans[t].param2].name  << endl; 
					break;
				case LOGNORM_DIST:
					cout << " Lognormally distributed with mean " << param[trans[t].param1].name 
							 << " and standard deviation " << param[trans[t].param2].name  << endl; 
					break;
			}
		}
		cout << endl;
	}
}

/// Gets the transition rate given "from" and "to" compartments
double MODEL::getrate(string from, string to)
{
	int tra;
	TRANS tr;
	
	for(tra = 0; tra < trans.size(); tra++){
		tr = trans[tra];
		if(comp[tr.from].name == from && comp[tr.to].name == to){
			return paramval[tr.param1];
		}
	}
	emsg("Could not find transition");
}

/// Adds a compartment to the model
void MODEL::addcomp(string name, double infectivity)
{
	COMP co;
	co.name = name;
	co.infectivity = infectivity;
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
void MODEL::addtrans(string from, string to, int type, string param1, string param2)
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
	
	tr.type = type;
	
	if(param1 != ""){
		p = 0; pmax = param.size(); while(p < pmax && param1 != param[p].name) p++;
		if(p == pmax) emsg("Cannot find parameter");	
		tr.param1 = p;
	}
	
	if(param2 != ""){
		p = 0; pmax = param.size(); while(p < pmax && param2 != param[p].name) p++;
		if(p == pmax) emsg("Cannot find parameter");	
		tr.param2 = p;
	}
	
	trans.push_back(tr);
}
	
// Converts the spline points to a finer timestep for use in simulations.
void MODEL::betaspline(double period)
{
  int p, s, n = nspline-1;
	double t,  fac, dt, a[n+1], b[n], c[n+1], d[n], h[n], alpha[n], l[n+1], mu[n+1], z[n+1];
	
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

/// This function simulates events as in individual passes through the compartmental model
void MODEL::simmodel(vector < vector <FEV> > &fev, int &tdnext, int &tdfnext, int i, double t, double period, vector <int> &N)
{
	int c, k, kmax, tra;
	double tnext, tnextmin, mean, sd, dt;
	TRANS tr;
	
	switch(modelsel){
	case MOD_IRISH:
		if(ran() < paramval[afracparam]) tra = 0;
		else tra = 1;
		break;
		
	case MOD_OLD:
		tra = 0;
		break;
	}

	addfev(fev,tdnext,tdfnext,t,tra,i,1,period,N);
	c = trans[tra].to;
	
	do{
		kmax = comp[c].trans.size();
		if(kmax == 0) break;
		
		tnextmin = period;
		for(k = 0; k < kmax; k++){
			tr = trans[comp[c].trans[k]];
			switch(tr.type){
				case EXP_DIST:
					tnext = t - log(ran())/paramval[tr.param1];
					break;
				
				case GAMMA_DIST:
					mean = paramval[tr.param1]; sd = paramval[tr.param2];
					dt = gammasamp(mean*mean/(sd*sd),mean/(sd*sd));
					tnext = t + dt; 
					break;
					
				case LOGNORM_DIST:
					mean = paramval[tr.param1]; sd = paramval[tr.param2];
					dt = exp(normal(mean,sd));
					tnext = t + dt; 
					break;
					
				default: emsg("MODEL: EC2b"); break;
			}
			
			if(tnext < tnextmin){ tnextmin = tnext; tra = comp[c].trans[k];}
		}
		
		if(tnextmin == period) break; 

		addfev(fev,tdnext,tdfnext,tnextmin,tra,i,0,period,N);

		c = trans[tra].to; t = tnextmin;
	}while(1 == 1);
}

/// Adds a future event to the timeline
void MODEL::addfev(vector < vector <FEV> > &fev, int &tdnext, int &tdfnext, double t, int tra, int i, int done, double period, vector <int> &N)
{
	int d, j, jmax;
	
	if(t >= period) return;
	
	FEV fe; fe.t = t; fe.trans = tra; fe.ind = i; fe.done = done;
	
	d = int((t/period)*fev.size());
	j = 0; jmax = fev[d].size();
	if(done == 0){ while(j < jmax && t >= fev[d][j].t) j++;}
	else{ while(j < jmax && t >= fev[d][j].t) j++;}
	
	if(j == jmax) fev[d].push_back(fe);
	else fev[d].insert(fev[d].begin()+j,fe);
	
	if(done == 0){
		if(d == tdnext){ if(j < tdfnext) tdfnext = j;}
		if(d < tdnext){ tdnext = d; tdfnext = j;}
	}
	else{
		TRANS tr = trans[tra];
		N[tr.from]--; if(N[tr.from] < 0) emsg("Part: EC12"); 
		N[tr.to]++;
		
		if(d == tdnext) tdfnext++;
	}
}
