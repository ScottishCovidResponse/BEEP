// This file defines the compartmental model and is used to simulate compartmental dynamics after an individual becomes infected

// All the code related to the compartmental model (i.e. what happens after exposure)

#include <string>
#include <sstream>
#include <iostream>

#include "math.h"

#include "consts.hh"
#include "utils.hh"
#include "model.hh"

using namespace std;

void MODEL::definemodel()
{
	// R0 determines how many individuals an infected individual on average infects
	// These five value represent how R0 changes over time in the simulation (this captures the effect of lockdown) 
	double R0sim[5] = {3,2.5,1.5,0.5,0.5};       
	
	double tEA = 2.5, tAI = 3, tAR = 3.5, tIR = 16.7;            // Estimates of transition times from literature
	double tIH = 3.1, tID = 12.9, tHD = 8, tHR = 13.1;
	
	double r, rAI, rAR, rIR, rIH, rID, tinfav;
	short p, c, t;
		
	addcomp("E",0); addcomp("A",0.2); addcomp("I",1);            // Different compartment in the model
	addcomp("H",0); addcomp("R",0); addcomp("D",0); 
	
	// This calculates the average integral of infectivity after an indvidual becomes infected
	rAI = 1.0/tAI; rAR = 1.0/tAR; rIH = 1.0/tIH; rIR = 1.0/tIR; rID = 1.0/tID;
	tinfav = 0.2/(rAI + rAR) + 1*(rAI/(rAI+rAR))*(1.0/(rIH + rIR + rID));
		
	nspline = 5;                                                 // 5 spline points represent time variarion in beta 
	for(p = 0; p < nspline; p++){
		splinet.push_back(double(p*tmax)/(nspline-1));
		stringstream ss; ss << "beta_" << p;
		r = R0sim[p]/tinfav; addparam(ss.str(),r,0,3*r);
	}		
	betaspline();

	addparam("tEA",tEA,tEA,tEA);                             // We define all the parameters in the model (with uniform priors)
	addparam("sdEA",tEA/2,tEA/2,tEA/2);
	r = 1.0/tEA; addparam("rEA",r,r,r);
	r = 1.0/tAI; addparam("rAI",r,r,r);
	r = 1.0/tAR; addparam("rAR",r,r,r);
	r = 1.0/tIR; addparam("rIR",r,r,r);
	r = 1.0/tIH; addparam("rIH",r,r,r);
	r = 1.0/tID; addparam("rID",r,r,r);
	r = 1.0/tHD; addparam("rHD",r,r,r);
	r = 1.0/tHR; addparam("rHR",r,r,r);

	addtrans("E","A",GAMMA_DIST,"tEA","sdEA");               // We define all the transition in the model
	addtrans("A","I",EXP_DIST,"rAI","");
	addtrans("A","R",EXP_DIST,"rAR","");
	addtrans("I","R",EXP_DIST,"rIR","");
	addtrans("I","H",EXP_DIST,"rIH","");
	addtrans("I","D",EXP_DIST,"rID","");
	addtrans("H","D",EXP_DIST,"rHD","");
	addtrans("H","R",EXP_DIST,"rHR","");
	
	cout << endl;                                            //Outputs a summary of the model
	cout << "Parameters:" << endl;
	for(p = 0; p < param.size(); p++){
		cout << param[p].name << " " << param[p].val << " (" << param[p].min << " - " << param[p].max << ")" << endl;
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
			case EXP_DIST:
				cout << " Exponentially distributed with rate " << param[trans[t].param1].name << endl;
				break;
			case GAMMA_DIST:
				cout << " Gamma distributed with mean " << param[trans[t].param1].name 
						<< " and standard deviation " << param[trans[t].param2].name  << endl; 
				break;
		}
	}
	cout << endl;
}

void MODEL::addcomp(string name, double infectivity)
{
	COMP co;
	co.name = name;
	co.infectivity = infectivity;
	comp.push_back(co);	
}

void MODEL::addparam(string name, double val, double min, double max)
{
	PARAM par;
	par.name = name; par.val = val; par.sim = val; par.min = min; par.max = max; par.jump = val/10; par.ntr = 0; par.nac = 0;

	param.push_back(par);
}

void MODEL::addtrans(string from, string to, short type, string param1, string param2)
{
	short c, cmax, p, pmax;
	TRANS tr;
	
	c = 0; cmax = comp.size(); while(c < cmax && from != comp[c].name) c++;
	if(c == cmax) emsg("Cannot find compartment");
	tr.from = c;	
	comp[c].trans.push_back(trans.size());
	
	c = 0; cmax = comp.size(); while(c < cmax && to != comp[c].name) c++;
	if(c == cmax) emsg("Cannot find compartment");	
	tr.to = c;
	
	tr.type = type;
	
	p = 0; pmax = param.size(); while(p < pmax && param1 != param[p].name) p++;
	if(p == pmax) emsg("Cannot find parameter");	
	tr.param1 = p;
	
	if(param2 !=""){
		p = 0; pmax = param.size(); while(p < pmax && param2 != param[p].name) p++;
		if(p == pmax) emsg("Cannot find parameter");	
		tr.param2 = p;
	}
	
	trans.push_back(tr);
}
	
// Converts the spline points to a finer timestep for use in simulations.
// At the moment the spline is just linear, but it will probably become cubic at some point.
void MODEL::betaspline()
{
  short p, s, n = nspline-1;
	double t,  fac, dt, a[n+1], b[n], c[n+1], d[n], h[n], alpha[n], l[n+1], mu[n+1], z[n+1];
	
	if(1 == 1){   // This uses a cubic spline
		for(p = 0; p <= n; p++) a[p] = log(param[p].val);
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
			settime[s] = double((s+1)*tmax)/nsettime;;
			
			t = double((s+0.5)*tmax)/nsettime;
			while(p < nspline-1 && t > splinet[p+1]) p++;
			
			dt = t-splinet[p];	
			beta[s] = exp(a[p]+ b[p]*dt + c[p]*dt*dt + d[p]*dt*dt*dt);
		}
	}
	else{  // This uses a linear spline
		p = 0;
		for(s = 0; s < nsettime; s++){
			settime[s] = double((s+1)*tmax)/nsettime;;
			
			t = double((s+0.5)*tmax)/nsettime;
			while(p < nspline-1 && t > splinet[p+1]) p++;
			
			fac = (t-splinet[p])/(splinet[p+1]-splinet[p]);
			beta[s] = param[p].val*(1-fac) + param[p+1].val*fac;
		}
	}
}
