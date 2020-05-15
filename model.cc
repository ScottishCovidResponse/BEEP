// This file defines the compartmental model and is used to simulate compartmental dynamics after an individual becomes infected

// All the code related to the compartmental model (i.e. what happens after exposure)

#include <string>
#include <sstream>
#include <iostream>

#include "math.h"

#include "var.hh"
#include "consts.hh"
#include "functions.hh"

using namespace std;

void addcomp(string name, double infectivity);
void addparam(string name, double val, double min, double max);
void addtrans(string from, string to, short type, string param1, string param2);
void betaspline();

void definemodel()
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
	
	cout << "\n";                                            //Outputs a summary of the model
	cout << "Parameters:\n";
	for(p = 0; p < param.size(); p++){
		cout << param[p].name << " " << param[p].val << " (" << param[p].min << " - " << param[p].max << ")\n";
	}
	cout << "\n";
		
	cout << "Compartments:\n"; 
	for(c = 0; c < comp.size(); c++){
		cout << comp[c].name << "  Infectivity: " << comp[c].infectivity << "\n"; 
	}
	cout << "\n";
	
	cout << "Transitions:\n"; 
	for(t = 0; t < trans.size(); t++){
		cout << "  From: " << comp[trans[t].from].name << "  To: " << comp[trans[t].to].name << "  ";
		switch(trans[t].type){
			case EXP_DIST:
				cout << " Exponentially distributed with rate " << param[trans[t].param1].name << "\n";
				break;
			case GAMMA_DIST:
				cout << " Gamma distributed with mean " << param[trans[t].param1].name 
						<< " and standard deviation " << param[trans[t].param2].name  << "\n"; 
				break;
		}
	}
	cout << "\n";
}

void addcomp(string name, double infectivity)
{
	COMP co;
	co.name = name;
	co.infectivity = infectivity;
	comp.push_back(co);	
}

void addparam(string name, double val, double min, double max)
{
	PARAM par;
	par.name = name; par.val = val; par.sim = val; par.min = min; par.max = max; par.jump = val/10; par.ntr = 0; par.nac = 0;

	param.push_back(par);
}

void addtrans(string from, string to, short type, string param1, string param2)
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
void betaspline()
{
	short s, p;
	double t, fac;
	
	p = 0;
	for(s = 0; s < nsettime; s++){
		settime[s] = double((s+1)*tmax)/nsettime;;
		
		t = double((s+0.5)*tmax)/nsettime;
		while(p < nspline-1 && t > splinet[p+1]) p++;
		
		fac = (t-splinet[p])/(splinet[p+1]-splinet[p]);
		beta[s] = param[p].val*(1-fac) + param[p+1].val*fac;
	}
}

double gammasamp(double a, double b)             // Draws a sample from the gamma distribution x^(a-1)*exp(-b*x)
{
  if(a < 0 || b < 0) emsg("Model: EC1");

  if(a < 1){
    double u = ran();
    return gammasamp(1.0 + a, b) * pow (u, 1.0 / a);
  }
  else{
    double x, v, u;
    double d = a - 1.0 / 3.0;
    double c = (1.0 / 3.0) / sqrt (d);
 
    while(1 == 1){
      do{
        x = sqrt(-2*log(ran()))*cos(2*M_PI*ran());
        v = 1.0 + c * x;
      }while (v < 0);

      v = v*v*v;
      u = ran();

      if (u < 1 - 0.0331*x*x*x*x) break;

      if (log(u) < 0.5*x*x + d*(1 - v + log(v))) break;
    }

    return d*v/b;
  }
}

// Adds an exposed indivdual on node c on the finest scale (i.e. level-1)
void PART::addinfc(long c, double t)
{
	long l, i, cc, k, kmax;
	double dR, sum;

	kmax = indinf[c].size();
	do{
		l = long(ran()*long(subpop[c].size()));
		i = subpop[c][l];
		for(k = 0; k < kmax; k++) if(indinf[c][k] == i) break;
	}while(k < kmax);
	indinf[c].push_back(i);
	
	l = level-1; cc = c;
	dR = -Rtot[l][c]/pop[l][cc];
	do{
		Rtot[l][cc] += dR;
		pop[l][cc]--;
		cc = lev[l].node[cc].parent; l--;
	}while(l >= 0);
	simmodel(i,0,t);
}

// Once an individual goes into the exposed class, this function simulates all the subsequent future events
void PART::simmodel(long i, short enter, double t)
{
	short c, k, kmax, tra;
	double tnext, tnextmin, mean, sd;
	TRANS tr;

	N[enter]++;
	double dt;
	c = enter;
	do{
		kmax = comp[c].trans.size();
		if(kmax == 0) break;
		
		tnextmin = large;
		for(k = 0; k < kmax; k++){
			tr = trans[comp[c].trans[k]];
			switch(tr.type){
				case EXP_DIST:
					tnext = t - log(ran())/param[tr.param1].val;
					break;
				
				case GAMMA_DIST:
					mean = param[tr.param1].val; sd = param[tr.param2].val;
					dt = gammasamp(mean*mean/(sd*sd),mean/(sd*sd));
					tnext = t + dt; 
					break;
			}
			
			if(tnext < tnextmin){ tnextmin = tnext; tra = comp[c].trans[k];}
		}
		if(tnextmin == large) emsg("Model: EC2");
		
		addfev(tnextmin,tra,i);

		c = trans[tra].to; t = tnextmin;
	}while(1 == 1);
}

// Adds a future event to the timeline
void PART::addfev(double t, long trans, long i)
{
	long d, j, jmax;
	
	if(t >= tmax) return;
	
	FEV fe; fe.t = t; fe.trans = trans; fe.ind = i; fe.done = 0;
	
	d = long((t/tmax)*fediv);
	j = 0; jmax = fev[d].size();
	while(j < jmax && t > fev[d][j].t) j++;
	if(j == jmax) fev[d].push_back(fe);
	else fev[d].insert(fev[d].begin()+j,fe);
	
	if(d == tdnext){ if(j < tdfnext) tdfnext = j;}
	if(d < tdnext){ tdnext = d; tdfnext = j;}
}
