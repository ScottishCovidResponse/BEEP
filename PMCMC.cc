// The particle MCMC algorithm

#include <fstream>
#include <iostream>
#include "math.h"
#include "assert.h"

using namespace std;

#include "var.hh"
#include "model.hh"
#include "functions.hh"
#include "PART.hh"

static void readdata();
static double sample();
static double bootstrap();

void PMCMC()
{
	long p, samp, burnin = nsamp/3;
	double Li, Lf, valst, al;

	readdata();                                                    // Reads in weekly case data 
	 
	npart = 100;  // This is the number of particles (which much be sufficiently large for the simulations to capture the data)
	for(p = 0; p < npart; p++){ part[p] = new PART; }
	
	ofstream trace("trace.txt");
	trace << "state"; for(p = 0; p < param.size(); p++) trace << "\t" << param[p].name; trace << "\n";
	
	Li = -large;
	for(samp = 0; samp < nsamp; samp++){
		if(samp%1 == 0) cout << "Sample: " << samp << "\n";

		// Each PMCMC step consists of making a change to a parameter 
		// This change is probablisitically either accepted or rejected based
		// on whether the observations agree better with new value or the old one.
		
		for(p = 0; p < param.size(); p++){
			if(param[p].min != param[p].max){
				valst = param[p].val;
				param[p].val += normal(0,param[p].jump); if(p < nspline) betaspline();
				if(param[p].val < param[p].min || param[p].val > param[p].max) al = 0;
				else{
					Lf = sample();
					al = exp(Lf-Li);
				}
				param[p].ntr++;
				if(ran() < al){
					param[p].nac++;
					
					Li = Lf;
					if(samp < burnin) param[p].jump *= 1.1;
				}
				else{
					param[p].val = valst; if(p < nspline) betaspline();
					if(samp < burnin) param[p].jump *= 0.95;
				}
			}
		}
	
		trace << samp; for(p = 0; p < param.size(); p++) trace << "\t" << param[p].val; trace << "\n";	
	}
	
	// This gives the acceptance rates for different MCMC proposals on different parameters
	cout << "MCMC diagnostics:\n";
	for(p = 0; p < param.size(); p++){
		cout << param[p].name << ": ";
		if(param[p].ntr == 0) cout << "Fixed\n";
		else cout << "Acceptance rate " << double(param[p].nac)/param[p].ntr << "\n";
	}
}

/// This samples from the model using particles and returns an overall measure of how well the observations agreed with it 
static double sample()
{
	short p, step = 7;
	double tt, ttnext, Liav;

	for(p = 0; p < npart; p++) part[p]->partinit(p);

	tt = 0; Liav = 0;
	do{                                // We step through one measurement at a time
		ttnext = tt+step; if(ttnext > tmax) tt = tmax;

		for(p = 0; p < npart; p++){
			timesim -= clock();
			assert(siminf == 0);

			part[p]->gillespie(tt,ttnext, 0 /* Inference */); // Simulates the particle
			part[p]->Lobs(tt,ttnext);      // Measures how well it agrees with the observations (weekly number of cases)
			timesim += clock();
		}
		timeboot -= clock();
		Liav += bootstrap();             // Culls or copies particles based on how well they represent observations
		timeboot += clock();
			
		tt = ttnext;
	}while(tt < tmax);	

	return Liav;
}

/// This step culls some particles (which are not agreeing well with the observations) and copies others 
static double bootstrap()
{
	long p, pp;
	double Limax, av, z, sum, sumst[npart], flag[npart];
	vector <long> copylist;
	
	Limax = -large;
	for(p = 0; p < npart; p++){
		if(part[p]->Li > Limax) Limax = part[p]->Li;
	}
	
	av = 0;	for(p = 0; p < npart; p++) av += exp(part[p]->Li - Limax);
	av /= npart;
	
	for(p = 0; p < npart; p++) flag[p] = 0;
	
	sum = 0;
	for(p = 0; p < npart; p++){
		sum += exp(part[p]->Li - Limax);
		sumst[p] = sum;
	}
	
	for(p = 0; p < npart; p++){
		z = sum*ran();
		pp = 0; while(pp < npart && z > sumst[pp]) pp++;
		if(pp == npart) emsg("PMCMC: EC1");
		
		if(flag[pp] == 0) flag[pp] = 1;
		else copylist.push_back(pp);
	}
	
	for(p = 0; p < npart; p++){
		if(flag[p] == 0){
			if(copylist.size() == 0) emsg("PMCMC: EC2");
			pp = copylist[long(copylist.size())-1];
			part[p]->copy(pp);
			
			copylist.pop_back();
		}
	}	
	if(copylist.size() != 0) emsg("PMCMC: EC3");
 
	return Limax + log(av);
}

/// Reads in simulated case data
static void readdata()
{
	long week, tt, r;
	
	ifstream regplot("Weekly case data.txt");
	for(week = 0; week < tmax/7; week++){
		regplot >> tt;
		for(r = 0; r < nregion; r++) regplot >> ncase[r][week];
	}
}
