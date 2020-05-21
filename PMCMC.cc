// The particle MCMC algorithm

#include <fstream>
#include <iostream>
#include "stdlib.h"
#include "math.h"
#include "assert.h"

#include <mpi.h>

using namespace std;

#include "timers.hh"
#include "model.hh"
#include "utils.hh"
#include "PART.hh"

static void readdata();
static double sample(short core, short ncore);
static double bootstrap();
static double bootstrap_parallel(short core, short ncore);
static double normal(float mu, float sd);

PART* part[partmax];                       // Pointers to each of the particles 
long npart;                                // The number of particles used within this process

long ncase[nregion][tmax/7+1];             // Number of cases in each region as a function of time

ofstream trace;

void PMCMC(MODEL &model, POPTREE &poptree, long nsamp, short core, short ncore)
{
	long p, samp, burnin = nsamp/3, ntr, nac, accept;
	double Li, Lf, valst, al;

	srand(core);
		
	vector <PARAM> &param(model.param);

	readdata();                                                    // Reads in weekly case data 
		
	npart = 5;                                                     // The number of particles per core
	
	for(p = 0; p < npart; p++){ part[p] = new PART(model,poptree); }

	if(core == 0){
		trace.open("trace.txt");		
		trace << "state"; for(p = 0; p < param.size(); p++) trace << "\t" << param[p].name; trace << "\tLi"; trace << endl;
	}
	
	Li = -large; ntr = 0; nac = 0;
	for(samp = 0; samp < nsamp; samp++){
		if(core == 0 && samp%1 == 0) cout << "Sample: " << samp << endl;

		// Each PMCMC step consists of making a change to a parameter 
		// This change is probablisitically either accepted or rejected based
		// on whether the observations agree better with new value or the old one.			
	
		for(p = 0; p < param.size(); p++){
			if(param[p].min != param[p].max){
				valst = param[p].val;
				
				if(core == 0) param[p].val += normal(0,param[p].jump); 
					
				MPI_Bcast(&param[p].val,1,MPI_DOUBLE,0,MPI_COMM_WORLD);

				if(param[p].betachange == 1) model.betaspline();
				if(param[p].suschange == 1) poptree.setsus(model);
				if(param[p].infchange == 1) poptree.setinf(model);
				
				if(param[p].val < param[p].min || param[p].val > param[p].max) al = 0;
				else{
					Lf = sample(core,ncore);
					if(core == 0) al = exp(Lf-Li);
				}
				if(core == 0){ if(ran() < al) accept = 1; else accept = 0;}
				
				MPI_Bcast(&accept,1,MPI_LONG,0,MPI_COMM_WORLD);
					
				param[p].ntr++;
				if(accept == 1){
					param[p].nac++;
					
					Li = Lf;
					if(samp < burnin) param[p].jump *= 1.1;
				}
				else{
					param[p].val = valst;
					if(param[p].betachange == 1) model.betaspline();
					if(param[p].suschange == 1) poptree.setsus(model);
					if(param[p].infchange == 1) poptree.setinf(model);
				
					if(samp < burnin) param[p].jump *= 0.95;
				}
			}
		}
	
		Lf = sample(core,ncore);
		if(core == 0){ al = exp(Lf-Li); ntr++; if(ran() < al){ nac++; Li = Lf;}}
	
		if(core == 0){
			trace << samp; for(p = 0; p < param.size(); p++) trace << "\t" << param[p].val; trace << "\t" << Li; trace << endl;
		}			
	}
		
	if(core == 0){      // This gives the acceptance rates for different MCMC proposals on different parameters
		cout << "MCMC diagnostics:" << endl;
		cout << "Base acceptance rate " << double(nac)/ntr << endl;
		for(p = 0; p < param.size(); p++){
			cout << param[p].name << ": ";
			if(param[p].ntr == 0) cout << "Fixed" << endl;
			else cout << "Acceptance rate " << double(param[p].nac)/param[p].ntr << endl;
		}
	}
}

/// This samples from the model using particles and returns an overall measure of how well the observations agreed with it 
static double sample(short core, short ncore)
{
	short p, step = 7;
	double tt, ttnext, Liav;

	for(p = 0; p < npart; p++) part[p]->partinit(p);

	tt = 0; Liav = 0;
	do{                                // We step through one measurement at a time
		ttnext = tt+step; if(ttnext > tmax) tt = tmax;

		for(p = 0; p < npart; p++){
			timers.timesim -= clock();
			part[p]->gillespie(tt,ttnext, 0 /* Inference */); // Simulates the particle
			part[p]->Lobs(tt,ttnext, ncase);      // Measures how well it agrees with the observations (weekly number of cases)
			timers.timesim += clock();
		}
		timers.timeboot -= clock();
		
		//Liav += bootstrap();             // Culls or copies particles based on how well they represent observations
		
		Liav += bootstrap_parallel(core,ncore);  // Culls or copies particles based on how well they represent observations
		timers.timeboot += clock();
			
		tt = ttnext;
	}while(tt < tmax);	
	
	return Liav;
}

/// This step culls some particles (which are not agreeing well with the observations) and copies others 
static double bootstrap_parallel(short core, short ncore)
{
	const long nparttot = npart*ncore;	
	long p, pp, ppp, k, kmax, j, co, si;
	short	copynum[nparttot], flag[nparttot], num, n;
	double res;
	vector <short> flaglist;
	vector < vector <short> > corelist;
	vector <short> colist;
	vector <double> pac; 

	if(core == 0){
		long pp, co;
		double Limax, av, z, Litot[nparttot], sum, sumst[nparttot], flag[nparttot];
		
		for(p = 0; p < npart; p++) Litot[p] = part[p]->Li;                          // Gathers together the Likelihoods
	
		for(co = 1; co < ncore; co++){
			MPI_Recv(Litot+co*npart,npart,MPI_DOUBLE,co,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
		}
		
		Limax = -large;
		for(p = 0; p < nparttot; p++){
			if(Litot[p] > Limax) Limax = Litot[p];
		}
		
		av = 0;	for(p = 0; p < nparttot; p++) av += exp(Litot[p] - Limax);
		av /= nparttot;

		sum = 0;
		for(p = 0; p < nparttot; p++){
			sum += exp(Litot[p] - Limax);
			sumst[p] = sum;
		}

		for(p = 0; p < nparttot; p++) copynum[p] = 0;
			
		for(p = 0; p < nparttot; p++){                                      // Calculates how many copies of particle to be kept
			z = sum*ran();
			pp = 0; while(pp < nparttot && z > sumst[pp]) pp++;
			if(pp == nparttot) emsg("PMCMC: EC1");	
			
			copynum[pp]++;
		}

		res = Limax + log(av);
	}
	else{
		double Li[npart];
		for(p = 0; p < npart; p++) Li[p] = part[p]->Li;
		MPI_Send(Li,npart,MPI_DOUBLE,0,0,MPI_COMM_WORLD);	
		res = 0;
	}
	
	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Bcast(copynum,nparttot,MPI_SHORT,0,MPI_COMM_WORLD);

	for(p = 0; p < nparttot; p++) flag[p] = 0;
	
	for(p = 0; p < nparttot; p++){                             // Particles just stay the same
		if(copynum[p] > 0){
			flag[p] = 1; copynum[p]--;
		}		
	}

	for(p = 0; p < nparttot; p++){                             // Performs internal copies
		if(copynum[p] > 0){            
			co = p/npart;
			for(pp = co*npart; pp < (co+1)*npart; pp++){  
				if(flag[pp] == 0){
					flag[pp] = 1;
					if(co == core){
						part[pp%npart]->copy(*part[p%npart]);
					}
					copynum[p]--; if(copynum[p] == 0) break;
				}
			}
		}
	}
	
	for(p = 0; p < nparttot; p++){ if(flag[p] == 0) flaglist.push_back(p);}

	corelist.resize(ncore);
	for(p = 0; p < nparttot; p++){                             // Performs copying between cores
		num = copynum[p];
		if(num > 0){      		
			colist.clear();
			for(n = 0; n < num; n++){
				if(flaglist.size() == 0) emsg("pmcmc: EC10");
				pp = flaglist[long(flaglist.size())-1]; flaglist.pop_back();
				co = pp/npart;
				if(corelist[co].size() == 0) colist.push_back(co);
				corelist[co].push_back(pp);
			}

			co = p/npart;
			if(core == co){
				part[p%npart]->partpack(pac);
				si = pac.size();
				for(k = 0; k < colist.size(); k++){
					MPI_Send(&si,1,MPI_LONG,colist[k],0,MPI_COMM_WORLD);	
					MPI_Send(&pac[0],si,MPI_DOUBLE,colist[k],0,MPI_COMM_WORLD);	
				}
			}
			else{
				if(corelist[core].size() > 0){
					MPI_Recv(&si,1,MPI_LONG,co,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
					pac.resize(si);
					MPI_Recv(&pac[0],si,MPI_DOUBLE,co,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
					pp = corelist[core][0];
					part[pp%npart]->partunpack(pac);
					for(k = 1; k < corelist[core].size(); k++){
						ppp = corelist[core][k];
						part[ppp%npart]->copy(*part[pp%npart]);                  // Copies multiple version of particle
					}
				}
			}
			for(k = 0; k < colist.size(); k++) corelist[colist[k]].clear();
		}
	}
	if(flaglist.size() != 0) emsg("pmcmc: EC11");

	return res;
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
			part[p]->copy(*part[pp]);
			
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

/// Draws a normally distributed number with mean mu and standard deviation sd
static double normal(float mu, float sd)
{
	return mu + sd*sqrt(-2*log(ran()))*cos(2*M_PI*ran());
}
