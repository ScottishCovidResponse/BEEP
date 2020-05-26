// The particle MCMC algorithm

#include <fstream>
#include <iostream>
#include <sstream>
#include "stdlib.h"
#include "math.h"
#include "assert.h"

using namespace std;

#include "timers.hh"
#include "model.hh"
#include "utils.hh"
#include "PART.hh"
#include "output.hh"

const int MAX_NUMBERS = 10000000;
double buffer[MAX_NUMBERS];

long timebo = 0; // temporary variables

static double sample(short core, short ncore, short npart, double tmax, double invT,	vector < vector <FEV> > &xinew);
static double bootstrap(short core, short ncore, short npart, short fedivmin, short w, short *back);
void geneventsample(short core, short ncore, short npart, short nweek, short *back, vector < vector <FEV> > &xif);
static double normal(float mu, float sd);

PART* part[partmax];                       // Pointers to each of the particles 

void PMCMC(DATA &data, MODEL &model, POPTREE &poptree, long nsamp, short core, short ncore, short npart)
{
	long p, samp, burnin = nsamp/3, accept, period=10;
	double Li, Lf, valst, al, ac = 0;
	double invT = 0.1;                 // The inverse temperature (used to relax the observation model)
	vector < vector <FEV> > xi, xif;   // Stores the current and proposed event sequences
	vector <SAMPLE> opsamp;            // Stores sample for generating satatistics later

	srand(core);
	srand(core);
	srand(core);
		
	vector <PARAM> &param(model.param);
                                                   // The number of particles per core	
	for(p = 0; p < npart; p++){ part[p] = new PART(data,model,poptree);}

	if(core == 0) outputinit(model,ncore*npart);
			
	MPI_Barrier(MPI_COMM_WORLD);
				
	Li = -large; 
	for(samp = 0; samp < nsamp; samp++){
		if(core == 0 && samp%1 == 0) cout << "Sample: " << samp << " " << invT << endl;

		if(core == 0 && samp < burnin){                      // Dynamically alters inverse temperature
			if(samp != 0 && samp%period == 0){
				if(ac/period < 0.6) invT *= 0.8;
				if(ac/period > 0.8) invT *= 1.2;
				ac = 0; Li = -large;
			}
		}
		
		Lf = sample(core,ncore,npart,data.tmax,invT,xif);

		if(core == 0){ 
			al = exp(Lf-Li); model.ntr++; 
			if(ran() < al){ xi = xif; model.nac++; ac++; Li = Lf;}
		}

		// Each PMCMC step consists of making a change to a parameter 
		// This change is probablisitically either accepted or rejected based
		// on whether the observations agree better with new value or the old one.			
		for(p = 0; p < long(param.size()); p++){
			if(param[p].min != param[p].max){
				valst = param[p].val;
				
				if(core == 0) param[p].val += normal(0,param[p].jump); 
				
				MPI_Bcast(&param[p].val,1,MPI_DOUBLE,0,MPI_COMM_WORLD);

				if(param[p].betachange == 1) model.betaspline(data.tmax);
				if(param[p].suschange == 1) poptree.setsus(model);
				if(param[p].infchange == 1) poptree.setinf(model);
				
				if(param[p].val < param[p].min || param[p].val > param[p].max) al = 0;
				else{
					Lf = sample(core,ncore,npart,data.tmax,invT,xif);
					if(core == 0){ al = exp(Lf-Li); cout << Lf << " " << Li << " " << al << "al\n";}
				}
				if(core == 0){ if(ran() < al) accept = 1; else accept = 0;}
				
				MPI_Bcast(&accept,1,MPI_LONG,0,MPI_COMM_WORLD);
					
				param[p].ntr++;
				if(accept == 1){
					if(core == 0){
						param[p].nac++;
						xi = xif;
						Li = Lf;
						if(samp < burnin) param[p].jump *= 1.1;
					}
				}
				else{
					param[p].val = valst;
					if(param[p].betachange == 1) model.betaspline(data.tmax);
					if(param[p].suschange == 1) poptree.setsus(model);
					if(param[p].infchange == 1) poptree.setinf(model);
				
					if(core == 0){
						if(samp < burnin) param[p].jump *= 0.95;
					}
				}
			}
		}
	
		if(core == 0) opsamp.push_back(outputsamp(invT,samp,Li,data,model,poptree,xi));	
	}
	
	if(core == 0) cout << double(timebo)/CLOCKS_PER_SEC << " Timebo" << endl;
	
	if(core == 0) outputresults(data,model,opsamp,0,ncore*npart);
	if(core == 0) cout << part[0]->Rtot[0][0] << " " << invT << " invT\n";
}

/// This samples from the model using particles and returns an overall measure of how well the observations agreed with it 
static double sample(short core, short ncore, short npart, double tmax, double invT, vector < vector <FEV> > &xif)
{
	short p, w, nweek = tmax/timestep, back[ncore*nweek];
	double Liav;
		
	for(p = 0; p < npart; p++) part[p]->partinit(p);

  Liav = 0;
	for(w = 0; w < nweek; w++){                                // We step through one measurement at a time
		for(p = 0; p < npart; p++){
			timers.timesim -= clock();
	
			part[p]->gillespie(w*timestep,(w+1)*timestep, 0 /* Inference */); // Simulates the particle
	
			part[p]->Lobs(w,invT);      // Measures how well it agrees with the observations (weekly number of cases)
			timers.timesim += clock();
		}
		
		timers.timeboot -= clock();
		
		// Culls or copies particles based on how well they represent observations
		Liav += bootstrap(core,ncore,npart,(fediv*(w+1))/nweek,w,back);
		timers.timeboot += clock();
	}
	
	geneventsample(core,ncore, npart,nweek,back,xif);
	
	return Liav;
}

 // Constructs the event sequence sample by gathering all the pieces from different particles
void geneventsample(short core, short ncore, short npart, short nweek, short *back, vector < vector <FEV> > &xif)	
{
	short p, w, co, d, fedivmin, fedivmax;
	long k;
	int siz;

	xif.resize(fediv);
	for(w = nweek-1; w >= 0; w--){
		fedivmin = (fediv*w)/nweek;	fedivmax = (fediv*(w+1))/nweek;
		if(w == nweek-1) p = back[ncore*w];  // Picks the first final particle
		else p = back[ncore*w + p/npart];
		
		co = p/npart;
		if(core == 0){  
			if(co == 0){
				for(d = fedivmin; d < fedivmax; d++) xif[d] = part[p%npart]->fev[d];
			}
			else{
				MPI_Status status;
				MPI_Recv(buffer,MAX_NUMBERS,MPI_DOUBLE,co,0,MPI_COMM_WORLD,&status);
				MPI_Get_count(&status, MPI_DOUBLE, &siz); if(siz >= MAX_NUMBERS) emsg("Buffer not big enough");
		
				k = 0;
				part[0]->unpack(k,buffer,xif,fedivmin,fedivmax);
				if(k != siz) emsg("PMCMC: EC10");
			}
		}
		else{
			if(co == core){
				siz = part[p%npart]->fevpack(buffer,fedivmin,fedivmax);
				MPI_Send(buffer,siz,MPI_DOUBLE,0,0,MPI_COMM_WORLD);	
			}
		}
	}
}

static double bootstrap(short core, short ncore, short npart, short fedivmin, short w, short *back)
{
	long p, psel, co, co2, j, k, nparttot = npart*ncore;
	double res; 
	int siz, flag;
	back += w*ncore;	
	
	if(core == 0){
		short sel[ncore];
		double Limax, av, z, Litot[nparttot], sumind, sumindst[npart], sum, sumst[ncore];
		vector <long> extra;
		
		for(p = 0; p < npart; p++) Litot[p] = part[p]->Li;                          // Gathers together the Likelihoods
	
		for(co = 1; co < ncore; co++){
			MPI_Recv(Litot+npart*co,npart,MPI_DOUBLE,co,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
		}
		
		Limax = -large;	for(p = 0; p < nparttot; p++){ if(Litot[p] > Limax) Limax = Litot[p];}
		
		av = 0;	for(p = 0; p < nparttot; p++) av += exp(Litot[p] - Limax); av /= ncore;

		res = Limax + log(av);
		
		sum = 0;
		for(co = 0; co < ncore; co++){
			sumind = 0;
			for(p = 0; p < npart; p++){
				sumind += exp(Litot[co*npart+p] - Limax);
				sumindst[p] = sumind;
			}
			z = sumind*ran(); p = 0; while(p < npart && z > sumindst[p]) p++;	if(p == npart) emsg("PMCMC: EC9");
			sel[co] = p; 
			
			sum += sumind;
			sumst[co] = sum;
		}

		for(co = 0; co < ncore; co++) back[co] = -1;
			
		timebo -= clock();
		for(co = 0; co < ncore; co++){                                      // Finds back for a core
			z = sum*ran();
			co2 = 0; while(co2 < ncore && z > sumst[co2]) co2++;
			if(co2 == ncore) emsg("PMCMC: EC1");	
			
			if(back[co2] == -1) back[co2] = co2*npart + sel[co2];
			else extra.push_back(co2);
		}

		for(co = 0; co < ncore; co++){  
			if(back[co] == -1){
				co2 = extra[extra.size()-1]; extra.pop_back();
				back[co] =  co2*npart + sel[co2];
			}
		}
		if(extra.size() != 0) emsg("PMCMC: EC9");
	}
	else{
		double Li[npart];
		for(p = 0; p < npart; p++) Li[p] = part[p]->Li;
		MPI_Send(Li,npart,MPI_DOUBLE,0,0,MPI_COMM_WORLD);	
		res = 0;
	}
	
	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Bcast(back,ncore,MPI_SHORT,0,MPI_COMM_WORLD);
	
	co = back[core]/npart;
	if(co == core){          // Send
		flag = 0;
		for(co = 0; co < ncore; co++){
			co2 = back[co]/npart;
			if(co2 == core){
				psel = back[co]%npart; 
				if(co2 == co){
					for(p = 0; p < npart; p++){
						if(p != psel) part[p]->copy(*part[psel],fedivmin);
					}
				}
				else{
					if(flag == 0){ siz = part[psel]->partpack(buffer,fedivmin); flag = 1;}
					MPI_Send(buffer,siz,MPI_DOUBLE,co,0,MPI_COMM_WORLD);	
				}
			}
		}
	}
	else{                  // Recieve
		MPI_Status status;
		MPI_Recv(buffer,MAX_NUMBERS,MPI_DOUBLE,co,0,MPI_COMM_WORLD,&status);
		MPI_Get_count(&status, MPI_DOUBLE, &siz); if(siz >= MAX_NUMBERS) emsg("Buffer not big enough");
		part[0]->partunpack(buffer,siz,fedivmin);
		for(p = 1; p < npart; p++) part[p]->copy(*part[0],fedivmin);
	}
	
	return res;
}

/// Draws a normally distributed number with mean mu and standard deviation sd
static double normal(float mu, float sd)
{
	return mu + sd*sqrt(-2*log(ran()))*cos(2*M_PI*ran());
}
