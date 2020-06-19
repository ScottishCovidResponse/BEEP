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
#include "pack.hh"
#include "consts.hh"
#include "obsmodel.hh"

MPI_Request reqs[partmax];                 // These are information used Isend and Irecv
MPI_Status stats[partmax];

vector < vector<double> > sendbuffer;
vector < vector<double> > recibuffer;
	
static double sample(DATA &data, MODEL &model, POPTREE &poptree, unsigned int core, unsigned int ncore, unsigned int npart, unsigned int period,	vector < vector <FEV> > &xinew);
static double bootstrap(unsigned int core, unsigned int ncore, unsigned int npart, unsigned int fedivmin, unsigned int w, unsigned int *backpart);
static void geneventsample(DATA &data, unsigned int core, unsigned int ncore, unsigned int npart, unsigned int nweek, unsigned int *backpart, vector < vector <FEV> > &xp);

PART* part[partmax];                       // Pointers to each of the particles 

void PMCMC(DATA &data, MODEL &model, POPTREE &poptree, unsigned int nsamp, unsigned int core, unsigned int ncore, unsigned int npart)
{
	const unsigned int fixinvT = 1;          // Detemines if the inverse temperature invT is fixed or is dynamically tuned 
	double invT = 0.1;                       // The inverse temperature (used to relax the observation model)
	vector < vector <FEV> > xi, xp;          // Stores the current and proposed event sequences
	vector <SAMPLE> opsamp;                  // Stores sample for generating statistics later
	vector <PARAM> &param(model.param);      // Copies parameters from the model
	unsigned int p, samp, burnin = nsamp/4, accept;
	double Li, Lf, valst, al;

	sendbuffer.resize(npart);                // Sets up the buffers for sending and receiving messages via MPI
	recibuffer.resize(npart);
	for(p = 0; p < npart; p++){
		sendbuffer[p].resize(BUFMAX);
		recibuffer[p].resize(BUFMAX);
	}

	srand(core);
	
	for(p = 0; p < npart; p++){ part[p] = new PART(data,model,poptree);}

	if(core == 0) outputinit(data,model);

	MPI_Bcast(&model.paramval[0],model.paramval.size(),MPI_DOUBLE,0,MPI_COMM_WORLD);
	model.betaspline(data.period);
	model.setsus(data);
 
	Li = -large; 
	for(samp = 0; samp < nsamp; samp++){
		if(core == 0 && samp%1 == 0) cout << "Sample: " << samp << " / " << nsamp << endl;
 
		model.parami = model.paramval; model.paramp = model.paramval;
		model.settransprob();
		
		Lf = sample(data,model,poptree,core,ncore,npart,data.period,xp);

		if(core == 0){ 
			al = exp(invT*(Lf-Li)); model.ntr++; //cout << Lf << " " << Li << " " << al << "al self acceptance rate" << endl;
			if(ran() < al){ 
				xi = xp; model.nac++; Li = Lf;
				if(fixinvT == 0 && samp < burnin) invT *= 1.05;
			}
			else{
				if(fixinvT == 0 && samp < burnin) invT *= 0.9;
			}
		}
 
		// Each PMCMC step consists of making a change to a parameter 
		// This change is probablisitically either accepted or rejected based
		// on whether the observations agree better with new value or the old one.			
		for(p = 0; p < param.size(); p++){
			if(param[p].min != param[p].max){
				valst = model.paramval[p];
				
				if(core == 0) model.paramval[p] += normal(0,param[p].jump); 
				
				MPI_Bcast(&model.paramval[p],1,MPI_DOUBLE,0,MPI_COMM_WORLD);

				if(param[p].betachange == 1) model.betaspline(data.period);
				if(param[p].suschange == 1) model.setsus(data);

				model.setsus(data);
				 
				if(model.paramval[p] < param[p].min || model.paramval[p] > param[p].max) al = 0;
				else{
					model.parami = model.paramval; model.paramp = model.paramval;
					if(model.settransprob() == 0) al = 0;
					else{		
						Lf = sample(data,model,poptree,core,ncore,npart,data.period,xp);
						if(core == 0){ 
							al = exp(invT*(Lf-Li)); //cout << Lf << " " << Li << " " << al << " al" << endl;
						}
					}
				}
				if(core == 0){ if(ran() < al) accept = 1; else accept = 0;}
				
				MPI_Bcast(&accept,1,MPI_UNSIGNED,0,MPI_COMM_WORLD);
					
				param[p].ntr++;
				if(accept == 1){
					if(core == 0){
						param[p].nac++;
						xi = xp;
						Li = Lf;
						if(samp < burnin) param[p].jump *= 1.1;
					}
				}
				else{
					model.paramval[p] = valst;
					if(param[p].betachange == 1) model.betaspline(data.period);
					if(param[p].suschange == 1) model.setsus(data);
				
					if(core == 0){
						if(samp < burnin) param[p].jump *= 0.95;
					}
				}
			}
		}
	
		if(core == 0) opsamp.push_back(outputsamp(invT,samp,Li,data,model,poptree,model.paramval,xi));	
	}
	
	if(core == 0) outputresults(data,model,opsamp);
}

/// This samples from the model using particles and returns an overall measure of how well the observations agreed with it 
static double sample(DATA &data, MODEL &model, POPTREE &poptree, unsigned int core, unsigned int ncore, unsigned int npart, unsigned int period, vector < vector <FEV> > &xp)
{
	unsigned int p, t, backpart[ncore*npart*period];
	double Liav;

	for(p = 0; p < npart; p++) part[p]->partinit(p);
	
  Liav = 0;
	for(t = 0; t < period; t++){                                      // We step through one measurement at a time
		for(p = 0; p < npart; p++){
			timers.timesim -= clock();
	
			part[p]->gillespie(t,t+1,0);                                  // Simulates the particle

			part[p]->Li = Lobs(data,model,poptree,t,part[p]->fev,1);      // Measures how well it agrees with the observations
			timers.timesim += clock();
		}
		
		timers.timeboot -= clock();

		// Culls or copies particles based on how well they represent observations
		Liav += bootstrap(core,ncore,npart,(data.fediv*(t+1))/period,t,backpart); 
		
		timers.timeboot += clock();
	}

	geneventsample(data,core,ncore,npart,period,backpart,xp);
	
	return Liav;
}

/// Constructs the event sequence sample by gathering all the pieces from different particles (from the bootstrap function)
void geneventsample(DATA &data, unsigned int core, unsigned int ncore, unsigned int npart, unsigned int period, unsigned int *backpart, vector < vector <FEV> > &xp)	
{
	unsigned int p=UNSET, co, d, fedivmin, fedivmax, nparttot = npart*ncore;
	int t, siz;
	
	if(checkon == 1){
		for(t = 0; t < int(period); t++){
			cout << t << ": "; for(p = 0; p < nparttot; p++) cout << backpart[t*nparttot+p] << ",";
			cout << " Backpart" << endl;  
		}
	}

	xp.resize(data.fediv);
	for(t = period-1; t >= 0; t--){
		fedivmin = (data.fediv*t)/period;	fedivmax = (data.fediv*(t+1))/period;
		if(t == int(period-1)) p = backpart[nparttot*t];  // Picks the first final particle
		else p = backpart[nparttot*t + p];
		
		co = p/npart;
		if(core == 0){  
			if(co == 0){
				for(d = fedivmin; d < fedivmax; d++) xp[d] = part[p%npart]->fev[d];
			}
			else{
				MPI_Status status;
				MPI_Recv(packbuffer(),MAX_NUMBERS,MPI_DOUBLE,co,0,MPI_COMM_WORLD,&status);
				MPI_Get_count(&status, MPI_DOUBLE, &siz); if(siz >= int(MAX_NUMBERS)) emsg("Buffer not big enough");
		
				packinit();
				unpack(xp,fedivmin,fedivmax);
				if(packsize() != siz) emsg("PMCMC: EC10");
			}
		}
		else{
			if(co == core){
				packinit();
				pack(part[p%npart]->fev,fedivmin,fedivmax);
				MPI_Send(packbuffer(),packsize(),MPI_DOUBLE,0,0,MPI_COMM_WORLD);	
			}
		}
	}
}

/// This function does the usual bootstrap step of randomly selecting particles based on their observation probability 
static double bootstrap(unsigned int core, unsigned int ncore, unsigned int npart, unsigned int fedivmin, unsigned int t, unsigned int *backpart)
{
	unsigned int p, pp, pmin, pmax, cor, j, jmax, k, kmax, nparttot = npart*ncore, partstep;
	unsigned int nreqs = 0, nsendbuf = 0, rec,	npreclist, ncorlist;
	double res, *buf = packbuffer();
	
	vector <unsigned int> corlist;
	vector <unsigned int> preclist;
	vector <vector <unsigned int> > preclistli;
	
	backpart += t*nparttot;	
	
	if(core == 0){
		double Limax, av, z, Litot[nparttot], sum, sumst[nparttot];
		vector <unsigned int> extra;
		
		for(p = 0; p < npart; p++) Litot[p] = part[p]->Li;                    // Gathers together the likelihoods from different cores 
	
		for(cor = 1; cor < ncore; cor++){
			MPI_Recv(Litot+npart*cor,npart,MPI_DOUBLE,cor,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
		}
		
		Limax = -large;	for(p = 0; p < nparttot; p++){ if(Litot[p] > Limax) Limax = Litot[p];}
		
		av = 0;	for(p = 0; p < nparttot; p++) av += exp(Litot[p] - Limax); av /= nparttot;

		res = Limax + log(av);
		
		sum = 0;
		for(p = 0; p < nparttot; p++){
			sum += exp(Litot[p] - Limax);
			sumst[p] = sum;
		}

		for(p = 0; p < nparttot; p++) backpart[p] = UNSET;
			
		partstep = nparttot/32; if(partstep < 2) partstep = 2;
		
		for(p = 0; p < nparttot; p++){                                // Finds "backpart", the particle which is being copied from
			z = sum*ran();
	
			pp = 0;
			while(pp < nparttot && z > sumst[pp]) pp += partstep;
	    if(pp > 0) pp -= partstep;
			
			while(pp < nparttot && z > sumst[pp]) pp++;
			if(pp == nparttot) emsg("PMCMC: EC1");	
			if(pp > 0){ if(z < sumst[pp-1]) emsg("PMCMC: EC1a");}
			
			if(backpart[pp] == UNSET) backpart[pp] = pp;
			else extra.push_back(pp);
		}
		
		j = 0; jmax = extra.size();
		while(j < jmax){                                             // Tries to copy to same core (to increase speed)
			pp = extra[j];
			p = pp - pp%npart; pmax = p+npart; while(p < pmax && backpart[p] != UNSET) p++;

			if(p < pmax){	backpart[p] = pp; jmax--; extra[j] = extra[jmax]; extra.pop_back();}
			else j++;
		}
		
		for(p = 0; p < nparttot; p++){  
			if(backpart[p] == UNSET){
				pp = extra[extra.size()-1]; extra.pop_back();
				backpart[p] = pp;
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
	
	MPI_Bcast(backpart,nparttot,MPI_UNSIGNED,0,MPI_COMM_WORLD);     // Broadcasts backpart to all cores
	
	pmin = core*npart;                                         // Sets up information to be recieved
	npreclist = 0;
	for(p = pmin; p < pmin+npart; p++){
		pp = backpart[p];
		if(p != pp){                   
			cor = pp/npart;
			if(cor != core){
				j = 0; while(j < npreclist && preclist[j] != pp) j++;
				if(j < npreclist) preclistli[j].push_back(p%npart);
				else{
					preclist.push_back(pp); 
					preclistli.push_back(vector <unsigned int> ()); 
					preclistli[npreclist].push_back(p%npart);
					
					MPI_Irecv(&recibuffer[npreclist][0],BUFMAX,MPI_DOUBLE,cor,pp,MPI_COMM_WORLD,&reqs[nreqs]); nreqs++;
					npreclist++; if(npreclist > npart) emsg("PMCMC: EC11"); 
				}
			}
		}
	}

	for(p = pmin; p < pmin+npart; p++){                          // Initiates information to be sent
		if(p == backpart[p]){ 
			corlist.clear(); ncorlist = 0;
			for(pp = 0; pp < nparttot; pp++){
				if(p == backpart[pp]){
					cor = pp/npart;
					if(cor == core){
						if(pp != p) part[pp%npart]->copy(*part[p%npart],fedivmin);
					}
					else{
						j = 0; while(j < ncorlist && corlist[j] != cor) j++;
						if(j == ncorlist){ corlist.push_back(cor); ncorlist++;}
					}
				}
			}
						
			if(ncorlist > 0){
				part[p%npart]->partpack(fedivmin);
				kmax = packsize(); 
				sendbuffer[nsendbuf][0] = kmax; if(kmax >= BUFMAX-1) emsg("PMCMC: EC20");
				for(k = 0; k < kmax; k++) sendbuffer[nsendbuf][k+1] = buf[k];
				
				for(j = 0; j < ncorlist; j++){
					MPI_Isend(&sendbuffer[nsendbuf][0],kmax+1,MPI_DOUBLE,corlist[j],p,MPI_COMM_WORLD,&reqs[nreqs]);
					nreqs++;		
				}				
				nsendbuf++; if(nsendbuf > npart) emsg("PMCMC: EC21");
			}
		}
	}
	
	if(nreqs > 0){
		if(MPI_Waitall(nreqs,reqs,stats) != MPI_SUCCESS) emsg("PMCMC: EC22");
	}
		
	for(rec = 0; rec < npreclist; rec++){                           	// Unpacks the recieved information
		kmax = recibuffer[rec][0]; for(k = 0; k < kmax; k++) buf[k] = recibuffer[rec][k+1];		
		p = preclistli[rec][0];
		part[p]->partunpack(fedivmin);
		
		for(j = 1; j < preclistli[rec].size(); j++){
			pp = preclistli[rec][j];
			part[pp]->copy(*part[p],fedivmin);
		}
	}
	
	MPI_Barrier(MPI_COMM_WORLD); 
		
	return res;
}

