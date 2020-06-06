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

MPI_Request reqs[partmax];               // These are information used Isend and Irecv
MPI_Status stats[partmax];

//const short SENDRECMAX = 10;        
//const long BUFMAX = 10000000;
double sendbuffer[SENDRECMAX][BUFMAX];
double recibuffer[SENDRECMAX][BUFMAX];
	
static double sample(short core, short ncore, short npart, double tmax, double invT,	vector < vector <FEV> > &xinew);
static double bootstrap(short core, short ncore, short npart, short fedivmin, short w, short *backpart);
static double bootstrap2(short core, short ncore, short npart, short fedivmin, short w, short *back);
void geneventsample(short core, short ncore, short npart, short nweek, short *backpart, vector < vector <FEV> > &xp);
void geneventsample2(short core, short ncore, short npart, short nweek, short *back, vector < vector <FEV> > &xp);

PART* part[partmax];                       // Pointers to each of the particles 

void PMCMC(DATA &data, MODEL &model, POPTREE &poptree, long nsamp, short core, short ncore, short npart)
{
	const short fixinvT = 1;           // Detemines if the inverse temperature invT is fixed or is dynamically tuned 
	long p, samp, burnin = nsamp/3, accept;
	double Li, Lf, valst, al;
	double invT = 1;                  // The inverse temperature (used to relax the observation model)
	double varfac = 1;                 // Used to change the variances in the observation model
	vector < vector <FEV> > xi, xp;    // Stores the current and proposed event sequences
	vector <SAMPLE> opsamp;            // Stores sample for generating satatistics later
	vector <PARAM> &param(model.param);// Copies parameters from the model

	srand(core);
	
	for(p = 0; p < npart; p++){ part[p] = new PART(data,model,poptree);}

	if(core == 0) outputinit(model,ncore*npart);

	Li = -large; 
	for(samp = 0; samp < nsamp; samp++){
		if(core == 0 && samp%1 == 0) cout << "Sample: " << samp << " " << invT << endl;
 
		Lf = sample(core,ncore,npart,data.tmax,varfac,xp);

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
		//for(p = 0; p < long(param.size()); p++){
		for(p = 0; p < 8; p++){
			if(param[p].min != param[p].max){
				valst = param[p].val;
				
				if(core == 0) param[p].val += normal(0,param[p].jump); 
				
				MPI_Bcast(&param[p].val,1,MPI_DOUBLE,0,MPI_COMM_WORLD);

				if(param[p].betachange == 1) model.betaspline(data.tmax);
				if(param[p].suschange == 1) poptree.setsus(model);
				if(param[p].infchange == 1) poptree.setinf(model);
				
				if(param[p].val < param[p].min || param[p].val > param[p].max) al = 0;
				else{
					Lf = sample(core,ncore,npart,data.tmax,varfac,xp);
					if(core == 0){ 
						al = exp(invT*(Lf-Li)); //cout << Lf << " " << Li << " " << al << " al" << endl;
					}
				}
				if(core == 0){ if(ran() < al) accept = 1; else accept = 0;}
				
				MPI_Bcast(&accept,1,MPI_LONG,0,MPI_COMM_WORLD);
					
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
	
	if(core == 0) outputresults(data,model,opsamp,0,ncore*npart);
}

/// This samples from the model using particles and returns an overall measure of how well the observations agreed with it 
static double sample(short core, short ncore, short npart, double tmax, double varfac, vector < vector <FEV> > &xp)
{
	short p, w, nweek = tmax/timestep, back[ncore*nweek], backpart[ncore*npart*nweek];
	double Liav;
	
	for(p = 0; p < npart; p++) part[p]->partinit(p);

  Liav = 0;
	for(w = 0; w < nweek; w++){                                // We step through one measurement at a time
		for(p = 0; p < npart; p++){
			timers.timesim -= clock();
	
			part[p]->gillespie(w*timestep,(w+1)*timestep, 0 /* Inference */); // Simulates the particle
	
			part[p]->Lobs(w,varfac);      // Measures how well it agrees with the observations (weekly number of cases)
			timers.timesim += clock();
		}
		
		timers.timeboot -= clock();
		
		// Culls or copies particles based on how well they represent observations
		Liav += bootstrap(core,ncore,npart,(fediv*(w+1))/nweek,w,backpart);
		//Liav += bootstrap2(core,ncore,npart,(fediv*(w+1))/nweek,w,back);
		timers.timeboot += clock();
	}

	geneventsample(core,ncore, npart,nweek,backpart,xp);
	//geneventsample2(core,ncore, npart,nweek,back,xp);
	
	return Liav;
}

/// Constructs the event sequence sample by gathering all the pieces from different particles (from the bootstrap function)
void geneventsample(short core, short ncore, short npart, short nweek, short *backpart, vector < vector <FEV> > &xp)	
{
	short p, w, co, d, fedivmin, fedivmax, nparttot = npart*ncore;
	int siz;
	
	if(1 == 0){
		for(w = 0; w < nweek; w++){
			cout << w << ": "; for(p = 0; p < nparttot; p++) cout << backpart[w*nparttot+p] << ",";
			cout << " Backpart" << endl;  
		}
	}

	xp.resize(fediv);
	for(w = nweek-1; w >= 0; w--){
		fedivmin = (fediv*w)/nweek;	fedivmax = (fediv*(w+1))/nweek;
		if(w == nweek-1) p = backpart[nparttot*w];  // Picks the first final particle
		else p = backpart[nparttot*w + p];
		
		co = p/npart;
		if(core == 0){  
			if(co == 0){
				for(d = fedivmin; d < fedivmax; d++) xp[d] = part[p%npart]->fev[d];
			}
			else{
				//cout << co << "recieve fro" << endl;
				MPI_Status status;
				MPI_Recv(packbuffer(),MAX_NUMBERS,MPI_DOUBLE,co,0,MPI_COMM_WORLD,&status);
				MPI_Get_count(&status, MPI_DOUBLE, &siz); if(siz >= MAX_NUMBERS) emsg("Buffer not big enough");
		
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

 // Constructs the event sequence sample by gathering all the pieces from different particles (from the bootstrap2 function)
void geneventsample2(short core, short ncore, short npart, short nweek, short *back, vector < vector <FEV> > &xp)	
{
	short p, w, co, d, fedivmin, fedivmax;
	int siz;

	xp.resize(fediv);
	for(w = nweek-1; w >= 0; w--){
		fedivmin = (fediv*w)/nweek;	fedivmax = (fediv*(w+1))/nweek;
		if(w == nweek-1) p = back[ncore*w];  // Picks the first final particle
		else p = back[ncore*w + p/npart];
		
		co = p/npart;
		if(core == 0){  
			if(co == 0){
				for(d = fedivmin; d < fedivmax; d++) xp[d] = part[p%npart]->fev[d];
			}
			else{
				MPI_Status status;
				MPI_Recv(packbuffer(),MAX_NUMBERS,MPI_DOUBLE,co,0,MPI_COMM_WORLD,&status);
				MPI_Get_count(&status, MPI_DOUBLE, &siz); if(siz >= MAX_NUMBERS) emsg("Buffer not big enough");
		
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
static double bootstrap(short core, short ncore, short npart, short fedivmin, short w, short *backpart)
{
	long p, pp, pmin, pmax, cor, j, jmax, k, kmax, nparttot = npart*ncore, partstep;
	double res; 
	int flag, siz;
	short nreqs = 0, nsendbuf = 0, rec,	npreclist, ncorlist;
	double *buf = packbuffer();
	
	vector <short> corlist;
	
	vector <short> preclist;
	vector <vector <short> > preclistli;
	
	backpart += w*nparttot;	
	
	if(core == 0){
		double Limax, av, z, Litot[nparttot], sum, sumst[nparttot];
		vector <short> extra;
		
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

		for(p = 0; p < nparttot; p++) backpart[p] = -1;
			
		partstep = nparttot/32; if(partstep < 2) partstep = 2;
		
		for(p = 0; p < nparttot; p++){                                            // Finds "backpart", the particle which is being copied from
			z = sum*ran();
	
			pp = 0;
			while(pp < nparttot && z > sumst[pp]) pp += partstep;
	    if(pp > 0) pp -= partstep;
			
			while(pp < nparttot && z > sumst[pp]) pp++;
			if(pp == nparttot) emsg("PMCMC: EC1");	
			if(pp > 0){ if(z < sumst[pp-1]) emsg("PMCMC: EC1a");}
			
			if(backpart[pp] == -1) backpart[pp] = pp;
			else extra.push_back(pp);
		}
		
		j = 0; jmax = extra.size();
		while(j < jmax){                                                            // Tries to copy to same core (to increase speed)
			pp = extra[j];
			p = pp - pp%npart; pmax = p+npart; while(p < pmax && backpart[p] != -1) p++;

			if(p < pmax){	backpart[p] = pp; jmax--; extra[j] = extra[jmax]; extra.pop_back();}
			else j++;
		}
		
		for(p = 0; p < nparttot; p++){  
			if(backpart[p] == -1){
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
	
	MPI_Bcast(backpart,nparttot,MPI_SHORT,0,MPI_COMM_WORLD);                                  // Broadcasts back part to all cores
	
	pmin = core*npart;                                                                        // Sets up information to be recieved
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
					preclistli.push_back(vector <short> ()); 
					preclistli[npreclist].push_back(p%npart);
					
					MPI_Irecv(recibuffer[npreclist],BUFMAX,MPI_DOUBLE,cor,pp,MPI_COMM_WORLD,&reqs[nreqs]); nreqs++;
					npreclist++; if(npreclist == SENDRECMAX) emsg("PMCMC: EC11"); 
				}
			}
		}
	}

	for(p = pmin; p < pmin+npart; p++){                                                    //Initiates information to be sent
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
					MPI_Isend(sendbuffer[nsendbuf],kmax+1,MPI_DOUBLE,corlist[j],p,MPI_COMM_WORLD,&reqs[nreqs]);
					nreqs++;		
				}				
				nsendbuf++; if(nsendbuf == SENDRECMAX) emsg("PMCMC: EC21");
			}
		}
	}
	
	if(nreqs > 0){
		if(MPI_Waitall(nreqs,reqs,stats) != MPI_SUCCESS) emsg("PMCMC: EC22");
	}
		
	for(rec = 0; rec < npreclist; rec++){                                              	// Unpacks the recieved information
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

/// This bootstrap proceeds in two stages: First a particle is selected on each core based on observed likeluhood. 
/// Then each core selects a code based on agregated likelihood
/// Finally the selected particle is copied across all particles within a core. 
/// This approach can reduces MPI time because each core only either sends or recieves.
static double bootstrap2(short core, short ncore, short npart, short fedivmin, short w, short *back)
{
	long p, psel, co, co2, j, k, nparttot = npart*ncore;
	double res; 
	int flag, siz;
	back += w*ncore;	
	
	if(core == 0){
		short sel[ncore];
		double Limax, av, z, Litot[nparttot], sumind, sumindst[npart], sum, sumst[ncore];
		vector <long> extra;
		
		for(p = 0; p < npart; p++) Litot[p] = part[p]->Li;                          // Gathers together the likelihoods from all cores
	
		for(co = 1; co < ncore; co++){
			MPI_Recv(Litot+npart*co,npart,MPI_DOUBLE,co,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
		}
		
		Limax = -large;	for(p = 0; p < nparttot; p++){ if(Litot[p] > Limax) Limax = Litot[p];}
		
		av = 0;	for(p = 0; p < nparttot; p++) av += exp(Litot[p] - Limax); av /= nparttot;

		res = Limax + log(av);
		
		sum = 0;
		for(co = 0; co < ncore; co++){                                       // Selects particles on each core
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
			
		for(co = 0; co < ncore; co++){                                      // Finds "back" for a core
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
	
	MPI_Bcast(back,ncore,MPI_SHORT,0,MPI_COMM_WORLD);
	
	co = back[core]/npart;
	if(co == core){          // Sends information
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
					if(flag == 0){
						part[psel]->partpack(fedivmin);
						flag = 1;
					}
					MPI_Send(packbuffer(),packsize(),MPI_DOUBLE,co,0,MPI_COMM_WORLD);	
				}
			}
		}
	}
	else{                  // Recieves information
		MPI_Status status;
		MPI_Recv(packbuffer(),MAX_NUMBERS,MPI_DOUBLE,co,0,MPI_COMM_WORLD,&status);
		MPI_Get_count(&status, MPI_DOUBLE, &siz); if(siz >= MAX_NUMBERS) emsg("Buffer not big enough");
		part[0]->partunpack(fedivmin); 
		if(siz != packsize()) emsg("PMCMC: EV12");
		for(p = 1; p < npart; p++) part[p]->copy(*part[0],fedivmin);
	}
	
	MPI_Barrier(MPI_COMM_WORLD);  
	
	return res;
}

