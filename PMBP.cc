// Performs particle MBPs

#include <math.h>
#include <iostream>

using namespace std;

#include "utils.hh"
#include "timers.hh"
#include "PART.hh"
#include "data.hh"
#include "output.hh"

static double mbpsample(short core, short ncore, short npart, double tmax, vector < vector <FEV> > &xi, vector < vector <FEV> > &xp);

long timegen = 0, timembp = 0;
	
MBPPART* mbppart[partmax]; 

void PMBP(DATA &data, MODEL &model, POPTREE &poptree, long nsamp, short core, short ncore, short npart)
{
	long p, pp, d, j, nparam, st, w, burnin = nsamp/3, samp, accept;
	double valst, al;
	vector <SAMPLE> opsamp;            // Stores sample for generating satatistics later

	vector <PARAM> &param(model.param);
		
	vector < vector <FEV> > xi, xp;
	
	npart = 1;
	for(p = 0; p < npart; p++){ mbppart[p] = new PART(data,model,poptree);}
	
	mbppart[0]->partinit(0);
	mbppart[0]->gillespie(0,data.tmax, 1 /* simulating */);

	xi = mbppart[0]->fev; // start

	for(st = 0; st < nsettime; st++){
		model.betai[st] = model.beta[st];
		model.betap[st] = model.beta[st];
	}
	
	if(core == 0) outputinit(model,ncore*npart);

	nparam = model.param.size();
	model.parami.resize(nparam); model.paramp.resize(nparam);
	 
	for(samp = 0; samp < nsamp; samp++){
		if(core == 0 && samp%1 == 0) cout << " Sample: " << samp << endl;
	
		for(p = 0; p < -8; p++){
		//for(p = 0; p < long(param.size()); p++){
			if(param[p].min != param[p].max){
				for(pp = 0; pp < nparam; pp++) model.parami[pp] = model.param[pp].val;
				for(st = 0; st < nsettime; st++) model.betai[st] = model.beta[st];
	
				valst = param[p].val;
				
				if(core == 0) param[p].val += normal(0,param[p].jump); 
				
				MPI_Bcast(&param[p].val,1,MPI_DOUBLE,0,MPI_COMM_WORLD);

				if(param[p].betachange == 1) model.betaspline(data.tmax);
				if(param[p].suschange == 1) poptree.setsus(model);
				if(param[p].infchange == 1) poptree.setinf(model);
				
				if(param[p].val < param[p].min || param[p].val > param[p].max) al = 0;
				else{
					for(st = 0; st < nsettime; st++) model.betap[st] = model.beta[st];
					for(pp = 0; pp < nparam; pp++) model.paramp[pp] = model.param[pp].val;
			
					al = mbpsample(core,ncore,npart,data.tmax,xi,xp);
				}
				if(core == 0){ if(ran() < al) accept = 1; else accept = 0;}
				
				MPI_Bcast(&accept,1,MPI_LONG,0,MPI_COMM_WORLD);
					
				param[p].ntr++;
				if(accept == 1){
					if(core == 0){
						param[p].nac++;
						xi = xp;
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
	
		if(core == 0) opsamp.push_back(outputsamp(1,samp,0,data,model,poptree,xi));	
	}
	
	if(core == 0) outputresults(data,model,opsamp,0,ncore*npart);
	
	cout << double(timembp)/CLOCKS_PER_SEC << " Timembp" << endl; // These are diagnostics I'm currently using
	cout << double(timegen)/CLOCKS_PER_SEC << " Timegen" << endl;
}

/// This samples from the model using particles and returns an overall measure of how well the observations agreed with it 
static double mbpsample(short core, short ncore, short npart, double tmax, vector < vector <FEV> > &xi, vector < vector <FEV> > &xp)
{
	short p, w, nweek = tmax/timestep, back[ncore*nweek];
	double Litot, Lftot, al;
			
	mbppart[0]->fev = xi;
			 
	Litot = 0;
	for(w = 0; w < nweek; w++){
		mbppart[0]->Lobs(w,1);
		Litot += mbppart[0]->Li; 
	}

	for(p = 0; p < npart; p++) mbppart[p]->mbpinit(p,xi);

	timembp -= clock();
	mbppart[0]->mbp(0,tmax,xi,timegen);
	
	Lftot = 0;
	for(w = 0; w < nweek; w++){
		mbppart[0]->Lobs(w,1);
		Lftot += mbppart[0]->Li; 
	}
	
	al = exp(Lftot-Litot); //cout << al << " " << Litot << " " << Lftot << " al" << endl;
	
	xp = mbppart[0]->fev;
	timembp += clock();
	
	return al;
}
