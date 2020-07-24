// Performs MBPs on multiple chains spanning from the posterior to the prior

#include <fstream>
#include <iostream>
#include <sstream>
#include "stdlib.h"
#include "math.h"

using namespace std;

#include "model.hh"
#include "MBPCHAIN.hh"
#include "output.hh"
#include "pack.hh"
#include "utils.hh"
#include "consts.hh"
#include "timers.hh"
#include "MBP.hh"

static void MBPoutput(DATA &data, MODEL &model, POPTREE &poptree, vector <SAMPLE> &opsamp, unsigned int core, unsigned int ncore, unsigned int nchain);
static void MBPdiagnostic(DATA &data, MODEL &model, unsigned int core, unsigned int ncore, unsigned int nchain);
static void swap(MODEL &model, unsigned int core, unsigned int ncore, unsigned int nchain);

static MBPCHAIN* mbpchain[chainmax]; 
static vector <vector <double> > Listore;
static vector <double> invTstore;
	
/// Performs the multi-temperature MBP algorithm	
void MBP(DATA &data, MODEL &model, POPTREE &poptree, unsigned int nsamp, unsigned int core, unsigned int ncore, unsigned int nchain, enum proposalsmethod propmethod)
{
	double invTmax = 1, invTmin = 0.1;

	unsigned int p, pp, th, nchaintot = ncore*nchain, co;
	unsigned int samp, burnin = nsamp/4, quench = nsamp/5;
	long time, timeprop=0, ntimeprop=0;
	double invT, timeloop=0.1, K;
	long timeproptot[ncore], ntimeproptot[ncore], timeproptotsum, ntimeproptotsum;
	vector <SAMPLE> opsamp;

	K = nchaintot-1; while(pow((K-(nchaintot-1))/K,5) < invTmin) K += 0.1;

	for(p = 0; p < nchain; p++){
		mbpchain[p] = new MBPCHAIN(data,model,poptree);
		
		pp = core*nchain+p;
		pp = core*nchain+p;
		if(nchaintot == 1 || duplicate == 1) invT = invTmax;
		else invT = pow((K-pp)/K,5);

		mbpchain[p]->init(data,model,poptree,invT,pp);
	}

	if(core == 0){
		Listore.resize(nchaintot); invTstore.resize(nchaintot);
		for(pp = 0; pp < nchaintot; pp++){
			if(nchaintot == 1 || duplicate == 1) invT = invTmax;
			else invT = pow((K-pp)/K,5);
			invTstore[pp] = invT;
		}
	}

	if(checkon == 1){
		for(th = 0; th < mbpchain[0]->paramval.size(); th++){
			cout << th << " " << model.param[th].name << " " << mbpchain[0]->paramval[th] <<  "param" << endl;		
		}
	}
	
	if(core == 0){ outputinit(data,model); outputLiinit(data,ncore*nchain);}

	for(samp = 0; samp < nsamp; samp++){	
		if(core == 0 && samp%1 == 0) cout << " Sample: " << samp << " / " << nsamp << endl; 
	
		for(p = 0; p < nchain; p++){
			if(samp < quench) mbpchain[p]->invT = (double(samp)/quench)*mbpchain[p]->invTtrue;
			else mbpchain[p]->invT = mbpchain[p]->invTtrue;
			//mbpchain[p]->invT = mbpchain[p]->invTtrue;
		}
			
		time = clock();

		for(p = 0; p < nchain; p++) mbpchain[p]->standard_prop(samp,burnin);
		
		propmethod = proposalsmethod::fixednum;
		
		timeprop -= clock();
		switch(propmethod){
		case proposalsmethod::allchainsallparams:
			for(p = 0; p < nchain; p++){
				for(th = 0; th < model.param.size(); th++){
					if(model.param[th].min != model.param[th].max){ mbpchain[p]->proposal(th,samp,burnin); ntimeprop++;}
				}
			}
			break;
			
		case proposalsmethod::fixednum:
			short lo; 
			for(lo = 0; lo < 10; lo++){
				p = int(ran()*nchain);
				th = (unsigned int)(ran()*model.param.size());
				if(model.param[th].min != model.param[th].max){ mbpchain[p]->proposal(th,samp,burnin); ntimeprop++;}
			}
			break;
			
		case proposalsmethod::fixedtime:
			do{                         // Does proposals for timeloop seconds (on average 20 proposals)
				p = int(ran()*nchain);
				//unsigned int thbeta = 3;
				//if(ran() < 0.7) th = (unsigned int)(ran()*thbeta);
				//else th = thbeta + (unsigned int)(ran()*(model.param.size()-thbeta));
				th = (unsigned int)(ran()*model.param.size());
				if(model.param[th].min != model.param[th].max){ mbpchain[p]->proposal(th,samp,burnin); ntimeprop++;}
			}while(double(clock()-time)/CLOCKS_PER_SEC < timeloop);
			break;
		}
		
		timeprop += clock();

		if(samp%10 == 0){ for(p = 0; p < nchain; p++) mbpchain[p]->setQmapi(0);}  // Recalcualtes Qmapi (numerical)
	
		timers.timewait -= clock();
		
		MPI_Barrier(MPI_COMM_WORLD); 
	
		timers.timewait += clock();

		// This part dynamically adjusts timeloop so that approximately 20 proposals are performed per swap
		MPI_Gather(&timeprop,1,MPI_LONG,timeproptot,1,MPI_LONG,0,MPI_COMM_WORLD); 
		MPI_Gather(&ntimeprop,1,MPI_LONG,ntimeproptot,1,MPI_LONG,0,MPI_COMM_WORLD);
		if(core == 0){
			timeproptotsum = 0; ntimeproptotsum = 0; 
			for(co = 0; co < ncore; co++){ timeproptotsum += timeproptot[co]; ntimeproptotsum += ntimeproptot[co];}
            
			// Update the time to run only if some proposals have run (otherwise it runs forever)
			if (ntimeproptotsum > 0 && timeproptotsum > 0) {
				timeloop = 20*double(timeproptotsum)/(ntimeproptotsum*CLOCKS_PER_SEC);
			}
		}
		MPI_Bcast(&timeloop,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
		
		if(duplicate == 0) swap(model,core,ncore,nchain);
		
		if(samp%1 == 0) MBPoutput(data,model,poptree,opsamp,core,ncore,nchain);
		
		if(samp == nsamp-1 || (samp != 0 && samp%100 == 0)){
		  //if(samp == nsamp-1){
			if(core == 0) outputresults(data,model,opsamp);
			MBPdiagnostic(data,model,core,ncore,nchain);
		}
	}
cout << ran() << "end ran\n";
	for(p = 0; p < nchain; p++){
		delete mbpchain[p];
	}
}

/// Stochastically swaps chains with similar inverse temperatures 
static void swap(MODEL &model, unsigned int core, unsigned int ncore, unsigned int nchain)
{
	unsigned int p, p1, p2, th, tempi;
	unsigned int nparam = model.param.size(), nchaintot = nchain*ncore, nchainparam = nchain*nparam, nparamtot = nchainparam*ncore;
	double temp, al;
	float tempf;
	long templ, loop, loopmax = nchaintot*nchaintot;
	double Li[nchain], Litot[nchaintot];
	double invT[nchain], invTtot[nchaintot];		
	double invTtrue[nchain], invTtruetot[nchaintot];
	unsigned int ch[nchain], chtot[nchaintot];
	long timeprop[nchain], timeproptot[nchaintot];
	unsigned int ntr[nchainparam], ntrtot[nparamtot];
	unsigned int nac[nchainparam], nactot[nparamtot];
	float paramjump[nchainparam], paramjumptot[nparamtot];
	unsigned int ntrxi[nchainparam], ntrxitot[nparamtot];
	unsigned int nacxi[nchainparam], nacxitot[nparamtot];
	float paramjumpxi[nchainparam], paramjumpxitot[nparamtot];
	float numaddrem[nchain], numaddremtot[nchaintot];
	unsigned int ntr_addrem[nchain], ntr_addremtot[nchaintot];
	unsigned int nac_addrem[nchain], nac_addremtot[nchaintot];

	timers.timeswap -= clock();
		
	for(p = 0; p < nchain; p++){ 
		Li[p] = mbpchain[p]->Li; 
		invT[p] = mbpchain[p]->invT; 
		invTtrue[p] = mbpchain[p]->invTtrue; 
		ch[p] = mbpchain[p]->ch;
		timeprop[p] = mbpchain[p]->timeprop;
		for(th = 0; th < nparam; th++){
			ntr[p*nparam+th] = mbpchain[p]->ntr[th];
			nac[p*nparam+th] = mbpchain[p]->nac[th];
			paramjump[p*nparam+th] = mbpchain[p]->paramjump[th];
			
			ntrxi[p*nparam+th] = mbpchain[p]->ntrxi[th];
			nacxi[p*nparam+th] = mbpchain[p]->nacxi[th];
			paramjumpxi[p*nparam+th] = mbpchain[p]->paramjumpxi[th];
		}
		numaddrem[p] = mbpchain[p]->numaddrem;
		ntr_addrem[p] = mbpchain[p]->ntr_addrem;
		nac_addrem[p] = mbpchain[p]->nac_addrem;
	}
		
	MPI_Gather(Li,nchain,MPI_DOUBLE,Litot,nchain,MPI_DOUBLE,0,MPI_COMM_WORLD);
	MPI_Gather(invT,nchain,MPI_DOUBLE,invTtot,nchain,MPI_DOUBLE,0,MPI_COMM_WORLD);
	MPI_Gather(invTtrue,nchain,MPI_DOUBLE,invTtruetot,nchain,MPI_DOUBLE,0,MPI_COMM_WORLD);
	MPI_Gather(ch,nchain,MPI_UNSIGNED,chtot,nchain,MPI_UNSIGNED,0,MPI_COMM_WORLD);
	MPI_Gather(timeprop,nchain,MPI_LONG,timeproptot,nchain,MPI_LONG,0,MPI_COMM_WORLD);
	MPI_Gather(ntr,nchainparam,MPI_UNSIGNED,ntrtot,nchainparam,MPI_UNSIGNED,0,MPI_COMM_WORLD);
	MPI_Gather(nac,nchainparam,MPI_UNSIGNED,nactot,nchainparam,MPI_UNSIGNED,0,MPI_COMM_WORLD);
	MPI_Gather(paramjump,nchainparam,MPI_FLOAT,paramjumptot,nchainparam,MPI_FLOAT,0,MPI_COMM_WORLD);
	MPI_Gather(ntrxi,nchainparam,MPI_UNSIGNED,ntrxitot,nchainparam,MPI_UNSIGNED,0,MPI_COMM_WORLD);
	MPI_Gather(nacxi,nchainparam,MPI_UNSIGNED,nacxitot,nchainparam,MPI_UNSIGNED,0,MPI_COMM_WORLD);
	MPI_Gather(paramjumpxi,nchainparam,MPI_FLOAT,paramjumpxitot,nchainparam,MPI_FLOAT,0,MPI_COMM_WORLD);
	MPI_Gather(numaddrem,nchain,MPI_FLOAT,numaddremtot,nchain,MPI_FLOAT,0,MPI_COMM_WORLD);
	MPI_Gather(ntr_addrem,nchain,MPI_UNSIGNED,ntr_addremtot,nchain,MPI_UNSIGNED,0,MPI_COMM_WORLD);
	MPI_Gather(nac_addrem,nchain,MPI_UNSIGNED,nac_addremtot,nchain,MPI_UNSIGNED,0,MPI_COMM_WORLD);

	if(core == 0){   // Swaps different chains
		for(loop = 0; loop < loopmax; loop++){
			p1 = (unsigned int)(ran()*nchaintot); p2 = (unsigned int)(ran()*nchaintot); 
			if(p1 != p2){
				al = exp((invTtot[p1]-invTtot[p2])*(Litot[p2]-Litot[p1]));
				
				if(ran() < al){
					temp = invTtot[p1]; invTtot[p1] = invTtot[p2]; invTtot[p2] = temp;
					temp = invTtruetot[p1]; invTtruetot[p1] = invTtruetot[p2]; invTtruetot[p2] = temp;
					tempi = chtot[p1]; chtot[p1] = chtot[p2]; chtot[p2] = tempi;
					templ = timeproptot[p1]; timeproptot[p1] = timeproptot[p2]; timeproptot[p2] = templ;
					for(th = 0; th < nparam; th++){
						tempi = ntrtot[p1*nparam+th]; ntrtot[p1*nparam+th] = ntrtot[p2*nparam+th]; ntrtot[p2*nparam+th] = tempi;
						tempi = nactot[p1*nparam+th]; nactot[p1*nparam+th] = nactot[p2*nparam+th]; nactot[p2*nparam+th] = tempi;
						tempf = paramjumptot[p1*nparam+th]; paramjumptot[p1*nparam+th] = paramjumptot[p2*nparam+th];
						paramjumptot[p2*nparam+th] = tempf;
						
						tempi = ntrxitot[p1*nparam+th]; ntrxitot[p1*nparam+th] = ntrxitot[p2*nparam+th]; ntrxitot[p2*nparam+th] = tempi;
						tempi = nacxitot[p1*nparam+th]; nacxitot[p1*nparam+th] = nacxitot[p2*nparam+th]; nacxitot[p2*nparam+th] = tempi;
						tempf = paramjumpxitot[p1*nparam+th]; paramjumpxitot[p1*nparam+th] = paramjumpxitot[p2*nparam+th];
						paramjumpxitot[p2*nparam+th] = tempf;
					}		
					tempf = numaddremtot[p1]; numaddremtot[p1] = numaddremtot[p2]; numaddremtot[p2] = tempf;
					tempi = ntr_addremtot[p1]; ntr_addremtot[p1] = ntr_addremtot[p2]; ntr_addremtot[p2] = tempi;
					tempi = nac_addremtot[p1]; nac_addremtot[p1] = nac_addremtot[p2]; nac_addremtot[p2] = tempi;
				}
			}
		}
	}
	
	MPI_Scatter(invTtot,nchain,MPI_DOUBLE,invT,nchain,MPI_DOUBLE,0,MPI_COMM_WORLD);
	MPI_Scatter(invTtruetot,nchain,MPI_DOUBLE,invTtrue,nchain,MPI_DOUBLE,0,MPI_COMM_WORLD);
	MPI_Scatter(chtot,nchain,MPI_UNSIGNED,ch,nchain,MPI_UNSIGNED,0,MPI_COMM_WORLD);
	MPI_Scatter(timeproptot,nchain,MPI_LONG,timeprop,nchain,MPI_LONG,0,MPI_COMM_WORLD);
	MPI_Scatter(ntrtot,nchainparam,MPI_UNSIGNED,ntr,nchainparam,MPI_UNSIGNED,0,MPI_COMM_WORLD);
	MPI_Scatter(nactot,nchainparam,MPI_UNSIGNED,nac,nchainparam,MPI_UNSIGNED,0,MPI_COMM_WORLD);
	MPI_Scatter(paramjumptot,nchainparam,MPI_FLOAT,paramjump,nchainparam,MPI_FLOAT,0,MPI_COMM_WORLD);
	MPI_Scatter(ntrxitot,nchainparam,MPI_UNSIGNED,ntrxi,nchainparam,MPI_UNSIGNED,0,MPI_COMM_WORLD);
	MPI_Scatter(nacxitot,nchainparam,MPI_UNSIGNED,nacxi,nchainparam,MPI_UNSIGNED,0,MPI_COMM_WORLD);
	MPI_Scatter(paramjumpxitot,nchainparam,MPI_FLOAT,paramjumpxi,nchainparam,MPI_FLOAT,0,MPI_COMM_WORLD);
	MPI_Scatter(numaddremtot,nchain,MPI_FLOAT,numaddrem,nchain,MPI_FLOAT,0,MPI_COMM_WORLD);
	MPI_Scatter(ntr_addremtot,nchain,MPI_UNSIGNED,ntr_addrem,nchain,MPI_UNSIGNED,0,MPI_COMM_WORLD);
	MPI_Scatter(nac_addremtot,nchain,MPI_UNSIGNED,nac_addrem,nchain,MPI_UNSIGNED,0,MPI_COMM_WORLD);

	for(p = 0; p < nchain; p++){ 
		mbpchain[p]->invT = invT[p]; 
		mbpchain[p]->invTtrue = invTtrue[p]; 
		mbpchain[p]->ch = ch[p];
		mbpchain[p]->timeprop = timeprop[p];
		for(th = 0; th < nparam; th++){
			mbpchain[p]->ntr[th] = ntr[p*nparam+th];
			mbpchain[p]->nac[th] = nac[p*nparam+th];
			mbpchain[p]->paramjump[th] = paramjump[p*nparam+th];
			
			mbpchain[p]->ntrxi[th] = ntrxi[p*nparam+th];
			mbpchain[p]->nacxi[th] = nacxi[p*nparam+th];
			mbpchain[p]->paramjumpxi[th] = paramjumpxi[p*nparam+th];
		}
		mbpchain[p]->numaddrem = numaddrem[p];
		mbpchain[p]->ntr_addrem = ntr_addrem[p];
		mbpchain[p]->nac_addrem = nac_addrem[p];
	}
	
	timers.timeswap += clock();
}

/// Calculates the model evidence
double calcME()
{
	unsigned int ch, j, jmin, jmax;
	double max, dinvT, ME, sum;
	
	ME = 0;
	for(ch = 1; ch < Listore.size(); ch++){
		dinvT = invTstore[ch-1]-invTstore[ch];
		jmax = Listore[ch].size(); jmin = jmax/3;
		max = -large; for(j = jmin; j < jmax; j++){ if(Listore[ch][j] > max) max = Listore[ch][j];}
		
		sum = 0; for(j = jmin; j < jmax; j++) sum += exp(dinvT*(Listore[ch][j]-max));
		ME += dinvT*max + log(sum/(jmax-jmin));
	}
	if(ME < -large) ME = -large;
	
	return ME;
}

/// Ouputs a sample from the MBP algorithm
void MBPoutput(DATA &data, MODEL &model, POPTREE &poptree, vector <SAMPLE> &opsamp, unsigned int core, unsigned int ncore, unsigned int nchain)
{
	unsigned int p, ppost, nchaintot = ncore*nchain, samp = long(opsamp.size());
	int siz;
	double Li[nchain], Litot[nchaintot], Liord[nchaintot];
	unsigned int ch[nchain], chtot[nchaintot];

	timers.timeoutput -= clock();
		
	for(p = 0; p < nchain; p++){ Li[p] = mbpchain[p]->Li; ch[p] = mbpchain[p]->ch;}
		
	MPI_Gather(Li,nchain,MPI_DOUBLE,Litot,nchain,MPI_DOUBLE,0,MPI_COMM_WORLD);
	MPI_Gather(ch,nchain,MPI_UNSIGNED,chtot,nchain,MPI_UNSIGNED,0,MPI_COMM_WORLD);

	if(core == 0){
		for(p = 0; p < nchaintot; p++){
			if(chtot[p] == 0) ppost = p;
			Listore[chtot[p]].push_back(Litot[p]);
			Liord[chtot[p]] = Litot[p];
		}
		outputLi(samp,nchaintot,Liord);
	}
	
	MPI_Bcast(&ppost,1,MPI_UNSIGNED,0,MPI_COMM_WORLD);

	if(core == 0){
		double L, Pr;
		vector < vector <EVREF> > trevplot;
		vector < vector <FEV> > indevplot;
		vector <double> paramplot;
		unsigned int ninfplot;

		if(ppost < nchain){
			trevplot = mbpchain[ppost]->trevi; indevplot = mbpchain[ppost]->indevi; 
			paramplot = mbpchain[ppost]->paramval; L = mbpchain[ppost]->Li; Pr = mbpchain[ppost]->Pri;
			ninfplot =  mbpchain[ppost]->xi.size();
		}
		else{
			MPI_Status status;
			MPI_Recv(packbuffer(),MAX_NUMBERS,MPI_DOUBLE,ppost/nchain,0,MPI_COMM_WORLD,&status);
			MPI_Get_count(&status, MPI_DOUBLE, &siz); if(siz >= int(MAX_NUMBERS)) emsg("Buffer not big enough");
		
			packinit();
			unpack(trevplot);
			unpack(indevplot,0,data.popsize);
			unpack(paramplot);
			unpack(L);
			unpack(Pr);
			unpack(ninfplot);
			if(packsize() != siz) emsg("MBP: EC10");
		}

		opsamp.push_back(outputsamp(calcME(),samp,L,Pr,data,model,poptree,paramplot,ninfplot,trevplot,indevplot));
	}
	else{
		if(core == ppost/nchain){
			p = ppost%nchain;
			packinit();
			pack(mbpchain[p]->trevi);
			pack(mbpchain[p]->indevi,0,data.popsize);
			pack(mbpchain[p]->paramval);
			pack(mbpchain[p]->Li);
			pack(mbpchain[p]->Pri);
			pack((unsigned int)mbpchain[p]->xi.size());
			MPI_Send(packbuffer(),packsize(),MPI_DOUBLE,0,0,MPI_COMM_WORLD);
		}
	}
	
	timers.timeoutput += clock();
}

/// Outputs MBP MCMC diagnoistic information
static void MBPdiagnostic(DATA &data, MODEL &model, unsigned int core, unsigned int ncore, unsigned int nchain)
{
	unsigned int p, c, cc, th, j, jmax;
	unsigned int nparam = model.param.size(), nchaintot = ncore*nchain, nchainparam = nchain*nparam, nparamtot = nchainparam*ncore;
	double av, ntrsum;
	double Li[nchain], Litot[nchaintot];
	unsigned int ch[nchain], chtot[nchaintot];
	long timeprop[nchain], timeproptot[nchaintot];
	unsigned int ntr[nchainparam], ntrtot[nparamtot];
	unsigned int nac[nchainparam], nactot[nparamtot];
	float paramjump[nchainparam], paramjumptot[nparamtot];
	unsigned int ntrxi[nchainparam], ntrxitot[nparamtot];
	unsigned int nacxi[nchainparam], nacxitot[nparamtot];
	float paramjumpxi[nchainparam], paramjumpxitot[nparamtot];
	float numaddrem[nchain], numaddremtot[nchaintot];
	unsigned int ntr_addrem[nchain], ntr_addremtot[nchaintot];
	unsigned int nac_addrem[nchain], nac_addremtot[nchaintot];

	for(p = 0; p < nchain; p++){ 
		Li[p] = mbpchain[p]->Li; 
		ch[p] = mbpchain[p]->ch;
		timeprop[p] = mbpchain[p]->timeprop;
		for(th = 0; th < nparam; th++){
			ntr[p*nparam+th] = mbpchain[p]->ntr[th];
			nac[p*nparam+th] = mbpchain[p]->nac[th];
			paramjump[p*nparam+th] = mbpchain[p]->paramjump[th];
			
			ntrxi[p*nparam+th] = mbpchain[p]->ntrxi[th];
			nacxi[p*nparam+th] = mbpchain[p]->nacxi[th];
			paramjumpxi[p*nparam+th] = mbpchain[p]->paramjumpxi[th];
		}
		numaddrem[p] = mbpchain[p]->numaddrem;
		ntr_addrem[p] = mbpchain[p]->ntr_addrem;
		nac_addrem[p] = mbpchain[p]->nac_addrem;
	}
		
	MPI_Gather(Li,nchain,MPI_DOUBLE,Litot,nchain,MPI_DOUBLE,0,MPI_COMM_WORLD);
	MPI_Gather(ch,nchain,MPI_UNSIGNED,chtot,nchain,MPI_UNSIGNED,0,MPI_COMM_WORLD);
	MPI_Gather(timeprop,nchain,MPI_LONG,timeproptot,nchain,MPI_LONG,0,MPI_COMM_WORLD);
	MPI_Gather(ntr,nchainparam,MPI_UNSIGNED,ntrtot,nchainparam,MPI_UNSIGNED,0,MPI_COMM_WORLD);
	MPI_Gather(nac,nchainparam,MPI_UNSIGNED,nactot,nchainparam,MPI_UNSIGNED,0,MPI_COMM_WORLD);
	MPI_Gather(paramjump,nchainparam,MPI_FLOAT,paramjumptot,nchainparam,MPI_FLOAT,0,MPI_COMM_WORLD);
	MPI_Gather(ntrxi,nchainparam,MPI_UNSIGNED,ntrxitot,nchainparam,MPI_UNSIGNED,0,MPI_COMM_WORLD);
	MPI_Gather(nacxi,nchainparam,MPI_UNSIGNED,nacxitot,nchainparam,MPI_UNSIGNED,0,MPI_COMM_WORLD);
	MPI_Gather(paramjumpxi,nchainparam,MPI_FLOAT,paramjumpxitot,nchainparam,MPI_FLOAT,0,MPI_COMM_WORLD);
	MPI_Gather(numaddrem,nchain,MPI_FLOAT,numaddremtot,nchain,MPI_FLOAT,0,MPI_COMM_WORLD);
	MPI_Gather(ntr_addrem,nchain,MPI_UNSIGNED,ntr_addremtot,nchain,MPI_UNSIGNED,0,MPI_COMM_WORLD);
	MPI_Gather(nac_addrem,nchain,MPI_UNSIGNED,nac_addremtot,nchain,MPI_UNSIGNED,0,MPI_COMM_WORLD);

	if(core == 0){
		stringstream ss; ss << data.outputdir << "/MCMCdiagnostic.txt";
		ofstream diag(ss.str().c_str()); 
		ofstream timings(data.outputdir+"/MCMCdiagnostic_timings.txt"); 

		for(c = 0; c < nchaintot; c++){
			cc = 0; while(cc < nchaintot && c != chtot[cc]) cc++;
			if(cc == nchaintot) emsg("MBP: EC9");			
		
			diag << "-------- Chain: " << c << " ----------" << endl << endl;
			diag << "invT: " << invTstore[c] << "  ";
			
			av = 0; jmax = Listore[c].size(); for(j = 0; j < jmax; j++) av += Listore[c][j];
			diag << "Liav: " << av/jmax << "  ";
			
			ntrsum = 0; for(th = 0; th < nparam; th++) ntrsum += ntrtot[cc*nparam+th];
			//diag << "Time: " << int(1000*timeproptot[cc]/(ntrsum*CLOCKS_PER_SEC));
			diag << endl;
			
			diag << "MBP Accept: " << endl;
			for(th = 0; th < nparam; th++){
				diag << model.param[th].name << ": ";
				if(ntrtot[cc*nparam+th] == 0) diag << "--- ";
				else diag << double(nactot[cc*nparam+th])/ntrtot[cc*nparam+th] << " jump:" << paramjumptot[cc*nparam+th] << ", ";
				diag << endl;
			}
			diag << endl;
					
			diag << "Standard Accept: " << endl;
			for(th = 0; th < nparam; th++){
				diag << model.param[th].name << ": ";
				if(ntrxitot[cc*nparam+th] == 0) diag << "--- ";
				else diag << double(nacxitot[cc*nparam+th])/ntrxitot[cc*nparam+th] << " jump:" << paramjumpxitot[cc*nparam+th] << ", ";
				diag << endl;
			}
			diag << endl;
	
			diag << "Add / remove infected: " << endl;
			if(ntr_addremtot[cc] == 0) diag << "--- "; 
			else diag << double(nac_addremtot[cc])/ntr_addremtot[cc] << " num:" << numaddremtot[cc] << ", ";
			diag << endl << endl;
		}
		
		timings << endl << "Timings for different parts of the algorithm:" << endl;
		timings << double(timers.timewait)/CLOCKS_PER_SEC << " MBP waiting time (seconds)" << endl;
		timings << double(timers.timembp)/CLOCKS_PER_SEC << " MBP time (seconds)" << endl;
		timings << double(timers.timembpinit)/CLOCKS_PER_SEC << " MBP init (seconds)" << endl;
		timings << double(timers.timembpQmap)/CLOCKS_PER_SEC << " MBP Qmap (seconds)" << endl;
		timings << double(timers.timembpconRtot)/CLOCKS_PER_SEC << " MBP conRtot (seconds)" << endl;
		timings << double(timers.timembpprop)/CLOCKS_PER_SEC << " MBP prop (seconds)" << endl;
		timings << double(timers.timembptemp)/CLOCKS_PER_SEC << " MBP temp (seconds)" << endl;
		timings << double(timers.timembptemp2)/CLOCKS_PER_SEC << " MBP temp2 (seconds)" << endl;
		timings << double(timers.timembptemp3)/CLOCKS_PER_SEC << " MBP temp3 (seconds)" << endl;
		timings << double(timers.timembptemp4)/CLOCKS_PER_SEC << " MBP temp4 (seconds)" << endl;
		timings << double(timers.timestandard)/CLOCKS_PER_SEC << " Standard (seconds)" << endl;			
		timings << double(timers.timeparam)/CLOCKS_PER_SEC << " Param (seconds)" << endl;			
		timings << double(timers.timebetaphiinit)/CLOCKS_PER_SEC << " Betaphiinit (seconds)" << endl;		
		timings << double(timers.timebetaphi)/CLOCKS_PER_SEC << " Betaphi (seconds)" << endl;	
		timings << double(timers.timecovarinit)/CLOCKS_PER_SEC << " Covarinit (seconds)" << endl;
		timings << double(timers.timecovar)/CLOCKS_PER_SEC << " Covar (seconds)" << endl;
		timings << double(timers.timecompparam)/CLOCKS_PER_SEC << " Compparam (seconds)" << endl;						
		timings << double(timers.timeaddrem)/CLOCKS_PER_SEC << " Add / rem (seconds)" << endl;	
		timings << double(timers.timeswap)/CLOCKS_PER_SEC << " Swap (seconds)" << endl;	
		timings << double(timers.timeoutput)/CLOCKS_PER_SEC << " Output (seconds)" << endl;			
	}
}
