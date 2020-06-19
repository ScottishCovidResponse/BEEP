// Performs MBPs on multiple chains spanning from the posterior to the prior

#include <fstream>
#include <iostream>
#include <sstream>
#include "stdlib.h"
#include "math.h"

using namespace std;

#include "model.hh"
#include "PART.hh"
#include "MBPCHAIN.hh"
#include "output.hh"
#include "pack.hh"
#include "utils.hh"
#include "consts.hh"
#include "timers.hh"

static void MBPoutput(DATA &data, MODEL &model, POPTREE &poptree, vector <SAMPLE> &opsamp, unsigned int core, unsigned int ncore, unsigned int nchain);
static void MBPdiagnostic(DATA &data, MODEL &model, unsigned int core, unsigned int ncore, unsigned int nchain);
static void swap(MODEL &model, unsigned int core, unsigned int ncore, unsigned int nchain);

static MBPCHAIN* mbpchain[chainmax]; 
static vector <vector <double> > Listore;
static vector <double> invTstore;
	
/// Performs the multi-temperature MBP algorithm	
void MBP(DATA &data, MODEL &model, POPTREE &poptree, unsigned int nsamp, unsigned int core, unsigned int ncore, unsigned int nchain)
{
	unsigned int p, pp, th, nchaintot = ncore*nchain, loop, loopmax=1000, co;
	unsigned int samp, burnin = nsamp/4;
	long time, timeprop=0, ntimeprop=0;
	double invT, timeloop=0.1;
	long timeproptot[ncore], ntimeproptot[ncore], timeproptotsum, ntimeproptotsum;
	PART *part;
	vector <SAMPLE> opsamp;

	srand(core);

	part = new PART(data,model,poptree);                                          // Initialises chains
	for(p = 0; p < nchain; p++){
		mbpchain[p] = new MBPCHAIN(data,model,poptree);
		
		loop = 0;
		do{
			do{
				model.priorsamp();                                                        // Randomly samples parameters from the prior	
			}while(model.settransprob() == 0);
			
			part->partinit(0);                                
			part->gillespie(0,data.period,0);                                         // Simulates from the model
	
			pp = core*nchain+p;
			if(nchaintot == 1) invT = 1;
			else invT = pow(double(nchaintot-1-pp)/(nchaintot-1),5);

			mbpchain[p]->init(data,model,poptree,invT,part->fev,part->indev,pp);
			loop++;
		}while(loop < loopmax && mbpchain[p]->ninftot >= INFMAX);                   // Checks not too many infected (based on prior)
		if(loop == loopmax) emsg("Cannot find initial state under INFMAX");
	}

	if(core == 0){
		Listore.resize(nchaintot); invTstore.resize(nchaintot);
		for(pp = 0; pp < nchaintot; pp++){
			if(nchaintot == 1) invT = 1;
			else invT = pow(double(nchaintot-1-pp)/(nchaintot-1),5);
			invTstore[pp] = invT;
		}
	}

	if(core == 0){ outputinit(data,model); outputLiinit(data,ncore*nchain);}
			
	for(samp = 0; samp < nsamp; samp++){	
		if(core == 0 && samp%1 == 0) cout << " Sample: " << samp << " / " << nsamp << endl; 

		time = clock();
		//short lo;
		//for(lo = 0; lo < 10; lo++){
		do{                         // Does proposals for timeloop seconds (on average 10 proposals)
			p = int(ran()*nchain);
			th = (unsigned int)(ran()*model.param.size());
			if(model.param[th].min != model.param[th].max){
				timeprop -= clock();
				mbpchain[p]->proposal(data,model,poptree,th,samp,burnin);
				timeprop += clock();
				ntimeprop++;
			}
		}while(double(clock()-time)/CLOCKS_PER_SEC < timeloop);
		//}
		
		timers.timewait -= clock();
		
		MPI_Barrier(MPI_COMM_WORLD); 
	
		timers.timewait += clock();
		
		// This part dynamically adjusts timeloop so that approximately 10 proposals are performed per swap
		MPI_Gather(&timeprop,1,MPI_LONG,timeproptot,1,MPI_LONG,0,MPI_COMM_WORLD); 
		MPI_Gather(&ntimeprop,1,MPI_LONG,ntimeproptot,1,MPI_LONG,0,MPI_COMM_WORLD);
		if(core == 0){
			timeproptotsum = 0; ntimeproptotsum = 0; 
			for(co = 0; co < ncore; co++){ timeproptotsum += timeproptot[co]; ntimeproptotsum += ntimeproptot[co];}
			timeloop = 10*double(timeproptotsum)/(ntimeproptotsum*CLOCKS_PER_SEC);
		}
		MPI_Bcast(&timeloop,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
		
		swap(model,core,ncore,nchain);
		
		MBPoutput(data,model,poptree,opsamp,core,ncore,nchain);
		
		if(samp != 0 && samp%100 == 0){
			if(core == 0) outputresults(data,model,opsamp);
			MBPdiagnostic(data,model,core,ncore,nchain);
		}
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
	unsigned int ch[nchain], chtot[nchaintot];
	long timeprop[nchain], timeproptot[nchaintot];
	unsigned int ntr[nchainparam], ntrtot[nparamtot];
	unsigned int nac[nchainparam], nactot[nparamtot];
	float paramjump[nchainparam], paramjumptot[nparamtot];
	
	for(p = 0; p < nchain; p++){ 
		Li[p] = mbpchain[p]->Li; 
		invT[p] = mbpchain[p]->invT; 
		ch[p] = mbpchain[p]->ch;
		timeprop[p] = mbpchain[p]->timeprop;
		for(th = 0; th < nparam; th++){
			ntr[p*nparam+th] = mbpchain[p]->ntr[th];
			nac[p*nparam+th] = mbpchain[p]->nac[th];
			paramjump[p*nparam+th] = mbpchain[p]->paramjump[th];
		}
	}
		
	MPI_Gather(Li,nchain,MPI_DOUBLE,Litot,nchain,MPI_DOUBLE,0,MPI_COMM_WORLD);
	MPI_Gather(invT,nchain,MPI_DOUBLE,invTtot,nchain,MPI_DOUBLE,0,MPI_COMM_WORLD);
	MPI_Gather(ch,nchain,MPI_UNSIGNED,chtot,nchain,MPI_UNSIGNED,0,MPI_COMM_WORLD);
	MPI_Gather(timeprop,nchain,MPI_LONG,timeproptot,nchain,MPI_LONG,0,MPI_COMM_WORLD);
	MPI_Gather(ntr,nchainparam,MPI_UNSIGNED,ntrtot,nchainparam,MPI_UNSIGNED,0,MPI_COMM_WORLD);
	MPI_Gather(nac,nchainparam,MPI_UNSIGNED,nactot,nchainparam,MPI_UNSIGNED,0,MPI_COMM_WORLD);
	MPI_Gather(paramjump,nchainparam,MPI_FLOAT,paramjumptot,nchainparam,MPI_FLOAT,0,MPI_COMM_WORLD);

	if(core == 0){   // Swaps different chains
		for(loop = 0; loop < loopmax; loop++){
			p1 = (unsigned int)(ran()*nchaintot); p2 = (unsigned int)(ran()*nchaintot); 
			if(p1 != p2){
				al = exp((invTtot[p1]-invTtot[p2])*(Litot[p2]-Litot[p1]));
				
				if(ran() < al){
					temp = invTtot[p1]; invTtot[p1] = invTtot[p2]; invTtot[p2] = temp;
					tempi = chtot[p1]; chtot[p1] = chtot[p2]; chtot[p2] = tempi;
					templ = timeproptot[p1]; timeproptot[p1] = timeproptot[p2]; timeproptot[p2] = templ;
					for(th = 0; th < nparam; th++){
						tempi = ntrtot[p1*nparam+th]; ntrtot[p1*nparam+th] = ntrtot[p2*nparam+th]; ntrtot[p2*nparam+th] = tempi;
						tempi = nactot[p1*nparam+th]; nactot[p1*nparam+th] = nactot[p2*nparam+th]; nactot[p2*nparam+th] = tempi;
						tempf = paramjumptot[p1*nparam+th]; paramjumptot[p1*nparam+th] = paramjumptot[p2*nparam+th];
						paramjumptot[p2*nparam+th] = tempf;
					}		
				}
			}
		}
	}
	
	MPI_Scatter(invTtot,nchain,MPI_DOUBLE,invT,nchain,MPI_DOUBLE,0,MPI_COMM_WORLD);
	MPI_Scatter(chtot,nchain,MPI_UNSIGNED,ch,nchain,MPI_UNSIGNED,0,MPI_COMM_WORLD);
	MPI_Scatter(timeproptot,nchain,MPI_LONG,timeprop,nchain,MPI_LONG,0,MPI_COMM_WORLD);
	MPI_Scatter(ntrtot,nchainparam,MPI_UNSIGNED,ntr,nchainparam,MPI_UNSIGNED,0,MPI_COMM_WORLD);
	MPI_Scatter(nactot,nchainparam,MPI_UNSIGNED,nac,nchainparam,MPI_UNSIGNED,0,MPI_COMM_WORLD);
	MPI_Scatter(paramjumptot,nchainparam,MPI_FLOAT,paramjump,nchainparam,MPI_FLOAT,0,MPI_COMM_WORLD);

	for(p = 0; p < nchain; p++){ 
		mbpchain[p]->invT = invT[p]; 
		mbpchain[p]->ch = ch[p];
		mbpchain[p]->timeprop = timeprop[p];
		for(th = 0; th < nparam; th++){
			mbpchain[p]->ntr[th] = ntr[p*nparam+th];
			mbpchain[p]->nac[th] = nac[p*nparam+th];
			mbpchain[p]->paramjump[th] = paramjump[p*nparam+th];
		}
	}
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
	
	return ME;
}

/// Ouputs a sample from the MBP algorithm
void MBPoutput(DATA &data, MODEL &model, POPTREE &poptree, vector <SAMPLE> &opsamp, unsigned int core, unsigned int ncore, unsigned int nchain)
{
	unsigned int p, ppost, nchaintot = ncore*nchain, samp = long(opsamp.size());
	int siz;
	double Li[nchain], Litot[nchaintot], Liord[nchaintot];
	unsigned int ch[nchain], chtot[nchaintot];

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
		double L;
		vector < vector <FEV> > xiplot;
		vector <double> paramplot;

		if(ppost < nchain){
			xiplot = mbpchain[ppost]->xi; paramplot = mbpchain[ppost]->paramval; L = mbpchain[ppost]->Li;
		}
		else{
			MPI_Status status;
			MPI_Recv(packbuffer(),MAX_NUMBERS,MPI_DOUBLE,ppost/nchain,0,MPI_COMM_WORLD,&status);
			MPI_Get_count(&status, MPI_DOUBLE, &siz); if(siz >= int(MAX_NUMBERS)) emsg("Buffer not big enough");
		
			packinit();
			unpack(xiplot,0,data.fediv);
			unpack(paramplot);
			unpack(L);
			if(packsize() != siz) emsg("PMBP: EC10");
		}
		
		opsamp.push_back(outputsamp(calcME(),samp,L,data,model,poptree,paramplot,xiplot));
	}
	else{
		if(core == ppost/nchain){
			p = ppost%nchain;
			packinit();
			pack(mbpchain[p]->xi,0,data.fediv);
			pack(mbpchain[p]->paramval);
			pack(mbpchain[p]->Li);
			MPI_Send(packbuffer(),packsize(),MPI_DOUBLE,0,0,MPI_COMM_WORLD);
		}
	}
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
	
	for(p = 0; p < nchain; p++){ 
		Li[p] = mbpchain[p]->Li; 
		ch[p] = mbpchain[p]->ch;
		timeprop[p] = mbpchain[p]->timeprop;
		for(th = 0; th < nparam; th++) ntr[p*nparam+th] =  mbpchain[p]->ntr[th];
		for(th = 0; th < nparam; th++) nac[p*nparam+th] =  mbpchain[p]->nac[th];
	}
		
	MPI_Gather(Li,nchain,MPI_DOUBLE,Litot,nchain,MPI_DOUBLE,0,MPI_COMM_WORLD);
	MPI_Gather(ch,nchain,MPI_UNSIGNED,chtot,nchain,MPI_UNSIGNED,0,MPI_COMM_WORLD);
	MPI_Gather(timeprop,nchain,MPI_LONG,timeproptot,nchain,MPI_LONG,0,MPI_COMM_WORLD);
	MPI_Gather(ntr,nchainparam,MPI_UNSIGNED,ntrtot,nchainparam,MPI_UNSIGNED,0,MPI_COMM_WORLD);
	MPI_Gather(nac,nchainparam,MPI_UNSIGNED,nactot,nchainparam,MPI_UNSIGNED,0,MPI_COMM_WORLD);

	if(core == 0){
		stringstream ss; ss << data.outputdir << "/MCMCdiagnostic.txt";
		ofstream diag(ss.str().c_str()); 
		for(c = 0; c < nchaintot; c++){
			cc = 0; while(cc < nchaintot && c != chtot[cc]) cc++;
			if(cc == nchaintot) emsg("MBP: EC9");			
		
			diag << "Chain: " << c << endl;
			diag << "invT: " << invTstore[c] << "  ";
			
			av = 0; jmax = Listore[c].size(); for(j = 0; j < jmax; j++) av += Listore[c][j];
			diag << "Liav: " << av/jmax << "  ";
			
			ntrsum = 0; for(th = 0; th < nparam; th++) ntrsum += ntrtot[cc*nparam+th];
			diag << "Time: " << int(1000*timeproptot[cc]/(ntrsum*CLOCKS_PER_SEC));
			diag << endl;
			
			diag << "Accept: ";
			for(th = 0; th < nparam; th++){
				if(ntrtot[cc*nparam+th] == 0) diag << "--- ";
				else diag << double(nactot[cc*nparam+th])/ntrtot[cc*nparam+th] << " ";
			}
			diag << endl << endl;
		}
	}
}
