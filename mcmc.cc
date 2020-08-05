// Performs MBPs on multiple chains spanning from the posterior to the prior

#include <fstream>
#include <iostream>
#include <sstream>
#include <cassert>
#include "stdlib.h"
#include "math.h"

using namespace std;

#include "model.hh"
#include "chain.hh"
#include "output.hh"
#include "pack.hh"
#include "utils.hh"
#include "consts.hh"
#include "timers.hh"
#include "mcmc.hh"
#include "obsmodel.hh"

ofstream quenchplot;

Mcmc::Mcmc(Details &details, DATA &data, MODEL &model, POPTREE &poptree, Mpi &mpi, Inputs &inputs, Mode mode, bool verbose) : details(details), data(data), model(model), poptree(poptree), mpi(mpi)
{
	nsamp = inputs.find_int("nsamp",UNSET);                                             // Sets the number of samples for inference
	if(nsamp == UNSET) emsgroot("The number of samples must be set");
	
	burnin = nsamp/4;
	quench = nsamp/4;	
	
	nchaintot = inputs.find_int("nchain",UNSET);                                        // Sets the total number of mcmc chains
	if(nchaintot == UNSET) emsgroot("The number of chains must be set");
	if(nchaintot%mpi.ncore != 0) emsgroot("The number of chains must be a multiple of the number of cores");

	nchain = nchaintot/mpi.ncore;                                                           // The number of chains per core
	
	invTmin = inputs.find_double("invTmin",0);                                             // The temperatures of the prior and posterior chains
	invTmax = inputs.find_double("invTmax",0.25);       
	
	assert(mpi.ncore > 0);
  assert(mpi.core < mpi.ncore);
  assert(nchain > 0);
	assert(quench >= 0);
	assert(burnin >= 0);
}

/// Runs the multi-temperature MCMC algorithm	
void Mcmc::run(enum proposalsmethod propmethod)
{
	unsigned int p, pp, th, nchaintot = mpi.ncore*nchain, co;
	unsigned int samp;
	long time, timeprop=0, ntimeprop=0;
	double timeloop=0.1, pmin, pmax, kappa, ppf, ppeff, fac;
	long timeproptot[mpi.ncore], ntimeproptot[mpi.ncore], timeproptotsum, ntimeproptotsum;
	vector <SAMPLE> opsamp;
	vector <PARAMSAMP> psamp;
	vector <vector <double> > Listore;
	vector <double> invTstore;
	vector <double> nac_swap;
	
	Output output(details,data,model);
	
	
	pmax = 1-pow(invTmin,1.0/Tpower);
	pmin = 1-pow(invTmax,1.0/Tpower);
	
	for(p = 0; p < nchain; p++){
		pp = mpi.core*nchain+p;
		Chain ch = Chain(details,data,model,poptree,pp);
		chain.push_back(ch);
	}

	if(mpi.core == 0){
		Listore.resize(nchaintot); invTstore.resize(nchaintot);
		for(pp = 0; pp < nchaintot; pp++){
			kappa = double(pp)/(nchaintot-1);
			if(nchaintot == 1 || duplicate == 1) invTstore[pp] = invTmax;
			else invTstore[pp] = pow(1-(pmin+kappa*(pmax-pmin)),Tpower);
		}
		
		nac_swap.resize(nchaintot);
		for(pp = 0; pp < nchaintot; pp++) nac_swap[pp] = 0;
	}

	if(checkon == 1){
		for(th = 0; th < chain[0].paramval.size(); th++){
			cout << th << " " << model.param[th].name << " " << chain[0].paramval[th] <<  "param" << endl;		
		}
	}
	
	if(mpi.core == 0){ output.init(); output.Liinit(mpi.ncore*nchain);}

	if(quenchpl == 1){ quenchplot.open((data.outputdir+"/quenchplot.txt").c_str());}
	
	for(samp = 0; samp < nsamp; samp++){	
		if(mpi.core == 0 && samp%1 == 0) cout << " Sample: " << samp << " / " << nsamp << endl; 
	
		for(p = 0; p < nchain; p++){
			pp = chain[p].ch;
		
			if(nchaintot == 1 || duplicate == 1) chain[p].invT = invTmax;
			else{
				if(samp < quench) fac = double(samp)/quench;
				else fac = 1;
			
				kappa = double(pp)/(nchaintot-1);
				ppf = pmin+kappa*(pmax-pmin);
				ppeff = 1-fac*(1-ppf);	
			
				chain[p].invT = pow(1-ppeff,Tpower);
			}

			if(quenchpl == 1) chain[p].invT = exp(-10+12*double(samp)/nsamp);
		}
			
		time = clock();

		for(p = 0; p < nchain; p++) chain[p].standard_prop(samp,burnin);
		
		timeprop -= clock();
		switch(propmethod){
		case proposalsmethod::allchainsallparams:
			for(p = 0; p < nchain; p++){
				for(th = 0; th < model.param.size(); th++){
					if(model.param[th].min != model.param[th].max){ chain[p].proposal(th,samp,burnin); ntimeprop++;}
				}
			}
			break;
			
		case proposalsmethod::fixednum:
			short lo; 
			for(lo = 0; lo < 10; lo++){
				p = int(ran()*nchain);
				th = (unsigned int)(ran()*model.param.size());
				if(model.param[th].min != model.param[th].max){ chain[p].proposal(th,samp,burnin); ntimeprop++;}
			}
			break;
			
		case proposalsmethod::fixedtime:
			do{                         // Does proposals for timeloop seconds (on average 20 proposals)
				p = int(ran()*nchain);
				th = (unsigned int)(ran()*model.param.size());
				if(model.param[th].min != model.param[th].max){ chain[p].proposal(th,samp,burnin); ntimeprop++;}
			}while(double(clock()-time)/CLOCKS_PER_SEC < timeloop);
			break;
		}
		
		timeprop += clock();

		if(samp%10 == 0){ for(p = 0; p < nchain; p++) chain[p].setQmapi(0);}  // Recalcualtes Qmapi (numerical)
	
		timers.timewait -= clock();
		
		MPI_Barrier(MPI_COMM_WORLD); 
	
		timers.timewait += clock();

		// This part dynamically swadjusts timeloop so that approximately 20 proposals are performed per swap
		MPI_Gather(&timeprop,1,MPI_LONG,timeproptot,1,MPI_LONG,0,MPI_COMM_WORLD); 
		MPI_Gather(&ntimeprop,1,MPI_LONG,ntimeproptot,1,MPI_LONG,0,MPI_COMM_WORLD);
		if(mpi.core == 0){
			timeproptotsum = 0; ntimeproptotsum = 0; 
			for(co = 0; co < mpi.ncore; co++){ timeproptotsum += timeproptot[co]; ntimeproptotsum += ntimeproptot[co];}
            
			// Update the time to run only if some proposals have run (otherwise it runs forever)
			if (ntimeproptotsum > 0 && timeproptotsum > 0) {
				timeloop = 20*double(timeproptotsum)/(ntimeproptotsum*CLOCKS_PER_SEC);
			}
		}
		MPI_Bcast(&timeloop,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
		
		if(duplicate == 0) swap(nac_swap,samp);
		
		if(samp%1 == 0) output_param(output,Listore,psamp);
		if(samp%5 == 0) output_meas(opsamp);
		
		if(mpi.core == 0 && quenchpl == 1){ 
			quenchplot << chain[0].invT;
			for(p = 0; p < chain[0].paramval.size(); p++) quenchplot << "\t" << chain[0].paramval[p]; 
			quenchplot << endl;
		}	
	
		if(samp == nsamp-1 || (samp != 0 && samp%1000 == 0)){
		  //if(samp == nsamp-1){
			if(mpi.core == 0){
				output.results(psamp,opsamp);
				cout << "Model evidence: " << calcME(Listore,invTstore) << endl;
			}
			diagnostic(Listore,invTstore,nac_swap);
		}
	}	
}

/// Stochastically swaps chains with similar inverse temperatures 
void Mcmc::swap(vector <double> &nac_swap, unsigned int samp)
{
  unsigned int p, p1, p2, th, tempi, chs;
	unsigned int nparam = model.param.size(), nchaintot = nchain*mpi.ncore, nchainparam = nchain*nparam, nparamtot = nchainparam*mpi.ncore;
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
	unsigned int ntrxi[nchainparam], ntrxitot[nparamtot];
	unsigned int nacxi[nchainparam], nacxitot[nparamtot];
	float paramjumpxi[nchainparam], paramjumpxitot[nparamtot];
	float numaddrem[nchain], numaddremtot[nchaintot];
	unsigned int ntr_addrem[nchain], ntr_addremtot[nchaintot];
	unsigned int nac_addrem[nchain], nac_addremtot[nchaintot];

	timers.timeswap -= clock();
		
	for(p = 0; p < nchain; p++){ 
		Li[p] = chain[p].Li; 
		invT[p] = chain[p].invT; 
		ch[p] = chain[p].ch;
		timeprop[p] = chain[p].timeprop;
		for(th = 0; th < nparam; th++){
			ntr[p*nparam+th] = chain[p].ntr[th];
			nac[p*nparam+th] = chain[p].nac[th];
			paramjump[p*nparam+th] = chain[p].paramjump[th];
			
			ntrxi[p*nparam+th] = chain[p].ntrxi[th];
			nacxi[p*nparam+th] = chain[p].nacxi[th];
			paramjumpxi[p*nparam+th] = chain[p].paramjumpxi[th];
		}
		numaddrem[p] = chain[p].numaddrem;
		ntr_addrem[p] = chain[p].ntr_addrem;
		nac_addrem[p] = chain[p].nac_addrem;
	}
		
	MPI_Gather(Li,nchain,MPI_DOUBLE,Litot,nchain,MPI_DOUBLE,0,MPI_COMM_WORLD);
	MPI_Gather(invT,nchain,MPI_DOUBLE,invTtot,nchain,MPI_DOUBLE,0,MPI_COMM_WORLD);
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

	if(mpi.core == 0){   // Swaps different chains
		for(loop = 0; loop < loopmax; loop++){
			p1 = (unsigned int)(ran()*nchaintot); p2 = (unsigned int)(ran()*nchaintot); 
			if(p1 != p2){
				al = exp((invTtot[p1]-invTtot[p2])*(Litot[p2]-Litot[p1]));
				
				if(ran() < al){
					if(samp >= burnin){
						chs = chtot[p1]; if(chtot[p2] < chs) chs = chtot[p2];
						nac_swap[chs]++;
					}
					
					temp = invTtot[p1]; invTtot[p1] = invTtot[p2]; invTtot[p2] = temp;
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
		chain[p].invT = invT[p]; 
		chain[p].ch = ch[p];
		chain[p].timeprop = timeprop[p];
		for(th = 0; th < nparam; th++){
			chain[p].ntr[th] = ntr[p*nparam+th];
			chain[p].nac[th] = nac[p*nparam+th];
			chain[p].paramjump[th] = paramjump[p*nparam+th];
			
			chain[p].ntrxi[th] = ntrxi[p*nparam+th];
			chain[p].nacxi[th] = nacxi[p*nparam+th];
			chain[p].paramjumpxi[th] = paramjumpxi[p*nparam+th];
		}
		chain[p].numaddrem = numaddrem[p];
		chain[p].ntr_addrem = ntr_addrem[p];
		chain[p].nac_addrem = nac_addrem[p];
	}
	
	timers.timeswap += clock();
}

/// Calculates the model evidence
double Mcmc::calcME(vector <vector <double> > &Listore,vector <double> &invTstore ) const
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

/// Ouputs a parameter sample from the MBP algorithm
void Mcmc::output_param(Output &output, vector <vector <double> > &Listore, vector <PARAMSAMP> &psamp) const
{
	unsigned int p, ppost, nchaintot = mpi.ncore*nchain, samp = psamp.size();
	int siz;
	double Li[nchain], Litot[nchaintot], Liord[nchaintot];
	unsigned int ch[nchain], chtot[nchaintot];
	PARAMSAMP paramsamp;

	timers.timeoutput -= clock();
		
	for(p = 0; p < nchain; p++){ Li[p] = chain[p].Li; ch[p] = chain[p].ch;}
		
	MPI_Gather(Li,nchain,MPI_DOUBLE,Litot,nchain,MPI_DOUBLE,0,MPI_COMM_WORLD);
	MPI_Gather(ch,nchain,MPI_UNSIGNED,chtot,nchain,MPI_UNSIGNED,0,MPI_COMM_WORLD);

	if(mpi.core == 0){
		for(p = 0; p < nchaintot; p++){
			if(chtot[p] == 0) ppost = p;
			Listore[chtot[p]].push_back(Litot[p]);
			Liord[chtot[p]] = Litot[p];
		}
		output.Li(samp,nchaintot,Liord);
	}
	
	MPI_Bcast(&ppost,1,MPI_UNSIGNED,0,MPI_COMM_WORLD);

	if(mpi.core == 0){
		double L, Pr;
		unsigned int ninfplot;
	
		if(ppost < nchain){
			paramsamp.paramval = chain[ppost].paramval;
			L = chain[ppost].Li; Pr = chain[ppost].Pri;
			ninfplot =  chain[ppost].xi.size();
		}
		else{
			MPI_Status status;
			MPI_Recv(packbuffer(),MAX_NUMBERS,MPI_DOUBLE,ppost/nchain,0,MPI_COMM_WORLD,&status);
			MPI_Get_count(&status, MPI_DOUBLE, &siz); if(siz >= int(MAX_NUMBERS)) emsg("The buffer is not big enough");
		
			packinit();
			unpack(paramsamp.paramval);
			unpack(L);
			unpack(Pr);
			unpack(ninfplot);
		}
	
		output.traceplot(samp,L,Pr,ninfplot,paramsamp.paramval);
		psamp.push_back(paramsamp);
	}
	else{
		if(mpi.core == ppost/nchain){
			p = ppost%nchain;
		
			packinit();
			pack(chain[p].paramval);
			pack(chain[p].Li);
			pack(chain[p].Pri);
			pack((unsigned int) chain[p].xi.size());
			MPI_Send(packbuffer(),packsize(),MPI_DOUBLE,0,0,MPI_COMM_WORLD);
		}
	}
	
	timers.timeoutput += clock();
}

/// Ouputs a measurement sample from the MBP algorithm
void Mcmc::output_meas(vector <SAMPLE> &opsamp) const
{
	unsigned int p,  ppost, nchaintot = mpi.ncore*nchain;
	int siz;
	unsigned int ch[nchain], chtot[nchaintot];
	MEAS meas;
	vector <double> R0;	 
	SAMPLE sample;

	timers.timeoutput -= clock();

	for(p = 0; p < nchain; p++) ch[p] = chain[p].ch;
			
	MPI_Gather(ch,nchain,MPI_UNSIGNED,chtot,nchain,MPI_UNSIGNED,0,MPI_COMM_WORLD);
	
	if(mpi.core == 0){
		for(p = 0; p < nchaintot; p++){ if(chtot[p] == 0){ ppost = p; break;}}
		if(p == nchaintot) emsgEC("MBP",100);
	}
	
	MPI_Bcast(&ppost,1,MPI_UNSIGNED,0,MPI_COMM_WORLD);

	if(mpi.core == 0){
		if(ppost < nchain){
			sample.meas = getmeas(data,model,chain[ppost].trevi,chain[ppost].indevi);
			model.setup(chain[ppost].paramval);
			sample.R0 = model.R0calc();
			sample.phi = model.phi; 
		}
		else{
			MPI_Status status;
			MPI_Recv(packbuffer(),MAX_NUMBERS,MPI_DOUBLE,ppost/nchain,0,MPI_COMM_WORLD,&status);
			MPI_Get_count(&status, MPI_DOUBLE, &siz); if(siz >= int(MAX_NUMBERS)) emsg("The buffer not big enough.");
		
			packinit();
			unpack(sample.meas.transnum);
			unpack(sample.meas.popnum);
			unpack(sample.meas.margnum);
			unpack(sample.R0);
			unpack(sample.phi);
		}
		opsamp.push_back(sample);
	}
	else{
		if(mpi.core == ppost/nchain){
			p = ppost%nchain;
			
			meas = getmeas(data,model,chain[p].trevi,chain[p].indevi);
			model.setup(chain[p].paramval);
			R0 = model.R0calc();
	
			packinit();
			pack(meas.transnum);
			pack(meas.popnum);
			pack(meas.margnum);
			pack(R0);
			pack(model.phi);
			MPI_Send(packbuffer(),packsize(),MPI_DOUBLE,0,0,MPI_COMM_WORLD);
		}
	}

	timers.timeoutput += clock();
}

/// Outputs MCMC diagnoistic information
void Mcmc::diagnostic(vector <vector <double> > &Listore, vector <double> &invTstore, vector <double> &nac_swap) const
{
	unsigned int p, c, cc, th, j, jmax;
	unsigned int nparam = model.param.size(), nchaintot = mpi.ncore*nchain, nchainparam = nchain*nparam, nparamtot = nchainparam*mpi.ncore;
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
		Li[p] = chain[p].Li; 
		ch[p] = chain[p].ch;
		timeprop[p] = chain[p].timeprop;
		for(th = 0; th < nparam; th++){
			ntr[p*nparam+th] = chain[p].ntr[th];
			nac[p*nparam+th] = chain[p].nac[th];
			paramjump[p*nparam+th] = chain[p].paramjump[th];
			
			ntrxi[p*nparam+th] = chain[p].ntrxi[th];
			nacxi[p*nparam+th] = chain[p].nacxi[th];
			paramjumpxi[p*nparam+th] = chain[p].paramjumpxi[th];
		}
		numaddrem[p] = chain[p].numaddrem;
		ntr_addrem[p] = chain[p].ntr_addrem;
		nac_addrem[p] = chain[p].nac_addrem;
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

	if(mpi.core == 0){
		stringstream ss; ss << data.outputdir << "/MCMCdiagnostic.txt";
		ofstream diag(ss.str().c_str()); 
		ofstream timings(data.outputdir+"/MCMCdiagnostic_timings.txt"); 

		for(c = 0; c < nchaintot; c++){
			cc = 0; while(cc < nchaintot && c != chtot[cc]) cc++;
			if(cc == nchaintot) emsgEC("MBP",101);			
		
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
		
		diag << "Swaping with a chain at lower inverse temperature" << endl;
		for(c = 0; c < nchaintot; c++){
			diag << "Chain " << c << ": " << nac_swap[c] << endl; 
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
