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
#include "model_evidence.hh"

ofstream quenchplot;

/// In initialises the inference algorithm
Mcmc::Mcmc(const Details &details, DATA &data, MODEL &model, const POPTREE &poptree, const Mpi &mpi, Inputs &inputs, Output &output, Obsmodel &obsmodel) : details(details), data(data), model(model), poptree(poptree), mpi(mpi), output(output), obsmodel(obsmodel)
{
	nsamp = inputs.find_int("nsamp",UNSET);                                                 // Sets the number of samples for inference
	if(nsamp == UNSET) emsgroot("The number of samples must be set");
	
	burnin = nsamp/4;                                                                       // Sets the burnin to be a quater of the samples
	
	nchaintot = inputs.find_int("nchain",UNSET);                                            // Sets the total number of mcmc chains
	if(nchaintot == UNSET) emsgroot("The number of chains must be set");
	if(nchaintot%mpi.ncore != 0) emsgroot("The number of chains must be a multiple of the number of cores");

	nchain = nchaintot/mpi.ncore;                                                           // The number of chains per core
	if(nchain == 0) emsgroot("'nchain' must be non-zero");
			
	string propsmethod_str = inputs.find_string("propsmethod","fixedtime");                 // Sets the proposal method
	if(propsmethod_str == "allchainsallparams") propsmethod = proposalsmethod::allchainsallparams;
	else{
		if(propsmethod_str == "fixednum") propsmethod = proposalsmethod::fixednum;
		else{
			if(propsmethod_str == "fixedtime") propsmethod = proposalsmethod::fixedtime;
			else emsg("Parameter propsmethod set to an unrecognised value \""+propsmethod_str+"\"");
		}
	}
		
	for(auto p = 0u; p < nchain; p++){                                                     // Initialises the chains
		Chain ch = Chain(details,data,model,poptree,obsmodel,mpi.core*nchain+p);
		chain.push_back(ch);
	}
	
	auto quench = burnin;                                                                  // Initialises model_evidence
	auto invTmin = inputs.find_double("invTmin",0);                 
	auto invTmax = inputs.find_double("invTmax",0.25);       
	model_evidence.init(nchaintot,quench,invTmin,invTmax);

	if(mpi.core == 0){ output.trace_plot_init(); output.Li_trace_plot_init(mpi.ncore*nchain);}
	if(quenchpl == 1){ quenchplot.open((details.outputdir+"/quenchplot.txt").c_str());}
	
	if(mpi.core == 0){ nac_swap.resize(nchaintot); for(auto pp = 0u; pp < nchaintot; pp++) nac_swap[pp] = 0;}

	assert(mpi.ncore > 0);
  assert(mpi.core < mpi.ncore);
  assert(nchain > 0);
}

/// Runs the multi-temperature MCMC algorithm	
void Mcmc::run()
{
	vector <SAMPLE> opsamp;        // Stores output samples
	vector <PARAMSAMP> psamp;      // Stores parameter samples
	
	long timeprop=0, ntimeprop=0;
	auto timeloop=0.1;
	vector <long> timeproptot(mpi.ncore), ntimeproptot(mpi.ncore);

	for(auto samp = 0u; samp < nsamp; samp++){	
		if(mpi.core == 0 && samp%1 == 0) cout << " Sample: " << samp << " / " << nsamp << endl; 
	
		model_evidence.set_invT(samp,chain);
		
		long time = clock();

		for(auto p = 0u; p < nchain; p++) chain[p].standard_prop(samp,burnin);
		
		timeprop -= clock();
		switch(propsmethod){
		case proposalsmethod::allchainsallparams:
			for(auto p = 0u; p < nchain; p++){
				for(auto th = 0u; th < model.param.size(); th++){
					if(model.param[th].min != model.param[th].max){ chain[p].proposal(th,samp,burnin); ntimeprop++;}
				}
			}
			break;
			
		case proposalsmethod::fixednum:
			for(auto lo = 0u; lo < 10; lo++){
				auto p = int(ran()*nchain);
				auto th = (unsigned int)(ran()*model.param.size());
				if(model.param[th].min != model.param[th].max){ chain[p].proposal(th,samp,burnin); ntimeprop++;}
			}
			break;
			
		case proposalsmethod::fixedtime:
			do{                         // Does proposals for timeloop seconds (on average 20 proposals)
				auto p = int(ran()*nchain);
				auto th = int(ran()*model.param.size());
				if(model.param[th].min != model.param[th].max){ chain[p].proposal(th,samp,burnin); ntimeprop++;}
			}while(double(clock()-time)/CLOCKS_PER_SEC < timeloop);
			break;
		}
		
		timeprop += clock();

		if(samp%10 == 0){ for(auto p = 0u; p < nchain; p++) chain[p].setQmapi(0);}  // Recalcualtes Qmapi (numerical)
	
		timers.timewait -= clock();
		
		MPI_Barrier(MPI_COMM_WORLD); 
	
		timers.timewait += clock();

		// This part dynamically swadjusts timeloop so that approximately 20 proposals are performed per swap
		MPI_Gather(&timeprop,1,MPI_LONG,&timeproptot[0],1,MPI_LONG,0,MPI_COMM_WORLD); 
		MPI_Gather(&ntimeprop,1,MPI_LONG,&ntimeproptot[0],1,MPI_LONG,0,MPI_COMM_WORLD);
		if(mpi.core == 0){
			auto timeproptotsum = 0.0, ntimeproptotsum = 0.0; 
			for(auto co = 0u; co < mpi.ncore; co++){ timeproptotsum += timeproptot[co]; ntimeproptotsum += ntimeproptot[co];}
            
			// Update the time to run only if some proposals have run (otherwise it runs forever)
			if(ntimeproptotsum > 0 && timeproptotsum > 0) {
				timeloop = 20*double(timeproptotsum)/(ntimeproptotsum*CLOCKS_PER_SEC);
			}
		}
		MPI_Bcast(&timeloop,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
		
		if(duplicate == 0) swap(nac_swap,samp,nchain);
		
		if(samp%1 == 0) model_evidence.store(mpi,chain);
		if(samp%1 == 0) output_param(output,psamp);
		if(samp%5 == 0) output_meas(opsamp,nchain);
		
		if(mpi.core == 0 && quenchpl == 1){ 
			quenchplot << chain[0].invT;
			for(auto pval : chain[0].paramval) quenchplot << "\t" << pval; 
			quenchplot << endl;
		}	
	
		if(samp == nsamp-1 || (samp != 0 && samp%1000 == 0)){
		  //if(samp == nsamp-1){
			if(mpi.core == 0){
				output.results(psamp,opsamp);
				cout << "Model evidence: " << model_evidence.calculate() << endl;
			}
			diagnostic(nac_swap);
		}
	}	
}

/// Stochastically swaps chains with similar inverse temperatures 
void Mcmc::swap(vector <double> &nac_swap, unsigned int samp, unsigned int nchain)
{
 	unsigned int nparam = model.param.size(), nchaintot = nchain*mpi.ncore, nchainparam = nchain*nparam, nparamtot = nchainparam*mpi.ncore;
	auto loopmax = nchaintot*nchaintot;
	vector <double> Li(nchain), Litot(nchaintot);
	vector <double> invT(nchain), invTtot(nchaintot);		
	vector <unsigned int> ch(nchain), chtot(nchaintot);
	vector <long> timeprop(nchain), timeproptot(nchaintot);
	vector <unsigned int> ntr(nchainparam), ntrtot(nparamtot);
	vector <unsigned int> nac(nchainparam), nactot(nparamtot);
	vector <float> paramjump(nchainparam), paramjumptot(nparamtot);
	vector <unsigned int> ntrxi(nchainparam), ntrxitot(nparamtot);
	vector <unsigned int> nacxi(nchainparam), nacxitot(nparamtot);
	vector <float> paramjumpxi(nchainparam), paramjumpxitot(nparamtot);
	vector <float> sigmajump(nchain), sigmajumptot(nchaintot);
	vector <float> numaddrem(nchain), numaddremtot(nchaintot);
	vector <unsigned int> ntr_addrem(nchain), ntr_addremtot(nchaintot);
	vector <unsigned int> nac_addrem(nchain), nac_addremtot(nchaintot);

	timers.timeswap -= clock();
		
	for(auto p = 0u; p < nchain; p++){ 
		Li[p] = chain[p].Li; 
		invT[p] = chain[p].invT; 
		ch[p] = chain[p].ch;
		timeprop[p] = chain[p].timeprop;
		for(auto th = 0u; th < nparam; th++){
			ntr[p*nparam+th] = chain[p].ntr[th];
			nac[p*nparam+th] = chain[p].nac[th];
			paramjump[p*nparam+th] = chain[p].paramjump[th];
			
			ntrxi[p*nparam+th] = chain[p].ntrxi[th];
			nacxi[p*nparam+th] = chain[p].nacxi[th];
			paramjumpxi[p*nparam+th] = chain[p].paramjumpxi[th];
		}
		sigmajump[p] = chain[p].sigmajump;
		numaddrem[p] = chain[p].numaddrem;
		ntr_addrem[p] = chain[p].ntr_addrem;
		nac_addrem[p] = chain[p].nac_addrem;
	}
		
	MPI_Gather(&Li[0],nchain,MPI_DOUBLE,&Litot[0],nchain,MPI_DOUBLE,0,MPI_COMM_WORLD);
	MPI_Gather(&invT[0],nchain,MPI_DOUBLE,&invTtot[0],nchain,MPI_DOUBLE,0,MPI_COMM_WORLD);
	MPI_Gather(&ch[0],nchain,MPI_UNSIGNED,&chtot[0],nchain,MPI_UNSIGNED,0,MPI_COMM_WORLD);
	MPI_Gather(&timeprop[0],nchain,MPI_LONG,&timeproptot[0],nchain,MPI_LONG,0,MPI_COMM_WORLD);
	MPI_Gather(&ntr[0],nchainparam,MPI_UNSIGNED,&ntrtot[0],nchainparam,MPI_UNSIGNED,0,MPI_COMM_WORLD);
	MPI_Gather(&nac[0],nchainparam,MPI_UNSIGNED,&nactot[0],nchainparam,MPI_UNSIGNED,0,MPI_COMM_WORLD);
	MPI_Gather(&paramjump[0],nchainparam,MPI_FLOAT,&paramjumptot[0],nchainparam,MPI_FLOAT,0,MPI_COMM_WORLD);
	MPI_Gather(&ntrxi[0],nchainparam,MPI_UNSIGNED,&ntrxitot[0],nchainparam,MPI_UNSIGNED,0,MPI_COMM_WORLD);
	MPI_Gather(&nacxi[0],nchainparam,MPI_UNSIGNED,&nacxitot[0],nchainparam,MPI_UNSIGNED,0,MPI_COMM_WORLD);
	MPI_Gather(&paramjumpxi[0],nchainparam,MPI_FLOAT,&paramjumpxitot[0],nchainparam,MPI_FLOAT,0,MPI_COMM_WORLD);
	MPI_Gather(&sigmajump[0],nchain,MPI_FLOAT,&sigmajumptot[0],nchain,MPI_FLOAT,0,MPI_COMM_WORLD);
	MPI_Gather(&numaddrem[0],nchain,MPI_FLOAT,&numaddremtot[0],nchain,MPI_FLOAT,0,MPI_COMM_WORLD);
	MPI_Gather(&ntr_addrem[0],nchain,MPI_UNSIGNED,&ntr_addremtot[0],nchain,MPI_UNSIGNED,0,MPI_COMM_WORLD);
	MPI_Gather(&nac_addrem[0],nchain,MPI_UNSIGNED,&nac_addremtot[0],nchain,MPI_UNSIGNED,0,MPI_COMM_WORLD);

	if(mpi.core == 0){   // Swaps different chains
		for(auto loop = 0u; loop < loopmax; loop++){
			auto p1 = (unsigned int)(ran()*nchaintot);
			auto p2 = (unsigned int)(ran()*nchaintot); 
			if(p1 != p2){
				auto al = exp((invTtot[p1]-invTtot[p2])*(Litot[p2]-Litot[p1]));
				
				if(ran() < al){
					if(samp >= burnin){
						auto chs = chtot[p1]; if(chtot[p2] < chs) chs = chtot[p2];
						nac_swap[chs]++;
					}
					
					double temp;
					unsigned tempi;
					float tempf;
					long templ;
					
					temp = invTtot[p1]; invTtot[p1] = invTtot[p2]; invTtot[p2] = temp;
					tempi = chtot[p1]; chtot[p1] = chtot[p2]; chtot[p2] = tempi;
					templ = timeproptot[p1]; timeproptot[p1] = timeproptot[p2]; timeproptot[p2] = templ;
					for(auto th = 0u; th < nparam; th++){
						tempi = ntrtot[p1*nparam+th]; ntrtot[p1*nparam+th] = ntrtot[p2*nparam+th]; ntrtot[p2*nparam+th] = tempi;
						tempi = nactot[p1*nparam+th]; nactot[p1*nparam+th] = nactot[p2*nparam+th]; nactot[p2*nparam+th] = tempi;
						tempf = paramjumptot[p1*nparam+th]; paramjumptot[p1*nparam+th] = paramjumptot[p2*nparam+th];
						paramjumptot[p2*nparam+th] = tempf;
						
						tempi = ntrxitot[p1*nparam+th]; ntrxitot[p1*nparam+th] = ntrxitot[p2*nparam+th]; ntrxitot[p2*nparam+th] = tempi;
						tempi = nacxitot[p1*nparam+th]; nacxitot[p1*nparam+th] = nacxitot[p2*nparam+th]; nacxitot[p2*nparam+th] = tempi;
						tempf = paramjumpxitot[p1*nparam+th]; paramjumpxitot[p1*nparam+th] = paramjumpxitot[p2*nparam+th];
						paramjumpxitot[p2*nparam+th] = tempf;
					}		
					tempf = sigmajumptot[p1]; sigmajumptot[p1] = sigmajumptot[p2]; sigmajumptot[p2] = tempf; 
					tempf = numaddremtot[p1]; numaddremtot[p1] = numaddremtot[p2]; numaddremtot[p2] = tempf;
					tempi = ntr_addremtot[p1]; ntr_addremtot[p1] = ntr_addremtot[p2]; ntr_addremtot[p2] = tempi;
					tempi = nac_addremtot[p1]; nac_addremtot[p1] = nac_addremtot[p2]; nac_addremtot[p2] = tempi;
				}
			}
		}
	}
	
	MPI_Scatter(&invTtot[0],nchain,MPI_DOUBLE,&invT[0],nchain,MPI_DOUBLE,0,MPI_COMM_WORLD);
	MPI_Scatter(&chtot[0],nchain,MPI_UNSIGNED,&ch[0],nchain,MPI_UNSIGNED,0,MPI_COMM_WORLD);
	MPI_Scatter(&timeproptot[0],nchain,MPI_LONG,&timeprop[0],nchain,MPI_LONG,0,MPI_COMM_WORLD);
	MPI_Scatter(&ntrtot[0],nchainparam,MPI_UNSIGNED,&ntr[0],nchainparam,MPI_UNSIGNED,0,MPI_COMM_WORLD);
	MPI_Scatter(&nactot[0],nchainparam,MPI_UNSIGNED,&nac[0],nchainparam,MPI_UNSIGNED,0,MPI_COMM_WORLD);
	MPI_Scatter(&paramjumptot[0],nchainparam,MPI_FLOAT,&paramjump[0],nchainparam,MPI_FLOAT,0,MPI_COMM_WORLD);
	MPI_Scatter(&ntrxitot[0],nchainparam,MPI_UNSIGNED,&ntrxi[0],nchainparam,MPI_UNSIGNED,0,MPI_COMM_WORLD);
	MPI_Scatter(&nacxitot[0],nchainparam,MPI_UNSIGNED,&nacxi[0],nchainparam,MPI_UNSIGNED,0,MPI_COMM_WORLD);
	MPI_Scatter(&paramjumpxitot[0],nchainparam,MPI_FLOAT,&paramjumpxi[0],nchainparam,MPI_FLOAT,0,MPI_COMM_WORLD);
	MPI_Scatter(&sigmajumptot[0],nchain,MPI_FLOAT,&sigmajump[0],nchain,MPI_FLOAT,0,MPI_COMM_WORLD);
	MPI_Scatter(&numaddremtot[0],nchain,MPI_FLOAT,&numaddrem[0],nchain,MPI_FLOAT,0,MPI_COMM_WORLD);
	MPI_Scatter(&ntr_addremtot[0],nchain,MPI_UNSIGNED,&ntr_addrem[0],nchain,MPI_UNSIGNED,0,MPI_COMM_WORLD);
	MPI_Scatter(&nac_addremtot[0],nchain,MPI_UNSIGNED,&nac_addrem[0],nchain,MPI_UNSIGNED,0,MPI_COMM_WORLD);

	for(auto p = 0u; p < nchain; p++){ 
		chain[p].invT = invT[p]; 
		chain[p].ch = ch[p];
		chain[p].timeprop = timeprop[p];
		for(auto th = 0u; th < nparam; th++){
			chain[p].ntr[th] = ntr[p*nparam+th];
			chain[p].nac[th] = nac[p*nparam+th];
			chain[p].paramjump[th] = paramjump[p*nparam+th];
			
			chain[p].ntrxi[th] = ntrxi[p*nparam+th];
			chain[p].nacxi[th] = nacxi[p*nparam+th];
			chain[p].paramjumpxi[th] = paramjumpxi[p*nparam+th];
		}
		chain[p].sigmajump = sigmajump[p]; 
		chain[p].numaddrem = numaddrem[p];
		chain[p].ntr_addrem = ntr_addrem[p];
		chain[p].nac_addrem = nac_addrem[p];
	}
	
	timers.timeswap += clock();
}

/// Ouputs a parameter sample from the MBP algorithm
void Mcmc::output_param(Output &output, vector <PARAMSAMP> &psamp) const
{
	unsigned int nchaintot = mpi.ncore*nchain, samp = psamp.size();
	 
	timers.timeoutput -= clock();
		
	vector <double> Li(nchain), Litot(nchaintot);
	vector <unsigned int> ch(nchain), chtot(nchaintot);
	for(auto p = 0u; p < nchain; p++){ Li[p] = chain[p].Li; ch[p] = chain[p].ch;}
		
	MPI_Gather(&Li[0],nchain,MPI_DOUBLE,&Litot[0],nchain,MPI_DOUBLE,0,MPI_COMM_WORLD);
	MPI_Gather(&ch[0],nchain,MPI_UNSIGNED,&chtot[0],nchain,MPI_UNSIGNED,0,MPI_COMM_WORLD);

	unsigned int ppost;
	if(mpi.core == 0){
		vector <double> Liord(nchaintot);
		for(auto p = 0u; p < nchaintot; p++){
			if(chtot[p] == 0) ppost = p;
			Liord[chtot[p]] = Litot[p];
		}
		output.Li_trace_plot(samp,nchaintot,Liord);
	}
	
	MPI_Bcast(&ppost,1,MPI_UNSIGNED,0,MPI_COMM_WORLD);

	if(mpi.core == 0){
		double L, Pr;
		unsigned int ninfplot;
		PARAMSAMP paramsamp;
	
		if(ppost < nchain){
			paramsamp.paramval = chain[ppost].paramval;
			L = chain[ppost].Li; Pr = chain[ppost].Pri;
			ninfplot =  chain[ppost].xi.size();
		}
		else{
			MPI_Status status;
			packinit(MAX_NUMBERS);
			MPI_Recv(packbuffer(),MAX_NUMBERS,MPI_DOUBLE,ppost/nchain,0,MPI_COMM_WORLD,&status);
			int siz;
			MPI_Get_count(&status, MPI_DOUBLE, &siz); if(siz >= int(MAX_NUMBERS)) emsg("The buffer is not big enough");
		
			unpack(paramsamp.paramval);
			unpack(L);
			unpack(Pr);
			unpack(ninfplot);
		}
	
		output.trace_plot(samp,L,Pr,ninfplot,paramsamp.paramval);
		psamp.push_back(paramsamp);
	}
	else{
		if(mpi.core == ppost/nchain){
			auto p = ppost%nchain;
		
			packinit(0);
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
void Mcmc::output_meas(vector <SAMPLE> &opsamp, unsigned int nchain) const
{
	unsigned int nchaintot = mpi.ncore*nchain;

	timers.timeoutput -= clock();

	vector <unsigned int> ch(nchain), chtot(nchaintot);
	for(auto p = 0u; p < nchain; p++) ch[p] = chain[p].ch;
			
	MPI_Gather(&ch[0],nchain,MPI_UNSIGNED,&chtot[0],nchain,MPI_UNSIGNED,0,MPI_COMM_WORLD);
	
	unsigned int ppost = UNSET;
	if(mpi.core == 0){
		for(auto p = 0u; p < nchaintot; p++){ if(chtot[p] == 0){ ppost = p; break;}}
		if(ppost == UNSET) emsgEC("MBP",100);
	}
	
	MPI_Bcast(&ppost,1,MPI_UNSIGNED,0,MPI_COMM_WORLD);

	if(mpi.core == 0){
		SAMPLE sample;
		if(ppost < nchain){
			sample.meas = obsmodel.getmeas(chain[ppost].trevi,chain[ppost].indevi);
			model.setup(chain[ppost].paramval);
			sample.R0 = model.R0calc(chain[ppost].paramval);
			sample.phi = model.phi; 
		}
		else{
			MPI_Status status;
			packinit(MAX_NUMBERS);
			MPI_Recv(packbuffer(),MAX_NUMBERS,MPI_DOUBLE,ppost/nchain,0,MPI_COMM_WORLD,&status);
			int siz;
			MPI_Get_count(&status, MPI_DOUBLE, &siz); if(siz >= int(MAX_NUMBERS)) emsg("The buffer not big enough.");
		
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
			auto p = ppost%nchain;
			
			MEAS meas = obsmodel.getmeas(chain[p].trevi,chain[p].indevi);
			model.setup(chain[p].paramval);
			vector <double> R0 = model.R0calc(chain[p].paramval);
	
			packinit(0);
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
void Mcmc::diagnostic(vector <double> &nac_swap) const
{
	unsigned int nparam = model.param.size(), nchaintot = mpi.ncore*nchain, nchainparam = nchain*nparam, nparamtot = nchainparam*mpi.ncore;
	vector <double> Li(nchain), Litot(nchaintot);
	vector <unsigned int> ch(nchain), chtot(nchaintot);
	vector <long> timeprop(nchain), timeproptot(nchaintot);
	vector <unsigned int> ntr(nchainparam), ntrtot(nparamtot);
	vector <unsigned int> nac(nchainparam), nactot(nparamtot);
	vector <float> paramjump(nchainparam), paramjumptot(nparamtot);
	vector <unsigned int> ntrxi(nchainparam), ntrxitot(nparamtot);
	vector <unsigned int> nacxi(nchainparam), nacxitot(nparamtot);
	vector <float> paramjumpxi(nchainparam), paramjumpxitot(nparamtot);
	vector <float> numaddrem(nchain), numaddremtot(nchaintot);
	vector <unsigned int> ntr_addrem(nchain), ntr_addremtot(nchaintot);
	vector <unsigned int> nac_addrem(nchain), nac_addremtot(nchaintot);

	for(auto p = 0u; p < nchain; p++){ 
		Li[p] = chain[p].Li; 
		ch[p] = chain[p].ch;
		timeprop[p] = chain[p].timeprop;
		for(auto th = 0u; th < nparam; th++){
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
		
	MPI_Gather(&Li[0],nchain,MPI_DOUBLE,&Litot[0],nchain,MPI_DOUBLE,0,MPI_COMM_WORLD);
	MPI_Gather(&ch[0],nchain,MPI_UNSIGNED,&chtot[0],nchain,MPI_UNSIGNED,0,MPI_COMM_WORLD);
	MPI_Gather(&timeprop[0],nchain,MPI_LONG,&timeproptot[0],nchain,MPI_LONG,0,MPI_COMM_WORLD);
	MPI_Gather(&ntr[0],nchainparam,MPI_UNSIGNED,&ntrtot[0],nchainparam,MPI_UNSIGNED,0,MPI_COMM_WORLD);
	MPI_Gather(&nac[0],nchainparam,MPI_UNSIGNED,&nactot[0],nchainparam,MPI_UNSIGNED,0,MPI_COMM_WORLD);
	MPI_Gather(&paramjump[0],nchainparam,MPI_FLOAT,&paramjumptot[0],nchainparam,MPI_FLOAT,0,MPI_COMM_WORLD);
	MPI_Gather(&ntrxi[0],nchainparam,MPI_UNSIGNED,&ntrxitot[0],nchainparam,MPI_UNSIGNED,0,MPI_COMM_WORLD);
	MPI_Gather(&nacxi[0],nchainparam,MPI_UNSIGNED,&nacxitot[0],nchainparam,MPI_UNSIGNED,0,MPI_COMM_WORLD);
	MPI_Gather(&paramjumpxi[0],nchainparam,MPI_FLOAT,&paramjumpxitot[0],nchainparam,MPI_FLOAT,0,MPI_COMM_WORLD);
	MPI_Gather(&numaddrem[0],nchain,MPI_FLOAT,&numaddremtot[0],nchain,MPI_FLOAT,0,MPI_COMM_WORLD);
	MPI_Gather(&ntr_addrem[0],nchain,MPI_UNSIGNED,&ntr_addremtot[0],nchain,MPI_UNSIGNED,0,MPI_COMM_WORLD);
	MPI_Gather(&nac_addrem[0],nchain,MPI_UNSIGNED,&nac_addremtot[0],nchain,MPI_UNSIGNED,0,MPI_COMM_WORLD);

	if(mpi.core == 0){
		stringstream ss; ss << details.outputdir << "/MCMCdiagnostic.txt";
		ofstream diag(ss.str().c_str()); 
		ofstream timings(details.outputdir+"/MCMCdiagnostic_timings.txt"); 

		for(auto c = 0u; c < nchaintot; c++){
			unsigned int cc = 0; while(cc < nchaintot && c != chtot[cc]) cc++;
			if(cc == nchaintot) emsgEC("MBP",101);			
		
			diag << "-------- Chain: " << c << " ----------" << endl << endl;
			diag << "invT: " << model_evidence.get_invT(c) << "  ";
			
			diag << "Li average: " <<  model_evidence.average_Li(c) << "  ";
			
			if (false) {
				double ntrsum = 0; for(auto th = 0u; th < nparam; th++) ntrsum += ntrtot[cc*nparam+th];
			}
			diag << endl;
			
			diag << "MBP Accept: " << endl;
			for(auto th = 0u; th < nparam; th++){
				diag << model.param[th].name << ": ";
				if(ntrtot[cc*nparam+th] == 0) diag << "--- ";
				else diag << double(nactot[cc*nparam+th])/ntrtot[cc*nparam+th] << " jump:" << paramjumptot[cc*nparam+th] << ", ";
				diag << endl;
			}
			diag << endl;
					
			diag << "Standard Accept: " << endl;
			for(auto th = 0u; th < nparam; th++){
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
		for(auto c = 0u; c < nchaintot; c++){
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
