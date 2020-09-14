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

/// In initialises the inference algorithm
Mcmc::Mcmc(const Details &details, const Data &data, const Model &model, const AreaTree &areatree, const Mpi &mpi, const Inputs &inputs, Output &output, const ObservationModel &obsmodel) : details(details), data(data), model(model), areatree(areatree), mpi(mpi), output(output), obsmodel(obsmodel)
{
	nsamp = inputs.find_integer("nsamp",UNSET);                               // Sets the number of samples for inference
	if(nsamp == UNSET) emsgroot("The number of samples must be set");
	
	burnin = inputs.find_integer("burnin",nsamp/4);                           // Sets the burnin period
	
	nchain_total = inputs.find_integer("nchain",UNSET);                       // Sets the total number of mcmc chains
	if(nchain_total == UNSET) emsgroot("The number of chains must be set");
	if(nchain_total%mpi.ncore != 0) emsgroot("The number of chains must be a multiple of the number of cores");

	nchain = nchain_total/mpi.ncore;                                          // The number of chains per core
	if(nchain == 0) emsgroot("'nchain' must be non-zero");
			
	string propsmethod_str = inputs.find_string("propsmethod","fixedtime");   // Sets the proposal method
	if(propsmethod_str == "allchainsallparams") propsmethod = proposal_method::allchainsallparams;
	else{
		if(propsmethod_str == "fixednum") propsmethod = proposal_method::fixednum;
		else{
			if(propsmethod_str == "fixedtime") propsmethod = proposal_method::fixedtime;
			else emsg("Parameter propsmethod set to an unrecognised value \""+propsmethod_str+"\"");
		}
	}
		
	for(auto p = 0u; p < nchain; p++){                                         // Initialises the chains
		chain.push_back(Chain(details,data,model,areatree,obsmodel,output,mpi.core*nchain+p));
	}
	
	auto quench = burnin;                                                      // Initialises model_evidence
	auto invTmin = inputs.find_double("invTmin",0);                 
	auto invTmax = inputs.find_double("invTmax",0.25);       
	model_evidence.init(nchain_total,quench,invTmin,invTmax);

	if(mpi.core == 0){                                                         // Initialises trace plots
		output.trace_plot_inititialise(); 
		output.L_trace_plot_inititialise(mpi.ncore*nchain);
	}
	
	if(mpi.core == 0){                                                         // Initialises diagnostics for swapping
 		nac_swap.resize(nchain_total); for(auto& nac_swa : nac_swap) nac_swa = 0;
	}

	timeprop = 0; ntimeprop = 0;                                               // Variable that inform time between swaps
	timeloop = 0.1;
		
	assert(mpi.ncore > 0);
  assert(mpi.core < mpi.ncore);
  assert(nchain > 0);
}


/// Runs the MC3 algorithm	
void Mcmc::run()
{
	vector <Sample> opsamp;                                                    // Stores output samples
	vector <ParamSample> psamp;                                                // Stores parameter samples
	
	for(auto samp = 0u; samp < nsamp; samp++){	
		if(mpi.core == 0 && samp%1 == 0) cout << " Sample: " << samp << " / " << nsamp << endl; 
		
		for(auto& cha : chain) cha.jump.setburnin(samp,burnin);                  // Updates information in jump 
			
		model_evidence.set_invT(samp,chain);                                     // Sets the inverse temperature of chains
		
		update();                                                                // Updates the MCMC chains
		
		 // Recalculates Qmapi (this done tp avoid accumulated numerical error)
		if(samp%10 == 0){ for(auto& cha : chain) cha.initial.set_Qmap(0);}      
	
		optimise_timeloop();                                                     // Calculate the time between swaps 
	
		swap_chains();                                                           // Swaps chains cased on likelihood
			
		if(samp%1 == 0) model_evidence.store(mpi,chain);                         // Outputs results
		if(samp%1 == 0) output_parameters(output,psamp);
		if(samp%5 == 0) output_measurements(opsamp,nchain);
	
		if(samp == nsamp-1 || (samp != 0 && samp%1000 == 0)){
			if(mpi.core == 0){
				output.final_results(psamp,opsamp);
				cout << "Model evidence: " << model_evidence.calculate() << endl;
			}

			diagnostics();                                                         // Outputs MCMC diagnostic information
		}
	}	
}


/// Calculate timeloop, the CPU time used between swaps are performed
/// This is dynamically adjusted so that approximately 20 proposals are performed per swap
void Mcmc::optimise_timeloop()
{
	mpi_barrier();                                                              // Calculate the time waiting for other cores

	auto sum_timeprop = mpi_sum(timeprop);                                      // The total time for proposals
	auto sum_ntimeprop = mpi_sum(ntimeprop);                                    // The total number of proposals
	
	if(mpi.core == 0 && sum_timeprop > 0 && sum_ntimeprop  > 0){				
		timeloop = 20*double(sum_timeprop)/(sum_ntimeprop*CLOCKS_PER_SEC);
	}
	
	mpi_bcast(timeloop);
}


/// Performs MCMC updates on all the chains
void Mcmc::update()
{
	long time = clock();
		
	for(auto& cha : chain) cha.standard_proposal();                            // Standard proposals
		
	timeprop -= clock();
	switch(propsmethod){                                                       // MBPs using different methods
	case proposal_method::allchainsallparams:
		for(auto& cha : chain){
			for(auto th = 0u; th < model.param.size(); th++){
				if(model.param[th].min != model.param[th].max){ cha.mbp_proposal(th); ntimeprop++;}
			}
		}
		break;
		
	case proposal_method::fixednum:
		for(auto lo = 0u; lo < 10; lo++){
			auto p = int(ran()*nchain);
			auto th = (unsigned int)(ran()*model.param.size());
			if(model.param[th].min != model.param[th].max){ chain[p].mbp_proposal(th); ntimeprop++;}
		}
		break;
		
	case proposal_method::fixedtime:
		do{                         // Does proposals for timeloop seconds (on average 20 proposals)
			auto p = int(ran()*nchain);
			auto th = int(ran()*model.param.size());
			if(model.param[th].min != model.param[th].max){ chain[p].mbp_proposal(th); ntimeprop++;}
		}while(double(clock()-time)/CLOCKS_PER_SEC < timeloop);
		break;
	}
	timeprop += clock();
}


/// Stochastically generates swap chain proposals 
void Mcmc::swap_chains()
{
	if(mpi.core == 0){                                                 // Gather chain infomation to core 0
		vector <ChainInfo> chaininfo(nchain*mpi.ncore);
		vector <double> Ltot(nchain*mpi.ncore);
		
		for(auto p = 0u; p < nchain; p++){
			Ltot[p] = chain[p].initial.L;
			chaininfo_set(chaininfo[p],chain[p]);
		}
		
		for(auto co = 1u; co < mpi.ncore; co++){
			pack_mpi_recv(co);
			for(auto p = 0u; p < nchain; p++){
				auto pp = co*nchain + p;
				unpack(Ltot[pp]);
				unpack(chaininfo[pp]);
			}
		}
		
		auto loopmax = nchain_total*nchain_total;
		for(auto loop = 0u; loop < loopmax; loop++){                      // Randomly selects chains and performs MH proposal
			auto p1 = (unsigned int)(ran()*nchain_total);
			auto p2 = (unsigned int)(ran()*nchain_total); 
			if(p1 != p2){
				auto al = exp((chaininfo[p1].invT-chaininfo[p2].invT)*(Ltot[p2]-Ltot[p1]));
				
				if(ran() < al){
					ChainInfo temp = chaininfo[p1]; 
					chaininfo[p1] = chaininfo[p2];
					chaininfo[p2] = temp; 					
				}
			}
		}
		
		for(auto p = 0u; p < nchain; p++) chain_set(chain[p],chaininfo[p]);
		
		for(auto co = 1u; co < mpi.ncore; co++){                          // Chain information is send to the other cores
			pack_initialise(0);
			for(auto p = 0u; p < nchain; p++) pack(chaininfo[co*nchain+p]);
			pack_mpi_send(co);
		}	
	}
	else{
		pack_initialise(0);
		for(auto p = 0u; p < nchain; p++){
			pack(chain[p].initial.L);
			ChainInfo chinfo;
			chaininfo_set(chinfo,chain[p]);
			pack(chinfo);
		}
		pack_mpi_send(0);
		
		pack_mpi_recv(0);
		for(auto p = 0u; p < nchain; p++){
			ChainInfo chinfo;
			unpack(chinfo);
			chain_set(chain[p],chinfo);
		}
	}
}


/// Sets chaininf properties from a chain
void Mcmc::chaininfo_set(ChainInfo &chinf, const Chain &cha) const
{
	chinf.invT = cha.invT;
	chinf.ch = cha.ch;
	chinf.jump = cha.jump;
}


/// Sets chain properties from chaininfo
void Mcmc::chain_set(Chain &cha, const ChainInfo &chinf) const 
{
	cha.invT = chinf.invT;
	cha.ch = chinf.ch;
	cha.jump = chinf.jump;
}		
	
	
/// Ouputs a parameter sample from the MBP algorithm
void Mcmc::output_parameters(Output &output, vector <ParamSample> &psamp) const
{
	unsigned int nchain_total = mpi.ncore*nchain, samp = psamp.size();
	 
	timers.timeoutput -= clock();
		
	vector <double> L(nchain);
	vector <unsigned int> ch(nchain);
	
	for(auto p = 0u; p < nchain; p++){ L[p] = chain[p].initial.L; ch[p] = chain[p].ch;}
	
	auto Ltot = mpi_gather(L);
	auto chtot = mpi_gather(ch);

	unsigned int ppost;
	if(mpi.core == 0){
		vector <double> Lord(nchain_total);
		for(auto p = 0u; p < nchain_total; p++){
			if(chtot[p] == 0) ppost = p;
			Lord[chtot[p]] = Ltot[p];
		}
		output.L_trace_plot(samp,Lord);
	}
	
	mpi_bcast(ppost);

	if(mpi.core == 0){
		double L, Pr;
		unsigned int ninfplot;
		ParamSample paramsamp;
	
		if(ppost < nchain){
			paramsamp.paramval = chain[ppost].initial.paramval;
			L = chain[ppost].initial.L; Pr = chain[ppost].initial.Pr;
			ninfplot = chain[ppost].initial.infev.size();
		}
		else{
			pack_mpi_recv(ppost/nchain);
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
		
			pack_initialise(0);
			pack(chain[p].initial.paramval);
			pack(chain[p].initial.L);
			pack(chain[p].initial.Pr);
			pack((unsigned int) chain[p].initial.infev.size());
			pack_mpi_send(0);
		}
	}
	
	timers.timeoutput += clock();
}

/// Ouputs a measurement sample from the MBP algorithm
void Mcmc::output_measurements(vector <Sample> &opsamp, unsigned int nchain) const
{
	unsigned int nchain_total = mpi.ncore*nchain;

	timers.timeoutput -= clock();

	vector <unsigned int> ch(nchain);
	for(auto p = 0u; p < nchain; p++) ch[p] = chain[p].ch;		
	auto chtot = mpi_gather(ch);
	
	unsigned int ppost = UNSET;
	if(mpi.core == 0){
		for(auto p = 0u; p < nchain_total; p++){ if(chtot[p] == 0){ ppost = p; break;}}
		if(ppost == UNSET) emsgEC("MBP",100);
	}
	
	mpi_bcast(ppost);

	if(mpi.core == 0){
		Sample sample;
		if(ppost < nchain){
			sample.meas = obsmodel.get_measured_quantities(chain[ppost].initial.transev,chain[ppost].initial.indev);
			sample.R0 = model.calculate_R_vs_time(chain[ppost].initial.paramval);
			sample.phi = model.create_disc_spline(model.phispline_ref,chain[ppost].initial.paramval); 
		}
		else{
			pack_mpi_recv(ppost/nchain);
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
			
			Measurements meas = obsmodel.get_measured_quantities(chain[p].initial.transev,chain[p].initial.indev);
			vector <double> R0 = model.calculate_R_vs_time(chain[p].initial.paramval);
	
			pack_initialise(0);
			pack(meas.transnum);
			pack(meas.popnum);
			pack(meas.margnum);
			pack(R0);
			pack(model.create_disc_spline(model.phispline_ref,chain[p].initial.paramval));
			pack_mpi_send(0);
		}
	}

	timers.timeoutput += clock();
}

/// Outputs MCMC diagnoistic information
void Mcmc::diagnostics() const
{
  if(mpi.core == 0){
		auto nchain_total = nchain*mpi.ncore;
		auto nparam = model.param.size();
		
		vector <ChainInfo> chaininfo(nchain_total);
		vector <double> Ltot(nchain_total);
		
		for(auto p = 0u; p < nchain; p++){
			Ltot[p] = chain[p].initial.L;
			chaininfo_set(chaininfo[p],chain[p]);
		}
		
		for(auto co = 1u; co < mpi.ncore; co++){
			pack_mpi_recv(co);
			for(auto p = 0u; p < nchain; p++){
				auto pp = co*nchain + p;
				unpack(Ltot[pp]);
				unpack(chaininfo[pp]);
			}
		}
		
		stringstream ss; ss << details.output_directory << "/MCMCdiagnostic.txt";
		ofstream diag(ss.str().c_str()); 
		ofstream timings(details.output_directory+"/MCMCdiagnostic_timings.txt"); 

		for(auto c = 0u; c < nchain_total; c++){
			unsigned int cc = 0; while(cc < nchain_total && c != chaininfo[cc].ch) cc++;
			if(cc == nchain_total) emsgEC("Mcmc",101);			
		
			diag << "-------- Chain: " << c << " ----------" << endl << endl;
			diag << "invT: " << model_evidence.get_invT(c) << "  ";
			
			diag << "Li average: " <<  model_evidence.average_L(c) << "  ";
			
			diag << endl;
			
			diag << "MBP Accept: " << endl;
			for(auto th = 0u; th < nparam; th++){
				diag << model.param[th].name << ": ";
				if(chaininfo[cc].jump.mbp_ntr[th] == 0) diag << "--- ";
				else{
					diag << double(chaininfo[cc].jump.mbp_nac[th])/chaininfo[cc].jump.mbp_ntr[th];
					diag << " jump:" << chaininfo[cc].jump.mbp[th] << ", ";
				}
				diag << endl;
			}
			diag << endl;
					
			diag << "Standard Accept: " << endl;
			for(auto th = 0u; th < nparam; th++){
				diag << model.param[th].name << ": ";
				if(chaininfo[cc].jump.stand_ntr[th] == 0) diag << "--- ";
				else{
					diag << double(chaininfo[cc].jump.stand_nac[th])/chaininfo[cc].jump.stand_ntr[th];
					diag << " jump:" << chaininfo[cc].jump.stand[th] << ", ";
				}
				diag << endl;
			}
			diag << endl;
	
			diag << "Add / remove infected: " << endl;
			if(chaininfo[cc].jump.standev_ntr == 0) diag << "--- "; 
			else{
				diag << double(chaininfo[cc].jump.standev_nac)/chaininfo[cc].jump.standev_ntr;
				diag << " num:" << chaininfo[cc].jump.naddrem << ", ";
			}
			diag << endl << endl;
		}
		
		diag << "Swaping with a chain at lower inverse temperature" << endl;
		for(auto c = 0u; c < nchain_total; c++){
			diag << "Chain " << c << ": " << nac_swap[c] << endl; 
		}
				
		timings << endl << "Timings for different parts of the algorithm:" << endl;
		timings << double(timers.timewait)/CLOCKS_PER_SEC << " MBP waiting time (seconds)" << endl;
		timings << double(timers.timembp)/CLOCKS_PER_SEC << " MBP time (seconds)" << endl;
		timings << double(timers.timembpinit)/CLOCKS_PER_SEC << " MBP init (seconds)" << endl;
		timings << double(timers.timembpQmap)/CLOCKS_PER_SEC << " MBP Qmap (seconds)" << endl;
		timings << double(timers.infection_sampler)/CLOCKS_PER_SEC << " MBP conRtot (seconds)" << endl;
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
	else{
		pack_initialise(0);
		for(auto p = 0u; p < nchain; p++){
			pack(chain[p].initial.L);
			ChainInfo chinfo;
			chaininfo_set(chinfo,chain[p]);
			pack(chinfo);
		}
		pack_mpi_send(0);
	}
}
