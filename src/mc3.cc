/// Implement the metropoluis-coupled MCMC algorithm (MC3)

#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <assert.h>
#include <math.h>

#include "mc3.hh"
#include "inputs.hh"
#include "output.hh"

bool start_from_sim = true;   // This is set to true when doing speed testing.

using namespace std;

/// Initilaises the MC3 class
MC3::MC3(const Details &details, const Data &data, const Model &model, Inputs &inputs, const Output &output, const ObservationModel &obsmodel, Mpi &mpi) : state(details,data,model,obsmodel), mbp(INVT,details,data,model,obsmodel,output,mpi), details(details), data(data), model(model), output(output), obsmodel(obsmodel),mpi(mpi)
{	
	inputs.find_nrun(nrun);
	switch(details.mode){
		case MC3_INF:
			inputs.find_nchain(nchain,Ntot,N,nrun,mpi.ncore);
			inputs.find_invT_start_invT_final(invT_start,invT_final);
			break;
			
		case MCMC_MBP:
			nchain = 1; Ntot = nrun; if(Ntot%mpi.ncore != 0) emsgroot("'nrun' must be a multiple of the number of cores");
			N = Ntot/mpi.ncore;
			invT_start = UNSET;
			inputs.find_invT(invT_final);
			break;
			
		default: emsgEC("MC3",1); break;
	}
	inputs.find_nsample_or_ESSmin_or_cpu_time(nsample,ESSmin,cpu_time);
	inputs.find_nburnin(nburnin,nsample);
	inputs.find_nquench(nquench,nburnin);
	inputs.find_Tpower(Tpower);
	inputs.find_nthin(thin,nsample);
	percentage = UNSET;
}


/// Runs the inference algorithm
void MC3::run()
{
	initialise();
	
	ofstream trace[N]; for(auto ch = 0u; ch < N; ch++) output.trace_plot_inititialise(chain[ch].name,trace[ch]);

	timer[TIME_ALG].start();
	auto samp = 0u;
	do{                                                       // Sequentially goes through MCMC samples
		update_burnin(samp);                                    // Updates the burnin procedure
		        
		set_invT(samp);                                         // Sets the inverce temperatures for the chains
																					
		for(auto ch = 0u; ch < N; ch++){                        // MBP-MCMC updates	on chains		                      
			chain[ch].nproposal = mbp.mc3_mcmc_updates(part[ch],chain[ch].param_samp,chain[ch].invT,pup,paramprop[ch]);	
			
			store_param_samp(ch);                                 // Stores the parameter sample
		
			output.trace_plot(samp,part[ch].EF,part[ch].paramval,trace[ch]);  // Outputs trace plot
		}
		
		swap_states();                                          // Swaps between neighbouring chains

		if(burnin == false && samp%thin == 0) store_sample();   // Stores samples for plotting later
		samp++;
	}while(!terminate(samp));

	timer[TIME_ALG].stop();

	if(Ntot == 64) find_nchain_optimum();  // TO DO
	else{
		auto alg_time_sum = mpi.sum(timer[TIME_ALG].val);
		auto wait_time_sum = mpi.sum(timer[TIME_WAIT].val);
		
		if(mpi.core == 0) cout << invT_final;
		find_EF_range();
		if(mpi.core == 0){
			auto samp_fac = double(samp-nburnin)/samp;
			cout << "," << samp_fac*double(alg_time_sum-wait_time_sum)/(60.0*CLOCKS_PER_SEC) << " RESULT" << endl;
			cout << endl;
		}
	}
	
	model_evidence();                                         // Calculates the model evidence

	output.generate_graphs(part_plot);	                      // Draws the final pdf

	diagnostics();                                            // Outputs diagnostic information
}
 

/// Initialises each of the chains
void MC3::initialise()
{
	for(auto ch = 0u; ch < N; ch++){
		auto param = model.sample_from_prior();                 // Samples parameters from the prior

		if(start_from_sim){                                     // Sets parameters to simulated values TURNOFF
			for(auto th = 0u; th < model.param.size(); th++) param[th] = model.param[th].value;
		}
		
		state.simulate(param);                                  // Simulates a state

		//state.load(details.output_directory+"/sample",details.ndivision); // TEMP
	
		auto n = mpi.core*N+ch;
		auto run = n/nchain;
		auto num = n%nchain;
		part.push_back(state.create_particle(run));             // Converts the initial state into a particle  
		
		string name;
		if(num == 0){ // These are posterior chains
			name = "Trace"; if(nrun > 1) name += "_Run"+to_string(run+1);
		}
		else{ 
			name = "Other Chains/Trace"+to_string(num+1); if(nrun > 1) name += "_Run"+to_string(run+1);
		}
		
		chain.push_back(Chain(name,num));                       // Creates the different chains
		paramprop.push_back(ParamProp(details,data,model,output,mpi));
	}
	
	for(auto loop = 0u; loop < initialise_param_samp; loop++){// Sets up the initial MVN proposal distributions
		ParamSample ps; ps.run = UNSET; ps.EF = UNSET;
		ps.paramval = model.sample_from_prior();                // Samples parameters from the prior
		for(auto ch = 0u; ch < N; ch++){
			chain[ch].param_samp.push_back(ps);                   // Stores the parameter sample
		}
	}
}

	
/// Updates the burnin procedure
void MC3::update_burnin(const unsigned int samp)
{
	if(samp < nburnin){
		burnin = true;	
		if(samp < nburnin/2) pup = FAST_UPDATE; else pup = SLOW_UPDATE;
	}
	else{
		burnin = false; pup = NO_UPDATE;
	}
}
 

/// Swaps state between neighbouring chains 
void MC3::swap_states()
{
	if(nchain == 1) return; 
		
	timer[TIME_SWAP].start();
	vector <unsigned int> run, order, ntr, nac;
	vector <double> invT, EF; 
	
	for(auto ch = 0u; ch < N; ch++){
		run.push_back(part[ch].run);
		order.push_back(mpi.core*N+ch);
		invT.push_back(chain[ch].invT);
		EF.push_back(part[ch].EF);	
		ntr.push_back(chain[ch].ntr);
		nac.push_back(chain[ch].nac);
	}
	
	auto run_tot = mpi.gather(run);
	auto order_tot = mpi.gather(order);
	auto invT_tot = mpi.gather(invT);
	auto EF_tot = mpi.gather(EF);
	auto ntr_tot = mpi.gather(ntr);
	auto nac_tot = mpi.gather(nac);
	
	if(mpi.core == 0){
		for(auto loop = 0u; loop < 5; loop++){
			for(auto i = 0u; i < Ntot-1; i++){
				if(run_tot[i] == run_tot[i+1]){
					auto al = exp(0.5*(invT_tot[i]-invT_tot[i+1])*(EF_tot[i]-EF_tot[i+1]));
					ntr_tot[i]++;
					if(ran() < al){
						nac_tot[i]++;
						auto temp = order_tot[i]; order_tot[i] = order_tot[i+1]; order_tot[i+1] = temp;
						auto tempf = EF_tot[i]; EF_tot[i] = EF_tot[i+1]; EF_tot[i+1] = tempf;
					}
				}
			}
		}
		
		if(false){
			for(auto i = 0u; i < Ntot; i++){
				cout << i << " " << invT_tot[i] << " " << order_tot[i] << " " << EF_tot[i] << " " << run_tot[i] << " ord" << endl;
			}
		}
	}
		
	ntr = mpi.scatter(ntr_tot);
	nac = mpi.scatter(nac_tot);
	
	for(auto ch = 0u; ch < N; ch++){
		chain[ch].ntr = ntr[ch];
		chain[ch].nac = nac[ch];
	}
	
	mpi.copy_particles(part,order_tot,N,Ntot);                              
	
	if(checkon == true){
		EF = mpi.scatter(EF_tot);
		for(auto ch = 0u; ch < N; ch++){
			if(EF[ch] != part[ch].EF) emsgEC("MC3",2);
		}
	}

	timer[TIME_SWAP].start();
}


/// Stores particle samples to be plotted later
void MC3::store_sample()
{
	for(auto ch = 0u; ch < N; ch++){
		if(chain[ch].num == 0) part_plot.push_back(part[ch]);
	}
}


/// Stores paramer samples to a) help tune parameter proposal kernels and b) to estimate the model evidence later 
void MC3::store_param_samp(unsigned int ch)
{
	if(burnin == true) chain[ch].param_samp.push_back(part[ch].create_param_samp());
	else chain[ch].EF_samp.push_back(part[ch].EF);
}		


/// Calculates an estimate for the model evidence
void MC3::model_evidence() const
{
	if(details.mode == MCMC_MBP){
		if(mpi.core == 0) cout << "The 'mcmcmbp' inference algorithm does not provide an estimate for model evidence." << endl;
		return;
	}
		
	vector <double> invT_total;

	auto EF_chain_sample = mpi.gather_EF_chain_sample(chain,part,N,nchain,nrun,invT_total);
	
	if(mpi.core == 0){		
		if(invT_total[nchain-1] != 0){
			cout << "The 'mc3' inference algorithm can only provide an estimate for model evidence if the hottest chain has an inverse temperature of zero." << endl;
			return;
		}
		
		if(false){
			for(auto ch = 0u; ch < nchain; ch++) cout << ch << " " <<  invT_total[ch] << " invT" << endl;

			for(auto run = 0u; run < nrun; run++){
				for(auto ch = 0u; ch < nchain; ch++){
					cout << run << " " << ch << " " << EF_chain_sample[run][ch].size() << " ' Samples" << endl;
				}
			}
		}

		vector <double> ME_list;
		for(auto run = 0u; run < nrun; run++){
			auto ME = 0.0;
			for(auto ch = nchain-1; ch > 0; ch--){
				auto DinvT = invT_total[ch-1] - invT_total[ch];
				auto mu = 0.0, n = 0.0;
				for(auto EF : EF_chain_sample[run][ch]){ mu += EF; n++;}
				mu /= n;
			
				auto sum = 0.0;
				for(auto EF : EF_chain_sample[run][ch]) sum += exp(-0.5*DinvT*(EF-mu));
				sum /= n;
				
				ME += -0.5*DinvT*mu + log(sum);
			}
			ME_list.push_back(ME);
		}

		output.final_model_evidence(ME_list,invT_final,UNSET);
	}
}


/// Outputs diagnostic information
void MC3::diagnostics()
{	
	for(auto ch = 0u; ch < N; ch++){
		auto chain_name = chain[ch].name+"_";
		if(details.mode == MCMC_MBP) chain_name = "";
		
		string filefull = details.output_directory+"/Diagnostics/" + chain_name+"MCMC_proposals.txt";
		ofstream dia(filefull);
		if(!dia) emsg("Cannot open the file '"+filefull+"'");
		
		paramprop[ch].set_ac_rate();
		dia << paramprop[ch].print_proposal_information(false);
	}	
	
	if(nchain > 1){
		vector <unsigned int> ntr, nac;
		vector <double> invT;
		for(auto ch = 0u; ch < N; ch++){
			ntr.push_back(chain[ch].ntr);
			nac.push_back(chain[ch].nac);
			invT.push_back(chain[ch].invT);
		}
		auto ntr_tot = mpi.gather(ntr);
		auto nac_tot = mpi.gather(nac);
		auto invT_tot = mpi.gather(invT);
			
		if(mpi.core == 0){
			string filefull = details.output_directory+"/Diagnostics/Chain_Swap.txt";
			ofstream dia(filefull);
			if(!dia) emsg("Cannot open the file '"+filefull+"'");
		
			dia << "This provides information about the rate at which adjacent chains swap state:" << endl;
			dia << endl;
			
			auto rate_min = LARGE, rate_av=0.0;
			for(auto ru = 0u; ru < nrun; ru++){
				if(nrun > 1) dia << "Run " << ru+1 << ":" << endl; 
				for(auto ch = 0u; ch < nchain-1; ch++){
					auto i = ru*nchain+ch;
					auto rate = double(nac_tot[i])/ntr_tot[i];
					rate_av += rate;
					if(rate < rate_min) rate_min = rate;
					dia << "Chain " << ch+1 << " <-> Chain " << ch+2 << "   " << per(rate) << endl;
				}
				dia << endl;
			}
			
			stringstream ss;
			ss << "The average acceptance rate is " << per(rate_av/(Ntot-1)) << " with a minimum of " << per(rate_min) << "." << endl;
			if(rate_min < 0.1) ss << "Consider making 'nchain' higher to increase this rate." << endl;
			else{
				if(rate_min > 0.5) ss << "Consider making 'nchain' smaller to improve algorithm speed." << endl;
				else ss << "This implies that 'nchain' is at a good level." << endl;
			}
			ss << endl;
			dia << ss.str();
			
			dia << "Here we give the inverse temperatures for the various chains:" << endl;
			for(auto ch = 0u; ch < nchain; ch++){
				dia << "Chain " << ch+1 << "   Inverse temperature: " << invT_tot[ch];
				if(ch == 0) dia << "             <<<< Posterior >>>>";
				if(ch == nchain-1) dia << "             <<<< Prior >>>>";
				dia << endl;
			}
		
			cout << "Chain swapping diagnostics can be found in 'Diagnostics/Chain_Swap.txt'" << endl;
			cout << ss.str();
		}
	}
}


/// Sets the inverse temperature of the chain
void MC3::set_invT(const unsigned int samp)
{
	auto pmax = 1-pow(invT_start,1.0/Tpower);
	auto pmin = 1-pow(invT_final,1.0/Tpower);
	
	for(auto ch = 0u; ch < N; ch++){
		auto num = chain[ch].num;
		
		auto fac = 1.0; if(samp < nquench && start_from_sim == false) fac = double(samp)/nquench;
	
		if(nchain == 1){
			chain[ch].invT = invT_final*pow(fac,Tpower);
		}
		else{
			auto kappa = double(num)/(nchain-1);
			auto ppf = pmin+kappa*(pmax-pmin);
			auto ppeff = 1-fac*(1-ppf);	
			chain[ch].invT = pow(1-ppeff,Tpower);
		}
	}
	
	if(false){
		for(auto ch = 0u; ch < N; ch++){
			cout << mpi.core << " "  << ch <<  "  num " << chain[ch].num;
			cout << "   invT" << chain[ch].invT << "  run" << part[ch].run << endl; 
		}
	}
}


/// Determines when to terminate the algorithm
bool MC3::terminate(const unsigned int samp)
{
	bool term = false;
	if(ESSmin != UNSET){
		if(samp%100 == 0){
			if(samp <= nburnin){
				if(mpi.core == 0) cout << "Burnin samples: " <<  samp << endl;	
			}
			else{
				auto psamp_tot = mpi.gather_psamp(part_plot);
				if(mpi.core == 0){
					auto ESS = output.get_effective_sample_size(psamp_tot);
					if(vec_min(ESS) > ESSmin && samp-nburnin >= 20) term = true; 
					cout << "Number of samples: " << samp-burnin;
					cout << "    Smallest ESS value:" << vec_min(ESS) << "     ESSmin: " << ESSmin << endl;
				}
			}
		}
	}
	
	if(cpu_time != UNSET){
		auto time_av = mpi.average((clock() - details.time_start)/(60.0*CLOCKS_PER_SEC));
		if(time_av > cpu_time){
			term = true;
			if(mpi.core == 0)	cout << "Maximum execution time reached." << endl << flush;
		}
	}
	
	if(nsample != UNSET){
		if(mpi.core == 0){
			if(samp == nsample) term = true;
			
			output.print_percentage(samp,nsample,percentage);       // Prints the percentage to show progress
		}
	}

	mpi.bcast(term);
	
	return term;
}


/// Finds the range in EF for the posterior chain
void MC3::find_EF_range()
{
	vector <double> invT_total;
	auto EF_chain_sample = mpi.gather_EF_chain_sample(chain,part,N,nchain,nrun,invT_total);
	
	if(mpi.core == 0){
		auto stat = output.get_statistic_75_percent(EF_chain_sample[0][0]);
		cout << "," << stat.CImax;
		stat = output.get_statistic(EF_chain_sample[0][0]);
		cout << "," << stat.CImax;
	}
}

/// Finds the optimum number of chains for a given inverse temperature
void MC3::find_nchain_optimum() const 
{
	vector <double> invT_total;

	auto EF_chain_sample = mpi.gather_EF_chain_sample(chain,part,N,nchain,nrun,invT_total);
	
	if(mpi.core == 0){		
		vector <Statistics> stat(nchain);
		for(auto ch = 0u; ch < nchain; ch++){
			stat[ch] = output.get_statistic_75_percent(EF_chain_sample[0][ch]);				
			//cout << ch << " " <<  invT_total[ch] << " " << stat.CImin << " " << stat.CImax << " invT" << endl;
		}
		
		for(auto ch = 0u; ch < nchain; ch++){
			cout << ch << " " <<  invT_total[ch] << " " << stat[ch].CImin << " " << stat[ch].CImax << " invT original" << endl;
		}
		
		for(auto ch = 0u; ch < nchain; ch+=2){
			auto kappa = double(ch)/(nchain-1);
			auto invT_fin = pow(1-kappa,Tpower)*invT_final;
		
			vector <Statistics> stat_new;
			
			unsigned int nch; 
			for(nch = 2; nch < 100; nch++){
				stat_new.resize(nch);
				for(auto c = 0u; c < nch; c++){
					auto kap = double(c)/(nch-1);
					auto invT = pow(1-kap,Tpower)*invT_fin;
					
					unsigned int chh;
					for(chh = 0; chh < nchain-1; chh++){
						auto kap1 = double(chh)/(nchain-1);
						auto invT1 =  pow(1-kap1,Tpower)*invT_final;
					
						auto kap2 = double(chh+1)/(nchain-1);
						auto invT2 = pow(1-kap2,Tpower)*invT_final;
						
						if(invT1 >= invT && invT2 <= invT){
							auto f = (invT1-invT)/(invT1-invT2);
							stat_new[c].CImin = stat[chh].CImin*(1-f) + stat[chh+1].CImin*f;
							stat_new[c].CImax = stat[chh].CImax*(1-f) + stat[chh+1].CImax*f;
							break;
						}
					}
					if(chh == nchain-1) emsgEC("Mc3",22);
				}					
				
				auto jmax = (unsigned int)(0.8*nch);
				auto j = 0u; while(j < jmax && stat_new[j].CImax > stat_new[j+1].CImin) j++;
				
				if(j == jmax) break;
			}
		
			cout << ch << " " << invT_fin << " " << nch << " Answer\n";
			/*
			for(auto c = 0u; c < nch; c++){
				auto kap = double(c)/(nch-1);
				auto invT = pow(1-kap,Tpower)*invT_fin;
				cout << c << " " <<  invT << " " << stat_new[c].CImin << " " << stat_new[c].CImax << " CI answer" << endl;
			}
			*/
			cout << "\n";
		}
	}
}


/// Initialises a chain
Chain::Chain(string name_, unsigned int num_)
{
	name = name_; num = num_; ntr = 0; nac = 0;
}
