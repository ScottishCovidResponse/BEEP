/// Implement the metropoluis-coupled MCMC algorithm (MC3)

#include <iostream>
#include <fstream>
#include <algorithm>
#include <assert.h>
#include <math.h>

#include "mc3.hh"
#include "inputs.hh"
#include "output.hh"

using namespace std;

/// Initilaises the MC3 class
MC3::MC3(const Details &details, const Data &data, const Model &model, const AreaTree &areatree, const Inputs &inputs, const Output &output, const ObservationModel &obsmodel, Mpi &mpi) : state(details,data,model,obsmodel), mbp(INVT,details,data,model,areatree,obsmodel,output,mpi), details(details), data(data), model(model), areatree(areatree), output(output), obsmodel(obsmodel),mpi(mpi)
{	
	inputs.find_nchain(Ntot,N,mpi.ncore);
	inputs.find_nsample(nsample);
	inputs.find_nburnin(nburnin,nsample);
	inputs.find_nquench(nquench,nburnin);
	inputs.find_invT_start_invT_final(invT_start,invT_final);
	inputs.find_nthin(thin,nsample);
	percentage = UNSET;
}


/// Runs the inference algorithm
void MC3::run()
{
	initialise();
	
	ofstream trace[N];
	for(auto ch = 0u; ch < N; ch++){
		string name = "Chain_"+to_string(mpi.core*N+ch);
		output.trace_plot_inititialise(name,trace[ch]);
	}
	
	timer[TIME_ALG].start();
	for(auto samp = 0u; samp < nsample; samp++){              // Sequentially goes through MCMC samples
		update_burnin(samp);                                    // Updates the burnin procedure
		        
		for(auto ch = 0u; ch < N; ch++){                        // Goes through all the chains
			set_invT(samp);
																														// Performs MBP-MCMC updates	
			chain[ch].nproposal = mbp.mc3_mcmc_updates(part[ch],chain[ch].param_samp,chain[ch].invT,pup,chain[ch].paramprop);	
			
			if(burnin == true) chain[ch].param_samp.push_back(part[ch].paramval); // Stores the parameter sample
			
			swap_states();                                        // Swaps between neighbouring chains
			
			if(burnin == false && samp%thin == 0) store_sample(); // Stores samples for plotting later
			
			output.trace_plot(samp,part[ch].EF,part[ch].paramval,trace[ch]);  // Outputs trace plot
		}
		
		output.print_percentage(samp,nsample,percentage);       // Prints the percentage to show progress
	}
	timer[TIME_ALG].stop();

	output.generate_graphs(part_plot);	                      // Draws the final pdf
	diagnostics();                                            // Outputs diagnostic information
}
 

/// Initialises each of the chains
void MC3::initialise()
{
	for(auto ch = 0u; ch < N; ch++){
		auto param = model.sample_from_prior();                 // Samples parameters from the prior

		if(false){                                              // Sets parameters to simulated values
			for(auto th = 0u; th < model.param.size(); th++) param[th] = model.param[th].value;
		}
				
		state.simulate(param);                                  // Simulates a state
	
		part.push_back(state.create_particle());                // Converts the initial state into a particle  
		
		chain.push_back(Chain (N*mpi.core+ch,details,data,model,areatree,output,mpi)); // Creates the different chains
	}

	for(auto loop = 0u; loop < 100; loop++){                  // Sets up the initial MVN proposal distributions
		auto param = model.sample_from_prior();                 // Samples parameters from the prior
		for(auto ch = 0u; ch < N; ch++){
			chain[ch].param_samp.push_back(param);                // Stores the parameter sample
		}
	}
}

	
/// Updates the burnin procedure
void MC3::update_burnin(unsigned int samp)
{
	if(samp < nburnin){
		burnin = true; if(samp < 200) pup = FAST_UPDATE; else pup = SLOW_UPDATE;
	}
	else{
		burnin = false; pup = NO_UPDATE;
	}
}
 

/// Swaps state between neighbouring chains 
void MC3::swap_states()
{
	timer[TIME_SWAP].start();
	vector <unsigned int> order, ntr, nac;
	vector <double> invT, EF; 
	
	for(auto ch = 0u; ch < N; ch++){
		order.push_back(chain[ch].num);
		invT.push_back(chain[ch].invT);
		EF.push_back(part[ch].EF);	
		ntr.push_back(chain[ch].ntr);
		nac.push_back(chain[ch].nac);
	}
	
	auto order_tot = mpi.gather(order);
	auto invT_tot = mpi.gather(invT);
	auto EF_tot = mpi.gather(EF);
	auto ntr_tot = mpi.gather(ntr);
	auto nac_tot = mpi.gather(nac);
	
	if(mpi.core == 0){
		for(auto loop = 0u; loop < 5; loop++){
			for(auto i = 0u; i < Ntot-1; i++){
				auto al = exp(0.5*(invT_tot[i]-invT_tot[i+1])*(EF_tot[i]-EF_tot[i+1]));
				ntr_tot[i]++;
				if(ran() < al){
					nac_tot[i]++;
					auto temp = order_tot[i]; order_tot[i] = order_tot[i+1]; order_tot[i+1] = temp;
					auto tempf = EF_tot[i]; EF_tot[i] = EF_tot[i+1]; EF_tot[i+1] = tempf;
				}
			}
		}
		
		if(false){
			for(auto i = 0u; i < Ntot; i++){
				cout << i << " " << invT_tot[i] << " " << order_tot[i] << " " << EF_tot[i] << " ord\n";
			}
			emsg("DO");
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
			if(EF[ch] != part[ch].EF) emsgEC("MC3",53);
		}
	}
	timer[TIME_SWAP].start();
}


/// Stores a particle sample to be plotted later
void MC3::store_sample()
{
	if(chain[0].num == 0) part_plot.push_back(part[0]);
}


/// Outputs diagnostic information
void MC3::diagnostics()
{	
	for(auto ch = 0u; ch < N; ch++){
		string filefull = details.output_directory+"/Diagnostics/Chain"+to_string(chain[ch].num)+"_MCMC_proposals.txt";
		ofstream dia(filefull);
		if(!dia) emsg("Cannot open the file '"+filefull+"'");
		
		dia << chain[ch].paramprop.print_proposal_information(false);
	}	
	
	vector <unsigned int> ntr, nac;
	for(auto ch = 0u; ch < N; ch++){
		ntr.push_back(chain[ch].ntr);
		nac.push_back(chain[ch].nac);
	}
	auto ntr_tot = mpi.gather(ntr);
	auto nac_tot = mpi.gather(nac);
	
	if(mpi.core == 0){
		string filefull = details.output_directory+"/Diagnostics/Chain_swap.txt";
		ofstream dia(filefull);
		if(!dia) emsg("Cannot open the file '"+filefull+"'");
	
		dia << "This provides information about the rate at which adjacent chains swap state:" << endl;
	
		auto rate_min = LARGE, rate_av=0.0;
		for(auto ch = 0u; ch < Ntot-1; ch++){
			auto rate = double(nac_tot[ch])/ntr_tot[ch];
			rate_av += rate;
			if(rate < rate_min) rate_min = rate;
			dia << "Chain " << ch << " <-> Chain " << ch+1 << "   " << per(rate) << endl;
		}
		
		cout << "The average acceptance rate is " << per(rate_av/(Ntot-1)) << " with a minimum of " << per(rate_min) << ". ";
		if(rate_min < 0.1) cout << "Consider making 'nchain' higher to increase this rate." << endl;
		else cout << "This implies that 'nchain' is at an acceptable level." << endl;
	}
}


/// Initialises a chain
Chain::Chain(unsigned int num_, const Details &details, const Data &data, const Model &model, const AreaTree &areatree, const Output &output, Mpi &mpi) : paramprop(details,data,model,areatree,output,mpi)
{
	num = num_; ntr = 0; nac = 0;
}


/// Sets the inverse temperature of the chain
void MC3::set_invT(unsigned int samp)
{
	auto pmax = 1-pow(invT_start,1.0/Tpower);
	auto pmin = 1-pow(invT_final,1.0/Tpower);
	
	for(auto ch = 0u; ch < N; ch++){
		auto num = chain[ch].num;
		
		if(Ntot == 1) chain[ch].invT = invT_final;
		else{
			auto fac=1.0;
			if(samp < nquench) fac = double(samp)/nquench;
		
			auto kappa = double(num)/(Ntot-1);
			auto ppf = pmin+kappa*(pmax-pmin);
			auto ppeff = 1-fac*(1-ppf);	
			chain[ch].invT = pow(1-ppeff,Tpower);
		}
	}
}
