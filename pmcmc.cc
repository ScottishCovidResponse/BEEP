/// Implements the PMCMC algorithm

#include <iostream>
#include <fstream>   
#include <cmath>

using namespace std;

#include "timers.hh"
#include "inputs.hh"
#include "state.hh"
#include "pmcmc.hh"
#include "mpi.hh"
#include "output.hh"

PMCMC::PMCMC(const Details &details, const Data &data, const Model &model, const AreaTree &areatree, const Inputs &inputs, Output &output, const ObservationModel &obsmodel, Mpi &mpi) : paramprop(details,data,model,areatree,output,mpi), details(details), data(data), model(model), areatree(areatree), output(output), obsmodel(obsmodel), mpi(mpi)
{
	inputs.find_nparticle(Ntot,N,mpi.ncore);
	inputs.find_nsample(nsample);
	inputs.find_nburnin(nburnin,nsample);
	inputs.find_nthin(thin,nsample);
	initialise_variables();
	percentage = UNSET;
}


/// Runs the PMCMC algorithm
void PMCMC::run()
{
	if(core == 0) output.trace_plot_inititialise("trace",trace);
	
	initialise();
	
	timer[TIME_ALG].start();
	for(auto samp = 0u; samp < nsample; samp++){                             // Sequentially goes through MCMC samples
		update_burnin(samp);                                                   // Updates the burnin procedure
		
		if(core == 0) output.trace_plot(samp,Li,Pi.paramval,trace);            // Outputs parameters to trace plot
		
		if(samp%10 == 0) get_proposals();                                      // Gets a list of proposals used as an "update"

		mcmc_updates();                                                        // Performs the MCMC updates

		if(core == 0 && burnin == true) param_samp.push_back(Pi.paramval);     // Samples adaptively improving proposals
		
		if(core == 0 && burnin == false && samp%thin == 0) particle_store.push_back(Pi); // Stores samples for plotting later
	
		output.print_percentage(samp,nsample,percentage);                      // Prints the percentage to show progress
	}
	timer[TIME_ALG].stop();
	
	output.generate_graphs(particle_store);		                               // Outputs the results
	
	if(core == 0) paramprop.set_ac_rate();                                   // Outputs diagnostic information regarding proposals
	paramprop.diagnostics();
}


/// Updates the burnin procedure
void PMCMC::update_burnin(unsigned int samp)
{
	if(samp < nburnin){
		burnin = true; if(samp < 200) pup = FAST_UPDATE; else pup = SLOW_UPDATE;
	}
	else{
		burnin = false; pup = NO_UPDATE;
	}
}

		
/// Initialises quantities in the class
void PMCMC::initialise()
{
	for(auto p = 0u; p < N; p++) particle.push_back(State(details,data,model,obsmodel));
	
	auto ninit_samp = 100u;
	Li = -LARGE;
	for(auto i = 0u; i < ninit_samp; i++){                             // Randomly samples parameters / states and picks the best
		if(core == 0){
			Pp.paramval = model.sample_from_prior();                       // Randomly samples from the prior
			if(true){                                                      // Sets parameters to those used in the simulation
				for(auto th = 0u; th < model.param.size(); th++){
					Pi.paramval[th] = model.param[th].value;
				}
			}
		}
		
		Lp = obs_prob(Pp.paramval);                                      // Calculates unbiased estimate of likelihood
		mpi.bcast(Lp);
		
		if(Lp > Li){
			Li = Lp;
			Pi = mpi.particle_sample(backpart,particle,obsmodel);
		}
	}
		
	if(core == 0){
		// This generates parameter samples near to the intial set (for the initial normal and MVN proposal distributions)
		for(auto loop = 0u; loop < 100; loop++){
			auto par = Pi.paramval;		
			
			auto par_samp = model.sample_from_prior();

			for(auto th = 0u; th < par.size(); th++) par[th] = 0.9*par[th] + 0.1*par_samp[th];
			param_samp.push_back(par);
		}
	}
	
	if(core == 0) Pri = model.prior(Pi.paramval);
}
	

/// Uses particle filtering to get an unbiased estimate for the observation probability for a given parameter set 
double PMCMC::obs_prob(vector <double> &paramv)
{
	timer[TIME_PMCMCLIKE].start();
	
	mpi.bcast(paramv);

	for(auto p = 0u; p < N; p++) particle[p].set_param(paramv);

	vector <double> L(N);
	
	auto obprob = 0.0; 
	for(auto sec = 0u; sec < obsmodel.nsection; sec++){
		if(sec > 0) mpi.end_sec_swap(particle,obsmodel.section_ti[sec],backpart[sec],buffersize);

		for(auto p = 0u; p < N; p++){ 	
			particle[p].simulate(obsmodel.section_ti[sec],obsmodel.section_tf[sec]);
		
			L[p] = obsmodel.calculate_section(&particle[p],sec);
		}
	
		obprob += bootstrap(sec,L);
		
		if(std::isnan(obprob) || std::isinf(obprob)) emsgEC("PMCMC",6);
	}
	
	timer[TIME_PMCMCLIKE].stop();
	
	return obprob;
}


/// Performs the bootstrap step which randomly selects particles based on their observation probability
double PMCMC::bootstrap(unsigned int sec, vector <double> &L)
{
	timer[TIME_BOOTSTRAP].start();
	
	vector <double> w(Ntot), wsum(Ntot);
	
	auto Ltot = mpi.gather(L);

	auto obprob = 0.0;
	auto sum = 0.0;
	if(core == 0){	
		double Lmax = -LARGE; 
		for(auto p = 0u; p < Ntot; p++){
			if(Ltot[p] > Lmax) Lmax = Ltot[p];
		}
	
		for(auto p = 0u; p < Ntot; p++){
			w[p] = exp(Ltot[p]-Lmax);
			sum += w[p];
			wsum[p] = sum;
		}
	
		obprob = Lmax + log(sum/Ntot);
	}
	
	auto &bp = backpart[sec+1];
		
	if(core == 0){
		for(auto p = 0u; p < Ntot; p++) bp[p] = UNSET;
		 
		vector <unsigned int> extra; 
	
		for(auto p = 0u; p < Ntot; p++){
			auto z = ran()*sum;
			auto pp = 0u; while(pp < Ntot && z > wsum[pp]) pp++;
			if(pp == Ntot) emsg("Zero probability");
			
			if(bp[pp] == UNSET) bp[pp] = pp;
			else extra.push_back(pp);	
		}
		
		for(auto p = 0u; p < Ntot; p++){  
			if(bp[p] == UNSET){
				auto pp = extra[extra.size()-1]; extra.pop_back();
				bp[p] = pp;
			}
		}
		if(extra.size() != 0) emsgEC("PMCMC",17);
	}

	mpi.bcast(bp);
	
	timer[TIME_BOOTSTRAP].stop();
	
	return obprob;	
}


/// Gets the proposals used for the next itermation of MCMC
void PMCMC::get_proposals()
{
	if(core == 0) prop_list = paramprop.get_proposal_list(param_samp); 
	mpi.bcast(prop_list);

	if(core == 0 && diagnotic_output == true) cout << "# Proposals " << prop_list.size() << endl;
}


/// Updates parameters using MCMC proposals
void PMCMC::mcmc_updates()
{
	timer[TIME_MCMCPROP].start();
	for(auto i = 0u; i < prop_list.size(); i++){
		auto& prop = prop_list[i];
		auto num = prop.num;
		
		switch(prop.type){
			case SELF_PROP:
				self_proposal(paramprop.self);    
				break;

			case MVN_PROP:
				{
					switch(paramprop.mvn[num].type){
						case SIGMA_PARAM: sigma_reff_proposal(paramprop.mvn[num]); break;
						default: mvn_proposal(paramprop.mvn[num]); break;
					}
				}
				break;
				
			case MEAN_TIME_PROP:
				mean_time_proposal(paramprop.mean_time[num]);    
				break;
				
			case NEIGHBOUR_PROP:
				neighbour_proposal(paramprop.neighbour[num]);    
				break;
				
			case JOINT_PROP:
				joint_proposal(paramprop.joint[num]);    
				break;
				
			case FIXEDTREE_PROP: case SLICETIME_PROP: emsgEC("PMCMC",32); break;
		}
		
	}
	timer[TIME_MCMCPROP].stop();
}


/// Returns the acceptance probability
double PMCMC::get_al()
{
	mpi.bcast(Pp.paramval);

	Lp = obs_prob(Pp.paramval);
	Pp = mpi.particle_sample(backpart,particle,obsmodel);
	if(core == 0) Prp = model.prior(Pp.paramval);
	
	double al = exp(invT*(Lp-Li) + Prp-Pri);
	if(core == 0 && false) cout << al << " " << invT << " " << Lp <<  " " << Li <<  " " << Prp << " " << Pri << " pr\n";
	
	if(std::isnan(al)) emsgEC("PMCMC",45);
				
	return al;
}


/// Copies the proposed state into the initial state
void PMCMC::copy_propose_state()
{
	Li = Lp;
	Pi = Pp;
	Pri = Prp;
}


/// Does a proposal without changing the parameter set
void PMCMC::self_proposal(Self &self)  
{
	timer[TIME_SELF].start();
	
	if(core == 0) Pp.paramval = Pi.paramval;
	
	auto al = get_al(); 
	if(core == 0){
		if(true) cout << al << " " << invT << " al self\n";
		ParamUpdate pu = pup; if(invT_dynamic == false) pu = NO_UPDATE;
		if(self.MH(al,invT,pu) == SUCCESS) copy_propose_state();
	}
	
	timer[TIME_SELF].stop();
}


/// Does a MVN proposals shifted by Langevin term to account for prior
void PMCMC::mvn_proposal(MVN &mvn)  
{
	timer[TIME_MVN].start();

	double probif;
	auto suc = 0u; if(core == 0){ if(mvn.propose_langevin(Pp.paramval,Pi.paramval,probif,model) == SUCCESS) suc = 1;}
	mpi.bcast(suc);
	
	if(suc == 1){
		auto al = get_al();
		if(core == 0){
			al *= exp(mvn.get_probfi(Pi.paramval,Pp.paramval,model) - probif);
			if(mvn.MH(al,pup) == SUCCESS) copy_propose_state();
		}
	}
	
	timer[TIME_MVN].stop();
}


/// Simulatenously changes the mean time of one transition and minus for the subsequent transitions 
void PMCMC::sigma_reff_proposal(MVN &mvn)
{
	timer[TIME_SIGMA].start();
	
	auto al = 0.0;
	auto suc = 0u; if(core == 0){ if(mvn.sigma_propose(Pp.paramval,Pi.paramval,model) == SUCCESS) suc = 1;}
	mpi.bcast(suc);
	
	if(suc == 1){
		al = get_al();
		if(core == 0){
			al *= exp(Pri - Prp);
			if(mvn.MH(al,pup) == SUCCESS) copy_propose_state();
		}
	}
	
	timer[TIME_SIGMA].stop();
}


/// Does joint proposals on splines
void PMCMC::joint_proposal(Joint &jn)  
{
	timer[TIME_JOINT].start();
	
	auto al = 0.0;
	auto suc = 0u; if(core == 0){ if(jn.propose(Pp.paramval,Pi.paramval,model) == SUCCESS) suc = 1;}
	mpi.bcast(suc);
	
	if(suc == 1){
		al = get_al();
		if(core == 0){
			if(jn.MH(al,pup) == SUCCESS) copy_propose_state();
		}
	}
	
	timer[TIME_JOINT].stop();
}


/// Changes a point on a spline and minus that change on the next point
void PMCMC::neighbour_proposal(Neighbour &nei)  
{
	timer[TIME_NEIGHBOUR].start();
	
	auto al = 0.0;
	auto suc = 0u; if(core == 0){ if(nei.propose(Pp.paramval,Pi.paramval,model) == SUCCESS) suc = 1;}
	mpi.bcast(suc);
	
	if(suc == 1){
		al = get_al();
		if(core == 0){
			if(nei.MH(al,pup) == SUCCESS) copy_propose_state();
		}
	}
	
	timer[TIME_NEIGHBOUR].stop();
}


/// Simulatenously changes the mean time of one transition and minus for the subsequent transitions 
void PMCMC::mean_time_proposal(MeanTime& mt)
{
	timer[TIME_MEANTIME].start();
	
	auto al = 0.0;
	auto suc = 0u; if(core == 0){ if(mt.propose(Pp.paramval,Pi.paramval,model) == SUCCESS) suc = 1;}
	mpi.bcast(suc);
	
	if(suc == 1){
		al = get_al();
		if(core == 0){
			if(mt.MH(al,pup) == SUCCESS) copy_propose_state();
		}
	}
	
	timer[TIME_MEANTIME].stop();
}


/// Initialises quantities before PMCMC can start
void PMCMC::initialise_variables()
{
	backpart.resize(obsmodel.nsection+1); for(auto sec = 0u; sec <= obsmodel.nsection; sec++) backpart[sec].resize(Ntot);
	
	buffersize = 1 + data.narea*(1+ model.comp.size()*(1+data.ndemocatpos));       // This stores population sizes
	buffersize += 1 + model.ninfection_trans*(1 + data.narage);                    // This stores Imap
	buffersize += 1 + model.ninfection_trans*(1 + data.narage);                    // This stores Idiag
	
	invT = 1; invT_dynamic = true;
	
	Pi.paramval.resize(model.param.size());
	Pp.paramval.resize(model.param.size());
	
	core = mpi.core;	
}
