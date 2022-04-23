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

PMCMC::PMCMC(const Details &details, const Data &data, const Model &model, Inputs &inputs, Output &output, const ObservationModel &obsmodel, Mpi &mpi) : paramprop(details,data,model,output,mpi), details(details), data(data), model(model), output(output), obsmodel(obsmodel), mpi(mpi)
{
	inputs.find_nrun(nrun);
	inputs.find_nparticle_pmcmc(Ntot,N,mpi.ncore);
	inputs.find_nsample_or_ESSmin_or_cpu_time(nsample,ESSmin,cpu_time);
	inputs.find_nburnin(nburnin,nsample);
	if(nburnin < pmcmc_init_samp) emsg("'nburnin' must be at least "+to_string(2*pmcmc_init_samp));
	inputs.find_nthin(thin,nsample);
	invT = inputs.find_double("invT",UNSET); 
	inputs.find_sd(sd_init);
	initialise_variables();
	
	percentage = UNSET;
}


/// Runs the PMCMC algorithm
void PMCMC::run()
{
	initialise();

	ofstream trace[nrun];
	if(core == 0){                                                           // Sets files for writing trace plots
		for(auto ru = 0u; ru < nrun; ru++){
			string name = "Trace"; if(nrun > 1) name += "_Run"+to_string(ru+1);
			output.trace_plot_inititialise(name,trace[ru]);
		}
	}

	if(pmcmc_start_param == true && true) calculate_start_invT();
	
	unsigned int get_prop_period = nburnin/20; if(get_prop_period == 0) get_prop_period = 1;
	
	timer[TIME_ALG].start();
	auto samp = 0u;
	do{     
		// Sequentially goes through MCMC samples
		
		update_burnin(samp);                                                   // Updates the burnin procedure
		
		if(core == 0){                                                         // Outputs parameters to trace plot
			for(auto ru = 0u; ru < nrun; ru++) output.trace_plot(samp,Li_run[ru],Pi_run[ru].paramval,trace[ru]);   
		}
		
		if(burnin == true && samp%get_prop_period == 0) get_proposals();       // Gets a list of proposals used as an "update"

		for(auto ru = 0u; ru < nrun; ru++){
			Li = Li_run[ru]; Pi = Pi_run[ru]; Pri = Pri_run[ru];                 // Loads stored values for run
			
			mcmc_updates();                                                      // Performs the MCMC updates

			if(core == 0 && burnin == true){                                     // Samples adaptively improving proposals
				ParamSample ps; ps.paramval = Pi.paramval; ps.run = UNSET; ps.EF = Li;
				param_samp.push_back(ps);
			}
		
			if(core == 0 && burnin == false && samp%thin == 0){                  // Stores samples for plotting later
				Pi.run = ru;
				particle_store.push_back(Pi);
			}
			
			Li_run[ru] = Li; Pi_run[ru] = Pi; Pri_run[ru] = Pri;                 // Saves stored values for run
		}
		
		if(samp%100000 == 0 && samp > 2*nburnin){
			output.generate_graphs(particle_store);
			if(core == 0) get_EF_dist();
			auto alg_time_sum = mpi.sum(timer[TIME_ALG].val+clock());
			auto wait_time_sum = mpi.sum(timer[TIME_WAIT].val);
			cout <<  "Algorithm time: " << double(alg_time_sum)/(60.0*CLOCKS_PER_SEC) << " minutes." << endl;
			cout << "Mpi wait time: "  <<double(wait_time_sum)/(60.0*CLOCKS_PER_SEC) << " minutes." << endl;
			auto samp_fac = double(samp-nburnin)/samp;
			cout <<  "Algorithm-wait time: " << samp_fac*double(alg_time_sum-wait_time_sum)/(60.0*CLOCKS_PER_SEC) << " minutes." << endl;
		}
		
		samp++;
	}while(!terminate(samp));
	timer[TIME_ALG].stop();
	
	output.generate_graphs(particle_store);		                               // Outputs the results
	
	paramprop.set_ac_rate();                                                 // Outputs diagnostic information about proposals
	paramprop.diagnostics();
	
	auto alg_time_sum = mpi.sum(timer[TIME_ALG].val);
	auto wait_time_sum = mpi.sum(timer[TIME_WAIT].val);
	if(core == 0){
		cout <<  "Algorithm time: " << double(alg_time_sum)/(60.0*CLOCKS_PER_SEC) << " minutes." << endl;
		cout << "Mpi wait time: "  <<double(wait_time_sum)/(60.0*CLOCKS_PER_SEC) << " minutes." << endl;
		auto samp_fac = double(samp-nburnin)/samp;
		cout << samp_fac << " samp fac" << endl;
		cout <<  "Algorithm-wait time-burn: " << samp_fac*double(alg_time_sum-wait_time_sum)/(60.0*CLOCKS_PER_SEC) << " minutes." << endl;
	
	
		cout << Ntot;
		get_EF_dist();
		cout << samp_fac*double(alg_time_sum-wait_time_sum)/(60.0*CLOCKS_PER_SEC) << "\n";
	}
}


/// Updates the burnin procedure
void PMCMC::update_burnin(const unsigned int samp)
{
	if(samp < nburnin){
		burnin = true; if(samp < 2*pmcmc_init_samp) pup = FAST_UPDATE; else pup = SLOW_UPDATE;
	}
	else{	
		if(samp == nburnin){
			if(core == 0 && invT_dynamic == true) cout << "'invT' is fixed to " << invT << "." << endl;
			burnin = false; pup = NO_UPDATE;
		}
	}
}

		
/// Initialises quantities in the class
void PMCMC::initialise()
{
	for(auto p = 0u; p < N; p++) particle.push_back(State(details,data,model,obsmodel));
		
	auto ninit_samp = 10u;
	
	for(auto ru = 0u; ru < nrun; ru++){                              // Goes over all the runs
		Li_run[ru] = -LARGE;
		for(auto i = 0u; i < ninit_samp; i++){                       // Randomly samples parameters / states and picks the best
			if(core == 0){
				Pi.paramval = model.sample_from_prior();             // Randomly samples from the prior
				if(pmcmc_start_param == true){                       // Sets parameters to those used in the simulation TURNOFF
					for(auto th = 0u; th < model.param.size(); th++) Pi.paramval[th] = model.param[th].value;
				}
			}
			
			Li = obs_prob(Pi.paramval);                              // Calculates unbiased estimate of likelihood
			mpi.bcast(Li);
			
			if(Li > Li_run[ru]){
				Li_run[ru] = Li;
				Pi_run[ru] = mpi.particle_sample(ru,backpart,particle,obsmodel);
			}
		}
			
		if(core == 0){
			// This generates parameter samples near to the intial set (for the initial normal and MVN proposal distributions)
			for(auto loop = 0u; loop < pmcmc_init_samp; loop++){
				auto par = Pi_run[ru].paramval;		
				auto par_samp = model.sample_from_prior();

				ParamSample ps; ps.paramval.resize(par.size());
				//for(auto th = 0u; th < par.size(); th++) ps.paramval[th] = 0.9*par[th] + 0.1*par_samp[th];
				//for(auto th = 0u; th < par.size(); th++) ps.paramval[th] = 0.5*par[th] + 0.5*par_samp[th];
				for(auto th = 0u; th < par.size(); th++) ps.paramval[th] = 0.95*par[th] + 0.05*par_samp[th];
				//for(auto th = 0u; th < par.size(); th++) ps.paramval[th] = par_samp[th];
				ps.run = UNSET; ps.EF = UNSET;
				param_samp.push_back(ps);
			}
		}
		
		if(core == 0) Pri_run[ru] = model.prior(Pi.paramval);
	}
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
		
		if(std::isnan(obprob) || std::isinf(obprob)) emsgEC("PMCMC",1);
	}
	timer[TIME_PMCMCLIKE].stop();
	
	return obprob;
}


/// Performs the bootstrap step which randomly selects particles based on their observation probability
double PMCMC::bootstrap(const unsigned int sec, vector <double> &L)
{
	timer[TIME_BOOTSTRAP].start();
	
	vector <double> w(Ntot), wsum(Ntot);
	
	mpi.barrier();
	
	auto Ltot = mpi.gather(L);

	auto obprob = 0.0;
	auto sum = 0.0;
	if(core == 0){	
		double Lmax = -LARGE; 
		for(auto p = 0u; p < Ntot; p++){
			if(Ltot[p] > Lmax) Lmax = Ltot[p];
		}
	
		for(auto p = 0u; p < Ntot; p++){
			w[p] = exp(invT*(Ltot[p]-Lmax));
			sum += w[p];
			wsum[p] = sum;
		}
	
		obprob = invT*Lmax + log(sum/Ntot);
	}
	
	auto &bp = backpart[sec+1];
		
	if(core == 0){
		for(auto p = 0u; p < Ntot; p++) bp[p] = UNSET;
		 
		vector <unsigned int> extra; 
	
		for(auto p = 0u; p < Ntot; p++){
			auto z = ran()*sum;
			auto pp = 0u; while(pp < Ntot && z > wsum[pp]) pp++;
			if(pp == Ntot) emsgEC("PMCMC",2);
			
			if(bp[pp] == UNSET) bp[pp] = pp;
			else extra.push_back(pp);	
		}
		
		for(auto p = 0u; p < Ntot; p++){  
			if(bp[p] == UNSET){
				auto pp = extra[extra.size()-1]; extra.pop_back();
				bp[p] = pp;
			}
		}
		if(extra.size() != 0) emsgEC("PMCMC",3);
	}

	mpi.bcast(bp);
	
	timer[TIME_BOOTSTRAP].stop();
	
	return obprob;	
}


/// Gets the proposals used for the next itermation of MCMC
void PMCMC::get_proposals()
{
	if(core == 0){
		vector <ParamSample> param_samp_top;  // Gets just the top half of the samples 
		auto N = param_samp.size();
		unsigned int step = N/2000; if(step == 0) step = 1;
		
		for(auto i = N/2; i < param_samp.size(); i += step){
			param_samp_top.push_back(param_samp[i]);
		}
		
		prop_list = paramprop.get_proposal_list(param_samp_top); 
		//paramprop.zero_ntr_nac();
	}
	mpi.bcast(prop_list);

	if(core == 0 && diagnotic_output == true) cout << "# Proposals " << prop_list.size() << endl;
}


/// Updates parameters using MCMC proposals
void PMCMC::mcmc_updates()
{
	timer[TIME_MCMCPROP].start();
	
	for(auto i = 0u; i < prop_list.size(); i++){
		auto &prop = prop_list[i];
		auto num = prop.num;

		switch(prop.type){
			case SELF_PROP:
				if(burnin == true) self_proposal(paramprop.self);    
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
				
			case COVAR_AREA_PROP:
				covar_area_proposal(paramprop.covar_area[num]);    
				break;
				
			case FIXEDTREE_PROP: case SLICETIME_PROP: 
				emsgEC("PMCMC",4);
				break;
		}
	}
	timer[TIME_MCMCPROP].stop();
}


/// Returns the acceptance probability
double PMCMC::get_al()
{
	mpi.bcast(Pp.paramval);

	Lp = obs_prob(Pp.paramval);
	Pp = mpi.particle_sample(UNSET,backpart,particle,obsmodel);
	if(core == 0) Prp = model.prior(Pp.paramval);
	
	auto al = exp(Lp-Li + Prp-Pri);
	
	if(core == 0 && false){
		cout << al << " " << invT << " " << Lp <<  " " << Li <<  " " << Prp << " " << Pri << " proposal prob" << endl;
	}

	if(std::isnan(al)) emsgEC("PMCMC",5);
				
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
		if(false) cout << al << " " << invT << " al self" << endl;
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
			if(mvn.MH_PMCMC(al,paramprop.self.av_ac,pup) == SUCCESS) copy_propose_state();
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
void PMCMC::mean_time_proposal(MeanTime &mt)
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


/// Looks at making proposals to change a covariate and area effects
void PMCMC::covar_area_proposal(CovarArea &ca)
{
	timer[TIME_COVAR_AREA].start();
		
	auto al = 0.0;
	auto suc = 0u; if(core == 0){ if(ca.propose(Pp.paramval,Pi.paramval,model,data) == SUCCESS) suc = 1;}
	mpi.bcast(suc);
	
	if(suc == 1){
		al = get_al();
		if(core == 0){
			if(ca.MH(al,pup) == SUCCESS) copy_propose_state();
		}
	}
	
	timer[TIME_COVAR_AREA].stop();
}


/// Initialises quantities before PMCMC can start
void PMCMC::initialise_variables()
{
	backpart.resize(obsmodel.nsection+1); for(auto sec = 0u; sec <= obsmodel.nsection; sec++) backpart[sec].resize(Ntot);
	
	buffersize = 1 + data.narea*(1+ model.comp.size()*(1+data.ndemocatpos)); // This stores population sizes
	buffersize += 1 + data.nstrain*(1 + data.narage);                        // This stores Imap
	buffersize += 1 + data.nstrain*(1 + data.narage);                        // This stores Idiag
	
	auto nparam = model.param.size();
	
	Li_run.resize(nrun);          
	Pi_run.resize(nrun); for(auto ru = 0u; ru < nrun; ru++) Pi_run[ru].paramval.resize(nparam);
	Pri_run.resize(nrun);

	Pi.paramval.resize(nparam); Pp.paramval.resize(nparam);
	
	core = mpi.core;	
	
	invT_dynamic = false; 
	if(invT == UNSET){
		invT = 1000;
		invT_dynamic = true; 
		if(core == 0){
			cout << "'invT' is dynamically altered during burnin for efficient operation." << endl << endl;
		}
	}
}


/// Determines when to terminate the algorithm
bool PMCMC::terminate(const unsigned int samp)
{
	bool term = false;
	if(ESSmin != UNSET){
		if(samp%100 == 0){
			if(samp <= nburnin){
				if(mpi.core == 0) cout << "Burnin samples: " <<  samp << endl;	
			}
			else{
				auto psamp_tot = mpi.gather_psamp(particle_store);
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


/// Calculates the start invT based on getting the variance in the unbiased predictor
void PMCMC::calculate_start_invT()
{
	if(core == 0) Pp.paramval = Pi.paramval;
	
	mpi.bcast(Pp.paramval);

	auto loopmax = 40u;
	
	auto invTav = 0.0;
	for(auto loop = 0u; loop < loopmax; loop++){
		auto N = 100u;

		auto av = 0.0, av2 = 0.0; 
		for(auto i = 0u; i < N; i++){
			Lp = obs_prob(Pp.paramval);
			av += Lp; av2 += Lp*Lp;
		}
		auto sd = sqrt(av2/N - (av/N)*(av/N));
		
		auto fac = sd_init-sd;
		if(fac < -2) fac = -2;
		if(fac > 2) fac = 2;
		
		invT *= exp(1*fac);
		if(core == 0) cout << loop << " " << invT << " " <<  sd << " standard deviation\n";
		if(loop >= loopmax/2) invTav += invT/(loopmax/2);
	}
	invT = invTav;
	if(core == 0) cout << invT << " invT estimate\n";
	
	for(auto ru = 0u; ru < nrun; ru++){                              // Goes over all the runs
		Li_run[ru] = obs_prob(Pi.paramval);       
		
		Pi_run[ru] = mpi.particle_sample(ru,backpart,particle,obsmodel);
	}		
}


/// Outputs the distribution in error function
void PMCMC::get_EF_dist() const 
{
	State s(details,data,model,obsmodel);
	
	vector <double> EF_sample;
	for(const auto &p : particle_store){
		s.initialise_from_particle(p);  
		s.set_EF();
		EF_sample.push_back(s.EF);
	}
	
	sort(EF_sample.begin(),EF_sample.end());
	
	auto n = EF_sample.size();
	
	/*
	cout << "Distribution in EF (using cumulative distribution): " << endl;
	auto i = (unsigned int)((n-1)*0.025); auto f = (n-1)*0.025 - i;
	cout << "2.5%: " << EF_sample[i]*(1-f) + EF_sample[i+1]*f << "  ";
				
	i = (unsigned int)((n-1)*0.25); f = (n-1)*0.25 - i;
	cout << "25%: " << EF_sample[i]*(1-f) + EF_sample[i+1]*f << "  ";
		
	i = (unsigned int)((n-1)*0.75); f = (n-1)*0.75 - i;
	cout << "75%: " << EF_sample[i]*(1-f) + EF_sample[i+1]*f << "  ";
			
	i = (unsigned int)((n-1)*0.975); f = (n-1)*0.975 - i;
	cout << "97.5%: " << EF_sample[i]*(1-f) + EF_sample[i+1]*f << endl;
	*/
	
	auto i = (unsigned int)((n-1)*0.75); auto f = (n-1)*0.75 - i;
	cout << "," << EF_sample[i]*(1-f) + EF_sample[i+1]*f << ",";
	i = (unsigned int)((n-1)*0.975); f = (n-1)*0.975 - i;
	cout << EF_sample[i]*(1-f) + EF_sample[i+1]*f << ",";
}
