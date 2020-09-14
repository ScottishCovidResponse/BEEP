// This file contains all the functions for running a MCMC er the MBP algorithm

#include <fstream>
#include <iostream>
#include <sstream>
#include <algorithm>  
#include "stdlib.h"
#include "math.h"
#include "assert.h"

using namespace std;

#include "timers.hh"
#include "model.hh"
#include "utils.hh"
#include "chain.hh"
#include "output.hh"
#include "pack.hh"
#include "obsmodel.hh"

/// Initialises a single mcmc chain
Chain::Chain(const Details &details, const Data &data, const Model &model, const AreaTree &areatree, const ObservationModel &obsmodel, const Output &output, unsigned int _ch) : initial(details,data,model,obsmodel), propose(details,data,model,obsmodel), comp(model.comp), lev(areatree.lev), trans(model.trans), details(details), data(data), model(model), areatree(areatree), obsmodel(obsmodel), output(output)
{
	ch = _ch;

	initialise_variables();                      // Initialises the variables in the class
	
	sample_state();                              // Samples the prior and simulates to generate an initial state 
	
	jump.init(initial.paramval,model);                 // Initialises the kernel used for jumping in parameter space
	
	if(details.mode != MCMCMC) return;
	
	initial.set_L_and_Pr();                      // Sets the initial observation likelihood and prior
}


/// Performs a "standard" set of proposals
void Chain::standard_proposal(double EFcut) 
{
	timers.timestandard -= clock();
	
	timers.timembptemp -= clock();
	initial.set_process_likelihood();              // First sets the latent process likelihood Lev
	timers.timembptemp += clock();
	
	timers.timeparam -= clock();
	initial.standard_parameter_prop(jump);         // Makes changes to parameters with fixed event sequence
	standard_event_prop(EFcut);                    // Makes changes to event sequences with fixed parameters
	timers.timeparam += clock();

	timers.timestandard += clock();
}

		
/// Randomly generates a state by sampling from the prior and simulating an event sequence
void Chain::sample_state()
{
	unsigned int loop, loopmax = 100;
	for(loop = 0; loop < loopmax; loop++){       // Keeps simulating until successful
		if(simulate(model.sample_from_prior()) == SUCCESS) break;
	}
	
	if(loop == loopmax){
		stringstream ss;
		ss << "After '"+to_string(loopmax)+"' random simulations, it was not possible to find an initial ";
		ss << "state with the number of infected individuals below the threshold 'maximum_infected' specified in the input TOML file.";
		emsg(ss.str());
	}
}


/// Performs an MBP proposal on parameter 'th'
void Chain::mbp_proposal(unsigned int th)  
{
	timers.timembpprop -= clock();

	auto paramv = jump.mbp_prop(initial.paramval,th); // Generates a proposed parameter set 

	double al = 0; 
	if(mbp(paramv) == SUCCESS){                       // Performs the MBP and if successful calculates acceptance probability
		propose.set_L_and_Pr();
		
		al = exp(propose.Pr-initial.Pr + invT*(propose.L-initial.L));		
		if(checkon == true) cout << al << " " << invT << " " << propose.L << " " << initial.L << " al" << endl;
	}

	if(ran() < al){
		initial.copy_state(propose);                    // If successful copies the proposed state into the initial state
		jump.mbp_accept(th);
	}
	else{
		jump.mbp_reject(th);
	}

	if(checkon == true) initial.check();
	
	timers.timembpprop += clock();
}


/// Performs a model-based proposal (MBP)
Status Chain::mbp(const vector<double> &paramv)
{
	timers.timembpinit -= clock();
	
	if(model.inbounds(paramv) == false) return FAIL;      // Checks parameters are within the prior bounds
		
	if(propose.set_param(paramv) == FAIL) return FAIL;    // Sets quantities derived from parameters (if not possible fails)
	
	mbp_initialise();                                     // Prepares for the proposal

	timers.timembpinit += clock();
	
	timers.timembp -= clock();
		
	auto n = 0u;                                          // Indexes infections in the initial state 
	double t = 0;                                         // Current time
	
	for(auto sett = 0u; sett < details.ndivision; sett++){// Goes across discrete timesteps
		if(details.mode == SIM) output.print_populations(sett,N);
		
		initial.set_beta_and_phi(sett);                     // Sets beta and phi variables used later
		propose.set_beta_and_phi(sett); 
		
		propose.set_Qmap_using_dQ(sett,initial,dQmap);      // Calculates proposed Qmap from initial Qmap plus dQ
	
		construct_infection_sampler(initial.Qmap[sett],propose.Qmap[sett]); // Sampler used to generate new infections 
	
		double tmax = details.division_time[sett+1];
		do{
			auto tini = initial.get_infection_time(n);        // Gets time of next infection event in initial state
			auto tinf = t + exp_sample(new_infection_rate[0][0]);// Gets time of new added infection event
			
			if(tini >= tmax && tinf >= tmax){ t = tmax; break;}
			
			if(tinf < tini){                                  // A new infection appears in the propsed state
				t = tinf;
				add_infection_in_area(area_of_next_infection(),t);	
			}
			else{                                             // An event on initial sequence is considered to be copies
				t = tini;
				copy_event_or_not(n);
				n++;
			}
	
			if(propose.infev.size() >= model.maximum_infected){// If the number of infections exceeds the prior limit then fails
				reset_susceptible_lists(); return FAIL;
			}
		}while(1 == 1);
		
		update_dQmap(initial.transev[sett],propose.transev[sett]);// Based on infections updates dQmap
		if(checkon == true) check(t,sett);
	}

	timers.timembp += clock();
		
	reset_susceptible_lists();
	
	return SUCCESS;
}

/// Decides to copy infection event n from the initial to proposed states or not
void Chain::copy_event_or_not(unsigned int n)
{
	auto i = initial.infev[n].ind;

	if(susceptible_status[i] == BOTH_SUSCEPTIBLE){    // The copy can only be done if individual in propose is suscepticle
		auto w = data.ind[i].area*data.ndemocatpos + data.ind[i].dp;

		auto al = propose.lambda[w]/initial.lambda[w];  // Acceptance probability
		if(ran() < al){                                
			change_susceptible_status(i,NOT_SUSCEPTIBLE,1);// Copies the infection event to propose
			
			if(do_mbp_event == true) mbp_compartmental_transitions(i);
			else propose.indev[i] = initial.indev[i];
			
			propose.add_indev(i);
		}
		else change_susceptible_status(i,ONLY_PROPOSE_SUSCEPTIBLE,1);      // Does not copy the infection event propose
	}
}


/// Based on compartmental trasitions in the initial state for i
/// this uses MBPs to generates transitions for the proposal
void Chain::mbp_compartmental_transitions(unsigned int i)
{
	const auto &evlisti = initial.indev[i];
	auto &evlistp = propose.indev[i];
	const auto &parami = initial.paramval, &paramp = propose.paramval;
	const auto &comptransprobi = initial.comptransprob, &comptransprobp = propose.comptransprob;
	
	evlistp.clear();
	
	auto ev = evlisti[0];
	auto tra = ev.trans;
	auto ti = ev.t; 
	auto timep = ev.timep;
	evlistp.push_back(ev);
	
	auto a = data.democatpos[data.ind[i].dp][0];
	auto c = trans[tra].to;
	
	auto tp = ti;
	unsigned int emax = evlisti.size();
	bool switch_branch = false;
	for(auto e = 1u; e < emax; e++){                              // Goes through each event in initial
		tra = evlisti[e].trans;
		if(trans[tra].istimep == false){		                        // Igonores transtions related to change in time period
			unsigned int kmax = comp[c].trans.size();
			if(kmax == 0) break;
			
			double dtnew;
			if(kmax > 1){                                             // Looks at switching to another branch
				auto k = 0u; while(k < kmax && tra != comp[c].trans[k]) k++;
				if(k == kmax) emsgEC("model",4);
				
				auto ratio = comptransprobp[c].prob[a][k]/comptransprobi[c].prob[a][k];
				if(ratio < 1){ 
					if(ran() < 1 - ratio){
						auto sum = 0.0;                                     // Selects new branch to go down
						vector <double> sumst;
						for(auto k = 0u; k < kmax; k++){
							auto dif = comptransprobp[c].prob[a][k] - comptransprobi[c].prob[a][k];
							if(dif > 0) sum += dif;
							sumst.push_back(sum);
						}
						
						auto z = ran()*sum; auto k = 0u; while(k < kmax && z > sumst[k]) k++;
						if(k == kmax) emsgEC("Model",5);
						tra = comp[c].trans[k];
						switch_branch = true;
						dtnew = propose.sample_duration(tra);     
					}
				}	
			}
			
			if(switch_branch == false){
				auto dt = evlisti[e].t - ti;
				ti = evlisti[e].t;	
			
				switch(trans[tra].type){
				case EXP_DIST:
					{
						auto p = trans[tra].param_mean;
						dtnew = dt*paramp[p]/parami[p];
					}
					break;
				
				case GAMMA_DIST:
					emsgEC("model",6);
					break;
					
				case LOGNORM_DIST:
					{
						auto p = trans[tra].param_mean, p2 = trans[tra].param_cv;
						
						auto mean_nsi = parami[p], cv_nsi = parami[p2];
						auto mean_nsp = paramp[p], cv_nsp = paramp[p2];
						
						if(mean_nsi == mean_nsp && cv_nsi == cv_nsp) dtnew = dt;
						else{
							auto sdi = sqrt(log(1+cv_nsi*cv_nsi)), meani = log(mean_nsi) - sdi*sdi/2;
							auto sdp = sqrt(log(1+cv_nsp*cv_nsp)), meanp = log(mean_nsp) - sdp*sdp/2;
							dtnew = exp(meanp + (log(dt) - meani)*sdp/sdi); 
						}
					}
					break;
					
				default:
					emsgEC("Model",7);
					break;
				}
		
				if(dtnew < TINY) dtnew = TINY;
			}
		
			tp += dtnew;
	
			while(timep < model.ntime_period-1 && model.time_period[timep].tend < tp){    // Adds in changes in time period
				ev.trans = comp[c].transtimep; ev.t = model.time_period[timep].tend;
				evlistp.push_back(ev);
				timep++;
				ev.timep = timep; 
			}
			
			ev.trans = tra; ev.t = tp;
			evlistp.push_back(ev);
	
			c = trans[tra].to; 
			if(switch_branch == true) break;
		}
	}
	
	
	if(switch_branch == true){  // If the branch is switched then need to simulate rest of sequence
		if(comp[c].trans.size() != 0)	propose.simulate_compartmental_transitions(i,c,tp);
	}
}


/// Constructs a fast Gillespie sampler for finding the area of the next infection to be added	
void Chain::construct_infection_sampler(const vector <double> &Qmi, const vector <double> &Qmp)
{
	timers.infection_sampler -= clock();
		
	int l = areatree.level-1;
	for(auto c = 0u; c < data.narea; c++){
		auto wmin = c*data.ndemocatpos; 
		auto wmax = wmin + data.ndemocatpos;
	
		auto faci = initial.beta*initial.areafactor[c];
		auto facp = propose.beta*propose.areafactor[c];
		
		double sum = 0; 
		auto dp = 0u; 
		auto v = c*data.nage; 
		for(auto w = wmin; w < wmax; w++){
			auto a = data.democatpos[dp][0];
			initial.lambda[w] = initial.susceptibility[dp]*(faci*Qmi[v+a] + initial.phi);
			propose.lambda[w] = propose.susceptibility[dp]*(facp*Qmp[v+a] + propose.phi);
			auto dlambda = nboth_sysceptible_list[w]*(propose.lambda[w] - initial.lambda[w]); if(dlambda < 0) dlambda = 0;
			if(std::isnan(dlambda)){ emsgEC("Chain",400);}
			sum += dlambda + npropose_only_susceptible_list[w]*propose.lambda[w];
			dp++;
		}
		if(std::isnan(sum)) emsgEC("Chain",4);
		new_infection_rate[l][c] = sum;
	}
	
	for(int l = areatree.level-2; l >= 0; l--){                                 
		auto cmax = lev[l].node.size();
		for(auto c = 0u; c < cmax; c++){
			double sum = 0; for(const auto& ch : lev[l].node[c].child) sum += new_infection_rate[l+1][ch];
			
			new_infection_rate[l][c] = sum;
		}
	}
	
	timers.infection_sampler += clock();
}
	

/// Sets up lists of individuals (those suscepticle in both initial and propose states, only propose, and not at all)
void Chain::setup_susceptible_lists()
{	
	BOTH_SUSCEPTIBLE_list.clear(); BOTH_SUSCEPTIBLE_list.resize(data.nardp); 
	propose_only_susceptible_list.clear(); propose_only_susceptible_list.resize(data.nardp); 
	not_susceptible_list.clear(); not_susceptible_list.resize(data.nardp);

	nboth_sysceptible_list.resize(data.nardp); 
	npropose_only_susceptible_list.resize(data.nardp); 
	nnot_susceptible_list.resize(data.nardp);
	
	susceptible_status.resize(data.popsize); susceptible_list_ref.resize(data.popsize); 

	for(auto c = 0u; c < data.narea; c++){
		for(auto dp = 0u; dp < data.ndemocatpos; dp++){
			auto w = c*data.ndemocatpos + dp;

			for(const auto& i : data.area[c].ind[dp]){
				susceptible_status[i] = BOTH_SUSCEPTIBLE;
				susceptible_list_ref[i] = BOTH_SUSCEPTIBLE_list[w].size();
				BOTH_SUSCEPTIBLE_list[w].push_back(i);
			}
			nboth_sysceptible_list[w] = data.area[c].ind[dp].size();
			npropose_only_susceptible_list[w] = 0;
			nnot_susceptible_list[w] = 0;
		}
	}	
}


/// Makes a change in the susceptibility status of an individual
void Chain::change_susceptible_status(unsigned int i, unsigned int st, unsigned int updateR)
{
	auto c = data.ind[i].area;
	auto w = c*data.ndemocatpos + data.ind[i].dp;
	
	double dval = 0;
	if(updateR == 1){		               // Updates the event sampler for new events 
		auto dlambda = nboth_sysceptible_list[w]*(propose.lambda[w] - initial.lambda[w]); if(dlambda < 0) dlambda = 0;
		dval = -(dlambda + npropose_only_susceptible_list[w]*propose.lambda[w]);
	}
	
	int l = susceptible_list_ref[i];   // Removes the einitial.infevsiting entry
	int n;
	switch(susceptible_status[i]){
  case BOTH_SUSCEPTIBLE:
		if(BOTH_SUSCEPTIBLE_list[w][l] != i) emsgEC("Chain",7);
		n = BOTH_SUSCEPTIBLE_list[w].size();
		if(l < n-1){
			BOTH_SUSCEPTIBLE_list[w][l] = BOTH_SUSCEPTIBLE_list[w][n-1];
			susceptible_list_ref[BOTH_SUSCEPTIBLE_list[w][l]] = l;
		}
		BOTH_SUSCEPTIBLE_list[w].pop_back();
		nboth_sysceptible_list[w]--;
		break;
		
	case ONLY_PROPOSE_SUSCEPTIBLE:
		if(propose_only_susceptible_list[w][l] != i) emsgEC("Chain",8);
		n = propose_only_susceptible_list[w].size();
		if(l < n-1){
			propose_only_susceptible_list[w][l] = propose_only_susceptible_list[w][n-1];
			susceptible_list_ref[propose_only_susceptible_list[w][l]] = l;
		}
		propose_only_susceptible_list[w].pop_back();
		npropose_only_susceptible_list[w]--;
		break;
		
	case NOT_SUSCEPTIBLE:
		if(not_susceptible_list[w][l] != i) emsgEC("Chain",9);
		n = not_susceptible_list[w].size();
		if(l < n-1){
			not_susceptible_list[w][l] = not_susceptible_list[w][n-1];
			susceptible_list_ref[not_susceptible_list[w][l]] = l;
		}
		not_susceptible_list[w].pop_back();
		nnot_susceptible_list[w]--;
		break;
	
	default: emsgEC("Chain",10); break;
	}

	susceptible_status[i] = st;
	switch(susceptible_status[i]){
	case ONLY_PROPOSE_SUSCEPTIBLE:
		susceptible_list_ref[i] = propose_only_susceptible_list[w].size();
		propose_only_susceptible_list[w].push_back(i);
		npropose_only_susceptible_list[w]++;
		break;
		
	case NOT_SUSCEPTIBLE:
		susceptible_list_ref[i] = not_susceptible_list[w].size();
		not_susceptible_list[w].push_back(i);
		nnot_susceptible_list[w]++;
		break;
	
	case BOTH_SUSCEPTIBLE:
		susceptible_list_ref[i] = BOTH_SUSCEPTIBLE_list[w].size();
		BOTH_SUSCEPTIBLE_list[w].push_back(i);
		nboth_sysceptible_list[w]++;
		break;
		
	default: emsgEC("Chain",11); break;
	}
	
	if(updateR == 1){                              // Updates the infection sampler
		auto dlambda = nboth_sysceptible_list[w]*(propose.lambda[w] - initial.lambda[w]); if(dlambda < 0) dlambda = 0;
		dval += dlambda + npropose_only_susceptible_list[w]*propose.lambda[w];
		
		if(dval != 0){
			int l = areatree.level-1;
			do{
				new_infection_rate[l][c] += dval;
				c = lev[l].node[c].parent; l--;
			}while(l >= 0);
		}
	}
}


/// Places all individuals back onto the "both susceptible" list 
void Chain::reset_susceptible_lists()
{
	for(auto w = 0u; w < data.nardp; w++){
		for(const auto& i : propose_only_susceptible_list[w]){
			susceptible_status[i] = BOTH_SUSCEPTIBLE;
			susceptible_list_ref[i] = BOTH_SUSCEPTIBLE_list[w].size();
			BOTH_SUSCEPTIBLE_list[w].push_back(i);
			nboth_sysceptible_list[w]++;
		}
		propose_only_susceptible_list[w].clear();
		npropose_only_susceptible_list[w] = 0;
		
		for(const auto& i : not_susceptible_list[w]){
			susceptible_status[i] = BOTH_SUSCEPTIBLE;
			susceptible_list_ref[i] = BOTH_SUSCEPTIBLE_list[w].size();
			BOTH_SUSCEPTIBLE_list[w].push_back(i);
			nboth_sysceptible_list[w]++;
		}
		not_susceptible_list[w].clear();
		nnot_susceptible_list[w] = 0;
	}
}


/// Updates dQmap based on events that occur in the initial and proposed states
void Chain::update_dQmap(const vector <EventRef> &trei, const vector <EventRef> &trep)
{
	timers.timembpQmap -= clock();
	
	auto nage = data.nage;

	for(const auto& tre : trei){
		auto i = tre.ind; 
		auto tra = initial.indev[i][tre.e].trans;
		indmap[i][tra] = 1;
	}
	
	for(const auto& tre : trep){
		auto i = tre.ind; 
		auto tra = propose.indev[i][tre.e].trans;	
		if(indmap[i][tra] == 0){
			auto v = data.ind[i].area*nage+data.democatpos[data.ind[i].dp][0];
			auto dq = trans[tra].DQ[propose.indev[i][tre.e].timep];

			if(dq != UNSET){
				for(auto loop = 0u; loop < 2; loop++){
					auto q = model.DQ[dq].q[loop];
					if(q != UNSET){
						if(dQbuf[v][q] == 0){ dQbuflistv.push_back(v); dQbuflistq.push_back(q);}
						dQbuf[v][q] += model.DQ[dq].fac[loop];
					}
				}
			}
		}
		else indmap[i][tra] = 0;
	}	
	
	for(const auto& tre : trei){
		auto i = tre.ind; 
		auto tra = initial.indev[i][tre.e].trans;
		if(indmap[i][tra] == 1){
			auto v = data.ind[i].area*nage+data.democatpos[data.ind[i].dp][0];
			auto dq = trans[tra].DQ[initial.indev[i][tre.e].timep];
			if(dq != UNSET){
				for(auto loop = 0u; loop < 2; loop++){
					auto q = model.DQ[dq].q[loop];
					if(q != UNSET){
						if(dQbuf[v][q] == 0){ dQbuflistv.push_back(v); dQbuflistq.push_back(q);}
						dQbuf[v][q] -= model.DQ[dq].fac[loop];
					}
				}
			}
			indmap[i][tra] = 0;
		}
	}
	
	if(details.mode == SIM){
		for(const auto& tre : trep){
			auto tra = propose.indev[tre.ind][tre.e].trans;
			N[trans[tra].from]--;
			N[trans[tra].to]++;
		}
	}

	nage = data.nage;
	auto jmax = dQbuflistv.size();
	for(auto j = 0u; j < jmax; j++){
		auto v = dQbuflistv[j]; 
		auto q = dQbuflistq[j]; 
		auto qt = data.Q[q].Qtenref;
		
		auto fac = dQbuf[v][q];
		if(fac < -VTINY || fac > VTINY){
			auto kmax = data.genQ.Qten[qt].ntof[v];
			
			auto& cref = data.genQ.Qten[qt].tof[v];
			auto& valref = data.genQ.Qten[qt].valf[v];
			if(nage == 1){
				for(auto k = 0u; k < kmax; k++){
					dQmap[cref[k]] += fac*valref[k][0];
				}
			}
			else{
				for(auto k = 0u; k < kmax; k++){
					auto vv = cref[k]*nage;	
					for(auto a = 0u; a < nage; a++){
						dQmap[vv] += fac*valref[k][a];
						vv++;
					}
				}
			}
		}
		dQbuf[v][q] = 0;
	}
	dQbuflistv.clear(); dQbuflistq.clear(); 
	
	timers.timembpQmap += clock();
}
	
	
/// This samples which area the next new infection occurs in 
/// Starting at the top level l=0 the algorithm proceeds to finer and finer scales
unsigned int Chain::area_of_next_infection()
{
	double sumst[4];
	
	auto l = 0u, c = 0u;                            
	auto lmax = areatree.level;
	while(l < lmax-1){
		auto& node = lev[l].node[c];
		
		double sum = 0;
		auto jmax = node.child.size();
		for(auto j = 0u; j < jmax; j++){
			sum += new_infection_rate[l+1][node.child[j]];	
			sumst[j] = sum;
		}
		
		double z = ran()*sum; auto j = 0u; while(j < jmax && z > sumst[j]) j++; if(j == jmax) emsgEC("Chain",12);
		
		c = node.child[j];
		l++;
	};
	
	return c;
}


/// Adds a individual infected at time t in area c
void Chain::add_infection_in_area(unsigned int c, double t)
{
	auto dpmax = data.ndemocatpos;
	
	vector <double> sumst;
	sumst.resize(dpmax);
	
	double sum = 0;                                      // Selects which demographic possibility from the area
	for(auto dp = 0u; dp < dpmax; dp++){
		auto w = c*dpmax + dp;
		double dlambda = nboth_sysceptible_list[w]*(propose.lambda[w] - initial.lambda[w]); if(dlambda < 0) dlambda = 0;
		sum += dlambda + npropose_only_susceptible_list[w]*propose.lambda[w];
		sumst[dp] = sum;
	}

	double z = ran()*sum; auto dp = 0u; while(dp < dpmax && z > sumst[dp]) dp++; if(dp == dpmax) emsgEC("Chain",13);
	
	auto w = c*dpmax + dp;                               // Next selects susceptible list type
	double dlambda = nboth_sysceptible_list[w]*(propose.lambda[w] - initial.lambda[w]); if(dlambda < 0) dlambda = 0;

	unsigned int i;
	if(ran() < dlambda/(dlambda + npropose_only_susceptible_list[w]*propose.lambda[w])){ // Both suscetible
		auto n = BOTH_SUSCEPTIBLE_list[w].size(); if(n == 0) emsgEC("Chain",14);
		i = BOTH_SUSCEPTIBLE_list[w][(unsigned int)(ran()*n)];
	}
	else{                                                // Only proposed state susceptible
		auto n = propose_only_susceptible_list[w].size(); if(n == 0) emsgEC("Chain",15);
		i = propose_only_susceptible_list[w][(unsigned int)(ran()*n)];
	}
	
	change_susceptible_status(i,NOT_SUSCEPTIBLE,1);      // Changes the susceptibility status
	
	propose.simulate_compartmental_transitions(i,0,t);   // Simulates compartmental transitions after infection

	propose.add_indev(i);
}


/// Simulates an event sequence given a parameter set (returns fail if it is found not to be possible)
Status Chain::simulate(const vector <double>& paramv)
{
	vector <double> zero(paramv.size());                  // Sets up a parameter set with all zeros
	for(auto& pval : zero) pval = 0;
	initial.set_param(zero);
	
	initial.clear();                                      // Clears the initial event sequence
	
	if(mbp(paramv) == FAIL) return FAIL;                  // Performs an MBP (this effectively simulates the system)
	
	initial.copy_state(propose);                          // Copies proposed state into initial state
	
	return SUCCESS;
}


/// Used for checking the code is running correctly
void Chain::check(double t, unsigned int sett) const
{
	for(auto i = 0u; i < data.popsize; i++){             // Checks susceptible lises are correct
	  auto w = data.ind[i].area*data.ndemocatpos + data.ind[i].dp;
		
		if(nboth_sysceptible_list[w] != BOTH_SUSCEPTIBLE_list[w].size()) emsgEC("Chain",24);
		if(npropose_only_susceptible_list[w] != propose_only_susceptible_list[w].size()) emsgEC("Chain",25);
		if(nnot_susceptible_list[w] != not_susceptible_list[w].size()) emsgEC("Chain",26);
		
		if((initial.indev[i].size() == 0 || t < initial.indev[i][0].t) && propose.indev[i].size() == 0){
			if(susceptible_status[i] != BOTH_SUSCEPTIBLE) emsgEC("Chain",27);
			if(BOTH_SUSCEPTIBLE_list[w][susceptible_list_ref[i]] != i) emsgEC("Chain",28);
		}
		else{
			if((initial.indev[i].size() != 0 && t >= initial.indev[i][0].t) && propose.indev[i].size() == 0){
				if(susceptible_status[i] != ONLY_PROPOSE_SUSCEPTIBLE) emsgEC("Chain",29);
				if(propose_only_susceptible_list[w][susceptible_list_ref[i]] != i) emsgEC("Chain",30);
			}
			else{
				if(susceptible_status[i] != NOT_SUSCEPTIBLE) emsgEC("Chain",31);
				if(not_susceptible_list[w][susceptible_list_ref[i]] != i) emsgEC("Chain",32);
			}
		}
	}
	
	auto l = areatree.level-1;                         // Checks infection sampler is correct
	for(auto c = 0u; c < data.narea; c++){
		auto wmin = c*data.ndemocatpos, wmax = wmin + data.ndemocatpos;
	
		double sum = 0; 
		auto dp = 0u; 
		auto v = c*data.nage; 
		for(auto w = wmin; w < wmax; w++){
			auto a = data.democatpos[dp][0];
			
			double dd;
			dd = initial.lambda[w] - initial.susceptibility[dp]*(initial.beta*initial.areafactor[c]*initial.Qmap[sett][v+a] 
			                                                    + initial.phi); 
			if(sqrt(dd*dd) > TINY) emsgEC("Chain",33);
			dd = propose.lambda[w] - propose.susceptibility[dp]*(propose.beta*propose.areafactor[c]*propose.Qmap[sett][v+a] 
			                                                    + propose.phi); 
			if(sqrt(dd*dd) > TINY) emsgEC("Chain",34);
	
			auto dlambda = nboth_sysceptible_list[w]*(propose.lambda[w] - initial.lambda[w]); if(dlambda < 0) dlambda = 0;
			sum += dlambda + npropose_only_susceptible_list[w]*propose.lambda[w];
			dp++;
		}
		auto dd = new_infection_rate[l][c] - sum; if(sqrt(dd*dd) > TINY){ emsgEC("Chain",35);}
	}
	
	for(auto i = 0u; i < data.popsize; i++){
		for(auto tra = 0u; tra < model.trans.size(); tra++){
			if(indmap[i][tra] != 0) emsgEC("Chain",36);
		}
	}
}


/// Calculates propose.Qmap based on the initial and final sequences
void Chain::calculate_Qmap_for_propose()
{	
	for(auto& dQma : dQmap) dQma = 0;
	
	for(auto sett = 0u; sett < details.ndivision; sett++){
		propose.set_Qmap_using_dQ(sett,initial,dQmap);
		update_dQmap(initial.transev[sett],propose.transev[sett]);	
	} 
}


/// Adds and removes infectious individuals
void Chain::standard_event_prop(double EFcut)
{	
	timers.timeaddrem -= clock();

	propose.copy_state(initial);                             // Copies initial state into proposed state
		
	for(const auto& x : propose.infev) change_susceptible_status(x.ind,NOT_SUSCEPTIBLE,0);
	 
	double dprob;
	if(ran() < 0.5) dprob = add_individuals();
	else dprob = remove_individuals();
	
	if(dprob == UNSET){ reset_susceptible_lists(); return;}
	
	timers.timembptemp4 -= clock();
	propose.set_process_likelihood();
	
	double al;
	switch(details.mode){
	case ABC_MBP:
	
		propose.EF = obsmodel.observation_likelihood(propose.transev,propose.indev);
		if(propose.EF > EFcut) al = 0;
		else al = exp(propose.Lev-initial.Lev + dprob);
		if(checkon == true) cout << al <<  " " << EFcut << " al" << endl;		
		break;
	
	case MCMCMC:
		propose.L = obsmodel.observation_likelihood(propose.transev,propose.indev);
		al = exp(invT*(propose.L-initial.L) + propose.Lev-initial.Lev + dprob);
		if(checkon == true) cout << al << " al" << endl;		
		break;
		
	default: 
		emsgEC("Chain",67);
		break;
	}
	timers.timembptemp4 += clock();
	
	if(ran() < al){
		initial.copy_state(propose);
		jump.standev_accept();
	}
	else{
		jump.standev_reject();
	}
	
	reset_susceptible_lists();

	if(checkon == true) initial.check();
	timers.timeaddrem += clock();
}


/// Adds jump.naddrem individuals into the propose state
double Chain::add_individuals()
{
	auto probif = 0.0, probfi = 0.0;
		
	timers.timembptemp2 -= clock();
	infection_sampler(initial.Qmap);
	timers.timembptemp2 += clock();
	
	for(auto j = 0u; j < jump.naddrem; j++){
		if(propose.infev.size() >= model.maximum_infected) return UNSET;
		
		auto sumtot = lambdasum[data.nsettardp-1]; if(sumtot == 0) emsgEC("Chain",56);
		auto z = ran()*sumtot;
		
		//k = 0; while(k < data.nsettardp && z > lambdasum[k]) k += 1;
		//if(k == data.nsettardp) emsg("pr"); 
		
		auto k = 0u; auto dk = data.nsettardp/10; if(dk == 0) dk = 1;
		do{
			while(k < data.nsettardp && z > lambdasum[k]) k += dk;
			if(dk == 1) break;
			if(k >= dk) k -= dk;
			dk /= 10; if(dk == 0) dk = 1;
		}while(1 == 1);
		if(k >= data.nsettardp) emsgEC("Chain",57);
	
		if(k > 0){
			if(k >= lambdasum.size()) emsgEC("Chain",58);
			if(!(z < lambdasum[k] && z > lambdasum[k-1])) emsgEC("Chain",59);
		}
		
		probif += log(lambda[k]/sumtot);

		auto sett = k/data.nardp, w = k%data.nardp;
		
		if(nboth_sysceptible_list[w] == 0) return UNSET;
		auto i = BOTH_SUSCEPTIBLE_list[w][int(ran()*nboth_sysceptible_list[w])];
		probif += log(1.0/nboth_sysceptible_list[w]);
	
		change_susceptible_status(i,NOT_SUSCEPTIBLE,0);
		
		auto dt = details.division_time[sett+1]-details.division_time[sett];
		auto t = details.division_time[sett] + ran()*dt;
		probif += log(1.0/dt);
		
		propose.simulate_compartmental_transitions(i,0,t);
	
		propose.add_indev(i);
		
		probfi += log(1.0/propose.infev.size());
	}
	
	timers.timembptemp3 -= clock();
	propose.sort_infev();
	calculate_Qmap_for_propose();
	timers.timembptemp3 += clock();
	
	return probfi - probif;
}


/// Removes jump.naddrem individuals from the propose state
double Chain::remove_individuals()
{
	auto probif = 0.0, probfi = 0.0;
	
	vector <int> kst;
	for(auto j = 0u; j < jump.naddrem; j++){
		if(propose.infev.size() == 0) return UNSET;
		
		auto l = int(ran()*propose.infev.size());
		probif += log(1.0/propose.infev.size());
		auto i = propose.infev[l].ind;
		auto sett = (unsigned int)(details.ndivision*initial.indev[i][propose.infev[l].e].t/details.period); 

		auto c = data.ind[i].area;
		auto w = c*data.ndemocatpos + data.ind[i].dp;

		auto dt = details.division_time[sett+1]-details.division_time[sett];

		probfi += log(1.0/dt);
		kst.push_back(sett*data.nardp + w);
		
		propose.indev[i].clear();
		
		propose.infev[l] = propose.infev[propose.infev.size()-1];
		propose.infev.pop_back();
		
		change_susceptible_status(i,BOTH_SUSCEPTIBLE,0);
		
		probfi += log(1.0/nboth_sysceptible_list[w]);
	}
	
	for(auto& transev : propose.transev){  // Removes events in propose.transev
		auto j = 0u;
		auto jmax = transev.size();
		while(j < jmax){
			if(propose.indev[transev[j].ind].size() == 0){
				jmax--;
				transev[j] = transev[jmax];
				transev.pop_back();
			}
			else j++;
		}
	}
	
	timers.timembptemp3 -= clock();
	propose.sort_infev();
	calculate_Qmap_for_propose();
	timers.timembptemp3 += clock();
	
	timers.timembptemp2 -= clock();
	infection_sampler(propose.Qmap);
	timers.timembptemp2 += clock();
	
	auto sumtot = lambdasum[data.nsettardp-1]; 
	for(const auto& ks : kst) probfi += log(lambda[ks]/sumtot);
	
	return probfi - probif;
}


/// Generates a sampler for adding infected individuals into the system based on the force of infection
void Chain::infection_sampler(const vector< vector<double> > &Qmap)
{
	auto sum = 0.0;
	for(auto sett = 0u; sett < details.ndivision; sett++){
		auto phi = initial.disc_spline[model.phispline_ref][sett]; 
		auto beta = initial.disc_spline[model.betaspline_ref][sett];
	
		for(auto c = 0u; c < data.narea; c++){
			auto fac = beta*initial.areafactor[c];
			for(auto dp = 0u; dp < data.ndemocatpos; dp++){
				auto w = c*data.ndemocatpos + dp;
				auto tot = sett*data.nardp + w;
				auto v = c*data.nage + data.democatpos[dp][0];
				
				auto val = nboth_sysceptible_list[w]*initial.susceptibility[dp]*(fac*Qmap[sett][v] + phi);
				sum += val;
				
				lambda[tot] = val;				
				lambdasum[tot] = sum;
			}
		}
	}
}

/// Makes a MBP to a defined parameter set
Status Chain::abcmbp_proposal(const vector <double> &paramv, double EFcut)  
{
	auto al = 0.0;
	if(mbp(paramv) == SUCCESS){
		propose.set_EF();
		if(propose.EF < EFcut){
			propose.Pr = model.prior(propose.paramval);
			al = exp(propose.Pr-initial.Pr);
			if(checkon == true) cout << al << " " << propose.Pr << " " << initial.Pr <<  "al\n";
		}
	}
	
	if(ran() < al){
		initial.copy_state(propose);
		return SUCCESS;
	}

	return FAIL;
}


/// Initialises the variables within Chain
void Chain::initialise_variables()
{
	setup_susceptible_lists();
	
	indmap.resize(data.popsize);
	for(auto i = 0u; i < data.popsize; i++){
		indmap[i].resize(model.trans.size());
		for(auto tra = 0u; tra < model.trans.size(); tra++) indmap[i][tra] = 0;
	}
	
	dQmap.resize(data.narage);                                                 // Initialises vectors
	
	dQbuf.resize(data.narage);
	for(auto v = 0u; v < data.narage; v++){
		dQbuf[v].resize(data.Q.size()); for(auto q = 0u; q < data.Q.size(); q++) dQbuf[v][q] = 0;
	}
	dQbuflistv.clear(); dQbuflistq.clear();
	
	
	new_infection_rate.resize(areatree.level); for(auto l = 0u; l < areatree.level; l++) new_infection_rate[l].resize(lev[l].node.size()); 
	N.resize(comp.size()); 
	                                   // Used for event based changes
	lambda.resize(data.nsettardp); lambdasum.resize(data.nsettardp);
}


/// Clears variables ready for a MBP
void Chain::mbp_initialise()
{
	for(auto c = 0u; c < comp.size(); c++) N[c] = 0;
	N[0] = data.popsize;
		
	propose.clear();
	
	for(auto v = 0u; v < data.narage; v++) dQmap[v] = 0;
	
	// Set to true if MBPs on compartmental transitions needed
	do_mbp_event = model.do_mbp_events(initial.paramval,propose.paramval); 
}	
	