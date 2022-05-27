// This file contains all the functions for running MCMC using MBPs

#include <fstream>
#include <iostream>
#include <sstream>
#include <algorithm>  
#include "stdlib.h"
#include "math.h"
#include "assert.h"

using namespace std;

#include "mbp.hh"
#include "param_prop.hh"
#include "state.hh"
#include "model.hh"
#include "mpi.hh"
#include "details.hh"
#include "mpi.hh"
#include "obsmodel.hh"

/// Initialises mbp update
Mbp::Mbp(ObsModelMode obsmodel_mode_, const Details &details, const Data &data, const Model &model, const ObservationModel &obsmodel, const Output &output, Mpi &mpi) : state1(details,data,model,obsmodel), state2(details,data,model,obsmodel), comp(model.comp), trans(model.trans), details(details), data(data), model(model), obsmodel(obsmodel), output(output), mpi(mpi)
{
	obsmodel_mode = obsmodel_mode_;                        // Determines the observation model mode (either INVT or CUTOFF)
	
	initial = &state1; propose = &state2;                  // Sets the initial and proposed states (these may be swapped)

	initialise_variables();                                // Initialises the variables in the class
}


/// Updates a vector of particles using MBP proposals
unsigned int Mbp::mcmc_updates(vector <Particle> &part, const vector <ParamSample> &param_samp, double EFcut_, double invT_, ParamUpdate pup_, ParamProp &paramprop, const double cor_max)
{
	EFcut = EFcut_;
	invT = invT_;
	pup = pup_;
	
	timer[TIME_UPDATE].start();    
	timer[TIME_MVNSETUP].start(); 
	for(auto &mv : paramprop.mvn) mv.setup(param_samp);
	timer[TIME_MVNSETUP].stop(); 
	
	auto nprop = 0u;   
	auto psamp_before = mpi.gather_psamp(part);                          // Stores parameter samples

	timer[TIME_MBPUPDATE].start();     
	vector <double> cor;                                                 // Stores the correlation with initial samples
	paramprop.zero_ntr_nac_state(); 
	
	paramprop.zero_ntr_nac();	
	do{
		//paramprop.zero_ntr_nac();	
			
		for(auto &pa : part) update_particle(pa,paramprop,nprop);

		mpi.barrier();
		
		timer[TIME_UPDATEPROP].start();
		if(pup == COMBINE_DYNAMIC_UPDATE){      // For the first generation dynamically update during proposal
			paramprop.combine_update_proposals();
			paramprop.zero_ntr_nac();	
		}
		timer[TIME_UPDATEPROP].stop();
		

		timer[TIME_WAIT].start();
		auto psamp_after = mpi.gather_psamp(part);
	
		cor = output.get_correlation(psamp_before,psamp_after);
		
		if(mpi.core == 0 && false){
			cout << vec_max(cor)  << " " << cor_max << " correlation" << endl;
			for(auto th = 0u; th < model.param.size(); th++){
				cout << model.param[th].name << " " << cor[th] << " " << cor_max << " correlation" << endl;
			}
		}
		
		timer[TIME_WAIT].stop();
	}while(nprop/part.size() < details.prop_max && vec_max(cor) > cor_max);  // Iterates until correlation less than threshold
	timer[TIME_MBPUPDATE].stop(); 
	
	timer[TIME_UPDATEPROP].start();
	if(pup == COMBINE_UPDATE) paramprop.combine_update_proposals();
	timer[TIME_UPDATEPROP].stop();
	
	timer[TIME_UPDATEPROP].start();
	if(pup != NO_UPDATE) paramprop.update_fixed_splice_proposals(pup);
	timer[TIME_UPDATEPROP].stop();
	
	timer[TIME_UPDATE].stop();  
			
	return nprop/part.size();
}


/// Makes a standard update
void Mbp::update_particle(Particle &pa, ParamProp &paramprop, unsigned int &nprop)
{
	initial->initialise_from_particle(pa);                           // Initialises mbp updates from a particle

	timer[TIME_MCMCPROP].start();                                    // Type I MBPs
	for(auto loop = 0u; loop < 3; loop++){ 
		for(auto &mv : paramprop.mvn){
			if(mv.type != PARAM_PROP){ mvn_proposal(mv); nprop++;}
		}
	}
	timer[TIME_MCMCPROP].stop();

	timer[TIME_PARAMPROP].start();                                   // Parameter proposals
	auto Li = initial->likelihood();
	copy_initial_to_propse();
	
	for(auto loop = 0u; loop < 3; loop++){ 
		for(auto &mv : paramprop.mvn){
			if(mv.type == PARAM_PROP){ parameter_proposal(mv,Li); nprop++;}
		}
	}
	//for(auto loop = 0u; loop < 3; loop++){ parameter_proposal(paramprop.mvn[1],Li); nprop++;}
	timer[TIME_PARAMPROP].stop();
	
	timer[TIME_STATEPROP].start();                                   // Type II MBPs
	mbp_slicetime(paramprop.slicetime[int(ran()*paramprop.slicetime.size())]); nprop++;

	if(paramprop.fixedtree.size() > 0){
		mbp_fixedtree(paramprop.fixedtree[int(ran()*paramprop.fixedtree.size())]); nprop++;
	}
	timer[TIME_STATEPROP].stop();     
	
	pa = initial->create_particle(pa.run);
	if(details.mode == ABC_MBP && pa.EF > EFcut) emsgEC("MBP",101);	
}


/// Updates a particle using MBP proposals
unsigned int Mbp::mc3_mcmc_updates(Particle &part, const vector <ParamSample> &param_samp, const double invT_, const ParamUpdate pup_, ParamProp &paramprop)
{
	EFcut = UNSET;
	invT = invT_;
	pup = pup_;
	
	auto nprop = 0u;   
	
	timer[TIME_MVNSETUP].start(); 
	for(auto &mv : paramprop.mvn) mv.setup(param_samp);
	timer[TIME_MVNSETUP].stop(); 
	
	for(auto loop = 0u; loop < model.param.size(); loop++) update_particle(part,paramprop,nprop);

	return nprop;
}


/*
/// Updates a particle using various types of MBP (specified in a proposal list)
void Mbp::update_particle(Particle &pa, const vector <Proposal> &prop_list, ParamProp &paramprop)
{
	initial->initialise_from_particle(pa);                           // Initialises mbp updates from a particle

	auto param_prop_flag = false;

	timer[TIME_MCMCPROP].start();
	for(auto i = 0u; i < prop_list.size(); i++){
		auto &prop = prop_list[i];
		auto num = prop.num;
		
		switch(prop.type){
			case MVN_PROP:
				{ 
					auto &mvn = paramprop.mvn[num];
					switch(mvn.type){
						case SIGMA_PARAM: sigma_reff_proposal(mvn); break;
						case DIST_R_JOINT: dist_R_joint_proposal(mvn); break;
						case PARAM_PROP: param_prop_flag = true; break;
						default: mvn_proposal(mvn); break;
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
			
			case FIXEDTREE_PROP:
				mbp_fixedtree(paramprop.fixedtree[num]);
				break;
				
			case SLICETIME_PROP:
				mbp_slicetime(paramprop.slicetime[num]);
				break;
			
			case SELF_PROP:
				emsgEC("MBP",1);
				break;
		}				
		
		if(checkon == true) initial->check(2);
	}
	timer[TIME_MCMCPROP].stop();
	
	timer[TIME_PARAMPROP].start();
	if(param_prop_flag == true){
		auto Li = initial->likelihood();
		copy_initial_to_propse();
		
		for(auto i = 0u; i < prop_list.size(); i++){
			auto &prop = prop_list[i];
			auto num = prop.num;
			
			if(prop.type == MVN_PROP){
				auto &mvn = paramprop.mvn[num];
				if(mvn.type == PARAM_PROP){
					parameter_proposal(mvn,Li); 
				}
				
				if(false){
					if(Li != initial->likelihood()) emsgEC("Mbp",43); 
				}
			}
		}
	}
	timer[TIME_PARAMPROP].stop();
		
	pa = initial->create_particle(pa.run);                       // Places final state back into particle
} 
*/


/// Swaps the initial and proposed states
void Mbp::swap_initial_propose_state()
{
	State* temp = initial;
	initial = propose;
	propose = temp;
}


/// Performs a model-based proposal (MBP)
Status Mbp::mbp(const vector<double> &paramv, const InfUpdate inf_update)
{	
	if(model.inbounds(paramv) == false) return FAIL;            // Checks parameters are within the prior bounds

	propose->set_param(paramv);                                 // Sets quantities derived from parameters (if not possible fails)

	timer[TIME_MBP].start();
	
	timer[TIME_MBPINIT].start();
	mbp_initialise();                                           // Prepares for the proposal
	timer[TIME_MBPINIT].stop();

	for(auto sett = 0u; sett < details.ndivision; sett++){      // Performs a pure MBPs or a combination of MBP and simulation
		propose->democat_change_pop_adjust(sett);
	
		switch(inf_update){                                       // Sets Imap  
			case INF_UPDATE: propose->set_Imap_sett(sett); break;
			case INF_DIF_UPDATE: propose->set_Imap_using_dI(sett,initial,dImap,dIdiag); break;
		}
	
		timer[TIME_TRANSNUM].start();
		double rate_i, rate_p;
		int num_i, num_p=0;
	
		for(auto c = 0u; c < data.narea; c++){                    // Performs simulation / MBPs on the transitions
			propose->set_transmean(sett,c);
			
			auto &init_tnum = initial->transnum[sett][c];
			auto &prop_tnum = propose->transnum[sett][c];
			auto &init_tmean = initial->transmean[sett][c];
			auto &prop_tmean = propose->transmean[sett][c];
			
			auto sorm = simu_or_mbp[sett][c];
			
			for(auto tr = 0u; tr < model.trans.size(); tr++){
				for(auto dp = 0u; dp < data.ndemocatpos; dp++){	 
					num_i = init_tnum[tr][dp];
					rate_p = prop_tmean[tr][dp];
					if(rate_p == 0) num_p = 0;
					else{
						switch(sorm){
							case SIMU:	
								num_p = poisson_sample(rate_p);								
								break;
						
							case MBP:
								rate_i = init_tmean[tr][dp];	
								if(rate_p > rate_i) num_p = num_i + poisson_sample(rate_p-rate_i);
								else{
									if(rate_p == rate_i) num_p = num_i;
									else num_p = binomial_sample(rate_p/rate_i,num_i);
								}
								break;
						}
					}
					
					auto num_diff = num_p - num_i;
					dtransnum[c][tr][dp] = num_diff;
					prop_tnum[tr][dp] = num_p;			
				}
			}
		}
		timer[TIME_TRANSNUM].stop();
			
		if(sett < details.ndivision-1){
			timer[TIME_UPDATEPOP].start();
			propose->update_pop(sett);
			timer[TIME_UPDATEPOP].stop();
			
			timer[TIME_UPDATEIMAP].start();
			if(inf_update == INF_DIF_UPDATE) propose->update_I_from_transnum(dImap,dIdiag,dtransnum);
			timer[TIME_UPDATEIMAP].stop();
		}
	}
	
	timer[TIME_MBP].stop();

	return SUCCESS;
}


/// Gets the acceptance probability
double Mbp::get_al()
{
	auto al = 0.0;
	
	propose->set_EF();
	propose->set_Pr();
		
	switch(obsmodel_mode){
		case CUTOFF: 
			if(propose->EF < EFcut) al = exp(propose->Pr - initial->Pr); 
			break;
		
		case INVT:  // Note EF = -2*log(obsmodelprob)
			al = exp(-invT*0.5*(propose->EF-initial->EF) + propose->Pr - initial->Pr);
			break;
	}		

	if(std::isnan(al)) emsgEC("Mbp",2);
				
	if(false){
		cout << al << " " << propose->EF << " " << initial->EF << " " <<  propose->Pr;
		cout << " " << initial->Pr <<  "mbp al" << endl;
	}
	
	if(al > 1) al = 1;

	return al;
}


/// Initialises the variables within Mbp
void Mbp::initialise_variables()
{
	dtransnum.resize(data.narea);
	for(auto c = 0u; c < data.narea; c++){
		dtransnum[c].resize(model.trans.size());
		for(auto tr = 0u; tr < model.trans.size(); tr++){
			dtransnum[c][tr].resize(data.ndemocatpos);
			for(auto dp = 0u; dp < data.ndemocatpos; dp++) dtransnum[c][tr][dp] = 0;
		}				
	}
	
	simu_or_mbp.resize(details.ndivision);
	for(auto sett = 0u; sett < details.ndivision; sett++) simu_or_mbp[sett].resize(data.narea);
	simu_or_mbp_reset();
	
	dImap.resize(data.nstrain); dIdiag.resize(data.nstrain);
	for(auto st = 0u; st < data.nstrain; st++){	
		dImap[st].resize(data.narage); dIdiag[st].resize(data.narage);      
	} 

	auto &tn = data.genQ.treenode;                                    // Sets mbp_sim for fixedtree proposals
	mbp_sim.resize(tn.size());
	for(auto n = 0u; n < tn.size(); n++){
		mbp_sim[n].resize(data.narea);
		for(auto c = 0u; c < data.narea; c++) mbp_sim[n][c] = MBP;
		for(auto c : tn[n].arearef) mbp_sim[n][c] = SIMU;
	}
	
	if(false){
		cout << mpi.core << "core" << endl;
		for(auto n = 0u; n < tn.size(); n++){  
			cout << "Node: " << n << "  "; 
			for(auto c = 0u; c < data.narea; c++) if(mbp_sim[n][c] == MBP) cout << "M"; else cout << "S";
			cout << endl;
		}
	}
}


/// Resets simu_or_mbp to just perform pure MBPs
void Mbp::simu_or_mbp_reset()
{
	for(auto sett = 0u; sett < details.ndivision; sett++){
		for(auto c = 0u; c < data.narea; c++) simu_or_mbp[sett][c] = MBP;
	}
}


/// Clears variables ready for a MBP
void Mbp::mbp_initialise()
{
	for(auto st = 0u; st < data.nstrain; st++){	
		for(auto v = 0u; v < data.narage; v++){
			dImap[st][v] = 0;
			dIdiag[st][v] = 0;
		} 
	}

	propose->pop_init();
}	


/// Below are all the possible proposals that can be carried out ///

/// Does a MVN proposals shifted by Langevin term to account for prior
void Mbp::mvn_proposal(MVN &mvn)  
{
	timer[TIME_MVN].start();

	InfUpdate inf_update = INF_DIF_UPDATE;
	if(mvn.type == INF_PARAM) inf_update = INF_UPDATE;
	
	auto al = 0.0;
	vector <double> param_prop;
	double probif;
	if(mvn.propose_langevin(param_prop,initial->paramval,probif,model) == SUCCESS){
		if(mbp(param_prop,inf_update) == SUCCESS){	
			al = get_al()*exp(mvn.get_probfi(initial->paramval,param_prop,model) - probif);
		}
	}
	
	
	if(mvn.MH(al,pup) == SUCCESS){
		if(mvn.type == R_NEIGH) mvn.add_prop_vec(initial->paramval,propose->paramval);
		swap_initial_propose_state();
	}
	else{
		if(mvn.type == R_NEIGH) mvn.add_prop_vec(initial->paramval,initial->paramval);
	}
	
	timer[TIME_MVN].stop();
}


/// Looks to make changes to parameters with a fixed state 
void Mbp::parameter_proposal(MVN &mvn, double &Li)
{
	if(checkon == true){
		if(Li != initial->likelihood()) emsgEC("MBP",102);
	}
	
	auto al = 0.0, Li_prop=0.0;
	vector <double> param_prop;
	double probif;
	if(mvn.propose_langevin(param_prop,initial->paramval,probif,model) == SUCCESS){
		propose->set_param(param_prop);   
	
		for(auto sett = 0u; sett < details.ndivision; sett++){   
			for(auto c = 0u; c < data.narea; c++){ 
				propose->set_transmean(sett,c);
			}
		}
		
		Li_prop = propose->likelihood();
		propose->set_Pr();
	
		al = exp(Li_prop -  Li + propose->Pr - initial->Pr)*exp(mvn.get_probfi(initial->paramval,param_prop,model) - probif);
	}

	if(mvn.MH(al,pup) == SUCCESS){
		propose->set_EF();
		if(propose->EF < EFcut){
			Li = Li_prop;
			swap_initial_propose_state();
		}
	}
}

	
/// Does a MVN proposals shifted by Langevin term to account for prior
void Mbp::mvn_proposal_rev(MVN &mvn)  
{
	timer[TIME_MVN].start();

	auto al = 0.0;
	vector <double> param_prop;
	double probif;
	
	//param_prop = initial->paramval;
	//mbp_rev(param_prop);
	
	if(mvn.propose_langevin(param_prop,initial->paramval,probif,model) == SUCCESS){
		if(mbp_rev(param_prop) == SUCCESS){	
			al = get_al()*exp(mvn.get_probfi(initial->paramval,param_prop,model) - probif);
		}
	}
	
	if(mvn.MH(al,pup) == SUCCESS) swap_initial_propose_state();
	
	timer[TIME_MVN].stop();
}


/// Does a joint proposal on means/ branching probs and R
void Mbp::dist_R_joint_proposal(MVN &mvn) 
{
	timer[TIME_MVN].start();
	
	InfUpdate inf_update = INF_DIF_UPDATE;
	if(mvn.type == INF_PARAM) inf_update = INF_UPDATE;
	
	auto al = 0.0;
	auto param_prop = mvn.propose(initial->paramval);
	if(model.inbounds(param_prop) == true){ 	
		auto fac = 1.0;                                      // This updates the value of R based on 
		propose->set_param(param_prop);
		
		auto genT_st = propose->genT;
		
		auto genT_fac = propose->genT/initial->genT;
		
		for(auto th : model.R_param){
			auto R = propose->paramval[th];
			propose->paramval[th] = pow(R,genT_fac);
			fac *= genT_fac*pow(R,genT_fac-1);
		}
		
		if(mbp(param_prop,inf_update) == SUCCESS){	
			if(propose->genT != genT_st) emsgEC("MBP",70);
		
			al = get_al()*fac;
		}
	}
	
	if(mvn.MH(al,pup) == SUCCESS) swap_initial_propose_state();
	
	timer[TIME_MVN].stop();
}


///  Simulatenously changes sigma and regional effects
void Mbp::sigma_reff_proposal(MVN &mvn)
{
	timer[TIME_SIGMA].start();
	
	auto al = 0.0;
	vector <double> param_prop;
	if(mvn.sigma_propose(param_prop,initial->paramval,model) == SUCCESS){
		if(mbp(param_prop,INF_DIF_UPDATE) == SUCCESS){ 
			al = get_al()*exp(initial->Pr-propose->Pr);
		}
	}
	else al = -1;
	
	if(mvn.MH(al,pup) == SUCCESS) swap_initial_propose_state();
	
	timer[TIME_SIGMA].stop();
}


/// Simulatenously changes the mean time of one transition and minus for the subsequent transitions 
void Mbp::mean_time_proposal(MeanTime &mt)
{
	timer[TIME_MEANTIME].start();
		
	auto al = 0.0;
	vector <double> param_prop;
	if(mt.propose(param_prop,initial->paramval,model) == SUCCESS){
		if(mbp(param_prop,INF_DIF_UPDATE) == SUCCESS){ 
			al = get_al();
		}
	}
	else al = -1;
	
	if(mt.MH(al,pup) == SUCCESS) swap_initial_propose_state();
	
	timer[TIME_MEANTIME].stop();
}


/// Looks at making proposals to between pair of adjacent points on a spline
void Mbp::neighbour_proposal(Neighbour &rn)
{
	timer[TIME_NEIGHBOUR].start();
		
	auto al = 0.0;
	vector <double> param_prop;
	if(rn.propose(param_prop,initial->paramval,model) == SUCCESS){
		if(mbp(param_prop,INF_DIF_UPDATE) == SUCCESS){ 
			al = get_al();
		}
	}
	else al = -1;
	
	if(rn.MH(al,pup) == SUCCESS) swap_initial_propose_state();
	
	timer[TIME_NEIGHBOUR].stop();
}


/// Looks at making proposals to between pair of adjacent points on a spline
void Mbp::joint_proposal(Joint &rn)
{
	timer[TIME_JOINT].start();
	
	auto al = 0.0;
	vector <double> param_prop;
	if(rn.propose(param_prop,initial->paramval,model) == SUCCESS){
		if(mbp(param_prop,INF_DIF_UPDATE) == SUCCESS){ 
			al = get_al();
		}
	}
	else al = -1;
	
	if(rn.MH(al,pup) == SUCCESS) swap_initial_propose_state();
	
	timer[TIME_JOINT].stop();
}


/// Looks at making proposals to change a covariate and area effects
void Mbp::covar_area_proposal(CovarArea &ca)
{
	timer[TIME_COVAR_AREA].start();
		
	auto al = 0.0;
	vector <double> param_prop;
	if(ca.propose(param_prop,initial->paramval,model,data) == SUCCESS){
		if(mbp(param_prop,INF_DIF_UPDATE) == SUCCESS){ 
			al = get_al();
		}
	}
	else al = -1;
			
	if(ca.MH(al,pup) == SUCCESS) swap_initial_propose_state();
	
	timer[TIME_COVAR_AREA].stop();
}


/// Keeps parameters fixed and mixs MBPs with simulation when performing the proposals on specific areas
void Mbp::mbp_fixedtree(FixedTree &ft)  
{
	timer[TIME_FIXEDTREE].start();
	
	auto frac = ft.sim_frac;                                              // Sets to simulate for specific areas
	for(auto c = 0u; c < data.narea; c++){             
		if(mbp_sim[ft.n][c] == SIMU){
			for(auto sett = 0u; sett < details.ndivision; sett++){ if(ran() < frac) simu_or_mbp[sett][c] = SIMU;}
		}
	}
	
	auto al = 0.0; 
	if(mbp(initial->paramval,INF_DIF_UPDATE) == SUCCESS) al = get_al();			

	if(details.mode == MC3_INF || details.mode == MCMC_MBP){                          // Dynamically updates simulation fraction
		if(al > 1) al = 1;
		if(pup == FAST_UPDATE) update(ft.sim_frac,al,eta_fast);
		if(pup == SLOW_UPDATE) update(ft.sim_frac,al,eta);
		if(ft.sim_frac > 1) ft.sim_frac = 1; 
	}
	
	ft.ntr++; ft.nac+= al;
	if(ran() < al) swap_initial_propose_state();
	
	simu_or_mbp_reset();
	
	timer[TIME_FIXEDTREE].stop();
}


/// Keeps parameters fixed and mixs MBPs with simulation when performing proposals on specific time periods
void Mbp::mbp_slicetime(SliceTime &st)  
{
	timer[TIME_SLICETIME].start();
		
	auto frac = st.sim_frac;                                          // Sets to simulate for a given time period
	for(auto sett = 0u; sett < details.ndivision; sett++){
		if(sett >= st.sett_i && sett < st.sett_f){
			for(auto c = 0u; c < data.narea; c++){ if(ran() < frac) simu_or_mbp[sett][c] = SIMU;}
		}
	}
			
	auto al = 0.0; 
	if(mbp(initial->paramval,INF_DIF_UPDATE) == SUCCESS){ 
		al = get_al();			
	}

	if(details.mode == MC3_INF || details.mode == MCMC_MBP){                          // Dynamically updates simulation fraction
		if(al > 1) al = 1;
		if(pup == FAST_UPDATE) update(st.sim_frac,al,eta_fast);
		if(pup == SLOW_UPDATE) update(st.sim_frac,al,eta);
		if(st.sim_frac > 1) st.sim_frac = 1; 
	}
	
	st.ntr++; st.nac += al;
	if(ran() < al) swap_initial_propose_state();
	
	simu_or_mbp_reset();
	
	timer[TIME_SLICETIME].stop();
}


/// Performs a reverse model-based proposal (MBP)
/// This was found not to work because it was difficult to get back to the initial starting state 
Status Mbp::mbp_rev(const vector<double> &paramv)
{	
	if(model.inbounds(paramv) == false) return FAIL;            // Checks parameters are within the prior bounds

	propose->set_param(paramv);                                 // Sets quantities derived from parameters (if not possible fails)

	timer[TIME_MBP].start();
	
	timer[TIME_MBPINIT].start();
	mbp_initialise();                                           // Prepares for the proposal
	timer[TIME_MBPINIT].stop();

	auto probif = 0.0, probfi = 0.0;	
	propose->Imap = initial->Imap;
	propose->Idiag = initial->Idiag;
	propose->transmean = initial->transmean;
	propose->transnum = initial->transnum;
	propose->pop = initial->pop;
	
	int sett_start = 100;
	
	for(auto sett = sett_start; sett >= 0; sett--){  
		cout << sett << " " << details.ndivision << " sett" << endl;
		//propose->democat_change_pop_adjust(sett);  TO DO
	
		propose->set_Imap_using_dI(sett,initial,dImap,dIdiag);
	
		timer[TIME_TRANSNUM].start();
		double rate_i, rate_p;
		int num_i, num_p=0;
	
		for(auto c = 0u; c < data.narea; c++){                    // Performs simulation / MBPs on the transitions
			propose->set_transmean(sett,c);
			
			auto &init_tnum = initial->transnum[sett][c];
			auto &prop_tnum = propose->transnum[sett][c];
			auto &init_tmean = initial->transmean[sett][c];
			auto &prop_tmean = propose->transmean[sett][c];
			
			for(auto tr = 0u; tr < model.trans.size(); tr++){
				for(auto dp = 0u; dp < data.ndemocatpos; dp++){	 
					num_i = init_tnum[tr][dp];
					rate_p = prop_tmean[tr][dp];
				
					rate_i = init_tmean[tr][dp];	
					
					if(rate_p > rate_i){
						num_p = num_i + poisson_sample(rate_p-rate_i);
						probif += poisson_probability(num_p-num_i,rate_p-rate_i);
						probfi += binomial_probability(num_i,rate_i/rate_p,num_p);
					} 
					else{
						if(rate_p == rate_i) num_p = num_i;
						else{
							num_p = binomial_sample(rate_p/rate_i,num_i);
							probif += binomial_probability(num_p,rate_p/rate_i,num_i);
							probfi += poisson_probability(num_i-num_p,rate_i-rate_p);			
						}		
					} 
								
					auto num_diff = num_p - num_i;
					dtransnum[c][tr][dp] = -num_diff;
					prop_tnum[tr][dp] = num_p;		
				}
			}
		}
		timer[TIME_TRANSNUM].stop();
			
		
		timer[TIME_UPDATEPOP].start();
		propose->update_pop_rev(sett);
		timer[TIME_UPDATEPOP].stop();
		
		timer[TIME_UPDATEIMAP].start();
		propose->update_I_from_transnum(dImap,dIdiag,dtransnum);
		timer[TIME_UPDATEIMAP].stop();
	}
	propose->EF = -2*obsmodel.calculate(propose);
	propose->Pr = model.prior(propose->paramval);
	 
	propose->check(0);
	
	//auto Li_initial = initial->likelihood();
	//auto Li_propose = propose->likelihood();
	
	timer[TIME_MBP].stop();

	return SUCCESS;
}


/// Copies the initial state to the proposed state (used in parameter proposals)
void Mbp::copy_initial_to_propse()
{
	propose->transnum = initial->transnum; propose->pop = initial->pop;
	propose->Imap = initial->Imap; propose->Idiag = initial->Idiag; 
	//propose->EF = initial->EF; propose->transmean = initial->transmean;
}
