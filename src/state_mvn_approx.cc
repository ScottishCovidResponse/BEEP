/// Functions related to making a MVN approximiation to the populations

#include <math.h>  

using namespace std;

#include "state.hh"
#include "matrix.hh"

/// Calculates likelihood based mon MVN approximation to state
double State::likelihood_approx(const vector <double> &paramval, vector <ObsSlice> &obs_slice, const double invT, const Accuracy ac)
{
	timer[TIME_LIKELIHOOD_APPROX].start();

	set_param(paramval);
	
	pop_init();
	
	pop_covar.resize(details.ndivision);
	
	vector < vector <double> > covar;
	covar.resize(ncomp_approx);
	for(auto c = 0u; c < ncomp_approx; c++){
		covar[c].resize(ncomp_approx);
		for(auto cc = 0u; cc < ncomp_approx; cc++){
			covar[c][cc] = 0;
		}
	}
	
	//test_transmean_gradients(100); return 0;

	auto L = 0.0;                         // Keeps track of the log-likelihood
 	
	auto o = 0u;
	for(auto sett = 0u; sett <= details.ndivision; sett++){
		//cout << sett << " sett" << endl;
	
		if(o < obs_slice.size()){
			auto &os = obs_slice[o];
			if(sett == os.sett && 1 == 1){ // Incorporate an observation
				timer[TIME_OBS_APPROX].start();
				if(sett > 0){
					switch(os.obs_type){
						case OBS_APPROX: L += apply_observation_fast(os,covar,invT,ac); break;	
						case OBS_EXACT: L += apply_exact_observation(os,covar); break;
					}
				}
				timer[TIME_OBS_APPROX].stop();
				o++;
			}
		}
		
		if(sett == details.ndivision) break;
		
		pop_covar[sett] = covar;
		
		set_Imap_sett(sett);
		
		for(auto c = 0u; c < data.narea; c++) set_transmean(sett,c);
	
		transnum[sett] = transmean[sett];
		
		timer[TIME_GRAD].start();
		set_transmean_gradients(sett);
		timer[TIME_GRAD].stop();
		
		auto mu = set_mu(sett);
		
		timer[TIME_COVAR].start();
		//covar_change(covar,pop_covar[sett],mu);
		covar_change_fast(covar,pop_covar[sett],mu,ac,true);	
		timer[TIME_COVAR].stop();
		
		timer[TIME_CORRECT].start();
		auto mean = set_p(sett);
		
		valid_matrix(covar,mean);
		
		timer[TIME_CORRECT].stop();
		
		update_pop(sett);
	}

	if(std::isnan(L)) emsgEC("State Approx",38);
	if(o != obs_slice.size()) emsgEC("State Approx",39);
		
	timer[TIME_LIKELIHOOD_APPROX].stop();

	return L;
}


/// Calculates the estimated future observation probability distribution based on MVN approximation
void State::future_obs_approx(vector <ObsSlice> &obs_slice, const double invT, const Accuracy ac)
{
	timer[TIME_FUTURE_OBS_APPROX].start();
	emsg("pro");
	//set_param(paramval);

	obs_mean.resize(details.ndivision);
	obs_covar.resize(details.ndivision);
	
	vector < vector <double> > covar;
	covar.resize(ncomp_approx);
	for(auto c = 0u; c < ncomp_approx; c++){
		covar[c].resize(ncomp_approx);
		for(auto cc = 0u; cc < ncomp_approx; cc++){
			if(c == cc) covar[c][cc] = 10000;
			else covar[c][cc] = 0;
		}
	}
	
	int o = obs_slice.size()-1;
	for(int sett = details.ndivision-1; sett >= 0; sett--){	
		set_I_from_pop(sett,false);
	
		if(o >= 0){
			auto &os = obs_slice[o];
			if(sett == int(os.sett)){ // Incorporate an observation
				timer[TIME_OBS_APPROX].start();
				if(sett > 0){
					switch(os.obs_type){
						case OBS_APPROX: apply_observation(os,covar,invT); break;
						case OBS_EXACT: apply_exact_observation(os,covar); break;
					}
					pop_positive(sett);
				}
				timer[TIME_OBS_APPROX].stop();
				o--;			
			}
		}
		
		obs_mean[sett] = set_p(sett);
		obs_covar[sett] = covar;
	
		for(auto c = 0u; c < data.narea; c++) set_transmean(sett,c);

		transnum[sett] = transmean[sett];
		
		timer[TIME_GRAD].start();
		set_transmean_gradients(sett);
		timer[TIME_GRAD].stop();
			
		auto mu = set_mu(sett);
	
		timer[TIME_COVAR].start();
		covar_change_fast(covar,pop_covar[sett],mu,ac,false);
		timer[TIME_COVAR].stop();
		
		timer[TIME_CORRECT].start();
		auto mean = set_p(sett);
		valid_matrix(covar,mean);
		//correct_mu_covar(mean,covar);
		timer[TIME_CORRECT].stop();
		
		if(sett > 0){ 
			update_pop_rev(sett);
			pop_positive(sett-1);
		}
	}
	timer[TIME_FUTURE_OBS_APPROX].stop();
}


/// Uses particles to estimate a sample from the posterior
Particle State::posterior_particle_sample(const vector <double> &paramval, vector <ObsSlice> &obs_slice, const unsigned int P, const double invT) 
{
	set_param(paramval);
	
	vector <ParticleApprox> part(P);
	for(auto p = 0u; p < P; p++){
		part[p].transnum.resize(details.ndivision);
		part[p].pop_store.resize(details.ndivision);
	}
	
	vector < vector <unsigned int> > back;
	back.resize(obs_slice.size());
	
	auto ti = 0u;
	for(auto o = 0u; o < obs_slice.size(); o++){
		const auto &os = obs_slice[o];
		
		auto tf = os.sett; 
		if(tf > details.ndivision) emsgroot("Observation out of range");

		if(tf != 0){
			if(os.obs_type == OBS_EXACT) emsgroot("Observations cannot be exact");
			
			for(auto p = 0u; p < P; p++){	
				if(ti > 0) pop[ti] = part[back[o-1][p]].pop_end;			
				
				simulate(ti,tf);
				for(auto t = ti; t < tf; t++){
					part[p].transnum[t] = transnum[t];
					part[p].pop_store[t] = pop[t];
				}
				
				part[p].pop_end_new = pop[tf];	
				
				// Calculates the likelihood of the obseration
				auto L = 0.0;
				for(const auto i : os.obs_ref){
					const auto &ob = data.obs[i];
						
					const auto &dt = data.datatable[ob.datatable];

					auto val = 0.0;
					auto var = ob.var_approx/invT;
					
					switch(dt.type){
						case POP:
							{
								for(auto ca : obs_clist[i]){
									const auto &cap = comp_approx[ca];
									val += pop[tf][cap.c][cap.co][cap.dp];
								}
							}
							break;
							
						case TRANS:
							{
								auto sett = ob.sett_i;
								
								auto pp = p, oo = o;
								
								while(oo > 0 && sett < obs_slice[oo-1].sett){
									//emsg("check working");
									oo--; pp = back[oo][pp];
								}
					
								for(auto ca : obs_clist[i]){
									const auto &cap = comp_approx[ca];
									val += pop[tf][cap.c][cap.co][cap.dp] - part[pp].pop_store[sett][cap.c][cap.co][cap.dp];
								}
							}
							break;
							
						case POPFRAC: emsgroot("'popfrac' not currently supported"); break;
						
						case MARGINAL: emsgroot("'marginal' cannot be used with maximum likelihood approaches."); break;	
					}
					
					auto value = ob.value;
					L += -0.5*(val-value)*(val-value)/var;
				}	
	
				part[p].L = L;
			}
			
			for(auto p = 0u; p < P; p++) part[p].pop_end = part[p].pop_end_new;
			
			auto Lav = 0.0; for(auto p = 0u; p < P; p++) Lav += part[p].L;
			Lav /= P;
			
			auto wsum = 0.0;
			for(auto p = 0u; p < P; p++){
				wsum += exp(part[p].L - Lav);
				part[p].wsum = wsum;
			}
			if(wsum == 0) emsgEC("State approx",54);
			for(auto p = 0u; p < P; p++) part[p].wsum /= wsum;
			
			for(auto p = 0u; p < P; p++){			
				back[o].resize(P);
		
				auto z = ran();
				auto pp = 0u; while(pp < P && z > part[pp].wsum) pp++;
				back[o][p] = pp;
			}
		}
		ti = tf;
	}
	
	auto p = 0u; 
	for(int o = obs_slice.size()-1; o >= 0; o--){
		if(back[o].size() > 0){
			p = back[o][p];
			auto ti = 0u; if(o > 0) ti = obs_slice[o-1].sett;
			auto tf = obs_slice[o].sett;
		
			for(auto t = ti; t < tf; t++) transnum[t] = part[p].transnum[t];
			
			if(o == int(obs_slice.size()-1)){ 
				pop[tf] = part[p].pop_end;			
				simulate(tf,details.ndivision);
			}
		}
	}
	
	return create_particle(UNSET);
}


/// Generates a MVN approximation to the observation model
Particle State::posterior_sample(const vector <double> &paramval, vector <ObsSlice> &obs_slice, const double invT, const Accuracy ac) 
{	
	likelihood_approx(paramval,obs_slice,invT,ac);
		
	future_obs_approx(obs_slice,invT,ac);
	
	pop_init();
	for(auto sett = 0u; sett < details.ndivision; sett++){
		democat_change_pop_adjust(sett);
		
		set_Imap_sett(sett);
	
		timer[TIME_TRANSNUM].start();
		
		for(auto c = 0u; c < data.narea; c++) set_transmean(sett,c);
		
		if(sett < details.ndivision-1){ 
			switch(1){
			case 0:  // This generates a forward prediction, multiplies by observation prob and estimates transition numbers
				{
					/* Works out the mean on the next time step */
					for(auto c = 0u; c < data.narea; c++){
						for(auto tr = 0u; tr < model.trans.size(); tr++){
							for(auto dp = 0u; dp < data.ndemocatpos; dp++){	 
								transnum[sett][c][tr][dp] = transmean[sett][c][tr][dp];							
							}
						}
					}
					
					update_pop(sett);
				
					auto mean1 = remove_S(set_p(sett+1));
							
					/* Works out the covariance of the next time step */		
					auto mu = set_mu(sett);
					
					vector < vector <double> > M1;
					M1.resize(ncomp_approx);
					for(auto c = 0u; c < ncomp_approx; c++){
						M1[c].resize(ncomp_approx);
						for(auto cc = 0u; cc < ncomp_approx; cc++) M1[c][cc] = 0;
					}
					
					for(auto ta = 0u; ta < ntrans_approx; ta++){
						const auto &tap = trans_approx[ta];
						auto cci = tap.ca_from;
						auto ccf = tap.ca_to;
						auto val = mu[ta];
						M1[ccf][ccf] += val;
						M1[cci][cci] += val;
						M1[cci][ccf] -= val;
						M1[ccf][cci] -= val;
					}
					valid_matrix(M1,mu);
					
					M1 = remove_S(M1);
					
					/* Distribution based on future observations */
					auto mean2 = remove_S(obs_mean[sett+1]);
					auto M2 = remove_S(obs_covar[sett+1]);
					valid_matrix(M2,mean2);
				
					/* Multiplies distributions to get posterior distribution */
					auto Minv1 = invert_matrix(M1);
					auto Minv2 = invert_matrix(M2);
					
					auto Minv = matrix_add(Minv1,Minv2);
					auto M = invert_matrix(Minv);
					
					auto Minv_mean1 = matrix_mult(Minv1,mean1);
					auto Minv_mean2 = matrix_mult(Minv2,mean2);
					auto Minv_mean = vec_add(Minv_mean1,Minv_mean2);
					auto mean = matrix_mult(M,Minv_mean);
					
					MVN mv("Distribution",vector <unsigned int>(),0,PARAM_PROP,MULTIPLE);
					//mv.sample_check(mean,M);
					auto p = mv.sample(mean,M);
					
					/* Adjusts transition number to ensure population p at the next time step */
					vector < vector <double> > T_new;
					for(auto j = 0u; j < ncomp_without_S; j++) T_new.push_back(T[comp_without_S[j]]);
								
					auto Tinv_new = invert_matrix(T_new);

					auto p_dif = vec_subtract(p,mean1);
				
					auto dtrans = matrix_mult(Tinv_new,p_dif);
				
					for(auto ta = 0u; ta < ntrans_approx; ta++){
						const auto &tap = trans_approx[ta];
						if(std::isnan(dtrans[ta])) emsg("NaN");
						transnum[sett][tap.c][tap.tr][tap.dp] += dtrans[ta];
					}
				}
				break;
			
			case 1:  // This adjusts rates according to observation model (steering)
				{
					for(auto c = 0u; c < data.narea; c++){
						for(auto tr = 0u; tr < model.trans.size(); tr++){
							for(auto dp = 0u; dp < data.ndemocatpos; dp++){	 
								transnum[sett][c][tr][dp] = transmean[sett][c][tr][dp];							
							}
						}
					}
					
					update_pop(sett);
					
					auto p = remove_S(set_p(sett+1));
					auto mean = remove_S(obs_mean[sett+1]);
					auto M = remove_S(obs_covar[sett+1]);
					auto Minv = invert_matrix(M);
					
					auto dif = vec_subtract(p,mean);
					auto dif_mag = vec_mult(dif,dif); 
					auto Minv_dif = matrix_mult(Minv,dif);
			
					vector <double> factor_rate_change(ntrans_approx);
					for(auto tr = 0u; tr < ntrans_approx; tr++){
						auto ta = trans_approx[tr];
						auto val = -Minv_dif[ta.ca_to-1];
						if(ta.ca_from > 0) val += Minv_dif[ta.ca_from-1];
						if(val < -5) val = -5; 
						if(val > 5) val = 5;
						factor_rate_change[tr] = val;		
					}	
					
					/* This loop ensures that the population aimed at is at least closer than the un steered version */
					auto limit = 1.0, difnew_mag = 0.0;
					do{
						for(auto c = 0u; c < data.narea; c++){
							for(auto tr = 0u; tr < model.trans.size(); tr++){
								for(auto dp = 0u; dp < data.ndemocatpos; dp++){	 
									auto ta = trans_approx_ref[c][tr][dp];
									transnum[sett][c][tr][dp] = exp(limit*factor_rate_change[ta])*transmean[sett][c][tr][dp];		
								}
							}
						}
						update_pop(sett);
						
						auto pnew = remove_S(set_p(sett+1));
						auto difnew = vec_subtract(pnew,mean);
						difnew_mag = vec_mult(difnew,difnew); 
						if(difnew_mag > dif_mag) limit *= 0.7;
					}while(difnew_mag > dif_mag);
					
					/* Finally samples transitions based on modified transition rates */
					for(auto c = 0u; c < data.narea; c++){
						for(auto tr = 0u; tr < model.trans.size(); tr++){
							for(auto dp = 0u; dp < data.ndemocatpos; dp++){	 
								auto ta = trans_approx_ref[c][tr][dp];
								transnum[sett][c][tr][dp] = poisson_sample(exp(limit*factor_rate_change[ta])*transmean[sett][c][tr][dp]);	
							}
						}
					}
				}
				break;
				
			case 2:    // This simulates multiple times and filters based on observation model
				{
					auto mean = remove_S(obs_mean[sett+1]);
					auto M = remove_S(obs_covar[sett+1]);
					auto Minv = invert_matrix(M);
				
					vector <double> val_st;   
					vector < vector < vector < vector <double> > > > transnum_store;
				
					auto valav = 0.0;
					const auto loopmax = 1000;
					for(auto loop = 0; loop < loopmax; loop++){
						for(auto c = 0u; c < data.narea; c++){
							for(auto tr = 0u; tr < model.trans.size(); tr++){
								for(auto dp = 0u; dp < data.ndemocatpos; dp++){	 
									transnum[sett][c][tr][dp] = poisson_sample(transmean[sett][c][tr][dp]);							
								}
							}
						}
					
						update_pop(sett);
						auto p = remove_S(set_p(sett+1));
						auto dif = vec_subtract(p,mean);
					
						auto Minv_dif = matrix_mult(Minv,dif);
					
						auto val = -0.5*vec_mult(dif,Minv_dif);
						valav += val;
						val_st.push_back(val);
					
						transnum_store.push_back(transnum[sett]);
					}					
					valav /= loopmax;
					
					auto wsum = 0.0;
					vector <double> wsum_st;  
					for(auto loop = 0; loop < loopmax; loop++){		
						auto dval = val_st[loop]-valav;
						if(dval < -10) dval = -10; 
						if(dval > 10) dval = 10;
						auto w = exp(dval);
						wsum += w;
						wsum_st.push_back(wsum);
					}
					
					auto z = ran()*wsum;
					auto k = 0u; while(k < loopmax && z > wsum_st[k]) k++;
					if(k == loopmax) emsgEC("State approx",10);
					transnum[sett] = transnum_store[k];
				}
				break;
				
			case 3: // Simply samples from the model
				{
					for(auto c = 0u; c < data.narea; c++){
						set_transmean(sett,c);
						
						for(auto tr = 0u; tr < model.trans.size(); tr++){
							for(auto dp = 0u; dp < data.ndemocatpos; dp++){	 
								transnum[sett][c][tr][dp] = poisson_sample(transmean[sett][c][tr][dp]);							
							}
						}
					}
				}
				break;
			}
		}
		
		timer[TIME_TRANSNUM].stop();
			
		if(sett < details.ndivision-1){ update_pop(sett);}
	}
	
	return create_particle(UNSET);
}


/// Applies an observation to the the system
double State::apply_observation(const ObsSlice &os, vector < vector <double> > &covar, const double invT)
{
	/* Loads exisiting state */
	auto sett = os.sett;          
	
	auto mu_start = set_p(sett);
	
	auto mu1 = remove_S(mu_start);
	auto M1 = remove_S(covar);
	
	vector < vector <double> > Minv1;
		
	bool pos_def;                               // This makes sure the matrix is positive definite
	auto loop = 0u, loopmax = 100u;
	do{
		Minv1 = invert_matrix(M1);
		
		pos_def = true;
		
		for(auto i = 0u; i < Minv1.size(); i++){
			if(Minv1[i][i] < 0) pos_def = false;
		}
		
		if(pos_def == false){
			pos_def = false;
			for(auto c = 0u; c < M1.size(); c++){
				for(auto cc = 0u; cc < M1[c].size(); cc++){
					if(c != cc) M1[c][cc] *= 0.8;
				}
			}
		}
		loop++;
	}while(loop < loopmax && pos_def == false);
	if(loop == loopmax) emsgEC("State approx",54);
	
	/*  Calculates the matrices for the observation */
	vector < vector <double> > Minv_single, Minv_obs;            // Calculates the matrices for the observation
	
	Minv_single.resize(ncomp_approx); Minv_obs.resize(ncomp_approx);
	for(auto cc = 0u; cc < ncomp_approx; cc++){
		Minv_single[cc].resize(ncomp_approx); Minv_obs[cc].resize(ncomp_approx);
		for(auto cc2 = 0u; cc2 < ncomp_approx; cc2++) Minv_obs[cc][cc2] = 0;
	}

	vector <double> Minv_mu_obs(ncomp_approx);
	for(auto cc = 0u; cc < ncomp_approx; cc++) Minv_mu_obs[cc] = 0;
	
	auto mu_Minv_mu_obs = 0.0;
	auto log_factor = 0.0;

	for(const auto ob_ref : os.obs_ref){
		const auto &ob = data.obs[ob_ref];
		const auto &dt = data.datatable[ob.datatable];

		vector <double> mu(ncomp_approx);
		for(auto cc = 0u; cc < ncomp_approx; cc++) mu[cc] = UNSET;
			
		for(auto cc = 0u; cc < ncomp_approx; cc++){
			for(auto cc2 = 0u; cc2 < ncomp_approx; cc2++){
				Minv_single[cc][cc2] = 0;
			}
		}
	
		double value = ob.value;
		
		if(dt.type == TRANS){  // In the case of transitions we add the total population as the start of the observation
			auto sett = ob.sett_i;
			for(auto ca : obs_clist[ob_ref]){
				const auto &cap = comp_approx[ca];
				value += pop[sett][cap.c][cap.co][cap.dp];
			}
		}
				
		const auto &clist = obs_clist[ob_ref];
					
		for(auto cc : clist) mu[cc] = value/clist.size();
		//print_vector("mu",mu);
		auto mu_shift = vec_subtract(mu,mu_start); //qq
		
		auto var = ob.var_approx/invT;
					
		auto val = 1.0/var;
				
		for(auto cc : clist){
			for(auto cc2 : clist){
				Minv_single[cc][cc2] = val;
			}
		}
	
		log_factor += log(val)-log(2*M_PI);
			
		auto Minv_mu = matrix_mult(Minv_single,mu_shift);
		auto mu_Minv_mu = vec_mult(mu_shift,Minv_mu);
		
		Minv_obs = matrix_add(Minv_obs,Minv_single);
		Minv_mu_obs = vec_add(Minv_mu_obs,Minv_mu);
		mu_Minv_mu_obs += mu_Minv_mu;
	}
	timer[TIME_TEMP2].stop();

	timer[TIME_TEMP3].start();
	/* Combines	existing distribution with observation distribution */
	
	auto Minv2 = remove_S(Minv_obs);

	auto Minv_mu2 = remove_S(Minv_mu_obs);

	auto Minv = matrix_add(Minv1,Minv2);
	auto Minv_mu = Minv_mu2;
	
	auto M = invert_matrix(Minv);
	timer[TIME_TEMP3].stop();
	
	timer[TIME_TEMP4].start();
	auto mu = matrix_mult(M,Minv_mu);
	
	auto mu_Minv_mu = vec_mult(mu,Minv_mu);
	
	const auto mu_Minv_mu1 = 0;
	const auto mu_Minv_mu2 = mu_Minv_mu_obs;
	
	auto prod = matrix_mult(Minv1,M);

	auto mu_final = vec_add(mu,mu1);

	for(auto &val : mu_final){ if(val < 0) val = 0;}
	
	auto vec = add_S(mu_final);
	
	covar = add_S(M);
	
	correct_mu_covar(vec,covar);
	valid_matrix(covar,vec);
	
	put_p(sett,vec);
	
	auto dL = 0.5*(log_factor + determinant_fast(prod) - (mu_Minv_mu1 + mu_Minv_mu2  - mu_Minv_mu));
	if(dL > 0) dL = 0;
	
	if(sett > 0){   // Works out what transitions numbers would have been need to get observed change in population
		auto dif = vec_subtract(vec,mu_start);
		
		auto dif_without_S = remove_S(dif);
		auto dtrans = matrix_mult(Tinv_without_S,dif_without_S);
	
		for(auto ta = 0u; ta < ntrans_approx; ta++){
			const auto &tap = trans_approx[ta];
			transnum[sett-1][tap.c][tap.tr][tap.dp] += dtrans[ta];
		}
	}
	timer[TIME_TEMP4].stop();
	
	return dL;
}


/// Applies an observation to the the system
double State::apply_observation_fast(const ObsSlice &os, vector < vector <double> > &covar, const double invT, const Accuracy ac)
{
	/* Loads exisiting state */
	auto sett = os.sett;          
	
	auto mu_start = set_p(sett);
	
	auto mu1 = remove_S(mu_start);
	auto M1 = remove_S(covar);
	
	//auto Minv = invert_matrix(M1);
	
	auto det1 = 0.0;
	auto Minv = invert_determinant_SIMD(M1,det1,ac);
	//auto Minv = invert_matrix_SIMD(M1,ac);

	/*
	vector < vector <double> > Minv;
		
	bool pos_def;                               // This makes sure the matrix is positive definite
	auto loop = 0u, loopmax = 100u;
	do{
		Minv = invert_matrix(M1);
		
		pos_def = true;
		
		for(auto i = 0u; i < Minv.size(); i++){
			if(Minv[i][i] < 0) pos_def = false;
		}
		
		if(pos_def == false){
			pos_def = false;
			for(auto c = 0u; c < M1.size(); c++){
				for(auto cc = 0u; cc < M1[c].size(); cc++){
					if(c != cc) M1[c][cc] *= 0.8;
				}
			}
		}
		loop++;
	}while(loop < loopmax && pos_def == false);
	cout << loop << " loop\n";
	if(loop == loopmax) emsgEC("State approx",54);
	*/
	
	/*  Calculates the matrices for the observation */

	vector <double> Minv_mu(ncomp_without_S,0);
	auto mu_Minv_mu2 = 0.0;
	
	auto log_factor = 0.0;

	for(const auto ob_ref : os.obs_ref){
		const auto &ob = data.obs[ob_ref];
		const auto &dt = data.datatable[ob.datatable];

		double value = ob.value;
		
		if(dt.type == TRANS){  // In the case of transitions we add the total population as the start of the observation
			auto sett = ob.sett_i;
			for(auto ca : obs_clist[ob_ref]){
				const auto &cap = comp_approx[ca];
				value += pop[sett][cap.c][cap.co][cap.dp];
			}
		}
				
		const auto &clist_without_S = obs_clist_without_S[ob_ref];
				
		auto var = ob.var_approx/invT;
					
		auto val = 1.0/var;
			
		auto value_scale = value/clist_without_S.size();
		for(auto cc : clist_without_S){
			for(auto cc2 : clist_without_S){
				Minv[cc][cc2] += val;
				Minv_mu[cc] += val*(value_scale-mu1[cc2]);
				mu_Minv_mu2 += (value_scale-mu1[cc])*val*(value_scale-mu1[cc2]);
			}
		}
	
		log_factor += log(val)-log(2*M_PI);
	}

	/* Combines	existing distribution with observation distribution */
	
	//auto M = invert_matrix(Minv);
	//auto M = invert_matrix_SIMD(Minv,DOUBLE);
	
	auto det = 0.0;
	auto M = invert_determinant_SIMD(Minv,det,ac);
	//auto M = invert_matrix_SIMD(Minv,ac);
	
	//auto M2 = invert_matrix_SIMD(Minv,DOUBLE);
	//cout << determinant_fast(M) << " " << determinant_fast(M2) << " j\n";
	//emsg("P");
	
	//invert_matrix_SIMD(const vector <vector <double> > &M_orig, Accuracy ac)   

	
	auto mu = matrix_mult(M,Minv_mu);
	
	auto mu_Minv_mu = vec_mult(mu,Minv_mu);
	
	const auto mu_Minv_mu1 = 0;

	auto mu_final = vec_add(mu,mu1);

	for(auto &val : mu_final){ if(val < 0) val = 0;}
	
	auto vec = add_S(mu_final);
	
	covar = add_S(M);
	
	correct_mu_covar(vec,covar);
	valid_matrix(covar,vec);
	
	put_p(sett,vec);
	
	/*
	cout << determinant_fast(M1) << "fast\n";
	cout << determinant_SIMD(M1) << "simd\n";
emsg("p");
	*/
	//auto det1 = determinant_SIMD(M1);
	//auto det2 = determinant_SIMD_float(M1);
	//cout << det1 << " " << det2 <<" " << det1-det2 << "j\n";
	// emsg("P");
	
	//auto dL = 0.5*(log_factor - determinant_SIMD(M1,ac) - determinant_SIMD(Minv,ac));
  //dL -= 0.5*(mu_Minv_mu1 + mu_Minv_mu2  - mu_Minv_mu);
	
	auto dL = 0.5*(log_factor - det1 - det);
  dL -= 0.5*(mu_Minv_mu1 + mu_Minv_mu2  - mu_Minv_mu);


	
	//cout << det << " " << determinant_SIMD(Minv,ac) << "K\n";
	//cout << det1 << " " <<  determinant_SIMD(M1,ac) << "K1\n";
//	cout << determinant_SIMD(M1,ac) << " " << det1 << " compar\n";
	
	//auto dL = 0.5*(log_factor - determinant_fast(M1) - determinant_fast(Minv)  - (mu_Minv_mu1 + mu_Minv_mu2  - mu_Minv_mu));
	if(dL > 0) dL = 0;
	
	if(sett > 0){   // Works out what transitions numbers would have been need to get observed change in population
		auto dif = vec_subtract(vec,mu_start);
		
		auto dif_without_S = remove_S(dif);
		auto dtrans = matrix_mult(Tinv_without_S,dif_without_S);
	
		for(auto ta = 0u; ta < ntrans_approx; ta++){
			const auto &tap = trans_approx[ta];
			transnum[sett-1][tap.c][tap.tr][tap.dp] += dtrans[ta];
		}
	}
	
	return dL;
}


/// Applies an observation to the the system
double State::apply_exact_observation(const ObsSlice &os, vector < vector <double> > &covar)
{
	auto sett = os.sett; 
	
	vector <double> mu_new(comp_approx.size());
	
	for(auto ca = 0u; ca < ncomp_approx; ca++) mu_new[ca] = UNSET;
	
	for(const auto i : os.obs_ref){
		const auto &ob = data.obs[i];
		const auto &dt = data.datatable[ob.datatable];
		
		if(dt.type != POP) emsgEC("State approx",13);
		
		vector <unsigned int> clist;   // Makes a list of all compartments 
		for(auto co : dt.complist){
			for(auto c : ob.area){
				for(auto dp : ob.dp_sel){
					auto cc = comp_approx_ref[c][co][dp];
					clist.push_back(cc);
				}
			}
		}	
		if(clist.size() != 1) emsg("Cannot have multiple state");
		
		mu_new[clist[0]] = ob.value;
	}
	
	vector <unsigned int> list1, list2;
	
	for(auto ca = 0u; ca < ncomp_approx; ca++){
		if(mu_new[ca] == UNSET) list1.push_back(ca);
		else list2.push_back(ca);
	}	
	auto nlist1 = list1.size(), nlist2 = list2.size();
	
	vector < vector <double> > M11, M12, M21, M22; 

	M11.resize(nlist1);
	for(auto j = 0u; j < nlist1; j++){
		M11[j].resize(nlist1);
		for(auto i = 0u; i < nlist1; i++){
			M11[j][i] = covar[list1[j]][list1[i]];
		}
	}
	
	M22.resize(nlist2);
	for(auto j = 0u; j < nlist2; j++){
		M22[j].resize(nlist2);
		for(auto i = 0u; i < nlist2; i++){
			M22[j][i] = covar[list2[j]][list2[i]];
		}
	}
	
	M12.resize(nlist1);
	for(auto j = 0u; j < nlist1; j++){
		M12[j].resize(nlist2);
		for(auto i = 0u; i < nlist2; i++){
			M12[j][i] = covar[list1[j]][list2[i]];
		}
	}
	
	M21.resize(nlist2);
	for(auto j = 0u; j < nlist2; j++){
		M21[j].resize(nlist1);
		for(auto i = 0u; i < nlist1; i++){
			M21[j][i] = covar[list2[j]][list1[i]];
		}
	}
	
	auto M22_inv = invert_matrix(M22);
	
	auto mu = set_p(sett);
		
	vector <double> vec_dif(nlist2);
	for(auto i = 0u; i < nlist2; i++){
		auto ca = list2[i];
		vec_dif[i] = mu_new[ca] - mu[ca];
	}
	
	auto Mmult = matrix_mult(M12,M22_inv);
	auto vec_new = matrix_mult(Mmult,vec_dif);
	
	for(auto i = 0u; i < nlist1; i++){
		auto ca = list1[i];
		mu_new[ca] = mu[ca] + vec_new[i];
	}
	
	auto M11_shift = matrix_mult(Mmult,M21);
	for(auto j = 0u; j < nlist1; j++){
		for(auto i = 0u; i < nlist1; i++){
			covar[list1[j]][list1[i]] -= M11_shift[j][i];
		}
	}
	
	for(auto j = 0u; j < nlist2; j++){
		for(auto i = 0u; i < nlist2; i++){
			covar[list2[j]][list2[i]] = 0;
		}
	}
	
	for(auto j = 0u; j < nlist1; j++){
		for(auto i = 0u; i < nlist2; i++){
			covar[list1[j]][list2[i]] = 0;
			covar[list2[i]][list1[j]] = 0;
		}
	}
	
	correct_mu_covar(mu_new,covar);
	valid_matrix(covar,mu_new);
	
	put_p(sett,mu_new);
	
	auto sum = 0.0;
	for(auto j = 0u; j < nlist2; j++){
		for(auto i = 0u; i < nlist2; i++){
			sum += vec_dif[j]*M22_inv[j][i]*vec_dif[i];
		}
	}
	
	return -0.5*nlist2*log(2*M_PI) - 0.5*determinant_fast(M22) - 0.5*sum;
}


/// Corrects mean and covariance to ensure correct 
void State::correct_mu_covar(vector <double> &mu, vector < vector <double> > &covar) const
{
	for(auto i = 0u; i < ncomp_with_S; i++){
		auto sum = 0.0; for(auto lj : comp_with_S_list[i]) sum += mu[lj];
		
		auto ca = comp_with_S[i];
		const auto &cap = comp_approx[ca];
		mu[ca] = total_pop[cap.c][cap.dp] - sum;
	}
	
	for(auto j = 0u; j < ncomp_with_S; j++){
		auto cj = comp_with_S[j];
		for(auto i = 0u; i < ncomp_with_S; i++){
			auto ci = comp_with_S[i];
					
			auto val = 0.0;
			for(auto lj : comp_with_S_list[j]){
				for(auto li : comp_with_S_list[i]) val += covar[lj][li];
			}
			covar[cj][ci] = val;
		}
		
		for(auto i = 0u; i < ncomp_without_S; i++){
			auto ci = comp_without_S[i];
			
			auto val = 0.0; for(auto lj : comp_with_S_list[j]) val += covar[lj][ci];
			
			covar[cj][ci] = -val; covar[ci][cj] = -val;
		}
	}
}


/// Removes the compartmental states related to susceptible compartment
vector <double> State::remove_S(const vector <double> &vec) const
{
	timer[TIME_ADD_REMOVE_S].start();
	vector <double> vec_new(ncomp_without_S);
	
	for(auto i = 0u; i < ncomp_without_S; i++) vec_new[i] = vec[comp_without_S[i]];
	timer[TIME_ADD_REMOVE_S].stop();
	return vec_new;
}


/// Adds the compartmental states related to susceptible compartment
vector <double> State::add_S(const vector <double> &vec) const
{
	timer[TIME_ADD_REMOVE_S].start();
	vector <double> vec_new(ncomp_approx);
	
	auto pop = total_pop;
	
	for(auto i = 0u; i < ncomp_without_S; i++){
		auto ca = comp_without_S[i];
		auto val = vec[i];
		vec_new[ca] = val;
		const auto &cap = comp_approx[ca]; 
		pop[cap.c][cap.dp] -= val;
	}
	
	for(auto i = 0u; i < ncomp_with_S; i++){
		auto ca = comp_with_S[i];
		const auto &cap = comp_approx[ca]; 
		vec_new[ca] = pop[cap.c][cap.dp];
	}
	
	timer[TIME_ADD_REMOVE_S].stop();
	
	return vec_new;
}


/// Removes states associated with S
vector < vector <double> > State::remove_S(const vector < vector <double> > &M) const
{
	timer[TIME_ADD_REMOVE_S].start();
	
	vector < vector <double> > M_new;
	M_new.resize(ncomp_without_S);
	for(auto j = 0u; j < ncomp_without_S; j++){
		M_new[j].resize(ncomp_without_S);
		for(auto i = 0u; i < ncomp_without_S; i++){
			M_new[j][i] = M[comp_without_S[j]][comp_without_S[i]];
		}
	}
	
	timer[TIME_ADD_REMOVE_S].stop();
	
	return M_new;
}


/// Adds states associated with S
vector < vector <double> > State::add_S(const vector < vector <double> > &M) const
{
	timer[TIME_ADD_REMOVE_S].start();
	
	vector < vector <double> > M_new;
	M_new.resize(ncomp_approx); for(auto j = 0u; j < ncomp_approx; j++) M_new[j].resize(ncomp_approx);
	
	for(auto j = 0u; j < ncomp_without_S; j++){
		for(auto i = 0u; i < ncomp_without_S; i++){
		  M_new[comp_without_S[j]][comp_without_S[i]] = M[j][i];
		}
	}
	
	for(auto j = 0u; j < ncomp_with_S; j++){
		auto cj = comp_with_S[j];
		for(auto i = 0u; i < ncomp_with_S; i++){
			auto ci = comp_with_S[i];
					
			auto val = 0.0;
			for(auto lj : comp_with_S_list[j]){
				for(auto li : comp_with_S_list[i]) val += M_new[lj][li];
			}
			M_new[cj][ci] = val;
		}
		
		for(auto i = 0u; i < ncomp_without_S; i++){
			auto ci = comp_without_S[i];
			
			auto val = 0.0; for(auto lj : comp_with_S_list[j]) val += M_new[lj][ci];
			
			M_new[cj][ci] = -val; M_new[ci][cj] = -val;
		}
	}
	
	timer[TIME_ADD_REMOVE_S].stop();
	
	return M_new;
}


/// Calculates the change in the covariance matrix
void State::covar_change_fast(vector < vector <double> > &covar_new, const vector < vector <double> > &covar, const vector <double> &mu, const Accuracy ac, const bool forward) const
{
	for(auto ta = 0u; ta < ntrans_approx; ta++){
		const auto &tap = trans_approx[ta];
		auto cci = tap.ca_from;
		auto ccf = tap.ca_to;
		auto val = mu[ta];
	
		covar_new[ccf][ccf] += val;
		covar_new[cci][cci] += val;
		covar_new[cci][ccf] -= val;
		covar_new[ccf][cci] -= val;
	}

	if(forward == true){
		if(true){
			set_mat(covar,ac);
			covar_matrix_mult_SIMD(dtransmean_dp,ncomp_approx,ac);
			covar_matrix_mult_SIMD(T,ncomp_approx,ac);
			covar_add_mat(covar_new,ac);
		}		
		else{
			auto prod = matrix_mult_sparse(dtransmean_dp,covar);
			auto prod2 = matrix_mult_sparse(T,prod);

			for(auto c = 0u; c < ncomp_approx; c++){
				for(auto c2 = 0u; c2 < ncomp_approx; c2++) covar_new[c][c2] += prod2[c][c2] + prod2[c2][c];
			}
		}
	}
}


/// Calculates the change in the covariance matrix
void State::covar_change(vector < vector <double> > &covar_new, const vector < vector <double> > &covar, const vector <double> &mu) const
{
	for(auto c = 0u; c < ncomp_approx; c++){
		for(auto cc = 0u; cc < ncomp_approx; cc++){
			for(auto j = 0u; j < ntrans_approx; j++){
				covar_new[c][cc] += T[c][j]*T[cc][j]*mu[j];
			}
			
			for(auto j = 0u; j < ntrans_approx; j++){
				for(auto d = 0u; d < ncomp_approx; d++){
					covar_new[c][cc] += (T[c][j]*covar[d][cc] + T[cc][j]*covar[d][c])*dtransmean_dp[j][d];
				}
			}
		}
	}
}

/// Ensures that a matrix is valid (diagonals are positive and odd diagonals contrained
void State::valid_matrix(vector < vector <double> > &M, const vector <double> mean) const
{
	auto diag_min = 0.001, cormax = 0.8;
	
	auto n = M.size(); if(mean.size() != n) emsgEC("State Approx",19);
	
	vector <double> sigma(n);
	for(auto c = 0u; c < n; c++){
		if(M[c][c] > mean[c]*mean[c]) M[c][c] = mean[c]*mean[c];
		if(M[c][c] < diag_min) M[c][c] = diag_min;
		sigma[c] = sqrt(M[c][c]);
	}
	
	for(auto c = 0u; c < n; c++){
		auto sig1 = sigma[c];
		for(auto cc = c+1; cc < n; cc++){
			auto sigmult = cormax*sig1*sigma[cc];
			auto val = M[c][cc];
			if(val > sigmult){ M[c][cc] = sigmult;  M[cc][c] = sigmult;}
			if(val < -sigmult){ M[c][cc] = -sigmult;  M[cc][c] = -sigmult;}
		}
	}
}
	
	
/// Set a vector giving the mean number of transitions
vector <double> State::set_mu(const unsigned int sett) const 
{
	vector <double> mu(ntrans_approx);
	
	for(auto ta = 0u; ta < ntrans_approx; ta++){
		const auto &tap = trans_approx[ta];
		mu[ta] = transmean[sett][tap.c][tap.tr][tap.dp];
		if(std::isnan(mu[ta])) emsgEC("State approx",14);
	}
	
	return mu;
}


/// Set a the vector giving the populations in different compartments
vector <double> State::set_p(const unsigned int sett) const 
{
	vector <double> p(ncomp_approx);
	for(auto ca = 0u; ca < ncomp_approx; ca++){   
		const auto &cap = comp_approx[ca];
		p[ca] = pop[sett][cap.c][cap.co][cap.dp];
	}
	
	return p;
}
	
	
/// Set a the vector giving the populations in different compartments
void State::put_p(const unsigned int sett, const vector <double> &p)
{
	for(auto ca = 0u; ca < ncomp_approx; ca++){   
		const auto &cap = comp_approx[ca];
		pop[sett][cap.c][cap.co][cap.dp] = p[ca];
	}
	
	if(sett < details.ndivision) set_I_from_pop(sett,false);
}
	

/// Sets the gradients in the number of transtions at a given time sett and within a given area c
void State::set_transmean_gradients(const unsigned int sett)
{
	auto dt = double(details.period)/details.ndivision;
	
	timer[TIME_TRANSMEAN].start();
	
	auto tr = model.infection_trans;
	auto from = model.trans[tr].from;
	
	auto st = 0; // Currently this only works for one strain
	
	for(auto c = 0u; c < data.narea; c++){
		auto NMI = get_NMI(sett,st,c);

		auto be = beta[st][c][sett];
		be *= areafactor[sett/details.division_per_time][c];
		
		auto efoi_info = model.efoispline_info[model.efoi_spl_ref[st][c]];
		auto et = disc_spline[efoi_info.spline_ref][sett];
		
		vector <double> eta_age(data.nage);
		const auto &agedist = efoi_info.efoi_agedist;

		for(auto a = 0u; a < data.nage; a++) eta_age[a] = et*agedist[a];
		
		for(auto dp = 0u; dp < data.ndemocatpos; dp++){	
			auto a = data.democatpos[dp][0];
			auto ca = comp_approx_ref[c][from][dp];
			if(ca != UNSET){
				dtransmean_dp[trans_approx_ref[c][tr][dp]][ca] = dt*susceptibility[dp]*(be*NMI[a] + eta_age[a]);
			}
		}
		
		unsigned int c2;
		double va;
		
		for(auto co = 0u; co < model.comp.size(); co++){
			auto inf = paramval[comp[co].infectivity_param];  // infectivity
			if(inf != 0){
				if(co == model.start_compartment) emsg("The start compartment cannot be infectious");
				
				auto diag = data.genQ.M.diag[c];
				auto &to = data.genQ.M.to[c];
				auto &val = data.genQ.M.val[c];
			
				auto jmax = to.size();
				for(auto j = 0u; j <= jmax; j++){
					if(j == jmax){ c2 = c; va = diag;}
					else{ c2 = to[j]; va = val[j];}
					
					auto be = beta[st][c2][sett];
					be *= areafactor[sett/details.division_per_time][c2];
		
					for(auto dp = 0u; dp < data.ndemocatpos; dp++){	
						auto a = data.democatpos[dp][0];
						for(auto dp2 = 0u; dp2 < data.ndemocatpos; dp2++){	
							auto a2 = data.democatpos[dp2][0];
					
							auto ta = trans_approx_ref[c2][tr][dp2];
							auto ca = comp_approx_ref[c][co][dp];
							if(ca == UNSET) emsgEC("State approx",11);
							
							auto valu = dt*susceptibility[dp2]*be*Ntime[sett][a][a2]*va*inf;
							
							dtransmean_dp[ta][ca] = valu*pop[sett][c2][from][dp2];
						}
					}
				}
			}
		}
	}

	for(auto tr = 0u; tr < model.trans.size(); tr++){                       // Non-infection transitions
		if(model.trans[tr].inf == TRANS_NOTINFECTION){
			//auto from = model.trans[tr].from;
			for(auto c = 0u; c < data.narea; c++){
				for(auto dp = 0u; dp < data.ndemocatpos; dp++){	
					auto ta = trans_approx_ref[c][tr][dp];
					dtransmean_dp[ta][trans_approx[ta].ca_from] = dt*transrate[tr][dp];
				}
			}
		}
	}
	
	timer[TIME_TRANSMEAN].stop();
}


/// Initialises parameters used when calculating the approximate likelihood
void State::initialise_approx()
{
	auto narea = data.narea;
	auto ndem = data.ndemocatpos;
	auto ncomp = model.comp.size();   
	auto ntrans = model.trans.size();   
		
	if(comp_approx.size() != 0) emsgEC("State approx",10);
	comp_approx_ref.resize(narea);                         // This maps from co c and dp to a unqiue compartment number 
	for(auto c = 0u; c < narea; c++){
		comp_approx_ref[c].resize(ncomp);
		for(auto co = 0u; co < ncomp; co++){
			comp_approx_ref[c][co].resize(ndem);
			for(auto dp = 0u; dp < ndem; dp++){	
				comp_approx_ref[c][co][dp] = comp_approx.size();
					
				CompApprox ca; ca.c = c; ca.co = co; ca.dp = dp;
				comp_approx.push_back(ca);
			}
		}
	}
	ncomp_approx = comp_approx.size();

	trans_approx_ref.resize(narea);                         // This maps from co c and dp to a unqiue compartment number 
	for(auto c = 0u; c < narea; c++){
		trans_approx_ref[c].resize(ntrans);
		for(auto tr = 0u; tr < ntrans; tr++){
			trans_approx_ref[c][tr].resize(ndem);
			for(auto dp = 0u; dp < ndem; dp++){	
				trans_approx_ref[c][tr][dp] = trans_approx.size();
				
				TransApprox ta; ta.c = c; ta.tr = tr; ta.dp = dp; 
				ta.ca_from = comp_approx_ref[c][model.trans[tr].from][dp];
				ta.ca_to = comp_approx_ref[c][model.trans[tr].to][dp];
				
				trans_approx.push_back(ta);
			}
		}
	}
	ntrans_approx = trans_approx.size();

	if(false){
		for(auto c = 0u; c < narea; c++){
			for(auto co = 0u; co < ncomp; co++){
				for(auto dp = 0u; dp < ndem; dp++){	
					cout << c << " "<< co << " " << dp << " " << comp_approx_ref[c][co][dp]  << "comp" << endl;
				}
			}
		}	
		
		for(auto c = 0u; c < narea; c++){
			for(auto tr = 0u; tr < ntrans; tr++){
				for(auto dp = 0u; dp < ndem; dp++){	
					cout << c << " "<< tr << " " << dp << " " << trans_approx_ref[c][tr][dp]  << "trans" << endl;
				}
			}
		}	
		emsg("test 1");
	}
	
	dtransmean_dp.resize(ntrans_approx);
	for(auto ta = 0u; ta < ntrans_approx; ta++){ 
		dtransmean_dp[ta].resize(ncomp_approx);
		for(auto ca = 0u; ca < ncomp_approx; ca++) dtransmean_dp[ta][ca] = 0;
	}
	
	T.resize(ncomp_approx);
	for(auto ca = 0u; ca < ncomp_approx; ca++){
		T[ca].resize(ntrans_approx);
		for(auto ta = 0u; ta < ntrans_approx; ta++){
			T[ca][ta] = 0;
		}
	}
	
	for(auto ta = 0u; ta < ntrans_approx; ta++){
		auto cci = trans_approx[ta].ca_from;
		if(cci != UNSET) T[cci][ta] = -1;
		
		auto ccf = trans_approx[ta].ca_to;
		if(ccf != UNSET) T[ccf][ta] = 1;
	}
	
	if(details.mode == ML_INF){
		vector < vector <double> > T_without_S;
		for(auto ca = 0u; ca < ncomp_approx; ca++){
			const auto &cap = comp_approx[ca];
			if(cap.co != model.start_compartment) T_without_S.push_back(T[ca]);
		}
		if(T_without_S.size() != T_without_S[0].size()) emsgEC("State approx",43); 
	
		Tinv_without_S = invert_matrix_permute(T_without_S);
	}
	
	if(false){
		print_matrix("T",T);
		print_matrix("Tinv_without_S",Tinv_without_S);
		emsg("Matrices");
	}

	for(auto ca = 0u; ca < ncomp_approx; ca++){
		const auto &cap = comp_approx[ca];
		if(cap.co != model.start_compartment) comp_without_S.push_back(ca);
		else comp_with_S.push_back(ca);
	}
	ncomp_without_S = comp_without_S.size();
	ncomp_with_S = comp_with_S.size();

	comp_with_S_list.resize(ncomp_with_S);
	for(auto i = 0u; i < ncomp_with_S; i++){
		auto cai = comp_with_S[i];
		const auto &capi = comp_approx[cai];
		for(auto j = 0u; j < ncomp_without_S; j++){
			auto caj = comp_without_S[j]; 
			const auto &capj = comp_approx[caj];
			if(capi.c == capj.c && capi.dp == capj.dp) comp_with_S_list[i].push_back(caj);
		}
	}

	if(false){
		for(auto i = 0u; i < ncomp_without_S; i++) cout << comp_without_S[i] << " woS" << endl;
		
		for(auto i = 0u; i < ncomp_with_S; i++){
			cout << comp_with_S[i] << ":";
			for(auto val : comp_with_S_list[i]) cout << val << ",";
			cout << "list" << endl;
		}
		emsg("woS");
	}
	
	total_pop.resize(narea);
	for(auto c = 0u; c < narea; c++){
		total_pop[c].resize(ndem);
		for(auto dp = 0u; dp < ndem; dp++){	
			auto sum = 0.0;
			for(auto co = 0u; co < ncomp; co++) sum += data.area[c].pop_init[co][dp];
			total_pop[c][dp] = sum;		
		}
	}
	
	/// Initialises the conpartment used for each of the observations
	for(auto o = 0u; o < data.nobs; o++){
		const auto &ob = data.obs[o];
		const auto &dt = data.datatable[ob.datatable];
		
		vector <unsigned int> complist;
		
		switch(dt.type){
			case POP: complist = dt.complist; break;
			
			case TRANS:
				for(auto tr : dt.translist){
					auto comp_after = model.trans_comp_divide(tr,1);
					for(auto c : comp_after){
						if(find_in(complist,c) == UNSET) complist.push_back(c);
						else{
							emsgroot("In datatable the observation '"+dt.observation+"' cannot contain multiple sequential transitions");
						}
					}
				}
				break;
				
			default: emsgroot("Only 'population' and 'transition' datatables can be used"); break;
		}
		
		vector <unsigned int> clist;   // Makes a list of all compartments 
		for(auto co : complist){
			if(co == model.start_compartment) emsgroot("Measurements can't be made on the susceptible compartment");
			for(auto c : ob.area){
				for(auto dp : ob.dp_sel){
					auto cc = comp_approx_ref[c][co][dp];
					clist.push_back(cc);
				}
			}
		}
		obs_clist.push_back(clist);	
		
		vector <unsigned int> clist_without_S;
		for(auto ca : clist){
			auto val = find_in(comp_without_S,ca);
			if(val != UNSET) clist_without_S.push_back(val);
		}
		obs_clist_without_S.push_back(clist_without_S);	
	}
	
	matrix_operation_initialise(ncomp_approx);
}
