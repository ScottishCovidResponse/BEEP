TAKEN FROM MBP


/// Performs MVN updates on type I MBPs and parameters and type II MBP with simple fraction
unsigned int Mbp::mvn_update(vector <Particle> &part, const vector <ParamSample> &param_samp, ParamProp &paramprop, const double cor_max)
{
	
}



unsigned int Mbp::mcmc_updates(vector <Particle> &part, const vector <ParamSample> &param_samp, double EFcut_, double invT_, ParamUpdate pup_, ParamProp &paramprop, const double cor_max)
{
	EFcut = EFcut_;
	invT = invT_;
	pup = pup_;

	timer[TIME_UPDATE].start();    
	switch(update_type){
		case MVN_UPDATE: return mvn_update(part,param_samp,paramprop,cor_max);
		case UNIVARIATE_UPDATE: return univariate_update(part,param_samp,paramprop,cor_max);	
		case COMPLEX_UPDATE: return complex_update(part,param_samp,paramprop,cor_max);
	}
	timer[TIME_UPDATE].stop();    
	return -1;
}




/// Performs univariate updates on type I MBPs and parameters and type II MBP with simple fraction
unsigned int Mbp::univariate_update(vector <Particle> &part, const vector <ParamSample> &param_samp, ParamProp &paramprop, const double cor_max)
{
	auto nprop = 0u;   
	
	auto psamp_before = mpi.gather_psamp(part);                          // Stores parameter samples

	vector <double> cor;                                                 // Stores the correlation with initial samples
	paramprop.zero_ntr_nac_state(); 
	auto N = paramprop.mvn.size()/2;
	
	auto slice_min = 0u;
	do{
		paramprop.zero_ntr_nac();	
	
		for(auto &mv : paramprop.mvn) mv.setup(param_samp);
			
		auto slice_max = slice_min+10;
		if(slice_max > paramprop.slicetime.size()) slice_max = paramprop.slicetime.size();
		 
		for(auto &pa : part){
			initial->initialise_from_particle(pa);                           // Initialises mbp updates from a particle

			if(checkon == true) initial->check(0);

			timer[TIME_MCMCPROP].start();                                    // Type I MBPs
			for(auto loop = 0u; loop < N; loop++){ mvn_proposal(paramprop.mvn[loop]); nprop++;}
			timer[TIME_MCMCPROP].stop();

			if(checkon == true) initial->check(1);
			
			timer[TIME_PARAMPROP].start();                                   // Parameter proposals
			auto Li = initial->likelihood();
			copy_initial_to_propse();
			
			for(auto loop = N; loop < 2*N; loop++){ parameter_proposal(paramprop.mvn[loop],Li); nprop++;}
			timer[TIME_PARAMPROP].stop();
			
			if(checkon == true) initial->check(2);
		
			timer[TIME_STATEPROP].start();                                   // Type II MBPs
			// TO DO
			
			//for(auto loop = 0u; loop < paramprop.fixedtree.size(); loop++){
				//mbp_fixedtree(paramprop.fixedtree[loop]); nprop++;
			//}
					
			//if(mpi.core == 0) cout <<  slice_min << " "<<  slice_max <<  " " <<  paramprop.slicetime.size() <<  " sl\n";
			 
			for(auto loop = slice_min; loop < slice_max; loop++){
				mbp_slicetime(paramprop.slicetime[loop]); nprop++;
			}
			
			//for(auto loop = 0u; loop < 3; loop++){ mbp_fixedtree(paramprop.fixedtree[0]); nprop++;}
			timer[TIME_STATEPROP].stop();     
			
			if(checkon == true) initial->check(3);
				
			pa = initial->create_particle(pa.run);   
			if(pa.EF > EFcut) emsgEC("MBP",101);
		}
		
		if(pup == COMBINE_UPDATE) paramprop.combine_update_proposals();
	
		//if(pup != NO_UPDATE) paramprop.update_fixed_splice_proposals(pup);
		
		slice_min = slice_max;
		if(slice_min == paramprop.slicetime.size()){
			if(pup != NO_UPDATE) paramprop.update_fixed_splice_proposals(pup);
			slice_min = 0;
		}
		
		timer[TIME_WAIT].start();
		auto psamp_after = mpi.gather_psamp(part);
		cor = output.get_correlation(psamp_before,psamp_after);
		timer[TIME_WAIT].stop();
	}while(vec_max(cor) > cor_max);                                       // Iterates until correlation less than threshold
	
	return nprop/part.size();
}

/// Performs a complex update
unsigned int Mbp::complex_update(vector <Particle> &part, const vector <ParamSample> &param_samp, ParamProp &paramprop, const double cor_max)
{
	auto nprop = 0u;   
	
	paramprop.zero_ntr_nac_state();                         // Updates state with fixed parameters
	
	auto state_prop_list = paramprop.get_state_proposal_list(); 
				
	for(auto &pa : part) update_particle(pa,state_prop_list,paramprop);
	
	if(pup != NO_UPDATE) paramprop.update_fixed_splice_proposals(pup);

	auto psamp_before = mpi.gather_psamp(part);             // Updates parameters
	
	for(auto &mv : paramprop.mvn){ mv.prop_vec.clear(); mv.after_prop_vec.clear();}
	
	auto ninteration = 1u;
	
	do{
		auto prop_list = paramprop.get_proposal_list(param_samp); 
		nprop += prop_list.size();
	
		paramprop.zero_ntr_nac();
	
		for(auto &pa : part) update_particle(pa,prop_list,paramprop);

		auto psamp_after = mpi.gather_psamp(part);
		
		auto cor = output.get_correlation(psamp_before,psamp_after);
	
		if(pup == COMBINE_UPDATE) paramprop.combine_update_proposals();
		
		if(vec_max(cor) < cor_max){
			paramprop.update_num_updates(cor,ninteration);
			break;
		}
		
		ninteration++;
	}while(true);
	
	return nprop;
}

PARAM PROP:

ParamProp::ParamProp(const Details &details, const Data &data, const Model &model, const Output &output, Mpi &mpi) : details(details), data(data), model(model), output(output), mpi(mpi)
{
	if(sim_only == true){
		add_single();
	}
	else{
	
		switch(update_type){
			// Performs MVN updates on type I MBPs and parameters and type II MBP with simple fraction
			case MVN_UPDATE: mvn_update_init(); break;
			
			// Performs univarte updates on type I MBPs and parameters and type II MBP with simple fraction
			case UNIVARIATE_UPDATE: univariate_update_init(); break;
			
			// Performs more complex updates
			case COMPLEX_UPDATE: complex_update_init(); break;
		}
	}
	
	zero_ntr_nac();
	zero_ntr_nac_state();
}



/// Initialises univariate particle updates
void ParamProp::univariate_update_init()
{
	for(auto th = 0u; th < model.param.size(); th++){
		auto type = model.param[th].type;
		if(model.param[th].priortype != FIXED_PRIOR){
			vector <unsigned int> vec; vec.push_back(th);
			MVN mv("Univariate: "+model.param[th].name,vec,1,type,SINGLE);
			mvn.push_back(mv);
		}
	}
	
	for(auto th = 0u; th < model.param.size(); th++){
		if(model.param[th].priortype != FIXED_PRIOR && (model.param[th].type != INF_PARAM)){
			vector <unsigned int> vec; vec.push_back(th);
			MVN mv("Param prop: "+model.param[th].name,vec,0.1,PARAM_PROP,SINGLE);
			mvn.push_back(mv);
		}
	}
	
	SliceTime st; st.sett_i = 0; st.sett_f = details.ndivision; st.sim_frac = 1;
	slicetime.push_back(st);
}


/// Initialises complex particle updates
void ParamProp::complex_update_init()
{
	if(details.mode == PMCMC_INF || false){
		MVN mv("MBP type I MVN",model.param_not_fixed,0.3,ALL_PARAM,MULTIPLE);
		mvn.push_back(mv);
	
		MVN mv_param("Parameter MVN",model.param_not_fixed,0.3,PARAM_PROP,MULTIPLE);
		mvn.push_back(mv_param);
	}
	else{				
		if(details.mcmc_update.full_mvn == true){
			MVN mv("Full MVN",model.param_not_fixed,0.3,ALL_PARAM,MULTIPLE);
			mvn.push_back(mv);
		}
		
		if(details.mcmc_update.mvn == true){
			add_mvn("Distribution value",DISTVAL_PARAM,0.5);
			add_mvn("Branching probability",BRANCHPROB_PARAM,0.5);
			add_mvn("Infectivity",INF_PARAM,0.5);
			add_mvn("Regional effect",RE_PARAM,0.5);
			add_mvn("Sigma",SIGMA_PARAM,0.3);
			add_mvn("Observation parameters",OBS_PARAM,0.5);
			add_mvn("Covariates",COVAR_PARAM,0.5);
			add_mvn("R / efoi",R_EFOI_PARAM,0.3);
			add_mvn("Geomix",GEOMIX_PARAM,0.5);
			add_mvn("Modify param",MODIFY_PARAM,0.5);
			add_mvn("Susceptibility",SUSCEPTIBILITY_PARAM,0.5);
			add_mvn("Demo Specific",DEMO_SPECIFIC_PARAM,0.5);
		}
		
		if(details.mcmc_update.single == true) add_single();

		if(details.mcmc_update.dist_R_joint == true) add_dist_R_joint();

		if(details.mode != PMCMC_INF) add_param_prop();

		//if(details.mcmc_update.demo_spec == true) add_demographic_specific();
		
		if(details.mcmc_update.mean_time == true) mean_time_init();
		
		if(details.mcmc_update.neighbour == true) neighbour_init(); 
		
		if(details.mcmc_update.joint == true) joint_init();
		
		covar_area_init();
	}
	
	if(details.mode == ABC_MBP || details.mode == MC3_INF || details.mode == MCMC_MBP || details.mode == PAS_INF){
		FixedTree ft; ft.n = 0; ft.sim_frac = 1;
		fixedtree.push_back(ft);
	
		SliceTime st; st.sett_i = 0; st.sett_f = details.ndivision; st.sim_frac = 1;
		slicetime.push_back(st);
	}
}
	
	
FROM PARAM_PROP

/*
/// Updates the number of different types of proposals using the correlation across generations
void ParamProp::update_num_updates(const vector <double> &cor, unsigned int ninteration)
{
	for(auto th = 0u; th < model.param.size(); th++){
		auto co = cor[th];
		if(co != UNSET){
			if(co < 0.1) co = 0.1;
			auto co_per_it = pow(co,1.0/ninteration);
		
			//auto fac = 0.95;
			//if(co_per_it > cor_update_num) fac = pow((1-cor_update_num)/(1-co_per_it),0.33);
			
			auto fac = pow((1-cor_update_num)/(1-co_per_it),0.33);
			if(fac < 0.9) fac = 0.9;
			
			
			for(auto &mv : mvn){
				if(find_in(mv.var,th) != UNSET) mv.num_updates *= fac;
			}
			
			for(auto &mt : mean_time){
				if(find_in(mt.param_mean,th) != UNSET || find_in(mt.param_mean_rev,th) != UNSET) mt.num_updates *= fac;
			}

			for(auto &ne : neighbour){
				if(ne.param1 == th || ne.param2 == th) ne.num_updates *= fac;
			}
			
			//for(auto &jo : joint){
			//if(find_in(jo.var_list,th) != UNSET) jo.num_updates *= fac;
			//}
			
			for(auto &ca : covar_area){ 
				if(model.covariate_param[ca.covar_ref] == th) ca.num_updates *= fac;
			}
		}
	}
	
	
	for(auto &mv : mvn){ 	                        // Increases if smaller than num_updates_min
		if(mv.num_updates < num_updates_min) mv.num_updates = num_updates_min;
		if(mv.num_updates > num_updates_max) mv.num_updates = num_updates_max;
	}
				
	for(auto &mt : mean_time){
		if(mt.num_updates < num_updates_min) mt.num_updates = num_updates_min;
		if(mt.num_updates > num_updates_max) mt.num_updates = num_updates_max;
	}

	for(auto &ne : neighbour){
		if(ne.num_updates < num_updates_min) ne.num_updates = num_updates_min;
		if(ne.num_updates > num_updates_max) ne.num_updates = num_updates_max;
	}
				
	for(auto &jo : joint){
		if(jo.num_updates < num_updates_min) jo.num_updates = num_updates_min;
		if(jo.num_updates > num_updates_max) jo.num_updates = num_updates_max;
	}
	
	for(auto &ca : covar_area){ 
		if(ca.num_updates < num_updates_min) ca.num_updates = num_updates_min;
		if(ca.num_updates > num_updates_max) ca.num_updates = num_updates_max;
	}
}
*/

