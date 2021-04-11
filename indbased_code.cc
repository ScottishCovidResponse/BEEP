/// This contains code for indvidual-based model-based proposals
/// These are not currently implemented into BEEPmbp

Status Mbp::mbp(const vector<double> &paramv, InfUpdate inf_update)
{	
	if(model.inbounds(paramv) == false) return FAIL;            // Checks parameters are within the prior bounds

	propose->set_param(paramv);                                 // Sets quantities derived from parameters (if not possible fails)
	
	timer[TIME_MBP].start();
	
	timer[TIME_MBPINIT].start();
	mbp_initialise();                                           // Prepares for the proposal
	timer[TIME_MBPINIT].stop();

	Status stat;
	switch(modeltype){
		case IND_MODEL:                                           // This performs a MBP for an individual-based model
			stat = mbp_ind_model(mbpsim_type,mbpsim_type_num,fixed_frac); 
			break; 
		case POP_MODEL:                                           // This performs a MBP for a population model
			stat = mbp_pop_model(mbpsim_type,mbpsim_type_num,fixed_frac,inf_update);
			break; 
	}
	
	timer[TIME_MBP].stop();

	return stat;
}


/// This performs a MBP for an individual-based model
Status Mbp::mbp_ind_model(MbpSimType mbpsim_type, unsigned int mbpsim_type_num, double fixed_frac)
{
	auto n = 0u;                                          // Indexes infections in the initial state 
	auto t = 0.0;                                         // Current time
	auto step = (unsigned int)(details.ndivision/10.0); if(step == 0) step = 1;
	for(auto sett = 0u; sett < details.ndivision; sett++){// Goes across discrete timesteps
		get_simu_or_mbp(sett,mbpsim_type,mbpsim_type_num,fixed_frac);		
		
		propose->set_Imap_using_dI(sett,initial,dImap,dIdiag);     // Calculates proposed Imap from initial Imap plus dI

		construct_infection_sampler(sett,initial->Imap[sett],initial->Idiag[sett],propose->Imap[sett],propose->Idiag[sett]); 

		double tmax = details.division_time[sett+1];
		do{
			auto tini = initial->get_infection_time(n);     	// Gets time of next infection event in initial state
			auto tinf = t + exp_sample(new_infection_rate[0][0]);// Gets time of new added infection event

			if(tini >= tmax && tinf >= tmax){ t = tmax; break;}
			
			if(tinf < tini){                                  // A new infection appears in the proposed state
				t = tinf;
				add_infection_in_area(area_of_next_infection(),t);	
			}
			else{                                             // An event on initial sequence is considered to be copied or not
				t = tini;
				copy_event_or_not(n);
				n++;
			}
	
			if(propose->infev.size() >= model.maximum_infected){// If the number of infections exceeds the prior limit then fails
				reset_susceptible_lists(); return FAIL;
			}
		}while(true);
		
		update_dImap(initial->transev[sett],propose->transev[sett]);// Based on infections updates dImap
		if(checkon == true) check(t,sett);
	}
	reset_susceptible_lists();
		
	return SUCCESS;
}

		
/// Decides to copy infection event n from the initial to proposed states or not
void Mbp::copy_event_or_not(unsigned int n)
{
	auto i = initial->infev[n].ind;

	if(susceptible_status[i] == BOTH_SUSCEPTIBLE){     // The copy can only be done if individual in propose is suscepticle
		auto c = data.ind[i].area;
		auto w = c*data.ndemocatpos + data.ind[i].dp;

		auto al = 0.0;	                                 // Acceptance probability
		if(simu_or_mbp[c] == MBP) al = propose->lambda[w]/initial->lambda[w];
		
		if(ran() < al){                                
			change_susceptible_status(i,NOT_SUSCEPTIBLE,1);// Copies the infection event to propose
			
			if(do_mbp_event == true) mbp_compartmental_transitions(i);
			else propose->indev[i] = initial->indev[i];
			
			propose->add_indev(i);
		}
		else change_susceptible_status(i,ONLY_PROPOSE_SUSCEPTIBLE,1);      // Does not copy the infection event propose
	}
}


/// Based on compartmental trasitions in the initial state for i
/// this uses MBPs to generates transitions for the proposal
void Mbp::mbp_compartmental_transitions(unsigned int i)
{
	const auto &evlisti = initial->indev[i];
	auto &evlistp = propose->indev[i];
	const auto &parami = initial->paramval, &paramp = propose->paramval;
	const auto &comptransprobi = initial->comptransprob, &comptransprobp = propose->comptransprob;
	
	evlistp.clear();
	
	auto ev = evlisti[0];
	auto tra = ev.trans;
	auto ti = ev.t; 
	evlistp.push_back(ev);
	
	auto dp = data.ind[i].dp;
	auto c = trans[tra].to;
	
	auto tp = ti;
	unsigned int emax = evlisti.size();
	bool switch_branch = false;
	for(auto e = 1u; e < emax; e++){                              // Goes through each event in initial
		tra = evlisti[e].trans;
		
		unsigned int kmax = comp[c].trans.size();
		if(kmax == 0) break;
		
		double dtnew;
		if(kmax > 1){                                             // Looks at switching to another branch
			auto k = 0u; while(k < kmax && tra != comp[c].trans[k]) k++;
			if(k == kmax) emsgEC("model",4);
			
			auto ratio = comptransprobp[c].prob[dp][k]/comptransprobi[c].prob[dp][k];
			if(ratio < 1){ 
				if(ran() < 1 - ratio){
					auto sum = 0.0;                                     // Selects new branch to go down
					vector <double> sumst;
					for(auto k = 0u; k < kmax; k++){
						auto dif = comptransprobp[c].prob[dp][k] - comptransprobi[c].prob[dp][k];
						if(dif > 0) sum += dif;
						sumst.push_back(sum);
					}
					
					auto z = ran()*sum; auto k = 0u; while(k < kmax && z > sumst[k]) k++;
					if(k == kmax) emsgEC("Model",5);
					tra = comp[c].trans[k];
					switch_branch = true;
					dtnew = propose->sample_duration(tra,dp);     
				}
			}	
		}
		
		if(switch_branch == false){
			auto dt = evlisti[e].t - ti;
			ti = evlisti[e].t;	
		
			switch(trans[tra].type){
			case EXP_DIST:
				{
					auto p = trans[tra].param_mean[dp]*trans[tra].mean_factor;
					dtnew = dt*paramp[p]/parami[p];
				}
				break;
			
			case GAMMA_DIST:
				emsgEC("model",6);
				break;
				
			case LOGNORM_DIST:
				{
					auto p = trans[tra].param_mean[dp], p2 = trans[tra].param_cv[dp];
					
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
	
		ev.trans = tra; ev.t = tp;
		evlistp.push_back(ev);

		c = trans[tra].to; 
		if(switch_branch == true) break;
	}
	
	
	if(switch_branch == true){  // If the branch is switched then need to simulate rest of sequence
		if(comp[c].trans.size() != 0)	propose->simulate_compartmental_transitions(i,c,tp);
	}
}


/// Constructs a fast Gillespie sampler for finding the area of the next infection to be added	
void Mbp::construct_infection_sampler(unsigned int sett, const vector < vector <double> > &Imi,  const vector< vector <double> > &Idi, const vector < vector <double> > &Imp, const vector < vector <double> > &Idp)
{
	emsg("TODO 1");
	// TO DO
cout << sett << Imi[0][0] << Idi[0][0] << Imp[0][0] << Idp[0][0];
/*
	int l = areatree.level-1;
	for(auto c = 0u; c < data.narea; c++){
		auto wmin = c*data.ndemocatpos; 
		auto wmax = wmin + data.ndemocatpos;
	
		auto faci = initial->beta[sett]*initial->areafactor[c];
		auto facp = propose->beta[sett]*propose->areafactor[c];
		
		auto sum = 0.0; 
		auto dp = 0u; 
		auto v = c*data.nage; 
		cout << faci << facp << v;
		for(auto w = wmin; w < wmax; w++){
			//auto a = data.democatpos[dp][0];
			//propose->lambda[w] = propose->susceptibility[dp]*(facp*Qmp[v+a] + propose->phi);
			
			switch(simu_or_mbp[c]){
			case SIMU: sum += (nboth_susceptible_list[w]+npropose_only_susceptible_list[w])*propose->lambda[w]; break;
			case MBP:
				//initial->lambda[w] = initial->susceptibility[dp]*(faci*Qmi[v+a] + initial->phi);
			
				auto dlambda = nboth_susceptible_list[w]*(propose->lambda[w] - initial->lambda[w]); if(dlambda < 0) dlambda = 0;
				if(std::isnan(dlambda)){ emsgEC("Mbp",400);}
				sum += dlambda + npropose_only_susceptible_list[w]*propose->lambda[w];
			}
			
			dp++;
		}
		if(std::isnan(sum)) emsgEC("Mbp",4);
		if(sum < 0) emsgEC("Mbp",4);
		new_infection_rate[l][c] = sum;
	}
	
	for(int l = areatree.level-2; l >= 0; l--){                                 
		auto cmax = lev[l].node.size();
		for(auto c = 0u; c < cmax; c++){
			auto sum = 0.0; for(const auto& ch : lev[l].node[c].child) sum += new_infection_rate[l+1][ch];
			
			new_infection_rate[l][c] = sum;
		}
	}
	*/
}
	

/// Sets up lists of individuals (those suscepticle in both initial and propose states, only propose, and not at all)
void Mbp::setup_susceptible_lists()
{	
	both_susceptible_list.clear(); both_susceptible_list.resize(data.nardp); 
	propose_only_susceptible_list.clear(); propose_only_susceptible_list.resize(data.nardp); 
	not_susceptible_list.clear(); not_susceptible_list.resize(data.nardp);

	nboth_susceptible_list.resize(data.nardp); 
	npropose_only_susceptible_list.resize(data.nardp); 
	nnot_susceptible_list.resize(data.nardp);
	
	susceptible_status.resize(data.popsize); susceptible_list_ref.resize(data.popsize); 

	for(auto c = 0u; c < data.narea; c++){
		for(auto dp = 0u; dp < data.ndemocatpos; dp++){
			auto w = c*data.ndemocatpos + dp;

			for(const auto& i : data.area[c].ind[dp]){
				susceptible_status[i] = BOTH_SUSCEPTIBLE;
				susceptible_list_ref[i] = both_susceptible_list[w].size();
				both_susceptible_list[w].push_back(i);
			}
			nboth_susceptible_list[w] = data.area[c].ind[dp].size();
			npropose_only_susceptible_list[w] = 0;
			nnot_susceptible_list[w] = 0;
		}
	}	
}


/// Makes a change in the susceptibility status of an individual
void Mbp::change_susceptible_status(unsigned int i, unsigned int st, unsigned int updateR)
{
	auto c = data.ind[i].area;
	auto w = c*data.ndemocatpos + data.ind[i].dp;
	
	auto dval = 0.0;
	if(updateR == 1){		               // Updates the event sampler for new events 
		switch(simu_or_mbp[c]){
		case SIMU:
			dval = -((nboth_susceptible_list[w]+npropose_only_susceptible_list[w])*propose->lambda[w]);
			break;
			
		case MBP:
			auto dlambda = nboth_susceptible_list[w]*(propose->lambda[w] - initial->lambda[w]); if(dlambda < 0) dlambda = 0;
			dval = -(dlambda + npropose_only_susceptible_list[w]*propose->lambda[w]);
			break;
		}
	}
	
	int l = susceptible_list_ref[i];   // Removes the einitial->infevsiting entry
	int n;
	switch(susceptible_status[i]){
  case BOTH_SUSCEPTIBLE:
		if(both_susceptible_list[w][l] != i) emsgEC("Mbp",7);
		n = both_susceptible_list[w].size();
		if(l < n-1){
			both_susceptible_list[w][l] = both_susceptible_list[w][n-1];
			susceptible_list_ref[both_susceptible_list[w][l]] = l;
		}
		both_susceptible_list[w].pop_back();
		nboth_susceptible_list[w]--;
		break;
		
	case ONLY_PROPOSE_SUSCEPTIBLE:
		if(propose_only_susceptible_list[w][l] != i) emsgEC("Mbp",8);
		n = propose_only_susceptible_list[w].size();
		if(l < n-1){
			propose_only_susceptible_list[w][l] = propose_only_susceptible_list[w][n-1];
			susceptible_list_ref[propose_only_susceptible_list[w][l]] = l;
		}
		propose_only_susceptible_list[w].pop_back();
		npropose_only_susceptible_list[w]--;
		break;
		
	case NOT_SUSCEPTIBLE:
		if(not_susceptible_list[w][l] != i) emsgEC("Mbp",9);
		n = not_susceptible_list[w].size();
		if(l < n-1){
			not_susceptible_list[w][l] = not_susceptible_list[w][n-1];
			susceptible_list_ref[not_susceptible_list[w][l]] = l;
		}
		not_susceptible_list[w].pop_back();
		nnot_susceptible_list[w]--;
		break;
	
	default: emsgEC("Mbp",10); break;
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
		susceptible_list_ref[i] = both_susceptible_list[w].size();
		both_susceptible_list[w].push_back(i);
		nboth_susceptible_list[w]++;
		break;
		
	default: emsgEC("Mbp",11); break;
	}
	
	if(updateR == 1){                              // Updates the infection sampler
		switch(simu_or_mbp[c]){
		case SIMU:
			dval += (nboth_susceptible_list[w]+npropose_only_susceptible_list[w])*propose->lambda[w];
			break;
		case MBP:
			auto dlambda = nboth_susceptible_list[w]*(propose->lambda[w] - initial->lambda[w]); if(dlambda < 0) dlambda = 0;
			dval += dlambda + npropose_only_susceptible_list[w]*propose->lambda[w];
			break;
		}
		
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
void Mbp::reset_susceptible_lists()
{
	for(auto w = 0u; w < data.nardp; w++){
		for(const auto& i : propose_only_susceptible_list[w]){
			susceptible_status[i] = BOTH_SUSCEPTIBLE;
			susceptible_list_ref[i] = both_susceptible_list[w].size();
			both_susceptible_list[w].push_back(i);
			nboth_susceptible_list[w]++;
		}
		propose_only_susceptible_list[w].clear();
		npropose_only_susceptible_list[w] = 0;
		
		for(const auto& i : not_susceptible_list[w]){
			susceptible_status[i] = BOTH_SUSCEPTIBLE;
			susceptible_list_ref[i] = both_susceptible_list[w].size();
			both_susceptible_list[w].push_back(i);
			nboth_susceptible_list[w]++;
		}
		not_susceptible_list[w].clear();
		nnot_susceptible_list[w] = 0;
	}
}


/// Updates dImap based on events that occur in the initial and proposed states
void Mbp::update_dImap(const vector <EventRef> &trei, const vector <EventRef> &trep)
{
	emsg("TODO 2");
	//timers.timembpImap -= clock();
		
	//TO DO
	auto nage = data.nage;

	for(const auto& tre : trei){
		auto i = tre.ind; 
		auto tra = initial->indev[i][tre.e].trans;
		indmap[i][tra] = 1;
	}
	
	for(const auto& tre : trep){
		auto i = tre.ind; 
		auto tra = propose->indev[i][tre.e].trans;	
		if(indmap[i][tra] == 0){
			
			auto v = data.ind[i].area*nage+data.democatpos[data.ind[i].dp][0];
			cout << v;
			/*
			auto dq = trans[tra].DQ[propose->indev[i][tre.e].timep];

			if(dq != UNSET){
				for(auto loop = 0u; loop < 2; loop++){
					auto q = model.DQ[dq].q[loop];
					if(q != UNSET){
						if(dQbuf[v][q] == 0){ dQbuflistv.push_back(v); dQbuflistq.push_back(q);}
						dQbuf[v][q] += model.DQ[dq].fac[loop];
					}
				}
			}
			*/
		}
		else indmap[i][tra] = 0;
	}	
	
	for(const auto& tre : trei){
		auto i = tre.ind; 
		auto tra = initial->indev[i][tre.e].trans;
		
		if(indmap[i][tra] == 1){
			auto v = data.ind[i].area*nage+data.democatpos[data.ind[i].dp][0];
			cout << v;
			/*
			auto dq = trans[tra].DQ[initial->indev[i][tre.e].timep];
			if(dq != UNSET){
				for(auto loop = 0u; loop < 2; loop++){
					auto q = model.DQ[dq].q[loop];
					if(q != UNSET){
						if(dQbuf[v][q] == 0){ dQbuflistv.push_back(v); dQbuflistq.push_back(q);}
						dQbuf[v][q] -= model.DQ[dq].fac[loop];
					}
				}
			}
			*/
			indmap[i][tra] = 0;
		}
		
	}
	
	if(details.mode == SIM){
		for(const auto& tre : trep){
			auto tra = propose->indev[tre.ind][tre.e].trans;
			N[trans[tra].from]--;
			N[trans[tra].to]++;
		}
	}

	nage = data.nage;
	/*
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
	*/
	
	//timers.timembpQmap += clock();
}
	
	
/// This samples which area the next new infection occurs in 
/// Starting at the top level l=0 the algorithm proceeds to finer and finer scales
unsigned int Mbp::area_of_next_infection()
{
	double sumst[4];
	
	auto l = 0u, c = 0u;                            
	auto lmax = areatree.level;
	while(l < lmax-1){
		auto& node = lev[l].node[c];
		
		auto sum = 0.0;
		auto jmax = node.child.size();
		for(auto j = 0u; j < jmax; j++){
			sum += new_infection_rate[l+1][node.child[j]];	
			sumst[j] = sum;
		}
		
		double z = ran()*sum; auto j = 0u; while(j < jmax && z > sumst[j]) j++; if(j == jmax) emsgEC("Mbp",12);
		
		c = node.child[j];
		l++;
	};

	return c;
}


/// Adds a individual infected at time t in area c
void Mbp::add_infection_in_area(unsigned int c, double t)
{
	auto dpmax = data.ndemocatpos;
	
	vector <double> sumst;
	sumst.resize(dpmax);
	
	auto i=0u;
	
	switch(simu_or_mbp[c]){
	case SIMU:
		{
			auto sum = 0.0;                                      // Selects which demographic possibility from the area
			for(auto dp = 0u; dp < dpmax; dp++){
				auto w = c*dpmax + dp;
				sum += (nboth_susceptible_list[w]+npropose_only_susceptible_list[w])*propose->lambda[w];
				sumst[dp] = sum;
			}

			double z = ran()*sum; auto dp = 0u; while(dp < dpmax && z > sumst[dp]) dp++;
			if(dp == dpmax) emsgEC("Mbp",13);
			
			auto w = c*dpmax + dp;                               // Next selects susceptible list type
			if(ran() < double(nboth_susceptible_list[w])/(nboth_susceptible_list[w] + npropose_only_susceptible_list[w])){
				auto n = both_susceptible_list[w].size(); if(n == 0) emsgEC("Mbp",14);
				i = both_susceptible_list[w][(unsigned int)(ran()*n)];
			}
			else{                                                // Only proposed state susceptible
				auto n = propose_only_susceptible_list[w].size(); if(n == 0) emsgEC("Mbp",15);
				i = propose_only_susceptible_list[w][(unsigned int)(ran()*n)];
			}
		}
		break;
	
	case MBP:
		{
			auto sum = 0.0;                                      // Selects which demographic possibility from the area
			for(auto dp = 0u; dp < dpmax; dp++){
				auto w = c*dpmax + dp;
				double dlambda = nboth_susceptible_list[w]*(propose->lambda[w] - initial->lambda[w]); if(dlambda < 0) dlambda = 0;
				sum += dlambda + npropose_only_susceptible_list[w]*propose->lambda[w];
				sumst[dp] = sum;
			}

			double z = ran()*sum; auto dp = 0u; while(dp < dpmax && z > sumst[dp]) dp++; 
			if(dp == dpmax) emsgEC("Mbp",13);
			
			auto w = c*dpmax + dp;                            // Next selects susceptible list type
			double dlambda = nboth_susceptible_list[w]*(propose->lambda[w] - initial->lambda[w]); if(dlambda < 0) dlambda = 0;

			if(ran() < dlambda/(dlambda + npropose_only_susceptible_list[w]*propose->lambda[w])){ // Both suscetible
				auto n = both_susceptible_list[w].size(); if(n == 0) emsgEC("Mbp",14);
				i = both_susceptible_list[w][(unsigned int)(ran()*n)];
			}
			else{                                             // Only proposed state susceptible
				auto n = propose_only_susceptible_list[w].size(); if(n == 0) emsgEC("Mbp",15);
				i = propose_only_susceptible_list[w][(unsigned int)(ran()*n)];
			}
		}
		break;
	}

	change_susceptible_status(i,NOT_SUSCEPTIBLE,1);       // Changes the susceptibility status
		
	propose->simulate_compartmental_transitions(i,0,t);   // Simulates compartmental transitions after infection

	propose->add_indev(i);
}

/// Clears variables ready for a MBP
void Mbp::mbp_initialise()
{
	for(auto c = 0u; c < comp.size(); c++) N[c] = 0;
	N[0] = data.popsize;
		
	for(auto inft = 0u; inft < model.ninfection_trans; inft++){	
		for(auto v = 0u; v < data.narage; v++){ dImap[inft][v] = 0; dIdiag[inft][v] = 0;} 
	}
	
	propose->clear();
	
	// Set to true if MBPs on compartmental transitions needed
	switch(modeltype){
		case IND_MODEL:
			do_mbp_event = model.do_mbp_events(initial->paramval,propose->paramval); 
			break;
			
		case POP_MODEL:
			break;
	}
}	

/// Generates a sampler for adding infected individuals into the system based on the force of infection
void Mbp::infection_sampler(const vector < vector< vector<double> > > &Imap, const vector < vector< vector<double> > > &Idiag)
{
	emsg("TODO 4");
	cout << Imap[0][0][0] << Idiag[0][0][0];
	/*
	//auto sum = 0.0;
	for(auto sett = 0u; sett < details.ndivision; sett++){
		auto phi = initial->disc_spline[model.phi_spline_ref][sett]; 
		cout << phi;
		auto beta = initial->beta[sett];
	
		for(auto c = 0u; c < data.narea; c++){
			auto fac = beta*initial->areafactor[c];
			cout << fac;
			for(auto dp = 0u; dp < data.ndemocatpos; dp++){
				auto w = c*data.ndemocatpos + dp;
				auto tot = sett*data.nardp + w;
				auto v = c*data.nage + data.democatpos[dp][0];
				cout << w << tot << v;
				
				//auto val = nboth_susceptible_list[w]*initial->susceptibility[dp]*(fac*Qmap[sett][v] + phi);
				//sum += val;
				
				//lambda[tot] = val;				
				//lambdasum[tot] = sum;
				
			}
		}
	}
	*/
}

/// Initialises the variables within Mbp
void Mbp::initialise_variables()
{
	switch(modeltype){
		case IND_MODEL:
			setup_susceptible_lists();
			
			indmap.resize(data.popsize);
			for(auto i = 0u; i < data.popsize; i++){
				indmap[i].resize(model.trans.size());
				for(auto tra = 0u; tra < model.trans.size(); tra++) indmap[i][tra] = 0;
			}
			
			emsg("TODO 5");
			//TO DO
			/*
			dQbuf.resize(data.narage);
			for(auto v = 0u; v < data.narage; v++){
				dQbuf[v].resize(data.Q.size()); for(auto q = 0u; q < data.Q.size(); q++) dQbuf[v][q] = 0;
			}
			dQbuflistv.clear(); dQbuflistq.clear();
			*/
			
			new_infection_rate.resize(areatree.level); 
			for(auto l = 0u; l < areatree.level; l++) new_infection_rate[l].resize(lev[l].node.size()); 
			
			lambda.resize(data.nsettardp); lambdasum.resize(data.nsettardp);	        // Used for event based changes
			break;
			
		case POP_MODEL:
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
			break;
	}

	dImap.resize(model.ninfection_trans); dIdiag.resize(model.ninfection_trans);                     // Initialises vectors
	for(auto inft = 0u; inft < model.ninfection_trans; inft++){	
		dImap[inft].resize(data.narage); dIdiag[inft].resize(data.narage);      
	}
		
	dIdiag.resize(data.narage);     

	N.resize(comp.size()); 
}


/// Shifts regional effect in response to changes in covaraites
void Mbp::shift_covariate_reff(vector <double>& param_prop)  
{
	for(auto j = 0u; j < model.covariate_param.size(); j++){
		auto th = model.covariate_param[j];
		auto dp = param_prop[th] - initial->paramval[th];
		
		for(auto r = 0u; r < data.region_effect.size(); r++){
			param_prop[model.region_effect_param[r]] += dp*model.region_effect_shift[j][r];
		}
	}
}



/// Calculates propose->Imap based on the initial and final sequences
void Mbp::calculate_Imap_for_propose()
{	
	for(auto inft = 0u; inft < model.ninfection_trans; inft++){
		for(auto& dIma : dImap[inft]) dIma = 0;
		for(auto& dIdia : dIdiag[inft]) dIdia = 0;
	}
	
	for(auto sett = 0u; sett < details.ndivision; sett++){
		propose->set_Imap_using_dI(sett,initial,dImap,dIdiag);
		update_dImap(initial->transev[sett],propose->transev[sett]);	
	} 
}

/// Clears variables ready for a MBP
void Mbp::mbp_initialise()
{
	for(auto c = 0u; c < comp.size(); c++) N[c] = 0;
	N[0] = data.popsize;
		
	for(auto inft = 0u; inft < model.ninfection_trans; inft++){	
		for(auto v = 0u; v < data.narage; v++){ dImap[inft][v] = 0; dIdiag[inft][v] = 0;} 
	}
	
	propose->clear();
}	



/// Used for checking the code is running correctly
void Mbp::check(double t, unsigned int sett) const
{
	cout << t << sett;
	// TO DO
	emsg("TODO 99");
	/*
	if(modeltype == POP_MODEL) emsgEC("Mbp",34);
	
	for(auto i = 0u; i < data.popsize; i++){             // Checks susceptible lises are correct
		auto w = data.ind[i].area*data.ndemocatpos + data.ind[i].dp;
		
		if(nboth_susceptible_list[w] != both_susceptible_list[w].size()) emsgEC("Mbp",24);
		if(npropose_only_susceptible_list[w] != propose_only_susceptible_list[w].size()) emsgEC("Mbp",25);
		if(nnot_susceptible_list[w] != not_susceptible_list[w].size()) emsgEC("Mbp",26);
		
		if((initial->indev[i].size() == 0 || t < initial->indev[i][0].t) && propose->indev[i].size() == 0){
			if(susceptible_status[i] != BOTH_SUSCEPTIBLE) emsgEC("Mbp",27);
			if(both_susceptible_list[w][susceptible_list_ref[i]] != i) emsgEC("Mbp",28);
		}
		else{
			if((initial->indev[i].size() != 0 && t >= initial->indev[i][0].t) && propose->indev[i].size() == 0){
				if(susceptible_status[i] != ONLY_PROPOSE_SUSCEPTIBLE) emsgEC("Mbp",29);
				if(propose_only_susceptible_list[w][susceptible_list_ref[i]] != i) emsgEC("Mbp",30);
			}
			else{
				if(susceptible_status[i] != NOT_SUSCEPTIBLE) emsgEC("Mbp",31);
				if(not_susceptible_list[w][susceptible_list_ref[i]] != i) emsgEC("Mbp",32);
			}
		}
	}
	
	auto l = areatree.level-1;                         // Checks infection sampler is correct
	for(auto c = 0u; c < data.narea; c++){
		auto wmin = c*data.ndemocatpos, wmax = wmin + data.ndemocatpos;
	
		auto sum = 0.0; 
		auto dp = 0u; 
		auto v = c*data.nage; 
		for(auto w = wmin; w < wmax; w++){
			auto a = data.democatpos[dp][0];
			
			double dd;
			switch(simu_or_mbp){
			case SIMU:
				dd = propose->lambda[w] - propose->susceptibility[dp]*(propose->beta*propose->areafactor[c]*propose->Qmap[sett][v+a] 
																													+ propose->phi); 
				if(sqrt(dd*dd) > TINY) emsgEC("Mbp",33);
				sum += (nboth_susceptible_list[w]+npropose_only_susceptible_list[w])*propose->lambda[w];
				break;
			
			case MBP:
				dd = initial->lambda[w] - initial->susceptibility[dp]*(initial->beta*initial->areafactor[c]*initial->Qmap[sett][v+a] 
																														+ initial->phi); 
				if(sqrt(dd*dd) > TINY) emsgEC("Mbp",33);
				dd = propose->lambda[w] - propose->susceptibility[dp]*(propose->beta*propose->areafactor[c]*propose->Qmap[sett][v+a] 
																														+ propose->phi); 
				if(sqrt(dd*dd) > TINY) emsgEC("Mbp",34);
		
				auto dlambda = nboth_susceptible_list[w]*(propose->lambda[w] - initial->lambda[w]); if(dlambda < 0) dlambda = 0;
				sum += dlambda + npropose_only_susceptible_list[w]*propose->lambda[w];
				break;
			}
			
			dp++;
		}
		auto dd = new_infection_rate[l][c] - sum; if(sqrt(dd*dd) > TINY){ emsgEC("Mbp",35);}
	}
	
	for(auto i = 0u; i < data.popsize; i++){
		for(auto tra = 0u; tra < model.trans.size(); tra++){
			if(indmap[i][tra] != 0) emsgEC("Mbp",36);
		}
	}
	*/
}



struct Event {                             // Stores information about a compartmental transition
  unsigned int trans;                      // References the transition type
	unsigned int ind;                        // The individual on which the transition happens
	double t;                                // The time of the transition
};

struct EventRef {                          // Used to reference an event
	unsigned int ind;                        // The individual
	unsigned int e;	                         // The event number
};


void Mpi::pack(const vector <Event> &vec)
{
	auto imax = vec.size(); buffer.push_back(imax); k++;
	for(auto i = 0u; i < imax; i++){
		buffer.push_back(vec[i].trans); k++;
		buffer.push_back(vec[i].ind); k++;
		buffer.push_back(vec[i].t); k++;
	}
}


/*
void Mpi::pack_item(const EventRef& ev)
{
	pack_item(ev.ind);
	pack_item(ev.e);
}
*/


/*
void Mpi::pack(const vector <vector <EventRef> > &vec)
{
	pack_item(vec);
}
*/


void Mpi::unpack(vector <Event> &vec)
{
	unsigned int imax = buffer[k]; k++; vec.resize(imax);
	for(auto i = 0u; i < imax; i++){
		vec[i].trans = buffer[k]; k++;
		vec[i].ind = buffer[k]; k++;
		vec[i].t = buffer[k]; k++;
	}
}


/*
void Mpi::unpack_item(EventRef& ev)
{
	unpack_item(ev.ind);
	unpack_item(ev.e);
}
*/

/*
void Mpi::unpack(vector <vector <EventRef> > &vec)
{
	unpack_item(vec);
}
*/

	void unpack(vector <Event> &vec);
	void unpack(vector< vector <Event> > &vec, unsigned int time_division_per_timesAmin, unsigned int time_division_per_timesAmax);
	
	
	
/// This simulates from the compartmental model and generates an event list for individual i
/// The individual states in state c at time t 
void State::simulate_compartmental_transitions(unsigned int i, unsigned int c, double t)
{
	vector <Event> &evlist = indev[i];

	Event ev;
	ev.ind = i;
		
	auto dp = data.ind[i].dp;
	
	unsigned int tra; 
	if(c == 0){                                            // If an infection then adds first 
		evlist.clear();	
		tra = 0;
		ev.trans = tra; ev.t = t; 
		evlist.push_back(ev);
		c = trans[tra].to;
	}
	
	while(comp[c].trans.size() > 0){                       // Generates the event sequence as a function of time
		auto tra = select_branch(c,dp);                       // Selects which branch to go down
		
		auto dt = sample_duration(tra,dp);                      // Samples the time for the transition
		
		t += dt;
		
		ev.trans = tra; ev.t = t;
		evlist.push_back(ev);

		c = trans[tra].to; 
	}
}


/// Selects which branch to follow
unsigned int State::select_branch(unsigned int c, unsigned int dp)
{
	auto kmax = comp[c].trans.size();
	if(kmax == 1) return comp[c].trans[0];
	else{                                                // Uses the branching probability to select branch
		auto z = ran(); auto k = 0u; while(k < kmax && z > comptransprob[c].probsum[dp][k]) k++;
		if(k == kmax) emsgEC("State",2);
		return comp[c].trans[k];
	}
}


		
/// Samples the time taken for a particular transition to happen
double State::sample_duration(unsigned int tra, unsigned int dp)
{
	double dt;
	
	switch(trans[tra].type){                             // Samples duration spent in compartment
	case EXP_DIST:
		dt = exp_sample_time(paramval[trans[tra].param_mean[dp]]*trans[tra].mean_factor);
		model.print_transition(tra); 
		break;
	
	case GAMMA_DIST:
		{
			auto mean = paramval[trans[tra].param_mean[dp]]; auto sd = paramval[trans[tra].param_cv[dp]]*mean;
			dt = gamma_sample(mean*mean/(sd*sd),mean/(sd*sd));
		}
		break;
		
	case LOGNORM_DIST:
		{
			auto mean_ns = paramval[trans[tra].param_mean[dp]], cv_ns = paramval[trans[tra].param_cv[dp]];
			auto sd = sqrt(log((1+cv_ns*cv_ns))), mean = log(mean_ns) - sd*sd/2;
			dt = lognormal_sample(mean,sd);
		}
		break;
		
	default: emsgEC("State",31); break;
	}

	if(dt < TINY) dt = TINY;
	
	return dt;
}


/// Outputs an event sequence (used for debugging)
void Model::print_events(const string& name, const vector <Event> &ev) const 
{
	cout << name << ":" << endl;
	for(auto e = 0u; e < ev.size(); e++){
		auto tra = ev[e].trans;
		cout << comp[trans[tra].from].name << "->" << comp[trans[tra].to].name << "  " << ev[e].t << endl;
	}
}

In class State
	vector < vector <Event> > indev;                     // The individual event sequences
		vector <EventRef> infev;                             // Ordered list of references to infection events 
		vector < vector <EventRef> > transev;                // Event references
	
	void save_event_sample(const vector < vector <Event> > &fev) const;
		

/// Outputs an event sample fev
void Output::save_event_sample(const vector < vector <Event> > &fev) const
{
	auto nind = data.ind.size();
	vector< vector <Event> > indev;
	indev.resize(nind);
	for(auto& fe: fev){
		for(auto& ev : fe) indev[ev.ind].push_back(ev);
	}
	
	auto file =  details.output_directory+"/events.txt";
	ofstream evsamp(file);
	if(!evsamp) emsg("Cannot output the file '"+file+"'");
}

State::State(const Details &details, const Data &data, const Model &model, const ObservationModel &obsmodel) : comp(model.comp), trans(model.trans), param(model.param), details(details), data(data), model(model), obsmodel(obsmodel)
{
	disc_spline.resize(model.spline.size());

	switch(modeltype){
	case IND_MODEL:
		transev.resize(details.ndivision); 
		indev.resize(data.popsize);                                   
		lambda.resize(data.nsettardp);
		break;
	
	case POP_MODEL:
		pop.resize(details.ndivision);
		for(auto sett = 0u; sett < details.ndivision; sett++){
			pop[sett].resize(data.narea); 
			for(auto c = 0u; c < data.narea; c++){
				pop[sett][c].resize(model.comp.size());
				for(auto co = 0u; co < model.comp.size(); co++) pop[sett][c][co].resize(data.ndemocatpos);
			}
		}
		
		transnum.resize(details.ndivision); transmean.resize(details.ndivision);
		for(auto sett = 0u; sett < details.ndivision; sett++){
		  transnum[sett].resize(data.narea); transmean[sett].resize(data.narea);
			for(auto c = 0u; c < data.narea; c++){
				transnum[sett][c].resize(model.trans.size()); transmean[sett][c].resize(model.trans.size());
				for(auto tr = 0u; tr < model.trans.size(); tr++){
					transnum[sett][c][tr].resize(data.ndemocatpos);
					transmean[sett][c][tr].resize(data.ndemocatpos);
				}		
			}
		}
		break;
	}
	
	Imap.resize(details.ndivision); Idiag.resize(details.ndivision);
	for(auto sett = 0u; sett < details.ndivision; sett++){
		Imap[sett].resize(model.ninfection_trans); Idiag[sett].resize(model.ninfection_trans);
		for(auto inft = 0u; inft < model.ninfection_trans; inft++){
			Imap[sett][inft].resize(data.narage); Idiag[sett][inft].resize(data.narage); 
			for(auto v = 0u; v < data.narage; v++){ 
				Imap[sett][inft][v] = 0; Idiag[sett][inft][v] = 0;
			}
		}
	}
	
	Ntime.resize(details.ndivision);
	for(auto sett = 0u; sett < details.ndivision; sett++){
		Ntime[sett].resize(data.nage);
		for(auto a = 0u; a < data.nage; a++) Ntime[sett][a].resize(data.nage);
	}
}
		
	
/// Clears all the infections from the state
void State::clear()
{
	switch(modeltype){
	  case IND_MODEL:
			for(const auto& i : infev) indev[i.ind].clear();
	
			infev.clear();
			transev.clear(); transev.resize(details.ndivision);
			break;
		
		case POP_MODEL:   // Sets the initial population
			for(auto c = 0u; c < data.narea; c++){
				for(auto co = 0u; co < model.comp.size(); co++){
					for(auto dp = 0u; dp < data.ndemocatpos; dp++){
						if(co == model.start_compartment) pop[0][c][co][dp] = data.area[c].pop[dp];
						else pop[0][c][co][dp] = 0;
					}
				}
			}
			break;
	}
}

void State::set_Imap(unsigned int check)
{
	vector < vector <double> > Ima(data.narage);
	vector < vector <double> > Idia(data.narage);
	
	Ima.resize(model.ninfection_trans); Idia.resize(model.ninfection_trans); 
	for(auto inft = 0u; inft < model.ninfection_trans; inft++){
		Ima[inft].resize(data.narage); Idia[inft].resize(data.narage); 
		for(auto v = 0u; v < data.narage; v++){ Ima[inft][v] = 0; Idia[inft][v] = 0;}
	}
	
	for(auto sett = 0u; sett < details.ndivision; sett++){
		for(auto inft = 0u; inft < model.ninfection_trans; inft++){
			auto& Imap_inft = Imap[sett][inft];
			auto& Idiag_inft = Idiag[sett][inft];
			auto& Ima_inft = Ima[inft];
			auto& Idia_inft = Idia[inft];
			for(auto v = 0u; v < data.narage; v++){
				auto val = Ima_inft[v];
				if(check == 1){
					if(false){ if(val < -TINY && modeltype == IND_MODEL) emsgEC("State",37);}
					if(val < Imap_inft[v]-SMALL || val > Imap_inft[v]+SMALL) emsgEC("State",3);
				}
				else{
					if(false){ if(modeltype == IND_MODEL){ if(val < 0){ val = 0; Ima_inft[v] = 0;}}}
			
					Imap_inft[v] = val;
				}
			
				val = Idia_inft[v];
				if(check == 1){
					if(false){ if(val < -TINY && modeltype == IND_MODEL) emsgEC("State",37);}
					if(val < Idiag_inft[v]-SMALL || val > Idiag_inft[v]+SMALL) emsgEC("State",3);
				}
				else{
					if(false){ if(modeltype == IND_MODEL){ if(val < 0){ val = 0; Idia_inft[v] = 0;}}}
			
					Idiag_inft[v] = val;
				}
			}
		}
		
		switch(modeltype){
			case IND_MODEL:
				emsg("TODO 9");
					/*
				for(const auto& tre : transev[sett]){
					//auto i = tre.ind;
					//Event fev = indev[i][tre.e];

					
					// TO DO
				
					auto v = data.ind[i].area*data.nage+data.democatpos[data.ind[i].dp][0];
					auto dq = trans[fev.trans].DQ[fev.timep];
					if(dq != UNSET){
						for(auto loop = 0u; loop < 2; loop++){
							auto q = model.DQ[dq].q[loop];
							if(q != UNSET){
								auto fac = model.DQ[dq].fac[loop];
								
								auto& genQten = data.genQ.Qten[data.Q[q].Qtenref];
								auto kmax = genQten.ntof[v];
								auto& cref = genQten.tof[v];
								auto& valref = genQten.valf[v];
								if(nage == 1){
									for(auto k = 0u; k < kmax; k++){
										Qma[cref[k]*nage] += fac*valref[k][0];
									}
								}
								else{
									for(auto k = 0u; k < kmax; k++){
										auto vv = cref[k]*nage;	
										for(auto a = 0u; a < nage; a++){
											Qma[vv] += fac*valref[k][a];
											vv++;
										}
									}
								}
							}
						}
					}
					
				}
				*/
				break;
			
			case POP_MODEL:
				update_I_from_transnum(Ima,Idia,transnum[sett]);
				break;
		}
	}
}

void State::initialise_from_particle(const Particle &part)
{
	paramval = part.paramval;
	EF = part.EF;
	Pr = model.prior(paramval);
	
	set_param(paramval);
	
	clear();                                               // Removes the existing initial sequence 

	switch(modeltype){
		case IND_MODEL:
			{
				vector <int> indlist;
				for(const auto& ev : part.ev){                         // Decompresses the events stored in the particles
					int i = ev.ind;
					if(indev[i].size() == 0) indlist.push_back(i);
					indev[i].push_back(ev);
				}	
				
				for(auto i : indlist) add_indev(i);
				
				sort_infev();
			}
			break;
	
		case POP_MODEL:
			transnum = part.transnum;
			for(auto sett = 0u; sett < details.ndivision-1; sett++){
				update_pop(sett);
				for(auto c = 0u; c < data.narea; c++) set_transmean(sett,c);
			}			
			break;
	}
	
	set_Imap(0);
	if(checkon == true) check();
}

Particle State::create_particle() const
{
	Particle part;
	
	part.EF = EF;
	part.paramval = paramval;
	
	switch(modeltype){
		case IND_MODEL:
			{
				vector <Event> store;
				for(const auto& inde : indev){              // Compresses the events to take up as little memory as possible 
					for(const auto& ev : inde) store.push_back(ev);
				}
				part.ev = store;
			}
			break;
			
		case POP_MODEL:
			part.transnum = transnum;
			break;
	}
	
	return part;
}


/// Gets the time of an infection event
double State::get_infection_time(unsigned int n) const
{
	if(n == infev.size()) return LARGE;
	return indev[infev[n].ind][infev[n].e].t;
}


/// Adds an individual event sequence defined in infev[i] to transev
void State::add_indev(unsigned int i)
{
	unsigned int e, emax, se;
	EventRef evref;
	
	auto iev = indev[i];
	emax = iev.size(); if(emax == 0) return;
	
	auto fac = double(details.ndivision)/details.period;
	
	evref.ind = i; evref.e = 0;
	infev.push_back(evref);
	for(e = 0; e < emax; e++){
		evref.e = e;
		se = (unsigned int)(fac*iev[e].t); 
		if(se < details.ndivision) transev[se].push_back(evref);
	}
}

		
/// Used for time ordering event references	
static bool compEventRefTime(EventRefTime lhs, EventRefTime rhs)
{
	return lhs.t < rhs.t;
};


/// Time orders infection events infev
void State::sort_infev()
{
	vector <EventRefTime> xt;
	for(const auto& iev : infev){
		EventRefTime evreft;
		evreft.ind = iev.ind; evreft.e = iev.e;
		if(indev[iev.ind].size() == 0) emsgEC("State",54);
		
		evreft.t = indev[iev.ind][iev.e].t;
		xt.push_back(evreft);	
	}
	sort(xt.begin(),xt.end(),compEventRefTime);

	for(auto i = 0u; i < infev.size(); i++){
		infev[i].ind = xt[i].ind;	infev[i].e = xt[i].e;
	}
}

void State::check()
{
	if(model.inbounds(paramval) == false) emsg("Parameters not in bounds");
		
	switch(modeltype){
		case IND_MODEL:
			{
				for(auto j = 0u; j < infev.size(); j++){                         // Checks bounds
					auto i = infev[j].ind, e = infev[j].e;
					if(i >= indev.size()) emsgEC("State_check",16);
					if(e >= indev[i].size()) emsgEC("State_check",17);
				}
				
				if(infev.size() > 0){
					for(auto j = 0u; j < infev.size()-1; j++){                       // Checks order of infection events
						if(get_infection_time(j) > get_infection_time(j+1)) emsgEC("State_check",18);
					}
				}
				
				for(const auto& iev : indev){                                    // Looks at individual event sequences 
					auto emax = iev.size();
					if(emax > 0){
						auto c = 0u; 
						auto t = 0.0;
						for(const auto& ev : iev){                         
							auto tt = ev.t; if(tt <= t) emsgEC("State_check",19);            // Checks time ordering
							auto tra = ev.trans;
							if(trans[tra].from != c) emsgEC("State_check",20);               // Checks for consistency
							c = trans[tra].to; t = tt;
						}
						if(comp[c].trans.size() != 0) emsgEC("State_check",21);
					}
				}
				
				vector < vector <int> > done;
				done.resize(indev.size());
				auto num = 0u;
				for(auto i = 0u; i < indev.size(); i++){                               // Check that infev correspond to indev
					if(indev[i].size() > 0){
						num++;
						done[i].resize(indev[i].size());
						for(auto e = 0u; e < indev[i].size(); e++) done[i][e] = 0;
					}
				}
				if(num != infev.size()) emsgEC("State_check",40);

				for(auto sett = 0u; sett < details.ndivision; sett++){                 // Chechs transev consistent with indev
					for(auto &trev : transev[sett]){
						auto i = trev.ind, e = trev.e;
						if(e >= indev[i].size()) emsgEC("State_check",41);
						
						auto se = (unsigned int)(details.ndivision*indev[i][e].t/details.period); 
						if(se != sett) emsgEC("State_check",42);
						if(done[i][e] != 0) emsgEC("State_check",43);
						done[i][e] = 1;
					}
				}
				
				for(auto i = 0u; i < indev.size(); i++){
					for(auto e = 0u; e < indev[i].size(); e++){
						if(indev[i][e].t < details.period){
							if(done[i][e] != 1) emsgEC("State_check",44);
						}
						else{
							if(done[i][e] != 0) emsgEC("State_check",45);
						}
					}
				}
			}
			break;
		
		case POP_MODEL:   // Checks transnum is consistent 
			for(auto sett = 0u; sett < details.ndivision; sett++){
				for(auto c = 0u; c < data.narea; c++){
					for(auto tr = 0u; tr < model.trans.size(); tr++){
						auto from = model.trans[tr].from;
						if(model.trans[tr].type == INFECTION_DIST || model.trans[tr].type == EXP_DIST){
							for(auto dp = 0u; dp < data.ndemocatpos; dp++){	 
								auto num = transnum[sett][c][tr][dp];
								if(num < 0) emsgEC("State_check",15);
								
								auto p = pop[sett][c][from][dp];
								if(p <= 0){
									if(num != 0) emsgEC("State_check",74);
								}
								else{
									if(transrate[tr][dp] == 0){
										if(num != 0) emsgEC("State_check",74);
									}
								}
							}
						}
					}
				}
			}		
					
			for(auto sett = 0u; sett < details.ndivision-1; sett++){  // Checks pop is correct
				auto popst = pop[sett+1];
				update_pop(sett);
				
				for(auto c = 0u; c < data.narea; c++){
					for(auto co = 0u; co < model.comp.size(); co++){
						for(auto dp = 0u; dp < data.ndemocatpos; dp++){
							if(popst[c][co][dp] != pop[sett+1][c][co][dp]) emsgEC("State_check",1);
						}
					}
				}
			}	
			break;
	}
	
	set_Imap(1);  // Checks to make sure Imap is OK
	
	if(details.mode == ABC_SMC || details.mode == ABC_MBP){
		auto EF_temp = -2*obsmodel.calculate(this);
		auto dd = EF - EF_temp;
		if(dd*dd > TINY) emsgEC("State_check",2);
	}
	
	auto dd = Pr - model.prior(paramval); if(sqrt(dd*dd) > TINY){ cout << Pr << " " << model.prior(paramval) << " h\n"; emsgEC("State_check",59);}
	
	auto af = model.create_areafactor(paramval);                      // Checks areafactor
	for(auto c = 0u; c < data.narea; c++){
		if(af[c] != areafactor[c]) emsgEC("State_check",88);
	}
	
	auto sus = model.create_susceptibility(paramval);                 // Checks susceptibility
	for(auto dp = 0u; dp < data.ndemocatpos; dp++){
		if(sus[dp] != susceptibility[dp]) emsgEC("State_check",89);
	}

	for(auto sp = 0u; sp < model.spline.size(); sp++){                // Checks disc_spline
		auto spl = model.create_disc_spline(sp,paramval);  
		for(auto se = 0u; se < details.ndivision; se++){
			if(spl[se] != disc_spline[sp][se]) emsgEC("State_check",90);
		}
	}
	
	auto ctprob = model.create_comptransprob(paramval);               // Checks comptransprob
	for(auto c = 0u; c < model.comp.size(); c++){
		for(auto j = 0u; j < comptransprob[c].prob.size(); j++){
			for(auto k = 0u; k < comptransprob[c].prob[j].size(); k++){
				if(ctprob[c].prob[j][k] != comptransprob[c].prob[j][k]) emsgEC("State_check",91);
			}
		}
	}
} 


/*
/// Saves a particle (used for diagnostics)
void State::save_particle()
{
	string file = details.output_directory+"/save.txt";
	ofstream save(file.c_str());
	
	for(auto th = 0u; th < model.param.size(); th++){
		save << paramval[th] << "\n";
	}
	
	for(auto sett = 0u; sett < details.ndivision; sett++){
		for(auto c = 0u; c < data.narea; c++){
			for(auto tr = 0u; tr < model.trans.size(); tr++){
				for(auto dp = 0u; dp < data.ndemocatpos; dp++){	
					save << transnum[sett][c][tr][dp] << "\n";
				}
			}
		}
	}
}

/// Loads a particle (used for diagnostics)
void State::load_particle()
{
	string file = details.output_directory+"/save.txt";
	ifstream load(file.c_str());
	
	clear();
	
	for(auto th = 0u; th < model.param.size(); th++){
		load >> paramval[th];
	}
	
	for(auto sett = 0u; sett < details.ndivision; sett++){
		for(auto c = 0u; c < data.narea; c++){
			for(auto tr = 0u; tr < model.trans.size(); tr++){
				for(auto dp = 0u; dp < data.ndemocatpos; dp++){	
					load >> transnum[sett][c][tr][dp];
				}
			}
		}
	}
	for(auto sett = 0u; sett < details.ndivision-1; sett++) update_pop(sett);
	
	set_param(paramval);
	
	set_Imap(0);
}
*/

In mbp.h
/*
		void shift_covariate_reff(vector <double>& param_propose);	
	
		void mbp_compartmental_transitions(unsigned int i);
		unsigned int area_of_next_infection();
		void add_infection_in_area(unsigned int c, double t);
		void check(double t, unsigned int sett) const;
		void update_dImap(const vector <EventRef> &trei, const vector <EventRef> &trep);
		void setup_susceptible_lists();
		void reset_susceptible_lists();
		void change_susceptible_status(unsigned int i, unsigned int st, unsigned int updateR);
		void construct_infection_sampler(unsigned int sett, const vector < vector <double> > &Imi,  const vector< vector <double> > &Idi, const vector < vector <double> > &Imp, const vector < vector <double> > &Idp);
		void infection_sampler(const vector < vector< vector<double> > > &Imap, const vector < vector< vector<double> > > &Idiag);
		void sortx(vector <EventRef> &x, vector <vector <Event> > &indev) const;
		void calculate_Imap_for_propose();
		*/

		//vector < vector <short> > indmap;									        	    // A map which is used for fast update in update_dImap 
		
		
		/*
		vector <double> lambda, lambdasum;                              // Used when adding and removing individuals
		
		vector < vector <unsigned int> > both_susceptible_list;         // List of individuals which are both susceptible 
		vector <unsigned int> nboth_susceptible_list;  
		vector < vector <unsigned int> > propose_only_susceptible_list; // List of individuals where proposed state is suscptivble
		vector <unsigned int> npropose_only_susceptible_list;       
		vector < vector <unsigned int> > not_susceptible_list;          // List of individuals not sus in either state
		vector <unsigned int> nnot_susceptible_list;
		vector <unsigned int> susceptible_list_ref;
		vector <unsigned int> susceptible_status;
		
		vector <vector <double> > new_infection_rate;                   // Tree giving rate of new infections
		
		vector <int> N;                                                 // The number of individuals in different compartments
		
		bool do_mbp_event;                                              // Set to true if MBPs on compartmental transitions needed
		*/
		
		
in struct Particle
//vector <Event> ev;                       // The event sequence for the particle (IND_MODEL)
	
	in state
		//vector <double> lambda;                              // Total force of infecion for an area
		
		