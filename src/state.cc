/// This class allows for manipution with the model state (i.e. this the latent state defined by transtions in the model)

#include <cmath>   
#include <iostream>
#include <sstream>

using namespace std;

#include "state.hh"
#include "output.hh"
#include "utils.hh"

/// Initialises the state class
State::State(const Details &details, const Data &data, const Model &model, const ObservationModel &obsmodel) : comp(model.comp), trans(model.trans), param(model.param), details(details), data(data), model(model), obsmodel(obsmodel)
{
	disc_spline.resize(model.spline.size());

	pop.resize(details.ndivision+1);
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
	
	Imap.resize(details.ndivision); Idiag.resize(details.ndivision);
	for(auto sett = 0u; sett < details.ndivision; sett++){
		Imap[sett].resize(data.nstrain); Idiag[sett].resize(data.nstrain);
		for(auto st = 0u; st < data.nstrain; st++){
			Imap[sett][st].resize(data.narage); Idiag[sett][st].resize(data.narage); 
			for(auto v = 0u; v < data.narage; v++){ 
				Imap[sett][st][v] = 0; Idiag[sett][st][v] = 0;
			}
		}
	}
	
	Ntime.resize(details.ndivision);
	for(auto sett = 0u; sett < details.ndivision; sett++){
		Ntime[sett].resize(data.nage);
		for(auto a = 0u; a < data.nage; a++) Ntime[sett][a].resize(data.nage);
	}
}
		
	
/// Sets up the initial popualtion
void State::pop_init()
{
	for(auto c = 0u; c < data.narea; c++) pop[0][c] = data.area[c].pop_init;
}


/// Sets up the state using a specified set of parameters
void State::set_param(const vector <double> &paramv)
{
	timer[TIME_SETPARAM].start();
	
	paramval = paramv;

	auto paramv_dir = model.dirichlet_correct(paramval);

	transrate = model.create_transrate(paramv_dir);
	
	disc_spline = model.create_disc_spline(paramv_dir);

	susceptibility = model.create_susceptibility(paramv_dir);   
	
	areafactor = model.create_areafactor(paramv_dir);  

	Ntime = model.create_Ntime(disc_spline);
 
	beta = model.calculate_beta_from_R(susceptibility,paramv_dir,Ntime,transrate,disc_spline);

	genT = model.calculate_generation_time(paramv_dir,susceptibility,Ntime[0],transrate,data.democatpos_dist,0);
	
	//model.eignevector_compare_models(susceptibility,paramv_dir,Ntime,transrate);
		
	timer[TIME_SETPARAM].stop();
}


/// Sets the error function (defined to be -2*log(observation model)
void State::set_EF()
{
	EF = -2*obsmodel.calculate(this);
}


/// Sets the prior
void State::set_Pr()
{
	Pr = model.prior(paramval);
}


/// For a given time sett, sets Imap by using the Imap from another state and the difference given by dImap
void State::set_Imap_using_dI(const unsigned int sett, const State *state, const vector < vector <double> > &dImap, const vector < vector <double> > &dIdiag)
{
	for(auto st = 0u; st < data.nstrain; st++){
		auto &stImap = state->Imap[sett][st];
		auto &neImap = Imap[sett][st];
		auto &stIdiag = state->Idiag[sett][st];
		auto &neIdiag = Idiag[sett][st];
		auto &dImap_inft = dImap[st];
		auto &dIdiag_inft = dIdiag[st];
		for(auto v = 0u; v < data.narage; v++){
			neImap[v] = stImap[v] + dImap_inft[v];
			neIdiag[v] = stIdiag[v] + dIdiag_inft[v];
		}
	}
}


/// Sets the infection map Imap and Idiag as a function of time (based on transitions in transnum)
void State::set_Imap(const unsigned int check)
{
	set_I_from_pop(0,false);
	auto Ima = Imap[0];
	auto Idia = Idiag[0];
	
	for(auto sett = 0u; sett < details.ndivision; sett++){
		for(auto st = 0u; st < data.nstrain; st++){
			auto &Imap_inft = Imap[sett][st];
			auto &Idiag_inft = Idiag[sett][st];
			auto &Ima_inft = Ima[st];
			auto &Idia_inft = Idia[st];
			for(auto v = 0u; v < data.narage; v++){
				auto val = Ima_inft[v];
				if(check == 1){
					if(val < Imap_inft[v]-SMALL || val > Imap_inft[v]+SMALL) emsgEC("State",4);
				}
				else{
					Imap_inft[v] = val;
				}
			
				val = Idia_inft[v];
				if(check == 1){
					if(val < Idiag_inft[v]-SMALL || val > Idiag_inft[v]+SMALL) emsgEC("State",5);
				}
				else{
					Idiag_inft[v] = val;
				}
			}
		}
		
		update_I_from_transnum(Ima,Idia,transnum[sett]);
	}
}


/// Sets the infectivity map at a particular time
void State::set_Imap_sett(const unsigned int sett)
{
	if(sett == 0) set_I_from_pop(sett,false);
	else{
		auto &Ima = Imap[sett];
		auto &Idia = Idiag[sett];
		Ima = Imap[sett-1];
		Idia = Idiag[sett-1];
		update_I_from_transnum(Ima,Idia,transnum[sett-1]);
	}
}


/// Changes a Imap in accordance with transitions in transnum
void State::update_I_from_transnum(vector < vector <double> > &Ima, vector< vector <double> > &Idia, const vector < vector < vector <double> > > &dtransnum) const
{	
	auto narea = data.narea;
	auto nage = data.nage;
	auto dpmax = data.ndemocatpos_per_strain;
	vector <double> dinf(nage);

	for(auto st = 0u; st < data.nstrain; st++){
		auto &Ima_inft = Ima[st];
		auto &Idia_inft = Idia[st];
		
		for(auto c = 0u; c < narea; c++){
			for(auto a = 0u; a < nage; a++) dinf[a] = 0;
			
			for(auto tr = 0u; tr < model.trans.size(); tr++){	
				auto di = model.get_infectivity_dif(tr,paramval); 
				if(di != 0){
					auto dpp = st*dpmax;
					for(auto dp = 0u; dp < dpmax; dp++){
						auto num = dtransnum[c][tr][dpp];
						if(num != 0) dinf[data.democatpos[dp][0]] += di*num; 
						dpp++;
					}
				}
			}
			
			auto diag = data.genQ.M.diag[c];
			auto &to = data.genQ.M.to[c];
			auto &val = data.genQ.M.val[c];
			
			auto jmax = to.size();
			auto v = c*nage;
			for(auto a = 0u; a < nage; a++){
				auto di = dinf[a];
				if(di != 0){
					Idia_inft[v+a] += di*diag;
					for(auto j = 0u; j < jmax; j++) Ima_inft[to[j]*nage+a] += di*val[j]; 
				}
			}
		}
	}
}


/// Sets the infectivty map based on the populations in compoartments
void State::set_I_from_pop(const unsigned int sett, const bool check)
{
	auto &Ima = Imap[sett];                        // Sets infectivity map to zero
	auto &Idia = Idiag[sett];

	vector< vector <double> > Ima_store;          // The infectivity map coming from other areas
	vector< vector <double> > Idia_store;  
	
	if(check == true){
		Ima_store = Ima; Idia_store = Idia;
	}
	
	for(auto st = 0u; st < data.nstrain; st++){
		auto &Ima_inft = Ima[st];
		auto &Idia_inft = Idia[st];
		for(auto v = 0u; v < data.narage; v++){ Ima_inft[v] = 0; Idia_inft[v] = 0;}
	}
	
	auto narea = data.narea;
	auto nage = data.nage;
	auto dpmax = data.ndemocatpos_per_strain;
	vector <double> dinf(nage);

	for(auto st = 0u; st < data.nstrain; st++){
		auto &Ima_inft = Ima[st];
		auto &Idia_inft = Idia[st];
		
		for(auto c = 0u; c < narea; c++){
			for(auto a = 0u; a < nage; a++) dinf[a] = 0;
			
			for(auto co = 0u; co < model.comp.size(); co++){
				auto inf = paramval[comp[co].infectivity_param];  // infectivity
			
				if(inf != 0){
					auto dpp = st*dpmax;
					for(auto dp = 0u; dp < dpmax; dp++){
						auto p = pop[sett][c][co][dpp]; 
						if(p != 0) dinf[data.democatpos[dp][0]] += inf*p;
						dpp++;
					}
				}
			}
			
			auto diag = data.genQ.M.diag[c];
			auto &to = data.genQ.M.to[c];
			auto &val = data.genQ.M.val[c];
			
			auto jmax = to.size();
			auto v = c*nage;
			for(auto a = 0u; a < nage; a++){
				auto di = dinf[a];
				if(di != 0){
					Idia_inft[v+a] += di*diag;
					for(auto j = 0u; j < jmax; j++) Ima_inft[to[j]*nage+a] += di*val[j]; 
				}
			}
		}
	}
	
	if(check == true){
		for(auto st = 0u; st < data.nstrain; st++){
			for(auto v = 0u; v < data.narage; v++){
				auto d = Ima[st][v] - Ima_store[st][v];
				if(d*d > SMALL) emsg("Problem with infectivty map");
			}
		}
	}
}


/// Initialises the state based on a particle
void State::initialise_from_particle(const Particle &part)
{
	timer[TIME_INITFROMPART].start();
	
	paramval = part.paramval;
	EF = part.EF;
	transnum = part.transnum;
		
	Pr = model.prior(paramval);
	if(paramval.size()  == 0) emsg("zero size");
	
	set_param(paramval);
	
	pop_init();
	for(auto sett = 0u; sett < details.ndivision; sett++){
		democat_change_pop_adjust(sett);
		if(sett < details.ndivision-1) update_pop(sett);
	}
	
	set_Imap(0);
		
	for(auto sett = 0u; sett < details.ndivision; sett++){
		for(auto c = 0u; c < data.narea; c++) set_transmean(sett,c);
	}			
	
	if(checkon == true) check(0);
	
	timer[TIME_INITFROMPART].stop();
}


/// Creates a particle from the state
Particle State::create_particle(const unsigned int run) const
{
	Particle part;
	part.run = run;
	part.EF = EF;
	part.paramval = paramval;
	part.transnum = transnum;
	
	return part;
}


/// Adjusts populations so to enforce democat_change 
void State::democat_change_pop_adjust(const unsigned int sett)
{
	for(const auto &dcc : data.democat_change){
		auto ncat = dcc.frac[sett].size();
	
		vector <double> pop_cat(ncat);
		for(auto f = 0u; f < ncat; f++) pop_cat[f] = 0;
			
		for(auto c : dcc.area){
			for(auto g = 0u; g < dcc.dp_group.size(); g++){
				for(auto co = 0u; co < model.comp.size(); co++){
					for(auto f = 0u; f < ncat; f++){
						pop_cat[f] += pop[sett][c][co][dcc.dp_group[g][f]];		
					}							
				}				
			}
		}
		auto pop_tot = 0.0;	for(auto f = 0u; f < ncat; f++)	pop_tot += pop_cat[f];
	
		vector <double> frac_dif(ncat);
		for(auto f = 0u; f < ncat; f++) frac_dif[f] = dcc.frac[sett][f] - pop_cat[f]/pop_tot;
		
		for(auto c : dcc.area){
			for(auto g = 0u; g < dcc.dp_group.size(); g++){
				for(auto co = 0u; co < model.comp.size(); co++){
					auto sum = 0.0;
					for(auto f = 0u; f < ncat; f++){
						sum += pop[sett][c][co][dcc.dp_group[g][f]];
					}
					
					for(auto f = 0u; f < ncat; f++){
						pop[sett][c][co][dcc.dp_group[g][f]] += frac_dif[f]*sum;
					}							
				}				
			}
		}
	}
}


/// Sets the mean number of transtions at a given time sett and within a given area c
void State::set_transmean(const unsigned int sett, const unsigned int c)
{
	timer[TIME_TRANSMEAN].start();
	
	auto dt = double(details.period)/details.ndivision;
				
	auto &tmean = transmean[sett][c];
	auto &p = pop[sett][c];
	
	auto dpmax = data.ndemocatpos_per_strain;
	
	auto tr = model.infection_trans;
	
	auto &sus_pop = p[model.trans[tr].from];
	
	if(data.nstrain > 1){
		for(auto dp = 0u; dp < dpmax; dp++){	                                           // Shifts susceptible populations 
			auto sum = 0.0;
			for(auto st = 1u; st < data.nstrain; st++){ 
				auto dpp = dpmax*st + dp; 
				sum += sus_pop[dpp]; sus_pop[dpp] = 0;
			}
			sus_pop[dp] += sum;
		}
	}
	
	for(auto st = 0u; st < data.nstrain; st++){                                        // Goes over all strains
		auto NMI = get_NMI(sett,st,c);

		auto be = beta[st][c][sett];
		be *= areafactor[sett/details.division_per_time][c];
		
		auto efoi_info = model.efoispline_info[model.efoi_spl_ref[st][c]];
		auto et = disc_spline[efoi_info.spline_ref][sett]/data.popsize;
		
		if(details.mode == PREDICTION) et *= model.modelmod.efoi_mult[sett][c][st];      // Potential model changes
		
		vector <double> eta_age(data.nage);
		const auto &agedist = efoi_info.efoi_agedist;
	
		for(auto a = 0u; a < data.nage; a++) eta_age[a] = et*agedist[a];
		
		auto dpp = st*dpmax;
		for(auto dp = 0u; dp < dpmax; dp++){	
			auto popu = sus_pop[dp];
			if(popu <= 0) tmean[tr][dpp] = 0;
			else{
				auto a = data.democatpos[dp][0];
				tmean[tr][dpp] = dt*popu*susceptibility[dpp]*(be*NMI[a] + eta_age[a]);
				if(std::isnan(tmean[tr][dpp])) emsgEC("State",4);
			}
			dpp++;
		}
	}
	
	for(auto tr = 0u; tr < model.trans.size(); tr++){                       // Non-infection transitions
		if(model.trans[tr].inf == TRANS_NOTINFECTION){
			auto from = model.trans[tr].from;
			for(auto dp = 0u; dp < data.ndemocatpos; dp++){	
				auto popu = p[from][dp];
				if(popu <= 0) tmean[tr][dp] = 0;
				else tmean[tr][dp] = dt*popu*transrate[tr][dp];
			}
		}
	}
	
	if(details.mode == PREDICTION){                                         // Incorporates model modification
		auto &tmean_mult = model.modelmod.transmean_mult[sett][c];
		for(auto tr = 0u; tr < model.trans.size(); tr++){   
			for(auto dp = 0u; dp < data.ndemocatpos; dp++){	
				tmean[tr][dp] *= tmean_mult[tr][dp];
			}
		}
	}
	
	timer[TIME_TRANSMEAN].stop();
}


/// Updates the population corresponding to transnum
void State::update_pop(const unsigned int sett)
{
	timer[TIME_UPDATEPOP].start();
	
	auto &popnext = pop[sett+1];
	popnext = pop[sett];
	for(auto c = 0u; c < data.narea; c++){
		auto &popnext_c = popnext[c];
		auto &tnum = transnum[sett][c];

		for(auto tr = 0u; tr < model.trans.size(); tr++){
			auto from = model.trans[tr].from;
			auto to = model.trans[tr].to;
			
			for(auto dp = 0u; dp < data.ndemocatpos; dp++){	
				auto num = tnum[tr][dp];				
				if(num != 0){
					popnext_c[from][dp] -= num;
					popnext_c[to][dp] += num; 
				}
			}		
		}
	}
	
	timer[TIME_UPDATEPOP].stop();
}


/// Updates the population corresponding to transnum going backwards in time (used for steering)
void State::update_pop_rev(const unsigned int sett)
{
	timer[TIME_UPDATEPOP].start();
	
	auto &popnext = pop[sett-1];
	popnext = pop[sett];

	for(auto c = 0u; c < data.narea; c++){
		auto &popnext_c = popnext[c];
		auto &tnum = transnum[sett][c];

		for(auto tr = 0u; tr < model.trans.size(); tr++){
			auto from = model.trans[tr].from;
			auto to = model.trans[tr].to;
			
			for(auto dp = 0u; dp < data.ndemocatpos; dp++){	
				int num = -tnum[tr][dp];				
				if(num != 0){
					popnext_c[from][dp] -= num;
					popnext_c[to][dp] += num; 
				}
			}		
		}
	}
	
	timer[TIME_UPDATEPOP].stop();
}


/// Makes sure that the populations are positive
void State::pop_positive(const unsigned int sett)
{
	timer[TIME_UPDATEPOP].start();
	
	for(auto c = 0u; c < data.narea; c++){
		for(auto dp = 0u; dp < data.ndemocatpos; dp++){	
			auto sum = 0.0;
			auto &p = pop[sett][c];
			for(auto co = 0u; co < model.comp.size(); co++){
				if(p[co][dp] < 0){ 
					//pop[sett][c][model.start_compartment][dp] += p[dp];
					p[co][dp] = 0; 
				}
				sum += p[co][dp];
			}

			auto fac = data.area[c].pop_dp[dp]/sum;
			if(fac < 0.999999 || fac > 1.000001){
				for(auto co = 0u; co < model.comp.size(); co++) p[co][dp] *= fac;
			}
		}
	}		
	
	timer[TIME_UPDATEPOP].stop();
}


/// Gets total infectivity as a function of age group at a given time sett for infection transition inft in an area c
vector <double> State::get_NMI(const unsigned int sett, const unsigned int st, const unsigned int c)
{ 
	auto nage = data.nage;
	vector <double> NMI(nage);			
	double d = disc_spline[model.geo_spline_ref][sett];
	auto omd = 1-d;
	auto fac = data.genQ.factor[c];

	if(nage == 1){   // Faster version when only 1 age group
		double val = Idiag[sett][st][c]*(d+fac*omd) + Imap[sett][st][c]*d;
		if(val < 0) val = 0;
		NMI[0] = val;
		return NMI;
	}
	
	vector <double> I(nage);	
	auto v = c*nage;
	auto &Idiag_st = Idiag[sett][st];
	auto &Imap_st = Imap[sett][st];
	for(auto a = 0u; a < nage; a++){
		I[a] = (Idiag_st[v+a] + Imap_st[v+a])*d + fac*Idiag_st[v+a]*omd;
	}
	
	for(auto a = 0u; a < nage; a++){
		auto sum = 0.0;
		for(auto aa = 0u; aa < nage; aa++) sum += Ntime[sett][a][aa]*I[aa];
		if(sum < 0) sum = 0;
		NMI[a] = sum;
	}
	
	return NMI;
}


/// Simulates the entire time for a set of parameter values
void State::simulate(const vector <double> &paramval)
{
	set_param(paramval);                                // Sets up the state with the parameter values                   
	simulate(0,details.ndivision);                      // Simulates from the model
	
	set_EF();                                           // Calculates the error function as -2*log(obsmodel)
	
	set_Pr();                                           // Calculates the prior
	
	if(checkon == true) check(1);
}


/// Simulates the state between two time points
void State::simulate(const unsigned int ti, const unsigned int tf)
{
	timer[TIME_SIMULATE].start();

	if(ti == 0) pop_init();
	
	auto step = (unsigned int)(details.ndivision/10.0); if(step == 0) step = 1;
	for(auto sett = ti; sett < tf; sett++){
		if(details.mode == SIM && (sett%step == 0 || sett == details.ndivision-1)){
			cout << print_populations(sett);   // Prints the compartmental populations to the terminal
		}
		democat_change_pop_adjust(sett);
	
		if(sett == ti) set_I_from_pop(sett,false);
		else set_Imap_sett(sett);

		timer[TIME_TRANSNUM].start();
		double mean;
		for(auto c = 0u; c < data.narea; c++){
			auto &prop_tnum = transnum[sett][c];
			auto &tmean = transmean[sett][c];
			
			set_transmean(sett,c);
			
			for(auto tr = 0u; tr < model.trans.size(); tr++){
				for(auto dp = 0u; dp < data.ndemocatpos; dp++){	 
					mean = tmean[tr][dp];
			
					if(mean == 0) prop_tnum[tr][dp] = 0;
					else{
						if(details.stochastic == true) prop_tnum[tr][dp] = poisson_sample(mean);		
						else prop_tnum[tr][dp] = mean;
					}						
				}
			}
		}
		timer[TIME_TRANSNUM].stop();
			
		if(sett < details.ndivision-1) update_pop(sett);
	}
	
	timer[TIME_SIMULATE].stop();
}


/// Prints the compartmental populations at a given time sett
string State::print_populations(const unsigned int sett) const
{
	stringstream ss;
	
	
	vector <double> N(model.comp.size());
	
	for(auto co = 0u; co < model.comp.size(); co++){
		N[co] = 0;
		for(auto c = 0u; c < data.narea; c++){
			for(auto dp = 0u; dp < data.ndemocatpos; dp++) N[co] += pop[sett][c][co][dp];
		}
	}

	if(sett == 0) ss << "  Start: ";
	else{
		if(sett == details.ndivision-1) ss << "  End:";
		else ss << "    t=" << details.division_time[sett];
	}
	
	auto len = ss.str().length();
	if(len < 15){ for(auto j = 0u; j < 15-len; j++) ss << " ";}
	
	for(const auto &cn : model.comp_name){
		auto sum = 0.0; for(auto c : cn.comp) sum += N[c];
		ss << "  " << cn.name << ":"	<< sum;
	}
	ss << endl;	

	return ss.str();
}


/// Creates a sample from the current state (these are used to create the final visualisations)
Sample State::create_sample() const 
{
	Sample sample;

	sample.graph_state = obsmodel.get_graph_state(this);

	auto paramv_dir = model.dirichlet_correct(paramval);
	
	sample.spline_output = model.get_spline_output(paramv_dir,pop);
	
	sample.derived_param = model.get_derived_param(paramv_dir,susceptibility,Ntime[0],transrate);
	
	sample.Rmap = model.get_Rmap(paramv_dir,disc_spline,areafactor,susceptibility,Ntime,transrate,pop);

	sample.Nsample.push_back(Ntime[0]);
		
	return sample;
}


/// Creates a parameter sample from the current state 
ParamSample State::create_param_sample(const unsigned int run) const
{
	ParamSample sample;
	sample.run = run;
	sample.paramval = paramval;
	sample.EF = EF;
	
	return sample;
}


/// Saves the state to a file
void State::save(const string file) const
{
	if(details.siminf == SIMULATE) return;
	
	auto filefull = file+"_parameter.csv";
	ofstream outp(filefull);
	if(!outp) emsg("Cannot open the file '"+filefull+"'");
	
	outp << "Parameter,Value" << endl;
	outp << "Error function," << EF << endl;
	for(auto th = 0u; th < param.size(); th++){
		outp << model.param[th].name << "," << paramval[th] << endl;
	}
	
	filefull = file+"_transition.csv"; 
	ofstream outp_trans(filefull);
	if(!outp_trans) emsg("Cannot open the file '"+filefull+"'");

	outp_trans << "Time,Transitions (ordered by area x trans x democatpos)" << endl;
	for(auto sett = 0u; sett < details.ndivision; sett++){
		outp_trans << details.division_time[sett];
		for(auto c = 0u; c < data.narea; c++){
			for(auto tr = 0u; tr < model.trans.size(); tr++){
				for(auto dp = 0u; dp < data.ndemocatpos; dp++){
					outp_trans << "," << transnum[sett][c][tr][dp];
				}
			}
		}
		outp_trans << endl;
	}
}


/// Loads the state from a file
void State::load(const string file, const unsigned int ndivision_post)
{
	string fr;
	Particle part;

	auto filefull = file+"_parameter.csv";
	ifstream inp(filefull);
	if(!inp) emsg("Cannot open the file '"+filefull+"'");
	
	string line;
	getline(inp,line);
	getline(inp,line);
	auto spl = split(line,','); if(spl.size() != 2) emsg("Problem loading the file '"+filefull+"'1");
	part.EF = get_double(spl[1],"In loading '"+filefull+"'");
	
	part.paramval.resize(param.size());
	for(auto th = 0u; th < param.size(); th++){
		getline(inp,line);
		auto spl = split(line,','); if(spl.size() != 2) emsg("Problem loading the file '"+filefull+"2'");
	
		if(spl[0] != model.param[th].name) emsg("Problem loading the file '"+filefull+"'3");
		part.paramval[th] = get_double(spl[1],"In file '"+filefull+"'");
	}
	inp >> fr; if(!inp.eof()) emsg("Problem loading the file '"+file+"'4");
	
	filefull = file+"_transition.csv";
	ifstream inp_trans(filefull);
	if(!inp_trans) emsg("Cannot open the file '"+filefull+"'");
	
	getline(inp_trans,line);
	part.transnum.resize(details.ndivision);
	for(auto sett = 0u; sett < ndivision_post; sett++){
		getline(inp_trans,line);
		auto spl = split(line,',');
	
		auto num = 1u;
		part.transnum[sett].resize(data.narea);
		for(auto c = 0u; c < data.narea; c++){
			part.transnum[sett][c].resize(model.trans.size());
			for(auto tr = 0u; tr < model.trans.size(); tr++){
				part.transnum[sett][c][tr].resize(data.ndemocatpos);
				for(auto dp = 0u; dp < data.ndemocatpos; dp++){
					part.transnum[sett][c][tr][dp] = get_double(spl[num],"Problem loading the file '"+filefull+"'");
					num++; if(num > spl.size()) emsg("Problem loading the file '"+filefull+"'");
				}
			}
		}
		if(num != spl.size()) emsg("Problem loading the file '"+filefull+"'");
	}
	inp_trans >> fr; if(!inp_trans.eof()) emsg("Problem loading the file '"+file+"'");
	
	for(auto sett = ndivision_post; sett < details.ndivision; sett++){
		part.transnum[sett].resize(data.narea);
		for(auto c = 0u; c < data.narea; c++){
			part.transnum[sett][c].resize(model.trans.size());
			for(auto tr = 0u; tr < model.trans.size(); tr++){
				part.transnum[sett][c][tr].resize(data.ndemocatpos);
				for(auto dp = 0u; dp < data.ndemocatpos; dp++){
					part.transnum[sett][c][tr][dp] = 0;
				}
			}
		}
	}
	
	initialise_from_particle(part);
}


/// Calculates the likelihood for the state
double State::likelihood() const
{
	auto Li = 0.0;
	for(auto sett = 0u; sett < details.ndivision; sett++){  
		for(auto c = 0u; c < data.narea; c++){      
			for(auto tr = 0u; tr < model.trans.size(); tr++){
				for(auto dp = 0u; dp < data.ndemocatpos; dp++){	 
					Li += poisson_probability(transnum[sett][c][tr][dp],transmean[sett][c][tr][dp]);	
				}
			}
		}
	}
	
	return Li;
}

