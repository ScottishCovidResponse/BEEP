/// This class allows for manipution with the model state (i.e. this the latent state defined by transtions in the model)

#include <cmath>   
#include <iostream>
#include <sstream>

using namespace std;

#include "state.hh"
#include "output.hh"

/// Initialises the state class
State::State(const Details &details, const Data &data, const Model &model, const ObservationModel &obsmodel) : comp(model.comp), trans(model.trans), param(model.param), details(details), data(data), model(model), obsmodel(obsmodel)
{
	disc_spline.resize(model.spline.size());

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
		
	
/// Sets up the initial popualtion
void State::pop_init()
{
	for(auto c = 0u; c < data.narea; c++){
		for(auto co = 0u; co < model.comp.size(); co++){
			for(auto dp = 0u; dp < data.ndemocatpos; dp++){
				if(co == model.start_compartment) pop[0][c][co][dp] = data.area[c].pop[dp];
				else pop[0][c][co][dp] = 0;
			}
		}
	}
}


/// Sets up the state using a specified set of parameters
void State::set_param(const vector <double> &paramv)
{
	timer[TIME_SETPARAM].start();
	
	paramval = paramv;

	comptransprob = model.create_comptransprob(paramval);
	
	transrate = model.create_transrate(comptransprob,paramval);
	
	for(auto sp = 0u; sp < model.spline.size(); sp++) disc_spline[sp] = model.create_disc_spline(sp,paramval);

	susceptibility = model.create_susceptibility(paramval);   
	
	areafactor = model.create_areafactor(paramval);  

	model.create_Ntime(Ntime,disc_spline);
 
	beta_R_ratio = model.create_beta_R_ratio(susceptibility,areafactor,paramval,Ntime,transrate);

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
void State::set_Imap_using_dI(unsigned int sett, const State *state, const vector < vector <double> > &dImap, const vector < vector <double> > &dIdiag)
{
	for(auto inft = 0u; inft < model.ninfection_trans; inft++){
		auto& stImap = state->Imap[sett][inft];
		auto& neImap = Imap[sett][inft];
		auto& stIdiag = state->Idiag[sett][inft];
		auto& neIdiag = Idiag[sett][inft];
		auto& dImap_inft = dImap[inft];
		auto& dIdiag_inft = dIdiag[inft];
		for(auto v = 0u; v < data.narage; v++){
			neImap[v] = stImap[v] + dImap_inft[v];
			neIdiag[v] = stIdiag[v] + dIdiag_inft[v];
		}
	}
}


/// Sets the infection map Imap and Idiag as a function of time (based on transitions in transnum)
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
					if(val < Imap_inft[v]-SMALL || val > Imap_inft[v]+SMALL) emsgEC("State",4);
				}
				else{
					if(false){ if(modeltype == IND_MODEL){ if(val < 0){ val = 0; Ima_inft[v] = 0;}}}
			
					Imap_inft[v] = val;
				}
			
				val = Idia_inft[v];
				if(check == 1){
					if(false){ if(val < -TINY && modeltype == IND_MODEL) emsgEC("State",37);}
					if(val < Idiag_inft[v]-SMALL || val > Idiag_inft[v]+SMALL) emsgEC("State",5);
				}
				else{
					if(false){ if(modeltype == IND_MODEL){ if(val < 0){ val = 0; Idia_inft[v] = 0;}}}
			
					Idiag_inft[v] = val;
				}
			}
		}
		
		update_I_from_transnum(Ima,Idia,transnum[sett]);
	}
}


/// Sets the infectivity map at a particular time
void State::set_Imap_sett(unsigned int sett)
{
	auto &Ima = Imap[sett];
	auto &Idia = Idiag[sett];
	
	if(sett == 0){
		for(auto inft = 0u; inft < model.ninfection_trans; inft++){
			auto &Ima_inft = Ima[inft];
			auto &Idia_inft = Idia[inft];
			for(auto v = 0u; v < data.narage; v++){ Ima_inft[v] = 0; Idia_inft[v] = 0;}
		}
	}
	else{
		Ima = Imap[sett-1];
		Idia = Idiag[sett-1];
		update_I_from_transnum(Ima,Idia,transnum[sett-1]);
	}
}


/// Changes a Imap in accordance with transitions in transnum
void State::update_I_from_transnum(vector < vector <double> > &Ima, vector< vector <double> > &Idia, const vector < vector < vector <int> > > &dtransnum) const
{	
	auto narea = data.narea;
	auto nage = data.nage;
	vector <double> dinf(nage);

	for(auto inft = 0u; inft < model.ninfection_trans; inft++){
		auto &Ima_inft = Ima[inft];
		auto &Idia_inft = Idia[inft];
		
		for(auto c = 0u; c < narea; c++){
			for(auto a = 0u; a < nage; a++) dinf[a] = 0;
			
			for(auto tr = 0u; tr < model.trans.size(); tr++){	
				auto di = model.get_infectivity_dif(inft,tr,paramval); 
				if(di != 0){
					for(auto dp = 0u; dp < data.ndemocatpos; dp++){
						auto num = dtransnum[c][tr][dp];
						if(num != 0) dinf[data.democatpos[dp][0]] += di*num; 
					}
				}
			}
			
			auto diag = data.genQ.M.diag[c];
			auto& to = data.genQ.M.to[c];
			auto& val = data.genQ.M.val[c];
			
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


/// Initialises the state based on a particle
void State::initialise_from_particle(const Particle &part)
{
	timer[TIME_INITFROMPART].start();
		
	paramval = part.paramval;
	EF = part.EF;
	transnum = part.transnum;
		
	Pr = model.prior(paramval);
	
	set_param(paramval);
	
	pop_init();	
	for(auto sett = 0u; sett < details.ndivision-1; sett++) update_pop(sett);

	set_Imap(0);
		
	for(auto sett = 0u; sett < details.ndivision; sett++){
		for(auto c = 0u; c < data.narea; c++) set_transmean(sett,c);
	}			
	
	if(checkon == true) check();
	
	timer[TIME_INITFROMPART].stop();
}


/// Creates a particle from the state
Particle State::create_particle() const
{
	Particle part;
	part.EF = EF;
	part.paramval = paramval;
	part.transnum = transnum;
	
	return part;
}


/// Sets the mean number of transtions at a given time sett and within a given area c
void State::set_transmean(unsigned int sett, unsigned int c)
{
	timer[TIME_TRANSMEAN].start();
	
	auto dt = double(details.period)/details.ndivision;
		
	auto &tmean = transmean[sett][c];
	auto &p = pop[sett][c];
		
	for(auto inft = 0u; inft < model.ninfection_trans; inft++){            // Goes over all infection transitions
		auto tr = model.infection_trans[inft];
		auto from = model.trans[tr].from;
		
		auto NMI = get_NMI(sett,inft,c);
		auto fac = beta_R_ratio[inft][sett]*disc_spline[model.R_spline_ref[inft][c]][sett]*areafactor[c];
		
		auto et = disc_spline[model.efoi_spline_ref[inft][c]][sett];
		
		if(details.mode == COUNTER){                                         // Potential counterfactual changes
			et *= model.countermod.efoi_mult[sett][c][inft];
		}
		
		vector <double> eta_age(data.nage);
		const auto& agedist = model.efoi_agedist[inft][c];
	
		for(auto a = 0u; a < data.nage; a++) eta_age[a] = et*agedist[a];
		
		for(auto dp = 0u; dp < data.ndemocatpos; dp++){	
			auto popu = p[from][dp];
			if(popu <= 0) tmean[tr][dp] = 0;
			else{
				auto a = data.democatpos[dp][0];
				tmean[tr][dp] = dt*popu*susceptibility[dp]*(fac*NMI[a] + eta_age[a]);
			}
		}
	}
	
	for(auto tr = 0u; tr < model.trans.size(); tr++){                       // Non-infection transitions
		if(model.trans[tr].type != INFECTION_DIST){
			auto from = model.trans[tr].from;
			for(auto dp = 0u; dp < data.ndemocatpos; dp++){	
				auto popu = p[from][dp];
				if(popu <= 0) tmean[tr][dp] = 0;
				else tmean[tr][dp] = dt*popu*transrate[tr][dp];
			}
		}
	}
	
	if(details.mode == COUNTER){                                            // Incorporates counterfactual modification
		auto &tmean_mult = model.countermod.transmean_mult[sett][c];
		for(auto tr = 0u; tr < model.trans.size(); tr++){   
			for(auto dp = 0u; dp < data.ndemocatpos; dp++){	
				tmean[tr][dp] *= tmean_mult[tr][dp];
			}
		}
	}
	
	timer[TIME_TRANSMEAN].stop();
}


/// Updates the population corresponding to transnum
void State::update_pop(unsigned int sett)
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


// TO DO
/// Makes changes to a spline used in the observation model
void State::obsmodel_prop(MVN &mvn, double EFcut, double invT, ObsModelMode obsmodel_mode)
{
	timer[TIME_OBSMODEL].start(); 
		emsg("PP");
	auto param_store = paramval;
	auto param_propose = mvn.propose(paramval);
	
	auto al = 0.0;
	
	vector <double> obs_value_prop;
	auto Pr_prop = Pr;
	auto EF_prop = EF;
	auto disc_spline_st = disc_spline;
		
	if(model.inbounds(param_propose) == false) mvn.nbo++;
	else{		
	// TO DO
		//for(auto spl = model.matrix_modify_ref+1; spl < model.spline.size(); spl++){
		for(auto spl = 0u; spl < model.spline.size(); spl++){
			disc_spline[spl] = model.create_disc_spline(spl,paramval);
		}
		
		EF_prop = -2*obsmodel.calculate(this);
		Pr_prop = model.prior(paramval);
		switch(obsmodel_mode){
			case CUTOFF: if(EF_prop < EFcut) al = exp(Pr_prop - Pr); break;
			case INVT: al = exp(-invT*(EF_prop - EF)); break;
		}	
	}
				
	mvn.ntr++;
	if(ran() < al){	
		mvn.nac++;
		EF = EF_prop; 
		Pr = Pr_prop;
	}
	else{
		disc_spline = disc_spline_st;
		paramval = param_store;
	}
	
	if(checkon == true) check();
	
	timer[TIME_OBSMODEL].stop(); 
}


/// Gets total infectivity as a function of age group at a given time sett for infection transition inft in an area c
vector <double> State::get_NMI(unsigned int sett, unsigned int inft, unsigned int c)
{ 
	auto nage = data.nage;
	vector <double> NMI(nage);			
	double d = disc_spline[model.geo_spline_ref][sett];
	auto omd = 1-d;
	auto fac = data.genQ.factor[c];

	if(nage == 1){   // Faster version when only 1 age group
		double val = Idiag[sett][inft][c]*(d+fac*omd) + Imap[sett][inft][c]*d;
		if(val < 0) val = 0;
	
		NMI[0] = val;
		return NMI;
	}
	
	vector <double> I(nage);	
	auto v = c*nage;
	auto &Idiag_inft = Idiag[sett][inft];
	auto &Imap_inft = Imap[sett][inft];
	for(auto a = 0u; a < nage; a++){
		I[a] = (Idiag_inft[v+a] + Imap_inft[v+a])*d + fac*Idiag_inft[v+a]*omd;
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
	
	if(checkon == true) check();
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
	
		set_Imap_sett(sett);

		timer[TIME_TRANSNUM].start();
		int num, p;
	
		for(auto c = 0u; c < data.narea; c++){
			auto &prop_pop = pop[sett][c];
			auto &prop_tnum = transnum[sett][c];
			auto &tmean = transmean[sett][c];
			
			set_transmean(sett,c);
			for(auto tr = 0u; tr < model.trans.size(); tr++){
				auto from = model.trans[tr].from;
				
				for(auto dp = 0u; dp < data.ndemocatpos; dp++){	 
					p = prop_pop[from][dp];
					if(p <= 0) num = 0;
					else{
						num = poisson_sample(tmean[tr][dp]);			
						if(num > p) num = p;
					}
					prop_tnum[tr][dp] = num;		
				}
			}
		}
		timer[TIME_TRANSNUM].stop();
			
		if(sett < details.ndivision-1) update_pop(sett);
	}
	
	timer[TIME_SIMULATE].stop();
}


/// Prints the compartmental populations at a given time sett
string State::print_populations(unsigned int sett) const
{
	vector <int> N(model.comp.size());
	
	for(auto co = 0u; co < model.comp.size(); co++){
		N[co] = 0;
		for(auto c = 0u; c < data.narea; c++){
			for(auto dp = 0u; dp < data.ndemocatpos; dp++) N[co] += pop[sett][c][co][dp];
		}
	}

	stringstream ss;
	if(sett == 0) ss << "  Start: ";
	else{
		if(sett == details.ndivision-1) ss << "  End:";
		else ss << "    t=" << details.division_time[sett];
	}
	
	auto len = ss.str().length();
	if(len < 15){ for(auto j = 0u; j < 15-len; j++) ss << " ";}
	
	for(auto c = 0u; c < model.comp.size(); c++) ss << "  " << model.comp[c].name << ":"	<< N[c];
	ss << endl;	
	return ss.str();
}


/// Gets the total number of infected individuals
unsigned int State::get_ninf_total()
{
	auto ninf_total = 0u;
	
	for(auto sett = 0u; sett < details.ndivision; sett++){
		for(auto c = 0u; c < data.narea; c++){   
			for(auto tr : model.infection_trans){
				for(auto dp = 0u; dp < data.ndemocatpos; dp++){	 
					ninf_total += transnum[sett][c][tr][dp];
				}
			}
		}
	}
	
	return ninf_total;
}


/// Creates a sample from the current state (these are used to create the final pdf)
Sample State::create_sample() const 
{
	Sample sample;

	sample.graph_state = obsmodel.get_graph_state(this);
	
	sample.spline_output = model.get_spline_output(paramval);
	
	sample.derived_param = model.get_derived_param(paramval,susceptibility,Ntime[0],transrate);
	
	return sample;
}


/// Creates a parameter sample from the current state 
ParamSample State::create_param_sample() const
{
	ParamSample sample;
	sample.paramval = paramval;
	
	return sample;
}


/// Saves the state to a file
void State::save(string file) const
{
	if(details.mode == COUNTER || details.mode == PPC) return;
	
	ofstream outp(file);
	if(!outp) emsg("Cound not open file '"+file+"'");
	
	outp << "EF: " << EF << endl;
	for(auto th = 0u; th < param.size(); th++){
		outp << model.param[th].name << ": " << paramval[th] << endl;
	}
	
	for(auto sett = 0u; sett < details.ndivision; sett++){
		for(auto c = 0u; c < data.narea; c++){
			for(auto tr = 0u; tr < model.trans.size(); tr++){
				for(auto dp = 0u; dp < data.ndemocatpos; dp++){
					outp << transnum[sett][c][tr][dp] << " ";
				}
			}
		}
		outp << endl;
	}
}


/// Loads the state from a file
void State::load(string file)
{
	string fr;
	Particle part;
	
	ifstream inp(file);
	if(!inp) emsg("Cound not open file '"+file+"'");
	
	inp >> fr; if(fr != "EF:") emsg("Problem loading file '"+file+"'");
	inp >> part.EF;
	
	part.paramval.resize(param.size());
	for(auto th = 0u; th < param.size(); th++){
		inp >> fr; if(fr != model.param[th].name+":") emsg("Problem loading file '"+file+"'");
		inp >> part.paramval[th];
	}
	
	part.transnum.resize(details.ndivision);
	for(auto sett = 0u; sett < details.ndivision; sett++){
		part.transnum[sett].resize(data.narea);
		for(auto c = 0u; c < data.narea; c++){
			part.transnum[sett][c].resize(model.trans.size());
			for(auto tr = 0u; tr < model.trans.size(); tr++){
				part.transnum[sett][c][tr].resize(data.ndemocatpos);
				for(auto dp = 0u; dp < data.ndemocatpos; dp++){
					if(inp.eof()) emsg("Problem loading file '"+file+"'");
					inp >> part.transnum[sett][c][tr][dp];
				}
			}
		}
	}
	inp >> fr; if(!inp.eof()) emsg("Problem loading file '"+file+"'");
	initialise_from_particle(part);
}
