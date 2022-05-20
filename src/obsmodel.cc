// This describes the observation model. This prescribes how likely the data is given a true event state.

#include <cmath>

using namespace std;

#include "obsmodel.hh"
#include "details.hh"
#include "data.hh"
#include "model.hh"
#include "state.hh"

ObservationModel::ObservationModel(const Details &details, Data &data, const Model &model) : details(details), data(data), model(model)
{
	initialise_obs_change();
	
	if(details.mode == PMCMC_INF) split_observations();
}


/// Measures how well the particle agrees with the observations 
double ObservationModel::calculate(const State *state) const
{
	timer[TIME_OBSPROB].start();
		
	auto obs_value = get_obs_value(state);

	auto L = 0.0; for(auto i = 0u; i < data.nobs; i++) L += obs_prob(obs_value[i],data.obs[i]);
	
	if(std::isnan(L)) emsgEC("ObsModel",1);

	timer[TIME_OBSPROB].stop();
	
	return L;
}


/// Measures how well the particle agrees with the observations within a given section
double ObservationModel::calculate_section(const State *state, unsigned int sec) const
{
	timer[TIME_OBSPROB].start();
	
	vector <double> obs_value(data.nobs);
	
	for(auto i : section_obs[sec])  obs_value[i] = 0;
	
	get_obs_value_section(state,obs_value,section_ti[sec],section_tf[sec]);
		
	auto L = 0.0;		

	for(auto i : section_obs[sec]) L += obs_prob(obs_value[i],data.obs[i]);

	if(std::isnan(L)) emsgEC("ObsModel",2);
	
	timer[TIME_OBSPROB].stop();
	
	return L;
}


/// Uses precalculated quantities to calculate measured quantities faster
vector <double> ObservationModel::get_obs_value(const State *state) const
{
	vector <double> obs_value(data.nobs);
	for(auto &ob : obs_value) ob = 0;
	
	get_obs_value_section(state,obs_value,0,details.ndivision);
		
	return obs_value;
}
	

/// Gets the observation values within a given time period
void ObservationModel::get_obs_value_section(const State *state,  vector <double> &obs_value, const unsigned int ti, const unsigned int tf) const
{
	for(auto sett = ti; sett < tf; sett++){
		for(auto c = 0u; c < data.narea; c++){
			for(auto tr = 0u; tr < model.trans.size(); tr++){
				for(auto dp = 0u; dp < data.ndemocatpos; dp++){
					double num;
					if(obsmodel_transmean == true && details.siminf == INFERENCE) num = state->transmean[sett][c][tr][dp];
					else num = state->transnum[sett][c][tr][dp];
					
					if(num != 0){
						for(auto ob : obs_trans[sett][c][tr][dp]){
							auto spl = data.obs[ob].factor_spline;
							if(spl == UNSET) obs_value[ob] += data.obs[ob].factor*num;
							else obs_value[ob] += data.obs[ob].factor*num*state->disc_spline[spl][sett];
						}
					}
				}
			}
			
			for(auto co = 0u; co < model.comp.size(); co++){
				for(auto dp = 0u; dp < data.ndemocatpos; dp++){
					if(obs_pop[sett][c][co][dp].size() > 0){
						auto num = state->pop[sett][c][co][dp];
						for(auto ob : obs_pop[sett][c][co][dp]){
							auto spl = data.obs[ob].factor_spline;
							if(spl != UNSET) emsgEC("ObsModel",3);
							if(spl == UNSET) obs_value[ob] += data.obs[ob].factor*num;
							else obs_value[ob] += data.obs[ob].factor*num*state->disc_spline[spl][sett];
						}
					}
				}
			}
		}
	}
}


/// Generates the underlying dynamics for each of the graphs
vector < vector <double> > ObservationModel::get_graph_state(const State *state) const
{
	vector < vector <double> > graph_state;

	auto obs_value = get_obs_value(state);

	graph_state.resize(data.graph.size());
	for(auto gr = 0u; gr < data.graph.size(); gr++){
		auto gra = data.graph[gr];
		auto dt = gra.datatable;
		auto spl = gra.factor_spline;
		
		switch(gra.type){
			case GRAPH_TIMESERIES:
				{
					auto graph_step = details.graph_step;
					auto fac_trans = gra.factor*double(details.division_per_time)/graph_step;
					auto fac_pop = gra.factor*double(1.0)/graph_step;

					for(auto sett = 0u; sett < details.ndivision - graph_step; sett += graph_step){
						auto sum = 0.0;
						for(auto step = 0u; step < graph_step; step++){
							auto sett2 = sett+step;
							
							for(auto c : gra.area){
								for(auto dp : gra.dp_sel){
									for(auto tr : data.datatable[dt].translist){	
										auto num = fac_trans*state->transnum[sett2][c][tr][dp];
										if(spl == UNSET) sum += num;
										else sum += num*state->disc_spline[spl][sett2];
									}
									
									for(auto co : data.datatable[dt].complist){
										auto num = fac_pop*state->pop[sett2][c][co][dp];
										if(spl == UNSET) sum += num;
										else sum += num*state->disc_spline[spl][sett2];
									}
								}
							}
						}

						graph_state[gr].push_back(sum);
					}
				}
				break;
					
			case GRAPH_MARGINAL:
				{
					for(const auto &p : gra.point){ 
						graph_state[gr].push_back(obs_value[p.obs]);
					}
				}
				break;
		}
	}
	
	return graph_state;
}
 
 
/// Initialises quantities related to how Measurements change when a transition changes
void ObservationModel::initialise_obs_change()
{
	obs_trans.resize(details.ndivision);
	for(auto sett = 0u; sett < details.ndivision; sett++){
		obs_trans[sett].resize(data.narea); obs_trans[sett].resize(data.narea);
		for(auto c = 0u; c < data.narea; c++){
			obs_trans[sett][c].resize(model.trans.size());
			for(auto tr = 0u; tr < model.trans.size(); tr++){
				obs_trans[sett][c][tr].resize(data.ndemocatpos);
			}
		}
	}
	
	obs_pop.resize(details.ndivision);
	for(auto sett = 0u; sett < details.ndivision; sett++){
		obs_pop[sett].resize(data.narea);
		for(auto c = 0u; c < data.narea; c++){
			obs_pop[sett][c].resize(model.comp.size());
			for(auto co = 0u; co < model.comp.size(); co++){
				obs_pop[sett][c][co].resize(data.ndemocatpos);
			}			
		}
	}
	
	for(auto i = 0u; i < data.nobs; i++){
		auto &ob = data.obs[i];
		
		auto dt = ob.datatable;
		
		for(auto sett = ob.sett_i; sett < ob.sett_f; sett++){
			for(auto c : ob.area){
				for(auto dp : ob.dp_sel){
					for(auto tr : data.datatable[dt].translist){					
						obs_trans[sett][c][tr][dp].push_back(i);
					}
					
					for(auto co : data.datatable[dt].complist){	
						obs_pop[sett][c][co][dp].push_back(i);
					}
				}
			}
		}
	}
}


/// Returns mean percentage error in observations
double ObservationModel::mean_percentage_error(const double invT) const 
{
	auto av = 0.0, nav = 0.0;
	for(const auto &gr : data.graph){
		if(gr.point.size() > 0){
			auto max = 0.0, per_max = 0.0;
			for(const auto &gp : gr.point){
				const auto &ob = data.obs[gp.obs];
				auto val = ob.value;
				if(val > max){
					max = val;
					per_max = sqrt(ob.var_approx/invT)/val;
				}
			}

			av += per_max; nav++;
		}
	}
	
	return 100*av/nav;
}


/// The error coming from a given observation
double ObservationModel::obs_prob(double value, const Observation& ob) const
{
	if(value < 0) value = 0;                                          // This accounts for small negative populations 
	
	auto val = ob.value;
	if(val == UNKNOWN) return 0;
	
	if(val == THRESH){
		auto thresh = data.datatable[ob.datatable].threshold;
		if(thresh == UNSET) emsgEC("ObsModel",4);
	
		switch(ob.obsmodel){
			case NORMAL_OBSMODEL: case NORMAL_PERCENT_OBSMODEL:
				{
					if(value < thresh) return normal_probability(thresh,thresh,ob.sd*ob.sd); 
					else return normal_probability(thresh,value,ob.sd*ob.sd); 
				}
				break;
				
			case POISSON_OBSMODEL:
				{
					auto lam = value; if(lam < 1) lam = 1;
					auto sum = 0.0; for(val = 0u; val <= thresh; val++) sum += exp(poisson_probability(val,lam));
					return log(sum);
				}
				break;
				
			case NEGBINO_OBSMODEL:
				{
					auto m = value; if(m < 1) m = 1;
					auto sum = 0.0; for(val = 0u; val <= thresh; val++) sum += exp(negative_binomial_probability(val,m,ob.shape));
					return log(sum);
				}
				break;
		}
	}
	else{
		switch(ob.obsmodel){
			case NORMAL_OBSMODEL: case NORMAL_PERCENT_OBSMODEL:
				{
					return normal_probability(val,value,ob.sd*ob.sd); 
				}
	
			case POISSON_OBSMODEL:
				{
					auto lam = value; if(lam < 1) lam = 1;
					return poisson_probability(int(val+0.5),lam);
				}
				
			case NEGBINO_OBSMODEL:
				{
					auto m = value; if(m < 1) m = 1;
					return negative_binomial_probability(int(val+0.5),m,ob.shape);
				}
		}
	}

	emsgEC("ObsModel",5);

	return 0;
}


/// Gets the error function split by data type	
vector <double> ObservationModel::get_EF_datatable(const State *state) const
{
	vector <double> EF_datatable(data.datatable.size());
	
	auto obs_value = get_obs_value(state);
		
	for(auto &val : EF_datatable) val = 0;
	
	for(auto i = 0u; i < data.nobs; i++){		
		EF_datatable[data.obs[i].datatable] += -2*obs_prob(obs_value[i],data.obs[i]);
	}
	
	return EF_datatable;
}


/// Splits up obsevations so that particle filtering can be performed
void ObservationModel::split_observations()
{
	nsection = (unsigned int)((details.ndivision+details.pmcmc_obs_period-1)/details.pmcmc_obs_period);

	section_obs.resize(nsection);
	for(auto s = 0u; s < nsection; s++){
		auto ti = s*details.pmcmc_obs_period;
		auto tf = ti + details.pmcmc_obs_period; if(tf > details.ndivision) tf = details.ndivision;
		
		section_ti.push_back(ti);
		section_tf.push_back(tf);
		
		for(auto i = 0u; i < data.nobs; i++){
			const auto &ob = data.obs[i];
			if(ob.sett_i < tf && ob.sett_f >= ti){ 
				if(ob.sett_i >= ti && ob.sett_f <= tf){
					section_obs[s].push_back(i);
				}
			}
		}
	}
	if(section_tf[nsection-1] != details.ndivision) emsgEC("Obsmodel",6);
}


// Used to order particles by EF
bool OS_ord (ObsSlice p1, ObsSlice p2)                      
{ return (p1.sett < p2.sett);};  


/// Converts the observations made on the system into a series of time slices
vector <ObsSlice> ObservationModel::generate_obs_slice() const
{
	vector <ObsSlice> obs_slice;

	for(auto o = 0u; o < data.nobs; o++){
		const auto &ob = data.obs[o];
		
		const auto &dt = data.datatable[ob.datatable];
		
		unsigned int sett = 0;
		
		switch(dt.type){
			case POP:	sett = ob.sett_i;	break;
			case TRANS:	sett = ob.sett_f; break;	
			case MARGINAL: sett = ob.sett_f; break;	
			case POPFRAC: sett = ob.sett_i; break;
		}
		
		auto i = 0u; while(i < obs_slice.size() && obs_slice[i].sett != sett) i++;
		if(i == obs_slice.size()){
			ObsSlice os; os.sett = sett;
			obs_slice.push_back(os);
		}
		obs_slice[i].obs_ref.push_back(o);
	}

	sort(obs_slice.begin(),obs_slice.end(),OS_ord);     
	
	return obs_slice;
}
