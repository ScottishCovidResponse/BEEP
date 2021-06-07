/// These function relate specifically to ABC-SMC

#include <algorithm>
#include <sstream> 

using namespace std;

#include "abcsmc.hh"
#include "timers.hh"
#include "output.hh"
#include "mpi.hh"


ABCSMC::ABCSMC(const Details &details, const Data &data, const Model &model, Inputs &inputs, const Output &output, const ObservationModel &obsmodel, Mpi &mpi) : state(details,data,model,obsmodel), details(details), model(model), output(output), obsmodel(obsmodel), mpi(mpi)
{
	inputs.find_nrun(nrun);
	inputs.find_nsample_GRmax(Ntot,GRmax,nrun);
	inputs.find_generation_or_cutoff_final(G,cutoff_final);
	inputs.find_cutoff_frac(cutoff_frac);
	inputs.find_propsize(propsize);
}


/// This is an implementation of an ABC-SMC algorithm, used to compare against the ABC-MBP approach 
void ABCSMC::run()
{	
	MVN mvn("ABCSMC",model.param_not_fixed,propsize,ALL_PARAM,MULTIPLE); // Used for sampling from MVN distribution
	
	for(auto g = 0u; g < G; g++){                                        // Sequentially goes through generations
		timer[TIME_ALG].start();
		
		auto acrate = 1.0;
		
		particle_store.clear();
		
		Generation gen;
		gen.time = clock();
		if(g == 0){                                                        // For the initial generation sample states
			do{          
				for(auto ru = 0u; ru < nrun; ru++){
					auto param = model.sample_from_prior();                      // Samples parameters from the prior
						
					state.simulate(param);                                       // Simulates the state
				
					store_sample(gen,g,ru,1);                                    // Stores the sample   
				}
			}while(!terminate(gen));
		}
		else{                                                              // Subsequent generations sample from particles
			Generation &gen_last = generation[g-1];
	
			mvn.setup(gen_last.param_samp,gen_last.w);                       // Sets up MVN using param samps from last gen
			
			setup_particle_sampler(gen_last);                                // Sets up sampler for particles
		
			auto ntr = 0u, nac = 0u;
			vector <double> wtot(nrun), wcut(nrun);                          // Sets up quantities for estimating model evidence
      for(auto ru = 0u; ru < nrun; ru++){ wtot[ru] = 0; wcut[ru] = 0;}      
			do{
				for(auto ru = 0u; ru < nrun; ru++){
					do{
						ntr++;
							
						auto p = particle_sampler(ru);                             // Samples from a particle in last generation
			
						auto param_prop = mvn.propose(gen_last.param_samp[p].paramval);// Proposes a new parameter set using MVN kernal
						
						
						if(model.inbounds(param_prop) == true){                    // Checks if parameters within bounds
							state.simulate(param_prop);                              // Simulates a new state

							auto w = calculate_particle_weight(param_prop,gen_last,ru,mvn); // Calculates the weight for the sample
	
							wtot[ru] += w;
							if(state.EF < gen_last.EFcut){                           // Checks if error function less than cutoff
								store_sample(gen,g,ru,w);
								nac++; wcut[ru] += w;
								break;
							}
						}
					}while(true);			
				}
			}while(!terminate(gen));
			
			for(auto ru = 0u; ru < nrun; ru++) gen.model_evidence.push_back(log(mpi.get_ratio(wcut[ru],wtot[ru])));
			acrate = mpi.get_acrate(nac,ntr);
		}
		
		implement_cutoff_frac(gen);                                        // Sets the cut-off based on cutoff_frac
		
		mpi.exchange_samples(gen);	                                       // Copies parameter samples and EF samples across cores
		
		normalise_particle_weights(gen);                                   // Normalises the particle weights
		
		generation.push_back(gen);                                         // Stores the current generation
		
		timer[TIME_ALG].stop();
		
		print_generation(generation,acrate);                               // Outputs statistics about generation
		
		auto time_av = mpi.average(timer[TIME_TOTAL].val+clock());
		if(false && mpi.core == 0) cout <<  "Total time: " << prec(double(time_av)/(60.0*CLOCKS_PER_SEC),3) << " minutes." << endl;

		if(gen.EFcut == cutoff_final) break;                               // Terminates if final EF is reached
	}
		
	if(mpi.core == 0) print_model_evidence();                            // Prints the final model evidence

	results();                                                           // Generates pdf of graphs
}


/// Stores parameter and state sample 
void ABCSMC::store_sample(Generation &gen, const unsigned int g, const unsigned int run, const double w)
{
	gen.param_samp.push_back(state.create_param_sample(run));            // Stores the parameter samples
	gen.EF_datatable.push_back(obsmodel.get_EF_datatable(&state));
	gen.w.push_back(w);
	if(g == G-1 || cutoff_final != UNSET){                               // In last generation stores particles for plotting  
		particle_store.push_back(state.create_particle(run));   
	}
}

							
/// Selects EFcut to remove samples at a certain acceptance rate
void ABCSMC::implement_cutoff_frac(Generation &gen)
{
	vector <double> EF_samp; for(const auto &ps : gen.param_samp) EF_samp.push_back(ps.EF);
	auto EFtot = mpi.combine(EF_samp);
	
	double EFcut;
	if(mpi.core == 0){
		sort(EFtot.begin(),EFtot.end());
		EFcut = EFtot[(unsigned int)(cutoff_frac*EFtot.size())]; 
		if(cutoff_final != UNSET && EFcut < cutoff_final) EFcut = cutoff_final;
	}
	
	mpi.bcast(EFcut);

	auto p = 0u;
	while(p < gen.param_samp.size()){
		if(gen.param_samp[p].EF >= EFcut){
			//wcut[gen.param_samp[p].run] -= gen.w[p];
			gen.EF_datatable.erase(gen.EF_datatable.begin()+p);
			gen.w.erase(gen.w.begin()+p);
			gen.param_samp.erase(gen.param_samp.begin()+p);
		}
		else p++;
	}

	p = 0;
	while(p < particle_store.size()){
		if(particle_store[p].EF >= EFcut){
			particle_store.erase(particle_store.begin()+p);
		}
		else p++;
	}
	 
	gen.EFcut = EFcut;
}


/// Determines when to terminate the algorithm
bool ABCSMC::terminate(const Generation &gen) const
{
	bool term = false;
	auto samptot = mpi.sum(long(gen.w.size()));
	
	if(GRmax != UNSET){
		auto psamp_tot = mpi.gather_psamp(gen.param_samp);
		auto w_tot = mpi.gather(gen.w);
		if(mpi.core == 0){
			auto GR = output.get_Gelman_Rubin_statistic(psamp_tot,w_tot,nrun);
			if(vec_max(GR) < GRmax && samptot >= nrun*20) term = true; 
			cout << "Number of samples: " << samptot << "    Largest GR value:" << vec_max(GR) << "     GRmax: " << GRmax << endl;
		}
	}
	else{
		if(mpi.core == 0){
			if(samptot >= nrun*Ntot/cutoff_frac) term = true;
		}
	}
	mpi.bcast(term);
	
	return term;
}


/// Generates a sampler to select a particle
void ABCSMC::setup_particle_sampler(const Generation &gen)
{
	auto NN = gen.w.size();
	wsum.resize(nrun);
	for(auto ru = 0u; ru < nrun; ru++){
		wsum[ru].resize(NN);
		auto sum = 0.0;
		auto num = 0u;
		for(auto i = 0u; i < NN; i++){ 
			if(gen.param_samp[i].run == ru){ sum += gen.w[i]; num++;}
			wsum[ru][i] = sum;
		}
		if(num == 0) emsgEC("ABCSMC",1);
	}
}
			

/// Samples the next particle
unsigned int ABCSMC::particle_sampler(const unsigned int run) const
{
	auto NN = wsum[run].size();
	double z = ran()*wsum[run][NN-1]; 
	auto p = 0u; while(p < NN && z > wsum[run][p]) p++;
	if(p == NN) emsgEC("ABCSMC",2);
	
	return p;
}
	
	
/// Finds the effective number of particles based on the particle weights
unsigned int ABCSMC::effective_particle_number(const vector <double> &w) const 
{
	auto sum = 0.0; for(auto i = 0u; i < w.size(); i++) sum += w[i];
	auto sum2 = 0.0; for(auto i = 0u; i < w.size(); i++) sum2 += (w[i]/sum)*(w[i]/sum);
	return (unsigned int)(TINY+1.0/sum2);
}


/// Calculates the weights for the different particles
double ABCSMC::calculate_particle_weight(const vector <double> param_prop, const Generation &gen_last, const unsigned int run, const MVN &mvn)
{	
	auto NN = gen_last.param_samp.size();
	
	auto sum = 0.0;
	for(auto j = 0u; j < NN; j++){
		const auto ps_last = gen_last.param_samp[j];
		if(ps_last.run == run){
			sum += gen_last.w[j]*exp(mvn.probability(param_prop,ps_last.paramval));
		}
	}
	
	return exp(model.prior(param_prop))/sum;
}	


/// Normalises the weights in a given generation 	
void ABCSMC::normalise_particle_weights(Generation &gen)
{
	auto NN = gen.param_samp.size();
		
	for(auto ru = 0u; ru < nrun; ru++){                                      // Normalises the results
		auto wsum = 0.0;
		auto num = 0u;
		for(auto i = 0u; i < NN; i++){ if(gen.param_samp[i].run == ru){ wsum += gen.w[i]; num++;}}
		wsum /= num;
		for(auto i = 0u; i < NN; i++){ if(gen.param_samp[i].run == ru) gen.w[i] /= wsum;}
	}
}

	
/// Generates pdf of graphs
void ABCSMC::results()
{
	output.generation_results(generation);   
	
	auto particle = mpi.gather_particle(particle_store);
	
	// Performs a bootstrap step which selects particles in proportion to their weight
	vector <Particle> particle_plot;
	if(mpi.core == 0){                      
		const Generation &gen_last = generation[generation.size()-1];
		
		setup_particle_sampler(gen_last);
		for(auto i = 0u; i < particle.size(); i++){
			auto p = particle_sampler(i%nrun);
			particle_plot.push_back(particle[p]);
		}
	}

	output.generate_graphs(particle_plot);					                     // Outputs diagnostic information
}


/// Outputs the model evidence
void ABCSMC::print_model_evidence()
{
	const auto &gen = generation[generation.size()-1];
	output.final_model_evidence(gen.model_evidence,UNSET,gen.EFcut);
}


/// Prints statistices from last generation
void ABCSMC::print_generation(const vector <Generation> &generation, const double acrate) const
{
	if(mpi.core == 0){
		auto g = generation.size()-1;
		const Generation &gen_last = generation[g]; 
		stringstream ss; ss << std::fixed; ss.precision(2);
		ss << "Generation " << g <<  " -     EF Cut-off: " << gen_last.EFcut << "     Acceptance rate: " << 100*acrate 
				 << "%     Effective number of particles: " << effective_particle_number(gen_last.w) << endl;
		cout << ss.str();
	}
}
