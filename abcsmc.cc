/// These function relate specifically to ABC-SMC

#include <algorithm>
#include <sstream> 

using namespace std;

#include "abcsmc.hh"
#include "timers.hh"
#include "output.hh"
#include "mpi.hh"


ABCSMC::ABCSMC(const Details &details, const Data &data, const Model &model, const Inputs &inputs, const Output &output, const ObservationModel &obsmodel, Mpi &mpi) : state(details,data,model,obsmodel), details(details), model(model), output(output), obsmodel(obsmodel), mpi(mpi)
{
	inputs.find_generation(G);
	inputs.find_nsample(Ntot);
	inputs.find_cutoff_frac(cutoff_frac);
	inputs.find_propsize(propsize);
}


/// This is an implementation of an ABC-SMC algorithm, used to compare against the ABC-MBP approach 
void ABCSMC::run()
{	
	MVN mvn("ABCSMC",model.param_not_fixed,propsize,ALL_PARAM,MULTIPLE); // Used for sampling from MVN distribution
	
	vector <Generation> generation; 
	
	for(auto g = 0u; g < G; g++){                                        // Sequentially goes through generations
		timer[TIME_ALG].start();
		
		auto acrate = 1.0;
		
		Generation gen;
		gen.time = clock();
		if(g == 0){                                                        // For the first generation sample states
			do{                                                             
				auto param = model.sample_from_prior();                        // Samples parameters from the prior
				
				state.simulate(param);                                         // Simulates the state
		
				store_sample(gen,g);                                           // Stores the sample                                    
			}while(!terminate(gen.EF_samp.size()));
		}
		else{                                                              // Subsequent generations sample from particles
			Generation &gen_last = generation[g-1];
		
			mvn.setup(gen_last.param_samp);                                  // Sets up MVN using covariance matrix from last gen
		
			setup_particle_sampler(gen_last);                                // Sets up sampler for particles
			
			auto ntr = 0u, nac = 0u;
			do{
				do{
					ntr++;
						
					auto p = particle_sampler();                                 // Samples from a particle in last generation
			
					auto param_prop = mvn.propose(gen_last.param_samp[p]);       // Proposes a new parameter set using MVN kernal
					
					if(model.inbounds(param_prop) == true){                      // Checks if parameters within bounds
						state.simulate(param_prop);                                // Simulates a new state
						
						if(state.EF < gen_last.EFcut){                             // Checks if error function less than last generation
							store_sample(gen,g);
							nac++;
							break;
						}
					}
				}while(true);			
			}while(!terminate(nac));
			acrate = mpi.get_acrate(nac,ntr);
		}
		
		implement_cutoff_frac(gen);                                        // Sets the cut-off based on cutoff_frac
		
		mpi.exchange_samples(gen,Ntot);	                                   // Copies parameter and EF samples across cores
		
		generation.push_back(gen);                                         // Stores the current generation
		
		calculate_particle_weight(generation,mvn);                         // Calculates the weights for the next generation
		
		timer[TIME_ALG].stop();
		
		print_generation(generation,acrate);                               // Outputs statistics about generation
		
		auto time_av = mpi.average(timer[TIME_TOTAL].val+clock());
		if(mpi.core == 0) cout <<  "Total time: " << prec(double(time_av)/(60.0*CLOCKS_PER_SEC),3) << " minutes." << endl;
	}
		
	output.generation_results(generation);                               // Generates pdf of graphs
	output.generate_graphs(particle_store);					                     // Outputs diagnostic information
}


/// Stores parameter and state sample 
void ABCSMC::store_sample(Generation &gen, unsigned int g)
{
	gen.param_samp.push_back(state.paramval);                            // Stores the parameter samples
	gen.EF_samp.push_back(state.EF);
	gen.EF_datatable.push_back(obsmodel.get_EF_datatable(&state));
	if(g == G-1) particle_store.push_back(state.create_particle());      // In last generation stores particles for plotting    
}

							
/// Selects EFcut to remove samples at a certain acceptance rate
void ABCSMC::implement_cutoff_frac(Generation &gen)
{
	auto EFtot = mpi.combine(gen.EF_samp);
	
	double EFcut;
	if(mpi.core == 0){
		sort(EFtot.begin(),EFtot.end());
		EFcut = EFtot[Ntot]; 
	}
	
	mpi.bcast(EFcut);

	auto p = 0u;
	while(p < gen.EF_samp.size()){
		if(gen.EF_samp[p] >= EFcut){
			gen.EF_samp.erase(gen.EF_samp.begin()+p);
			gen.EF_datatable.erase(gen.EF_datatable.begin()+p);
			gen.param_samp.erase(gen.param_samp.begin()+p);
		}
		else p++;
	}

	gen.EFcut = EFcut;
}


/// Determines when to terminate the algorithm
bool ABCSMC::terminate(long nac) const
{
	bool term = false;
	auto samptot = mpi.sum(nac);
	if(samptot >= Ntot/cutoff_frac) term = true;
	mpi.bcast(term);
	
	return term;
}


/// Generates a sampler to select a particle
void ABCSMC::setup_particle_sampler(const Generation &gen)
{
	wsum.resize(Ntot);
	auto sum = 0.0;
	for(auto i = 0u; i < Ntot; i++){ sum += gen.w[i]; wsum[i] = sum;}
}
			

/// Samples the next particle
unsigned int ABCSMC::particle_sampler() const
{
	double z = ran()*wsum[Ntot-1]; 
	auto p = 0u; while(p < Ntot && z > wsum[p]) p++;
	if(p == Ntot) emsgEC("ABCSMC",14);
	
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
void ABCSMC::calculate_particle_weight(vector <Generation> &generation, const MVN &mvn)
{
	unsigned int g = generation.size()-1;
	Generation &gen = generation[g];

	gen.w.resize(Ntot);		
		
	if(g == 0){
		for(auto i = 0u; i < Ntot; i++){
			 if(gen.EF_samp[i] >= gen.EFcut) emsg("Prob");
			gen.w[i] = 1;
		}
	}
	else{
		Generation &gen_last = generation[g-1];

		for(auto i = 0u; i < Ntot; i++){
			 if(gen.EF_samp[i] >= gen.EFcut) emsg("Prob");
		
			auto sum = 0.0;
			for(auto j = 0u; j < Ntot; j++){
				sum += gen_last.w[j]*exp(mvn.probability(gen.param_samp[i],gen_last.param_samp[j]));
			}
				
			gen.w[i] = exp(model.prior(gen.param_samp[i]))/sum;
		}	
	}
	
	auto wsum = 0.0; for(auto i = 0u; i < Ntot; i++) wsum += gen.w[i];
	wsum /= Ntot;
	for(auto i = 0u; i < Ntot; i++) gen.w[i] /= wsum;
}


/// Prints statistices from last generation
void ABCSMC::print_generation(const vector <Generation> &generation, double acrate) const
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
