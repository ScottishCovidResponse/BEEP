/// Implements a particle version of annealed importance sampling (PAIS)

#include <iostream>
#include <fstream>
#include <algorithm>
#include <assert.h>
#include <math.h>
#include <sstream>

using namespace std;

#include "pais.hh"
#include "inputs.hh"
#include "output.hh"

/// Initilaises the PAIS class
PAIS::PAIS(const Details &details, const Data &data, const Model &model, const AreaTree &areatree, const Inputs &inputs, const Output &output, const ObservationModel &obsmodel, Mpi &mpi) : state(details,data,model,obsmodel), mbp(INVT,details,data,model,areatree,obsmodel,output,mpi), paramprop(details,data,model,areatree,output,mpi), details(details), data(data), model(model), areatree(areatree), output(output), obsmodel(obsmodel),mpi(mpi)
{	
	inputs.find_generation_or_invT_final(G,invT_final);
	inputs.find_nparticle(Ntot,N,mpi.ncore);
	inputs.find_quench_factor(quench_factor);
	part.resize(N);
	partcopy.resize(Ntot);	
}


/// Runs the inference algorithm
void PAIS::run()
{
	for(auto g = 0u; g < G; g++){
		timer[TIME_ALG].start();
		
		Generation gen;
		gen.time = clock();
		
		if(g == 0){                                             // For the first generation sample states from the prior
			gen.invT = 0;                                         // Starts at zero inverse temperature (i.e. the prior)
			
			if(mpi.core == 0) cout << "Samping particles for the first generation..." << endl;

			for(auto p = 0u; p < part.size(); p++){
				auto param = model.sample_from_prior();             // Samples parameters from the prior

				state.simulate(param);                              // Simulates a state
		
				part[p] = state.create_particle();                  // Converts the state into a particle    
			}
		}
		else{                                                   // In subsequent generations filter particles and perform MCMC
			Generation &gen_last = generation[generation.size()-1];
			       
			bootstrap(gen,part,partcopy,gen_last.invT);           // Performs a bootstrap step 

			                                                      // Performs MBP-MCMC updates
			nproposal = mbp.mcmc_updates(part,gen_last.param_samp,gen.EFcut,gen.invT,NO_UPDATE,paramprop); 
		}
		
		timer[TIME_GEN].start();
		store_sample(gen);                                      // Stores a sample    
	
		mpi.exchange_samples(gen,Ntot);	                        // Exchanges parameter samples across MPI cores

		generation.push_back(gen);                              // Adds the new generation
		timer[TIME_GEN].stop();
		
		timer[TIME_ALG].stop();
		
		print_generation(gen,g);                                // Outputs statistics about generation
		
		if(gen.invT == invT_final) break;                       // Terminates if final invT reached
	}
	
	output.generation_results(generation);                    // Generates pdf of graphs
	output.generate_graphs(part);	
	paramprop.diagnostics();                                  // Outputs diagnostic information
}
 
 
/// Stores a sample from each of the particles 
void PAIS::store_sample(Generation &gen)
{
	for(const auto &pa : part){
		gen.param_samp.push_back(pa.paramval);
		gen.EF_samp.push_back(pa.EF);
		state.initialise_from_particle(pa);
		gen.EF_datatable.push_back(obsmodel.get_EF_datatable(&state));
	}
}


// Used to order particles by EF
bool PartEF_ord2 (PartEF p1,PartEF p2)                      
{ return (p1.EF < p2.EF); };  


/// Works out which particles should be copied and sets new inverse temperature
void PAIS::bootstrap(Generation &gen, vector<Particle> &part, vector <unsigned int> &partcopy, double invT)
{
	auto N = part.size();

	vector <double> EF(N);                                    // Gathers EFs from all particles
	for(auto i = 0u; i < N; i++) EF[i] = part[i].EF; 
	auto EFtot = mpi.gather(EF);

	double EFmin, EFmax;
	double DinvT;
	if(mpi.core == 0){
		vector <PartEF> partef(Ntot);
		for(auto i = 0u; i < Ntot; i++){ partef[i].i = i; partef[i].EF = EFtot[i];}
		
		sort(partef.begin(),partef.end(),PartEF_ord2);           // Sorts EFs
		EFmin = partef[int(0.025*Ntot)].EF; EFmax = partef[int(0.975*Ntot)].EF;
		
		auto av = 0.0, av2 = 0.0;                               // Calculate the mean and variance of EFtot
		for(auto val : EFtot){ av += val; av2 += val*val;}
		av /= Ntot; av2 /= Ntot;
		auto var = av2-av*av;
		
		DinvT = quench_factor/sqrt(var);                        // Sets the increase in inverse temperature
		if(invT_final != UNSET){
			if(invT+DinvT > invT_final) DinvT = invT_final-invT;  // Limits invT to invT_final
		}
		
		vector <double> wsum(Ntot);                             // Weights particles
		auto sum = 0.0;
		for(auto i = 0u; i < Ntot; i++){ 
			sum += exp(-0.5*DinvT*(EFtot[i]-av));
			wsum[i] = sum;
		}
		
		vector <unsigned int> num(Ntot);                        // Performs a bootstrap step
		for(auto i = 0u; i < Ntot; i++) num[i] = 0;
		
		for(auto i = 0u; i < Ntot; i++){                        // Samples in proportion to particle weights
			auto z = ran()*sum;
			auto sel=0u; while(sel < Ntot && z > wsum[sel]) sel++;
			if(sel == Ntot) emsg("problem");
			num[sel]++;
		}
		
		vector <unsigned int> list;                             // Works out which particles to copy and which to discard
		for(auto i = 0u; i < Ntot; i++){ 
			if(num[i] > 1){ for(auto j = 1u; j < num[i]; j++) list.push_back(i);}
		}
		
		if(false) cout << Ntot -list.size() << " / "<< Ntot << " Particles kept\n";
			 
		for(auto i = 0u; i < Ntot; i++){
			if(num[i] > 0) partcopy[i] = UNSET;
			else{
				if(list.size() == 0) emsgEC("PAIS",34);
				partcopy[i] = list[list.size()-1];
				list.pop_back();
			}
		}
		if(list.size() != 0) emsgEC("PAIS",35);
		
		if(false){
			cout << invT << " invT\n";
			for(auto i = 0u; i < Ntot; i++){ 
				cout << i << " " << exp(-DinvT*EFtot[i]) << " " << num[i] << " " << partcopy[i] << " weight\n";
			}
		}
	}

	mpi.bcast(DinvT); mpi.bcast(EFmin); mpi.bcast(EFmax);

	mpi.copy_particles(part,partcopy,N,Ntot);                // Copies particles between mpi processes

	gen.invT = invT+DinvT; gen.EFmin = EFmin; gen.EFmax = EFmax;
}


/// Prints statistics from a generation
void PAIS::print_generation(const Generation &gen, unsigned int g) const
{
	double timetaken = timer[TIME_ALG].val/(60.0*CLOCKS_PER_SEC);	
	mpi.bcast(timetaken);
	
	if(mpi.core == 0 && g > 0){
		stringstream ss; 
		ss << std::scientific; ss.precision(2);
		ss << "Generation " << g << " -   invT: " << prec(gen.invT,3) << "    ";
		ss << "EF 95% range (" << prec(gen.EFmin,3) << " - " <<  prec(gen.EFmax,3) << ")    ";
		ss << "Proposals: " << nproposal << "    ";
		ss << "Time taken: " << prec(timetaken,3) << " minutes.";
		cout << ss.str() << endl;
	}
}
	