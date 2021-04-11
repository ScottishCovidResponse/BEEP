/// Implements a version of ABC which uses model-based proposals in MCMC

#include <iostream>
#include <fstream>
#include <algorithm>
#include <assert.h>
#include <math.h>

#include "abcmbp.hh"
#include "inputs.hh"
#include "output.hh"

using namespace std;

/// Initilaises the ABCMBP class
ABCMBP::ABCMBP(const Details &details, const Data &data, const Model &model, const AreaTree &areatree, const Inputs &inputs, const Output &output, const ObservationModel &obsmodel, Mpi &mpi) : state(details,data,model,obsmodel), mbp(CUTOFF,details,data,model,areatree,obsmodel,output,mpi), paramprop(details,data,model,areatree,output,mpi), details(details), data(data), model(model), areatree(areatree), output(output), obsmodel(obsmodel),mpi(mpi)
{	
	inputs.find_generation(G);
	inputs.find_nparticle(Ntot,N,mpi.ncore);
	part.resize(N);
	partcopy.resize(Ntot);	
}


/// Runs the inference algorithm
void ABCMBP::run()
{				
	for(auto g = 0u; g < G; g++){
		timer[TIME_ALG].start();
		
		Generation gen; gen.time = clock();
		
		if(g == 0){                                             // For the first generation sample states from the prior
			gen.EFcut = LARGE; 
			if(mpi.core == 0) cout << "Samping particles for the first generation..." << endl;

			for(auto p = 0u; p < part.size(); p++){
				auto param = model.sample_from_prior();             // Samples parameters from the prior

				state.simulate(param);                              // Simulates a state
		
				part[p] = state.create_particle();                  // Converts the state into a particle    
			}
		}
		else{                                                   // In subsequent generations perform MBPs 
			EF_cutoff(gen,part,partcopy);                         // Sets EF cutoff so half the particles are cut and half copies
			nproposal = mbp.mcmc_updates(part,generation[generation.size()-1].param_samp,gen.EFcut,UNSET,NO_UPDATE,paramprop);  
		}
		
		timer[TIME_GEN].start();
		store_sample(gen);                                      // Stores samples    
	
		mpi.exchange_samples(gen,Ntot);	                        // Exchanges parameter samples across MPI cores

		generation.push_back(gen);                              // Adds the new generation
		timer[TIME_GEN].stop();
		timer[TIME_ALG].stop();
		
		print_generation(gen,g);                                // Outputs statistics about generation
	}
	
	if(mpi.core == 0) cout << "Generating extra samples in the last generation..." << endl;
	auto part_plot = part;
	for(auto loop = 0; loop < 5; loop++){                     // Performs extra samples in the last generation
		mbp.mcmc_updates(part,generation[G-1].param_samp,generation[G-1].EFcut,UNSET,NO_UPDATE,paramprop); 
		for(const auto &pa : part) part_plot.push_back(pa);
	}
	
	output.generation_results(generation);                    // Generates pdf of graphs
	output.generate_graphs(part_plot);	
	paramprop.diagnostics();                                  // Outputs diagnostic information
}
 
 
/// Stores a sample from each of the particles 
void ABCMBP::store_sample(Generation &gen)
{
	for(const auto &pa : part){
		gen.param_samp.push_back(pa.paramval);
		gen.EF_samp.push_back(pa.EF);
		state.initialise_from_particle(pa);
		gen.EF_datatable.push_back(obsmodel.get_EF_datatable(&state));
	}
}


// Used to order particles by EF
bool PartEF_ord (PartEF p1,PartEF p2)                      
{ return (p1.EF < p2.EF); };  


/// Sets EF cut-off and works out which particles should be copied and which culled
void ABCMBP::EF_cutoff(Generation &gen, vector <Particle> &part, vector <unsigned int> &partcopy)
{	
	auto N = part.size();
	
	vector <double> EF(N);                                    // Gathers EFs from all particles
	for(auto i = 0u; i < N; i++) EF[i] = part[i].EF; 
	auto EFtot = mpi.gather(EF);

	double EFcut, EFmin, EFmax;
	if(mpi.core == 0){
		vector <PartEF> partef(Ntot);
		for(auto i = 0u; i < Ntot; i++){ partef[i].i = i; partef[i].EF = EFtot[i];}
		
		sort(partef.begin(),partef.end(),PartEF_ord);           // Sorts EFs
		
		auto mid = Ntot/2;
		EFcut = 0.5*(partef[mid-1].EF + partef[mid].EF);        // Sets cut-off so half particles copied and half discarded
		EFmin = partef[int(0.025*Ntot)].EF; EFmax = partef[int(0.975*Ntot)].EF;
		
		for(auto j = 0u; j < Ntot; j++){
			auto i = partef[j].i;
			if(j < mid) partcopy[i] = UNSET;
			else partcopy[i] = partef[j-mid].i;
		}
	}
	mpi.bcast(EFcut); mpi.bcast(EFmin); mpi.bcast(EFmax);

	mpi.copy_particles(part,partcopy,N,Ntot);                 // Copies particles to be duplicated
	
	gen.EFcut = EFcut; gen.EFmin = EFmin; gen.EFmax = EFmax;
}


/// Prints statistics from a generation
void ABCMBP::print_generation(const Generation &gen, unsigned int g) const
{
	double timetaken = timer[TIME_ALG].val/(60.0*CLOCKS_PER_SEC);	
	mpi.bcast(timetaken);
	
	if(mpi.core == 0 && g > 0){
		cout << "Generation " << g << " -     EF cut-off: " << prec(gen.EFcut,3);
		cout << " (" << prec(gen.EFmin,3) << " - " << prec(gen.EFmax,3) << ")";
		cout << "   Proposals: " << nproposal << "   Time taken: " << prec(timetaken,3) << " minutes." << endl;
	}
}
	