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
ABCMBP::ABCMBP(const Details &details, const Data &data, const Model &model, Inputs &inputs, const Output &output, const ObservationModel &obsmodel, Mpi &mpi) : state(details,data,model,obsmodel), mbp(CUTOFF,details,data,model,obsmodel,output,mpi), paramprop(details,data,model,output,mpi), details(details), data(data), model(model), output(output), obsmodel(obsmodel),mpi(mpi)
{	
	inputs.find_generation_or_cutoff_final(G,cutoff_final);
	inputs.find_nrun(nrun);
	inputs.find_GRmax_nupdate(GRmax,nrun,nupdate);
	inputs.find_nparticle(npart,Ntot,N,nrun,mpi.ncore);
	part.resize(N);
	partcopy.resize(Ntot);	
}


/// Runs the inference algorithm
void ABCMBP::run()
{				
	//model.check_prior();
	
	for(auto g = 0u; g < G; g++){
		timer[TIME_ALG].start();
		
		Generation gen; gen.time = clock();
		
		if(g == 0){                                             // For the initial generation sample states from the prior
			gen.EFcut = LARGE; 
			if(mpi.core == 0) cout << "Samping particles for the initial generation..." << endl;

			for(auto p = 0u; p < N; p++){
				auto param = model.sample_from_prior();             // Samples parameters from the prior
		
				state.simulate(param);                              // Simulates a state
		
				part[p] = state.create_particle((mpi.core*N+p)/npart);// Converts the state into a particle    
			}		
		}
		else{                                                   // In subsequent generations performs MBPs 
			EF_cutoff(gen,part,partcopy);                         // Sets EF cutoff so half the particles are cut and half copied

			param_GR_clear();
			do{
				nproposal += mbp.mcmc_updates(part,generation[generation.size()-1].param_samp,gen.EFcut,UNSET,NO_UPDATE,paramprop);  
			}while(!terminate());
		}

		timer[TIME_GEN].start();
	
		store_sample(gen);                                      // Stores samples    

		mpi.exchange_samples(gen);	                            // Exchanges parameter samples across MPI cores

		generation.push_back(gen);                              // Adds the new generation
		timer[TIME_GEN].stop();
		timer[TIME_ALG].stop();
		
		print_generation(gen,g);                                // Outputs statistics about generation

		if(gen.EFcut == cutoff_final) break;                    // Terminates if final EF is reached
	}
	
	if(mpi.core == 0) model_evidence(generation);             // Calculates the model evidence
	
	output.generation_results(generation);                    // Generates pdf of graphs
	output.generate_graphs(part);	
	paramprop.diagnostics();                                  // Outputs diagnostic information
}
 
 
/// Stores a sample from each of the particles 
void ABCMBP::store_sample(Generation &gen)
{
	for(const auto &pa : part){
		gen.param_samp.push_back(pa.create_param_samp());
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

		if(cutoff_final != UNSET && EFcut < cutoff_final){      // Limits EFcut to cutoff_final
			EFcut = cutoff_final;
		}		

		EFmin = partef[int(0.025*Ntot)].EF; EFmax = partef[int(0.975*Ntot)].EF;
		
		for(auto ru = 0u; ru < nrun; ru++){
			vector <unsigned int> list;
			for(auto j = 0u; j < npart; j++){
				auto i = ru*npart + j;
				if(EFtot[i] < EFcut){ 
					partcopy[i] = UNSET; list.push_back(i);
				}
			}
			if(list.size() == 0) emsgEC("ABCMBP",1);
			
			for(auto j = 0u; j < npart; j++){
				auto i = ru*npart + j;
				if(EFtot[i] >= EFcut){ 
					partcopy[i] = list[(unsigned int)(ran()*list.size())]; 
				}
			}
		}
	}
	mpi.bcast(EFcut); mpi.bcast(EFmin); mpi.bcast(EFmax);

	mpi.copy_particles(part,partcopy,N,Ntot);                 // Copies particles to be duplicated
	
	gen.EFcut = EFcut; gen.EFmin = EFmin; gen.EFmax = EFmax;
}


/// Clears the parameter samples
void ABCMBP::param_GR_clear()
{
	nproposal = 0; loop = 0; psamp_GR.clear();
}


/// Calculates when to terminate algorithm
bool ABCMBP::terminate()
{
	auto term = false;
	loop++;
	
	if(GRmax != UNSET){                                      // Uses a Gelman Rubin statistics to determine termination
		auto psamp_tot = mpi.gather_psamp(part);
	
		if(mpi.core == 0){
			for(const auto &ps : psamp_tot) psamp_GR.push_back(ps);

			auto GR = output.get_Gelman_Rubin_statistic(psamp_GR);
			if(vec_max(GR) < GRmax) term = true; 
		
			if(false){
				for(auto th = 0u; th < model.param.size(); th++) cout << model.param[th].name << " " << GR[th] << endl;
			}
		}
	}
	else{
		if(nupdate != UNSET){                                  // Terminates after a specified number of updates   
			if(loop == nupdate) term = true; 
		}
	}
	
	mpi.bcast(term);

	return term;
}


/// Calculates the log of the model evidence for each generation
void ABCMBP::model_evidence(vector <Generation> &generation)
{
	vector <double> ME_final;
	for(auto run = 0u; run < nrun; run++){
		auto ME = 0.0;
		for(auto g = 0u; g < generation.size(); g++){
			auto &gen = generation[g];
	
			if(g > 0){
				auto num = 0.0, num_below = 0.0;
				for(auto gg = 0u; gg < g; gg++){
					for(auto ps : generation[gg].param_samp){
						if(ps.run == run){
							if(ps.EF < generation[g-1].EFcut){
								num++;
								if(ps.EF < gen.EFcut) num_below++;
							}
						}
					}
				}
				ME += log(num_below/num);
			}
			gen.model_evidence.push_back(ME);
		}
		ME_final.push_back(ME);
	}

	if(cutoff_final != UNSET) output.final_model_evidence(ME_final,UNSET,cutoff_final);
}


/// Prints statistics from a generation
void ABCMBP::print_generation(const Generation &gen, const unsigned int g) const
{
	double timetaken = timer[TIME_ALG].val/(60.0*CLOCKS_PER_SEC);	
	mpi.bcast(timetaken);
	
	if(mpi.core == 0 && g > 0){
		cout << "Generation " << g << " -     EF cut-off: " << prec(gen.EFcut,3);
		cout << " (" << prec(gen.EFmin,3) << " - " << prec(gen.EFmax,3) << ")";
		cout << "   Proposals: " << nproposal << endl;
	}
}
