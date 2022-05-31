/// Implements a particle version of annealed importance sampling (PAS)

#include <iostream>
#include <fstream>
#include <algorithm>
#include <assert.h>
#include <math.h>
#include <sstream>

using namespace std;

#include "pas.hh"
#include "inputs.hh"
#include "output.hh"

/// Initilaises the PAS class
PAS::PAS(const Details &details, const Data &data, const Model &model, Inputs &inputs, const Output &output, const ObservationModel &obsmodel, Mpi &mpi) : state(details,data,model,obsmodel), mbp(INVT,details,data,model,obsmodel,output,mpi), paramprop(details,data,model,output,mpi), details(details), data(data), model(model), output(output), obsmodel(obsmodel),mpi(mpi)
{	
	inputs.find_generation_or_invT_final_or_cpu_time(G,invT_final,cpu_time);
	inputs.find_nrun(nrun);
	inputs.find_nparticle(npart,Ntot,N,nrun,mpi.ncore);
	inputs.find_quench_factor(quench_factor);
	inputs.find_cor_max(cor_max,cor_max_last);
	part.resize(N);
	partcopy.resize(Ntot);	
	
	//if(mpi.core == 0) invT_vs_EF.open(details.output_directory+"/invT_vs_EF.csv");
}


/// Runs the inference algorithm
void PAS::run()
{
	timer[TIME_ALG].start();
		
	auto g = 0u;
	do{	
		Generation gen;
		if(g == 0){                                             // For the first generation sample states from the prior
			gen.invT = 0;                                         // Starts at zero inverse temperature (i.e. the prior)
			
			if(mpi.core == 0) cout << "Samping particles for the first generation..." << endl;

			for(auto p = 0u; p < N; p++){
				auto param = model.sample_from_prior();             // Samples parameters from the prior

				state.simulate(param);                              // Simulates a state
		
				part[p] = state.create_particle((mpi.core*N+p)/npart);// Converts the state into a particle    
			}
			nproposal = 0;
		}
		else{                                                   // In subsequent generations filter particles and perform MCMC
			Generation &gen_last = generation[generation.size()-1];
			       
			bootstrap(gen,part,partcopy,gen_last.invT);           // Performs a bootstrap step 

			auto cmax = cor_max; if(g == G-1) cmax = cor_max_last;
			     
			auto pup = COMBINE_UPDATE; if(g == 1) pup = COMBINE_DYNAMIC_UPDATE;
			nproposal = mbp.mcmc_updates(part,gen_last.param_samp,gen.EFcut,gen.invT,pup,paramprop,cmax); 
		}
		
		timer[TIME_GEN].start();
		store_sample(gen);                                      // Stores a sample    
	
		mpi.exchange_samples(gen);	                            // Exchanges parameter samples across MPI cores

		output.set_generation_time(gen);                        // Sets the CPU time for the generation
	
		generation.push_back(gen);                              // Adds the new generation
		timer[TIME_GEN].stop();
		
		print_generation(gen,g);                                // Outputs statistics about generation
		
		g++;
	}while(terminate_generation(g,generation[g-1].invT) == false);
	timer[TIME_ALG].stop();
		
	if(mpi.core == 0) model_evidence(generation);             // Calculates the model evidence 

	output.generation_results(generation);                    // Generates visualisation
	output.generate_graphs(part,generation[g-1].invT);	
	paramprop.diagnostics();                                  // Outputs diagnostic information
}
 

/// Determines when generations are stopped
bool PAS::terminate_generation(unsigned int g, double invT) const 
{
	if(invT == invT_final) return true;
	
	if(g == G-1) return true;
	
	auto time_av = mpi.average((clock() - details.time_start)/(60.0*CLOCKS_PER_SEC));
	if(time_av > cpu_time){
		if(mpi.core == 0)	cout << "Maximum execution time reached." << endl << flush;
		return true;
	}
	
	if(nproposal >= details.prop_max){
		if(mpi.core == 0) cout << "Maximum proposal number reached!" << endl << flush;
		return true;
	}
	
	return false;
}

 
/// Stores a sample from each of the particles 
void PAS::store_sample(Generation &gen)
{
	for(const auto &pa : part){
		gen.param_samp.push_back(pa.create_param_samp());
		state.initialise_from_particle(pa);
		gen.EF_datatable.push_back(obsmodel.get_EF_datatable(&state));
	}
}


// Used to order particles by EF
bool PartEF_ord2 (PartEF p1,PartEF p2)                      
{ return (p1.EF < p2.EF); };  


/// Works out which particles should be copied and sets new inverse temperature
void PAS::bootstrap(Generation &gen, vector<Particle> &part, vector <unsigned int> &partcopy, const double invT)
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
		
		sort(partef.begin(),partef.end(),PartEF_ord2);          // Sorts EFs
		//EFmin = partef[int(0.025*Ntot)].EF; EFmax = partef[int(0.975*Ntot)].EF;
		EFmin = partef[int(0.25*Ntot)].EF; EFmax = partef[int(0.75*Ntot)].EF;
		
		//invT_vs_EF << invT << "," << partef[int(0.75*Ntot)].EF << "," << partef[int(0.975*Ntot)].EF << endl;
		
		auto EF_range = partef[int(0.5*Ntot)].EF - partef[int(0.25*Ntot)].EF;
		
		auto muav = 0.0;                                        // Calculates the average in EF across runs
		for(auto val : EFtot) muav += val;
		muav /= EFtot.size();
		
		//DinvT = 0.5*quench_factor/EF_range;                     // Sets the increase in inverse temperature
	  DinvT = 0.25*quench_factor/EF_range;                     // Sets the increase in inverse temperature
	                    
		if(invT_final != UNSET){
			if(invT+DinvT > invT_final) DinvT = invT_final-invT;   // Limits invT to invT_final
		}
		
		vector <double> wsum(npart);                             // Weights particles
		vector <unsigned int> num(npart);   
		for(auto ru = 0u; ru < nrun; ru++){
			auto sum = 0.0;
			for(auto j = 0u; j < npart; j++){
				sum += exp(-0.5*DinvT*(EFtot[ru*npart + j]-muav));
				wsum[j] = sum;
			}
		
			for(auto j = 0u; j < npart; j++) num[j] = 0;            // Performs a bootstrap step

			for(auto j = 0u; j < npart; j++){                       // Samples in proportion to particle weights
				auto z = ran()*sum;
				auto sel=0u; while(sel < npart && z > wsum[sel]) sel++;
				if(sel == npart) emsgEC("PAS",1);
				num[sel]++;
			}
		
			vector <unsigned int> list;                             // Works out which particles to copy and which to discard
			for(auto j = 0u; j < npart; j++){     
				if(num[j] > 1){ for(auto k = 1u; k < num[j]; k++) list.push_back(ru*npart+j);}
			}
		
			if(false) cout << npart - list.size() << " / "<< npart << " Particles kept" << endl;
			 
			for(auto j = 0u; j < npart; j++){  
				if(num[j] > 0) partcopy[ru*npart + j] = UNSET;
				else{
					if(list.size() == 0) emsgEC("PAS",2);
					partcopy[ru*npart + j] = list[list.size()-1];
					list.pop_back();
				}
			}
			if(list.size() != 0) emsgEC("PAS",3);
		}
		
		if(false){
			cout << invT << " invT" << endl;
			for(auto i = 0u; i < Ntot; i++){ 
				cout << i << " " << exp(-DinvT*EFtot[i]) << " " << num[i] << " " << partcopy[i] << " weight" << endl;
			}
		}
	}

	mpi.bcast(DinvT); mpi.bcast(EFmin); mpi.bcast(EFmax);

	mpi.copy_particles(part,partcopy,N,Ntot);                // Copies particles between mpi processes

	gen.invT = invT+DinvT; gen.EFmin = EFmin; gen.EFmax = EFmax;
}


/// Calculates the log of the model evidence for each generation
void PAS::model_evidence(vector <Generation> &generation)
{
	vector <double> ME_final;
	for(auto run = 0u; run < nrun; run++){
		auto ME = 0.0;
		for(auto g = 0u; g < generation.size(); g++){
			auto &gen = generation[g];

			gen.model_evidence.push_back(ME);
		
			if(g < generation.size()-1){
				auto DinvT = generation[g+1].invT - gen.invT;
				
				auto mu = 0.0, n = 0.0;
				for(auto ps : gen.param_samp){
					if(ps.run == run){ mu += ps.EF; n++;}
				}
				mu /= n;
			
				auto sum = 0.0;
				for(auto ps : gen.param_samp){
					if(ps.run == run) sum += exp(-0.5*DinvT*(ps.EF-mu));
				}
				sum /= n;
				
				ME += -0.5*DinvT*mu + log(sum);
			}
		}
		ME_final.push_back(ME);
	}

	if(invT_final != UNSET) output.final_model_evidence(ME_final,invT_final,UNSET);
}


/// Prints statistics from a generation
void PAS::print_generation(const Generation &gen, const unsigned int g) const
{
	double timetaken = timer[TIME_ALG].val/(60.0*CLOCKS_PER_SEC);	
	mpi.bcast(timetaken);
	
	if(mpi.core == 0 && g > 0){
		stringstream ss; 
		ss << std::scientific; ss.precision(2);
		ss << "Generation " << g << " -   invT: " << prec(gen.invT,3) << "    ";
		ss << "EF 25% to 75% range (" << prec(gen.EFmin,3) << " - " <<  prec(gen.EFmax,3) << ")    ";
		ss << "Proposals: " << nproposal << "    ";
		//ss << "Time taken: " << prec(timetaken,3) << " minutes.";
		cout << ss.str() << endl;
	}
}
	