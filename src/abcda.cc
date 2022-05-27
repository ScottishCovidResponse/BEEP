/// Implements a version of ABC which uses the diffusion approximation

#include <iostream>
#include <fstream>
#include <algorithm>
#include <assert.h>
#include <math.h>

#include "abcda.hh"
#include "inputs.hh"
#include "output.hh"

using namespace std;

/// Initilaises the ABCDA class
ABCDA::ABCDA(const Details &details, const Data &data, const Model &model, Inputs &inputs, const Output &output, const ObservationModel &obsmodel, Mpi &mpi) : state(details,data,model,obsmodel), paramprop(details,data,model,output,mpi), details(details), data(data), model(model), output(output), obsmodel(obsmodel),mpi(mpi)
{	
	inputs.find_generation_or_cutoff_final_or_cpu_time(G,cutoff_final,cpu_time);
	
	inputs.find_nrun(nrun);
	inputs.find_nparticle(npart,Ntot,N,nrun,mpi.ncore);
	inputs.find_cor_max(cor_max,cor_max_last);
	part.resize(N);
	partcopy.resize(Ntot);	
	
	obs_slice = obsmodel.generate_obs_slice();
	
	state.initialise_approx();
}


/// Runs the inference algorithm
void ABCDA::run()
{	
	timer[TIME_ALG].start();
	auto param = model.sample_from_prior();  
	
	//cout << "\nSTART\n";
	//for(auto th = 0u; th < param.size(); th++) cout << model.param[th].name << "  par approx\n";
  //emsg("P");
	
	//state.simulate(param);
	//state.test_transmean_gradients(100);
	//emsg("gg");

	//param[3] = 4;// 6.3;
		//state.EF = -0.5*state.likelihood_approx(param,obs_slice,invT);
	
		//cout << param[3] << " " << state.EF  << "\n";
	
	//return;
	/*
	param[2] = 1;
		state.EF = -0.5*state.likelihood_approx(param,obs_slice,invT);
		cout << state.EF << " result\n";
		
	return;
	*/
	//param[3] = 4;
	//state.EF = -0.5*state.likelihood_approx(param,obs_slice,invT);
	
	//return;
	
	/*
	auto filesim = details.output_directory+"/op.txt";
	ofstream opsim(filesim);
	for(auto val = 1.0; val < 9; val+=0.1){
	
		param[3] = val;
		state.EF = -0.5*state.likelihood_approx(param,obs_slice,invT);
	
		opsim << param[3] << " " << state.EF  << "\n";
		
		cout << param[3] << " " << state.EF  << "\n";
	}
	timer[TIME_ALG].stop();
	return;
*/

	/*
		auto param = model.sample_from_prior();             // Samples parameters from the prior
		
		for(auto val = 1.0; val < 10; val++){
			param[3] = val;
			auto EF = -0.5*state.likelihood_approx(param,obs_slice);
			cout << val << " " << EF << " EF\n";
		}
		*/
			//param[3] = 6;
			//auto EF = -0.5*state.likelihood_approx(param,obs_slice);
			//cout << EF << " an\n";
		//emsg("don");
	/*
	auto param = model.sample_from_prior();  
	state.simulate(param);
	state.EF = -0.5*state.likelihood_approx(param,obs_slice);
	return;
	*/
	//for(auto th = 0; th < model.param.size(); th++) cout << th << " " << model.param[th].name << " \n";
	
	auto g = 0u;
	do{
		Generation gen;
		
		if(g == 0){                                             // For the initial generation sample states from the prior
			gen.EFcut = LARGE; 
			if(mpi.core == 0) cout << "Samping particles for the initial generation..." << endl;

			for(auto p = 0u; p < N; p++){
				auto param = model.sample_from_prior();             // Samples parameters from the prior
		
				state.EF = -0.5*state.likelihood_approx(param,obs_slice,1,DOUBLE);
				part[p] = state.create_particle((mpi.core*N+p)/npart);// Converts the state into a particle    
			}		
			nproposal = 0;
		}
		else{                                                   // In subsequent generations performs MBPs 
			timer[TIME_CUTOFF].start();		
			EF_cutoff(gen,part,partcopy);       
			
			timer[TIME_CUTOFF].stop();
		
			auto cmax = cor_max; if(g == G-1) cmax = cor_max_last;
			
			timer[TIME_PROP].start();	
			auto pup = COMBINE_UPDATE; if(g == 1) pup = COMBINE_DYNAMIC_UPDATE;
			nproposal = mcmc_updates(generation[generation.size()-1].param_samp,gen.EFcut,pup,cmax);
			timer[TIME_PROP].stop();
			
			if(mpi.core == 0) cout << paramprop.print_proposal_information(false) << endl; 
		}
	
		timer[TIME_GEN].start();
		store_sample(gen);                                      // Stores samples    
	
		mpi.exchange_samples(gen);	                            // Exchanges parameter samples across MPI cores
	
		output.set_generation_time(gen);                        // Sets the CPU time for the generation

		generation.push_back(gen);                              // Adds the new generation
		timer[TIME_GEN].stop();
	
		print_generation(gen,g);                        
	
		g++;
	}while(terminate_generation(g,generation[g-1].EFcut) == false);

	timer[TIME_ALG].stop();

	if(mpi.core == 0) model_evidence(generation);             // Calculates the model evidence
	
	paramprop.output_prop_vec();

	output.generation_results(generation);                    // Generates pdf of graphs

	output.generate_graphs(part,UNSET);	
	
	paramprop.diagnostics();                                  // Outputs diagnostic information
}


/// Performs the mcmc proposals to change the parameter values
unsigned int ABCDA::mcmc_updates(const vector <ParamSample> &param_samp, double EFcut, ParamUpdate pup, double cor_max)
{
	timer[TIME_UPDATE].start();    
	timer[TIME_MVNSETUP].start(); 
	for(auto &mv : paramprop.mvn) mv.setup(param_samp);
	timer[TIME_MVNSETUP].stop(); 
	
	auto nprop = 0u;   
	auto psamp_before = mpi.gather_psamp(part);                          // Stores parameter samples

	timer[TIME_MBPUPDATE].start();     
	vector <double> cor;                                                 // Stores the correlation with initial samples
	
	paramprop.zero_ntr_nac();	
	do{	
		for(auto p = 0u; p < N; p++){
			auto param = part[p].paramval;
			auto Pr = model.prior(param);
			
			for(auto &mv : paramprop.mvn){	
				auto param_prop = mv.propose(param);
				if(model.inbounds(param_prop) == true){
					auto Pr_prop = model.prior(param_prop);
					state.EF = -0.5*state.likelihood_approx(param_prop,obs_slice,1,DOUBLE);
				
					auto al = 0.0;
					if(state.EF < EFcut) al = exp(Pr_prop - Pr);
			
					if(mv.MH(al,pup) == SUCCESS){
						param = param_prop;
						Pr = Pr_prop;
						part[p] = state.create_particle((mpi.core*N+p)/npart);// Converts the state into a particle    
					}
				}
				nprop++;
			}
		}
	
		mpi.barrier();
		
		timer[TIME_UPDATEPROP].start();
		if(pup == COMBINE_DYNAMIC_UPDATE){      // For the first generation dynamically update during proposal
			paramprop.combine_update_proposals();
			paramprop.zero_ntr_nac();	
		}
		timer[TIME_UPDATEPROP].stop();
		
		timer[TIME_WAIT].start();
		auto psamp_after = mpi.gather_psamp(part);
	
		cor = output.get_correlation(psamp_before,psamp_after);
		
		if(mpi.core == 0){
			cout << vec_max(cor)  << " " << cor_max << "cor\n";
			//for(auto th = 0u; th < model.param.size(); th++){
			//	cout << model.param[th].name << " " << cor[th] << " " << cor_max << " cor\n";
			//}
		}
		
		timer[TIME_WAIT].stop();
	}while(nprop/part.size() < details.prop_max && vec_max(cor) > cor_max);  // Iterates until correlation less than threshold
	timer[TIME_MBPUPDATE].stop(); 
	
	timer[TIME_UPDATEPROP].start();
	if(pup == COMBINE_UPDATE) paramprop.combine_update_proposals();
	timer[TIME_UPDATEPROP].stop();

	timer[TIME_UPDATE].stop();  
			
	return nprop/part.size();
}

  
/// Determines when generations are stopped
bool ABCDA::terminate_generation(unsigned int g, double EFcut) const 
{
	if(cutoff_final != UNSET && EFcut <= cutoff_final) return true;
	
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
void ABCDA::store_sample(Generation &gen)
{
	for(const auto &pa : part){
		gen.param_samp.push_back(pa.create_param_samp());
		state.initialise_from_particle(pa);
		gen.EF_datatable.push_back(obsmodel.get_EF_datatable(&state));
	}
}


// Used to order particles by EF
bool PartEF_ord4 (PartEF p1,PartEF p2)                      
{ return (p1.EF < p2.EF); };  


/// Sets EF cut-off and works out which particles should be copied and which culled
void ABCDA::EF_cutoff(Generation &gen, vector <Particle> &part, vector <unsigned int> &partcopy)
{	
	auto N = part.size();
	
	vector <double> EF(N);                                    // Gathers EFs from all particles
	for(auto i = 0u; i < N; i++) EF[i] = part[i].EF; 
	auto EFtot = mpi.gather(EF);

	double EFcut, EFmin, EFmax;
	if(mpi.core == 0){
		vector <PartEF> partef(Ntot);
		for(auto i = 0u; i < Ntot; i++){ partef[i].i = i; partef[i].EF = EFtot[i];}
		
		sort(partef.begin(),partef.end(),PartEF_ord4);           // Sorts EFs
		
		auto mid = Ntot/2;
		EFcut = 0.5*(partef[mid-1].EF + partef[mid].EF);        // Sets cut-off so half particles copied and half discarded

		//if(cutoff_final != UNSET && EFcut < cutoff_final){      // Limits EFcut to cutoff_final
			//EFcut = cutoff_final;
		//}		

		EFmin = partef[int(0.025*Ntot)].EF; EFmax = partef[int(0.975*Ntot)].EF;
	
	/*
		for(i = 0; i < npart; i++) cout << EFtot[i] << ","; cout << "EFlist" << endl << flush;
		
		cout << EFcut << " EF cut\n";
		auto num = 0u;
		for(i = 0; i < npart; i++){
				if(EFtot[i] < EFcut) num++;
		}
		cout << npart << " " << num << "  num\n";
	*/
	
		if(true){   // Copies and particles
			for(auto j = 0u; j < Ntot/2; j++){ 
				auto i = partef[j].i;
				partcopy[i] = UNSET;
				auto ii = partef[j+Ntot/2].i;
				partcopy[ii] = i;
			}
			
			/*
			for(auto ru = 0u; ru < nrun; ru++){	
				for(auto j = 0u; j < npart; j++){
					auto i = ru*npart + j;
					if(EFtot[i] >= EFcut){ 
						if(list[ru].size() == 0) emsgEC("ABCDA",10);
						
						partcopy[i] = list[ru][list[ru].size()-1]; 
						list[ru].pop_back();
					}
				}
				if(list[ru].size() != 0) emsgEC("ABCDA",11);
			}
			*/
		}
		else{
			vector < vector <unsigned int> > list(nrun);
		
			for(auto ru = 0u; ru < nrun; ru++){
				for(auto j = 0u; j < npart; j++){
					auto i = ru*npart + j;
					if(EFtot[i] < EFcut){ 
						partcopy[i] = UNSET; list[ru].push_back(i);
					}
				}
			}
		
			for(auto ru = 0u; ru < nrun; ru++){
				auto ru_sel = ru;
				if(list[ru_sel].size() == 0){  // No valid state in this run so sample from another
					cout << "Sample from another run" << endl;
					do{
						ru_sel = (unsigned int)(ran()*nrun);
					}while(list[ru_sel].size() == 0);
				}
				
				for(auto j = 0u; j < npart; j++){
					auto i = ru*npart + j;
					if(EFtot[i] >= EFcut){ 
						partcopy[i] = list[ru_sel][(unsigned int)(ran()*list[ru_sel].size())]; 
					}
				}
			}
		}
	}
	mpi.bcast(EFcut); mpi.bcast(EFmin); mpi.bcast(EFmax);

	mpi.copy_particles(part,partcopy,N,Ntot);                 // Copies particles to be duplicated
	
	gen.EFcut = EFcut; gen.EFmin = EFmin; gen.EFmax = EFmax;
}


/// Calculates the log of the model evidence for each generation
void ABCDA::model_evidence(vector <Generation> &generation)
{
	vector <double> ME_final;
	for(auto run = 0u; run < nrun; run++){
		auto ME = 0.0;
		for(auto g = 0u; g < generation.size(); g++){
			auto &gen = generation[g];
			if(g > 0){
				auto num = 0.0, num_below = 0.0;
				for(auto gg = 0u; gg < g; gg++){
					for(const auto &ps : generation[gg].param_samp){
						if(ps.run == run){
							if(ps.EF < generation[g-1].EFcut){
								num++;
								if(ps.EF < gen.EFcut) num_below++;
							}
						}
					}
				}
				if(num_below == 0 || num == 0) ME += -1000;
				else ME += log(num_below/num);
			}
			gen.model_evidence.push_back(ME);
		}
		ME_final.push_back(ME);
	}
	
	if(cutoff_final != UNSET) output.final_model_evidence(ME_final,UNSET,cutoff_final);
}

	
/// Prints statistics from a generation
void ABCDA::print_generation(const Generation &gen, const unsigned int g) const
{
	double timetaken = timer[TIME_ALG].val/(60.0*CLOCKS_PER_SEC);	
	mpi.bcast(timetaken);
	
	if(mpi.core == 0){
		if(g == 0) cout << "Starting Generations..." << endl;
		else{
			cout << "Generation " << g;
		
			if(G != UNSET) cout << " / " << G-1;
			cout << " -     EF cut-off: " << prec(gen.EFcut,3);
			cout << " (" << prec(gen.EFmin,3) << " - " << prec(gen.EFmax,3) << ")";
			cout << "   Proposals: " << nproposal << endl;
		}
	}
}
