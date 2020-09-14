/// Implements different types of ABC algorithms

#include <iostream>
#include <fstream>
#include <algorithm>

#include <assert.h>
#include <math.h>

#include "abc.hh"

#include "utils.hh"
#include "timers.hh"
#include "chain.hh"
#include "data.hh"
#include "output.hh"
#include "obsmodel.hh"
#include "pack.hh"
#include "timers.hh"

using namespace std;

struct PartEF                                                              // Structure used to order particle EFs
{
	unsigned int i;                                                          // The number of the particle
	double EF;                                                               // The error function
};

bool PartEF_ord (PartEF p1,PartEF p2) { return (p1.EF < p2.EF); }          // Used to order by EF

/// Initilaises the ABC class
ABC::ABC(const Details &details, const Data &data, const Model &model, const AreaTree &areatree, const Mpi &mpi, const Inputs &inputs, const Output &output, const ObservationModel &obsmodel) : chain(details,data,model,areatree,obsmodel,output,0), jump(chain.jump), details(details), data(data), model(model), areatree(areatree), mpi(mpi), output(output), obsmodel(obsmodel)
{	
	total_time = inputs.find_integer("cputime",UNSET);  
	
	G = inputs.find_integer("ngeneration",UNSET);                  
	
	if(total_time == UNSET && G == UNSET) emsgroot("The algorithm must be limited by either 'cputime' or 'generation'");

	Ntot = inputs.find_integer("nparticle",UNSET);                                            // Sets the total number of mcmc chains
	if(Ntot == UNSET) emsgroot("The number of particles must be set");
	if(Ntot%mpi.ncore != 0) emsgroot("The number of particles must be a multiple of the number of cores");

	N = Ntot/mpi.ncore;                                                           // The number of particles per core
	if(N == 0) emsgroot("'nparticle' must be non-zero");
}

/// Implements a version of abc which uses model-based proposals in MCMC
void ABC::mbp()
{
	chain.jump.setburnin(0,1);
	 
	vector <Particle> part;					
	part.resize(N);
	
	vector <unsigned int> partcopy(Ntot);

	vector <Generation> generation; 
	
	for(auto g = 0u; g < G; g++){
		timers.timeabc -= clock();

		Generation gen;
		gen.time = clock();
		
		if(g == 0){
			gen.EFcut = LARGE;
			
			for(auto& pa : part){ 
				chain.sample_state();
				
				chain.initial.EF = obsmodel.observation_likelihood(chain.initial.transev,chain.initial.indev);
				chain.initial.generate_particle(pa);
		
				gen.param_samp.push_back(chain.initial.paramval);
				gen.EF_samp.push_back(chain.initial.EF);
			}
		}
		else{
			gen.EFcut = next_generation_mpi(part,partcopy);
			
			mcmc_updates(gen,part,chain);

			if(mpi.core == 0) cout << "Generation " << g << ": EFcut " << gen.EFcut << endl;
		}

		generation.push_back(gen);

		exchange_samples_mpi(generation[g]);	
		jump.calculate_cholesky_matrix(generation[g].param_samp);
			
		timers.timeabc += clock();
		
		if(g%5 == 0){
			if(mpi.core == 0){
				string file = "generation_mbp.txt";
				output.generation_plot(file,generation);
			}
		}
		
		if(g%5 == 0){
			results_mpi(generation,part,chain);
			if(mpi.core == 0){
				cout << int((100.0*timers.timeabcprop)/timers.timeabc) << "% CPU time on proposals\n";
				cout << int((100.0*timers.timestandard)/timers.timeabc) << "% CPU time on standard\n";
			}
		}

		if(g%5 == 0 && mpi.core == 0){
			jump.output_M(details.output_directory+"/M.txt");
		}
		
		double timetaken = timers.timeabc/(60.0*CLOCKS_PER_SEC);	
		mpi_bcast(timetaken);
	
		if(timetaken > total_time) break;
	}
	
	if(mpi.core == 0) cout << int((100.0*timers.timeabcprop)/timers.timeabc) << "% CPU time on proposals\n";

	results_mpi(generation,part,chain);
	
	if(mpi.core == 0){
		string file = "model_evidence_mbp.txt";
		output.model_evidence_plot(file,generation);
	}	
}
 
/// This is an implementation of an ABC-SMC algorithm, which is used to compare against the MBP-MCMC approach 
void ABC::smc()
{	
	Chain chain(details,data,model,areatree,obsmodel,output,0);
	
	const double jumpsize = 1;
	
	vector <Generation> generation; 
	
	auto nparam = model.param.size();

	for(auto g = 0u; g < G; g++){
		timers.timeabc -= clock();
			
		Generation gen;
		gen.time = clock();
		if(g == 0){
			for(auto i = 0u; i < N; i++){     // For the first generation 
				chain.sample_state();
				
				gen.param_samp.push_back(chain.initial.paramval);
				gen.EF_samp.push_back(obsmodel.observation_likelihood(chain.initial.transev,chain.initial.indev));
			}
		}
		else{
			Generation &gen_last = generation[g-1];
			
			jump.calculate_cholesky_matrix(gen_last.param_samp);      // Generated matrix for sampling MVN samples
		
			double EFcut = gen_last.EFcut;

			double sumst[Ntot];	           // Generate particle sampler
			double sum = 0; for(auto i = 0u; i < Ntot; i++){ sum += gen_last.w[i]; sumst[i] = sum;}
			
			double EF=UNSET;
			auto ntr = 0u;
			auto nac = 0u;
			for(auto i = 0u; i < N; i++){ 
				unsigned int fl;
				do{
					fl = 0;
					double z = ran()*sum; unsigned int k = 0; while(k < Ntot && z > sumst[k]) k++;
					if(k == Ntot) emsg("Problem");
			
					vector <double> param_propose = gen_last.param_samp[k];
					jump.mvn_propose(param_propose,jumpsize);
					
					for(auto th = 0u; th < nparam; th++){
						if(param_propose[th] < model.param[th].min || param_propose[th] > model.param[th].max) fl = 1;
					}
				
					if(fl == 0){
						fl = chain.simulate(param_propose);
						if(fl == 0){
							EF = obsmodel.observation_likelihood(chain.propose.transev,chain.propose.indev);
							if(EF > EFcut) fl = 1;
						}
					}
					ntr++;
				}while(fl == 1);
				nac++;
				
				gen.param_samp.push_back(chain.initial.paramval);
				gen.EF_samp.push_back(EF);
			}
			
			if(mpi.core == 0) cout << "Generation " << g <<  ":  Cuttoff " << EFcut << "  Acceptance " << double(nac)/ntr << "  Neff " << effective_particle_number(gen_last.w) << endl;
		}
		
		exchange_samples_mpi(gen);	   // Copies parameter and EF samples across cores
		
		generation.push_back(gen);
		
		calculate_particle_weight(generation,jumpsize);  // Calculates the weights for the next generation
		
		timers.timeabc += clock();
			
		if(mpi.core == 0){
			string file = "generation_smc.txt";
			output.generation_plot(file,generation);
		}
		
		if(timers.timeabc/(60.0*CLOCKS_PER_SEC) > total_time) break;
	}

	if(mpi.core == 0){
		string file = "Posterior_parameters.txt";
		output.plot_distribution(file,generation[generation.size()-1]);
	}	
}

/// Updates particles using MBPs
void ABC::mcmc_updates(Generation &gen, vector <Particle> &part, Chain &chain)
{	
	const auto nvar = jump.nvar;
	const auto p_mbp = 0.9;
	const auto sampstep = 5u;
	const double beta = 2;
	const double facup = 1.2, facdown = 0.8;

	unsigned int ntr_v[nvar], nac_v[nvar];
	for(auto v = 0u; v < nvar; v++){ ntr_v[v] = 0; nac_v[v] = 0;}
	
	double EFcut = gen.EFcut;
	
	auto jmax = 0u;
	for(auto& ju : jump.jumpv){
		auto num = 1.0/(ju*ju);
		if(num < 1) num = 1;
		jmax += beta*num;
	}
	
	if(jmax < 1) jmax = 1;
	if(mpi.core == 0) cout << jmax << "jmax\n";
	
	for(auto& pa : part){
		chain.initial.initialise_from_particle(pa);

		for(auto j = 0u; j < jmax; j++){
			if(j%sampstep == 0){
				gen.param_samp.push_back(chain.initial.paramval);				
				gen.EF_samp.push_back(chain.initial.EF);
			}
			
			if(ran() < p_mbp){
				auto v = (unsigned int)(ran()*nvar);
				auto th = jump.param_not_fixed[v];
					
				vector <double> param_propose = chain.initial.paramval;
				param_propose[th] += normal_sample(0,jump.jumpv[v]*sqrt(jump.M[v][v]));
						
				timers.timeabcprop -= clock();
				ntr_v[v]++;
				if(chain.abcmbp_proposal(param_propose,EFcut) == 1) nac_v[v]++;
				timers.timeabcprop += clock();
			}
			else{
				chain.standard_proposal(EFcut);
			}
		}
			
		chain.initial.generate_particle(pa);
	}	
			
	for(auto v = 0u; v < nvar; v++){
		double ac_rate = acceptance(double(nac_v[v])/(double(ntr_v[v])+0.01));
		
		if(ac_rate > 0.4){ jump.jumpv[v] *= facup; if(jump.jumpv[v] > 2) jump.jumpv[v] = 2;}
		else{
			if(ac_rate < 0.3) jump.jumpv[v] *= facdown; 
		}
	}
}

/// Finds the effective number of particles
unsigned int ABC::effective_particle_number(vector <double> w)
{
	sort(w.begin(),w.end());

	double sum = 0; for(auto i = 0u; i < w.size(); i++) sum += w[i];
	double sum2 = 0; auto i = 0u; while(i < w.size() && sum2 < sum/2){ sum2 += w[i]; i++;}
	return 2*(1+w.size()-i); 
}

/// Calculates the weights for the different particles
void ABC::calculate_particle_weight(vector <Generation> &generation, double jumpsize)
{
	unsigned int g = generation.size()-1;
	Generation &gen = generation[g];

	auto Ntot = gen.param_samp.size();
	auto Ncut = int(0.5*Ntot);
		
	double ef[Ntot];
	for(auto i = 0u; i < Ntot; i++) ef[i] = gen.EF_samp[i];
	sort(ef,ef+Ntot);
	double EFcut = ef[Ncut];
	
	gen.EFcut = EFcut;
	gen.w.resize(Ntot);		
		
	if(g == 0){
		for(auto i = 0u; i < Ntot; i++){
			double w;
			if(gen.EF_samp[i] >= EFcut) w = 0; else w = 1;
			gen.w[i] = w;
		}
	}
	else{
		Generation &gen_last = generation[g-1];

		for(auto i = 0u; i < Ntot; i++){
			double w;
			if(gen.EF_samp[i] >= EFcut) w = 0;
			else{
				
				double sum = 0;
				for(auto j = 0u; j < Ntot; j++)	sum += gen_last.w[j]*jump.mvn_prob(gen.param_samp[i],gen_last.param_samp[j],jumpsize);
			
				w = exp(model.prior(gen.param_samp[i]))/sum;
			}
			gen.w[i] = w;
		}	
	}
	
	double wsum = 0;
	for(auto i = 0u; i < Ntot; i++) wsum += gen.w[i];
	wsum /= Ncut;
	for(auto i = 0u; i < Ntot; i++) gen.w[i] /= wsum;
}

/// Calculate a measure of how well a generation is mixed (by comparing the similarity between the two sets of copied particles)
double ABC::calculate_mixing(const vector <Particle> &part, vector <unsigned int> &partcopy) const
{
	auto nparam = model.param.size();
	vector <double> paramval(nparam*N);

	for(auto p = 0u; p < N; p++){ 
		for(auto th = 0u; th < nparam; th++){
			paramval[p*nparam+th] = part[p].paramval[th];
		}
	}
	auto paramvaltot = mpi_gather(paramval);

	double mix = 0;
	if(mpi.core == 0){
		vector < vector <double> > between;
		for(auto i = 0u; i < Ntot; i++){
			auto p = partcopy[i];
			if(p != UNSET){
				vector <double> vec; for(auto th = 0u; th < nparam; th++) vec.push_back(paramvaltot[p*nparam+th] - paramvaltot[i*nparam+th]);
				between.push_back(vec);
			}
		}
		vector <double> var_between	= jump.variance_vector(between);
	
		vector < vector <double> > across;
		for(auto i = 0u; i < 10*Ntot; i++){
			unsigned int p1, p2;
			p1 = (unsigned int)(Ntot*ran());
			do{	p2 = (unsigned int)(Ntot*ran());}while(p2 == p1);
			vector <double> vec; for(auto th = 0u; th < nparam; th++) vec.push_back(paramvaltot[p1*nparam+th] - paramvaltot[p2*nparam+th]);
			across.push_back(vec);
		}
		vector <double> var_across = jump.variance_vector(across);
		
		for(auto i = 0u; i < jump.nvar; i++) mix += sqrt(var_between[i]/var_across[i]);
		mix /= jump.nvar;
	}

	return mix;
}

/// Calculates the acceptance rate averaged across cores
double ABC::acceptance(double rate) const 
{
	auto rate_av = mpi_sum(rate);
	if(mpi.core == 0) rate_av /=  mpi.ncore;
	mpi_bcast(rate);

	return rate;	
}

/// Gets a sample from a particle
Sample ABC::get_sample(const Particle &part, Chain &chain) const
{
	Sample sample;
	chain.initial.initialise_from_particle(part);
	sample.meas = obsmodel.get_measured_quantities(chain.initial.transev,chain.initial.indev);;
	sample.R0 = model.calculate_R_vs_time(chain.initial.paramval);
	sample.phi = model.create_disc_spline(model.phispline_ref,chain.initial.paramval); 
	
	return sample;
}

/// Gathers together results from all the cores
void ABC::results_mpi(const vector <Generation> &generation, const vector <Particle> &part, Chain &chain) const
{
	auto g = generation.size()-1;
	const Generation &gen = generation[g];
	
	auto N = part.size();
	
	if(mpi.core == 0){
		vector <ParamSample> psamp;      // Stores parameter samples
		for(auto& samp : gen.param_samp){
			ParamSample psa; 
			psa.paramval = samp;
			psamp.push_back(psa);
		}
		
		vector <Sample> opsamp;        // Stores output samples
		for(const auto& pa : part) opsamp.push_back(get_sample(pa,chain));
		
		for(auto co = 1u; co < mpi.ncore; co++){
			pack_mpi_recv(co);
			for(auto i = 0u; i < N; i++){
				Particle part;
				unpack(part);
				opsamp.push_back(get_sample(part,chain));
			}
		}
		
		output.final_results(psamp,opsamp);
	}
	else{
		pack_initialise(0);
		for(auto i = 0u; i < N; i++) pack(part[i]);
		pack_mpi_send(0);
	}
}

/// Uses mpi to swap parameter and EF samples across all chains
void ABC::exchange_samples_mpi(Generation &gen)
{
	if(mpi.core == 0){
		for(auto co = 1u; co < mpi.ncore; co++){
			pack_mpi_recv(co);		
			vector< vector <double> > paramsa;
			unpack(paramsa);
			for(auto i = 0u; i < paramsa.size(); i++) gen.param_samp.push_back(paramsa[i]);
			
			vector <double> EFsa;
			unpack(EFsa);
			for(auto i = 0u; i < EFsa.size(); i++) gen.EF_samp.push_back(EFsa[i]);
		}
	}
	else{
		pack_initialise(0);
		pack(gen.param_samp);
		pack(gen.EF_samp);
		pack_mpi_send(0);
	}
	
	if(mpi.core == 0){
		pack_initialise(0);
		pack(gen.param_samp);
		pack(gen.EF_samp);
	}
	
	pack_mpi_bcast();
	
	if(mpi.core != 0){ unpack(gen.param_samp); unpack(gen.EF_samp);}
}

/// Uses mpi to get the next generation of particles. It returns the EF cutoff used
double ABC::next_generation_mpi(vector<Particle> &part, vector <unsigned int> &partcopy)
{
	auto N = part.size();
	
	vector <double> EF(N);
	for(auto i = 0u; i < N; i++) EF[i] = part[i].EF; 
	auto EFtot = mpi_gather(EF);

	double EFcut;
	if(mpi.core == 0){
		PartEF partef[Ntot];
		for(auto i = 0u; i < Ntot; i++){ partef[i].i = i; partef[i].EF = EFtot[i];}
		
		sort(partef,partef+Ntot,PartEF_ord);
		
		if(Ntot%2 != 0) emsg("The number of particles must be even.");
		
		auto mid = Ntot/2;
		EFcut = 0.5*(partef[mid-1].EF + partef[mid].EF);
		
		for(auto j = 0u; j < Ntot; j++){
			auto i = partef[j].i;
			if(j < mid) partcopy[i] = UNSET;
			else partcopy[i] = partef[j-mid].i;
		}
	}
	mpi_bcast(EFcut);
	mpi_bcast(partcopy);

	vector< vector<double> > sendbuffer;
	vector< vector<double> > recibuffer;
	vector <unsigned int> sendbuffer_to;
	vector <unsigned int> recibuffer_from;
	vector <unsigned int> recibuffer_to;
		
	vector <unsigned int> buffersize(N);

	for(auto i = 0u; i < N; i++){
		buffersize[i] = 0;
		
		auto itot = mpi.core*N+i;
		if(partcopy[itot] == UNSET){
			for(auto iitot = 0u; iitot < Ntot; iitot++){
				if(partcopy[iitot] == itot && iitot/N != itot/N){
					pack_initialise(0);
					pack(part[i]);
		
					sendbuffer.push_back(copybuffer()); 
					sendbuffer_to.push_back(iitot);
				
					buffersize[i] = packsize();
					//cout << "send " << itot << " -> " << iitot << " \n";
				}				
			}
		}
		else{
			auto iitot = partcopy[itot];
			if(iitot/N == itot/N){
				part[i] = part[iitot%N]; //cout << "internal copy " <<iitot << " -> " << itot << "\n";
			}
			else{
				recibuffer.push_back(vector <double> ());
				recibuffer_from.push_back(iitot);
				recibuffer_to.push_back(itot);
				//cout << "receive " << iitot << " -> " << itot << " \n"; 
			}
		}
	}
	
	auto buffersizetot = mpi_gather(buffersize);
	mpi_bcast(buffersizetot);
	
	auto buftot = recibuffer.size()+sendbuffer.size();
	if(buftot > 0){
		MPI_Request reqs[buftot];                 // These are information used Isend and Irecv
		MPI_Status stats[buftot];
		unsigned int nreqs = 0;

		for(auto j = 0u; j < recibuffer.size(); j++){
			auto from = recibuffer_from[j];
			auto to = recibuffer_to[j];
			recibuffer[j].resize(buffersizetot[from]);
			MPI_Irecv(&recibuffer[j][0],recibuffer[j].size(),MPI_DOUBLE,from/N,to,MPI_COMM_WORLD,&reqs[nreqs]); nreqs++;
		}
		
		for(auto j = 0u; j < sendbuffer.size(); j++){
			auto to = sendbuffer_to[j];
			MPI_Isend(&sendbuffer[j][0],sendbuffer[j].size(),MPI_DOUBLE,to/N,to,MPI_COMM_WORLD,&reqs[nreqs]); nreqs++;		
		}	
	
		if(MPI_Waitall(nreqs,reqs,stats) != MPI_SUCCESS) emsgEC("ABC",1);
		
		for(auto j = 0u; j < recibuffer.size(); j++){
			setbuffer(recibuffer[j]);
			unpack(part[recibuffer_to[j]%N]);
		}
	}
	
	return EFcut;
}
