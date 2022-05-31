/// Information and routines for transferring data between cores using MPI

#include <sstream>

using namespace std;

#include "mpi.hh"
#include "param_prop.hh"
#include "state.hh"

Mpi::Mpi(const Details &details): details(details)
{
	#ifdef USE_MPI
	int num;
	MPI_Comm_size(MPI_COMM_WORLD,&num); ncore = (unsigned int) num;
  MPI_Comm_rank(MPI_COMM_WORLD,&num); core = (unsigned int) num;
	#endif
	
	#ifndef USE_MPI
	ncore = 1;
	core = 0;
	#endif
}

/// Copies data from core zero to all the others
void Mpi::copy_data(unsigned int &narea, vector <Area> &area, unsigned int &nobs, vector <Observation> &obs, GenerateQ &genQ, vector <Modification> &modification, vector <Covariate> &covar, LevelEffect &level_effect, vector <DemocatChange> &democat_change)
{
	if(core == 0){                                     				        // Copies the above information to all the other cores
		pack_initialise(0);
		pack(narea);
		pack(area);
		pack(nobs);	for(const auto &ob : obs)	pack(ob);
		pack(genQ);
		
		pack((unsigned int) modification.size()); 
		for(const auto &cf : modification) pack(cf);
		
		for(const auto &cov : covar) pack(cov.value);
		
		pack(level_effect.param_map);
		pack(level_effect.frac);
		
		pack((unsigned int) democat_change.size()); 
		for(const auto &dcc	: democat_change){ 
			pack(dcc.d);
			pack(dcc.area);
			pack(dcc.frac);
			pack(dcc.dp_group);
		}
	}

	pack_bcast();

	if(core != 0){
		unpack(narea);
		unpack(area);
		unpack(nobs); obs.resize(nobs); for(auto &ob : obs) unpack(ob);
		unpack(genQ);
		
		unsigned int nmodification; unpack(nmodification); modification.resize(nmodification); 
		for(auto &cf : modification)	unpack(cf);
		
		for(auto &cov : covar) unpack(cov.value);
		
		unpack(level_effect.param_map);
		unpack(level_effect.frac);
		
		unsigned int ndemocat_change;	unpack(ndemocat_change); democat_change.resize(ndemocat_change);
		for(auto &dcc	: democat_change){ 
			unpack(dcc.d);
			unpack(dcc.area);
			unpack(dcc.frac);
			unpack(dcc.dp_group);
		}
		unpack_check();
	}
}


/// Copies particles across MPI processes (ABC-MBP, ABC-MBP, PAS)
void Mpi::copy_particles(vector<Particle> &part, vector <unsigned int> &partcopy, const unsigned int N, const unsigned int Ntot)
{
	bcast(partcopy);
	
	vector<Particle> part_internal_copy(N);
	
	vector< vector<double> > sendbuffer;
	vector< vector<double> > recibuffer;
	vector <unsigned int> sendbuffer_to;
	vector <unsigned int> recibuffer_from;
	vector <unsigned int> recibuffer_to;
		
	vector <unsigned int> buffersize(N);

	for(auto i = 0u; i < N; i++){
		buffersize[i] = 0;
		
		auto itot = core*N+i;
		for(auto iitot = 0u; iitot < Ntot; iitot++){
			if(partcopy[iitot] == itot){
				if(iitot/N == itot/N){
					if(iitot != itot) part_internal_copy[i] = part[i];
				}
				else{
					pack_initialise(0);
					pack(part[i]);
		
					sendbuffer.push_back(copybuffer()); 
					sendbuffer_to.push_back(iitot);
				
					buffersize[i] = packsize();
					if(false) cout << "send " << itot << " -> " << iitot << endl;
				}
			}				
		}
	}
	
	for(auto i = 0u; i < N; i++){
		auto itot = core*N+i;
		auto iitot = partcopy[itot];
		if(iitot != UNSET){
			if(iitot/N == itot/N){
				if(i != iitot%N){
					part[i] = part_internal_copy[iitot%N];
					if(false) cout << "internal copy " <<iitot << " -> " << itot << endl;
				}
			}
			else{
				recibuffer.push_back(vector <double> ());
				recibuffer_from.push_back(iitot);
				recibuffer_to.push_back(itot);
				if(false) cout << "receive " << iitot << " -> " << itot << endl; 
			}
		}
	}
	
	auto buffersizetot = gather(buffersize);
	bcast(buffersizetot);
	
	auto buftot = recibuffer.size()+sendbuffer.size();
	if(buftot > 0){
		MPI_Request reqs[buftot];                                       // These store information used in Isend and Irecv
		MPI_Status stats[buftot];
		auto nreqs = 0u;

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
	
		if(MPI_Waitall(nreqs,reqs,stats) != MPI_SUCCESS) emsgEC("Mpi",1);
		
		for(auto j = 0u; j < recibuffer.size(); j++){
			setbuffer(recibuffer[j]);
			unpack(part[recibuffer_to[j]%N]);
		}
	}
}


/// Exchanges samples within a generation (used in ABC-SMC, ABC-MBP and PAS)
void Mpi::exchange_samples(Generation &gen)
{
	if(core == 0){
		for(auto co = 1u; co < ncore; co++){
			pack_recv(co);		
			vector<ParamSample> paramsa;
			unpack(paramsa);
			for(auto i = 0u; i < paramsa.size(); i++) gen.param_samp.push_back(paramsa[i]);
			
			vector <double> wsa;
			unpack(wsa);
			for(auto i = 0u; i < wsa.size(); i++) gen.w.push_back(wsa[i]);
			
			vector < vector <double> > EFDT;
			unpack(EFDT);
			for(auto i = 0u; i < EFDT.size(); i++) gen.EF_datatable.push_back(EFDT[i]);
			unpack_check();
		}
	}
	else{
		pack_initialise(0);
		pack(gen.param_samp);
		pack(gen.w);
		pack(gen.EF_datatable);
		pack_send(0);
	}
	
	if(ncore > 1){
		if(core == 0){
			pack_initialise(0);
			pack(gen.param_samp);
			pack(gen.w);
		}
	
		pack_bcast();
	
		if(core != 0){ 
			unpack(gen.param_samp); 
			unpack(gen.w);
			unpack_check();
		}
	}
}


/// At the end of an observed section, this swaps information between different particles (PMCMC)
void Mpi::end_sec_swap(vector <State> &particle, const unsigned int sett, const vector <unsigned int> &backpart, const unsigned int buffersize)
{
	timer[TIME_PMCMCSWAP].start();
		
	auto N = particle.size();
	auto Ntot = backpart.size();
	
	MPI_Request reqs[Ntot];                                        // These store information used in Isend and Irecv
	MPI_Status stats[Ntot];
	
	vector < vector<double> > sendbuffer;                          // Buffers used for Isend and Irecv
	vector < vector<double> > recibuffer;
	sendbuffer.resize(N); recibuffer.resize(N);

	if(false && core == 0){
		for(auto p = 0u; p < Ntot; p++) cout << p << " " << backpart[p] << "  backpart" << endl;
	}
	
	auto npreclist = 0u;
	vector <unsigned int> preclist;
	vector < vector <unsigned int> > preclistli;
	
	if(false){ for(auto p = 0u; p < N; p++){ particle[p].pop[sett][0][0][0] = core*N+p;}} // For testing
	
	auto nreqs = 0u;
	auto pmin = core*N;                                             // Sets up information to be recieved
	for(auto p = pmin; p < pmin+N; p++){
		auto pp = backpart[p];
		if(p != pp){                   
			auto cor = pp/N;
			if(cor != core){
				auto j = 0u; while(j < npreclist && preclist[j] != pp) j++;
				if(j < npreclist) preclistli[j].push_back(p%N);
				else{
					preclist.push_back(pp); 
					preclistli.push_back(vector <unsigned int> ()); 
					preclistli[npreclist].push_back(p%N);
					
					recibuffer[npreclist].resize(buffersize);
					MPI_Irecv(&recibuffer[npreclist][0],buffersize,MPI_DOUBLE,cor,pp,MPI_COMM_WORLD,&reqs[nreqs]);
					nreqs++; if(nreqs > Ntot) emsgEC("Mpi",2); 
					npreclist++; if(npreclist > N) emsgEC("Mpi",3); 
				}
			}
		}
	}

	auto nsendbuf = 0u;
	for(auto p = pmin; p < pmin+N; p++){                           // Initiates information to be sent
		if(p == backpart[p]){ 
			auto ncorlist = 0u;
			vector <unsigned int> corlist;
			for(auto pp = 0u; pp < Ntot; pp++){
				if(p == backpart[pp]){
					auto cor = pp/N;
					if(cor == core){
						if(pp != p){
							particle[pp%N].pop[sett] = particle[p%N].pop[sett];
							particle[pp%N].Imap[sett-1] = particle[p%N].Imap[sett-1];
							particle[pp%N].Idiag[sett-1] = particle[p%N].Idiag[sett-1];
						}
					}
					else{
						auto j = 0u; while(j < ncorlist && corlist[j] != cor) j++;
						if(j == ncorlist){ corlist.push_back(cor); ncorlist++;}
					}
				}
			}
						
			if(ncorlist > 0){
				pack_initialise(0);
				pack(particle[p%N].pop[sett]); 
				pack(particle[p%N].Imap[sett-1]);
				pack(particle[p%N].Idiag[sett-1]);
				sendbuffer[nsendbuf] = copybuffer(); if(sendbuffer[nsendbuf].size() != buffersize) emsgEC("Mpi",4);
				
				for(auto j = 0u; j < ncorlist; j++){
					MPI_Isend(&sendbuffer[nsendbuf][0],buffersize,MPI_DOUBLE,corlist[j],p,MPI_COMM_WORLD,&reqs[nreqs]);
					nreqs++; if(nreqs > Ntot) emsgEC("Mpi",5); 
				}				
				nsendbuf++; if(nsendbuf > N) emsgEC("Mpi",6);
			}
		}
	}
	
	if(nreqs > 0){
		if(MPI_Waitall(nreqs,reqs,stats) != MPI_SUCCESS) emsgEC("Mpi",7);
	}
		
	for(auto rec = 0u; rec < npreclist; rec++){                       	// Unpacks the recieved information
		auto p = preclistli[rec][0];
		
		setbuffer(recibuffer[rec]);
		unpack(particle[p].pop[sett]);
		unpack(particle[p].Imap[sett-1]);
		unpack(particle[p].Idiag[sett-1]);
	
		for(auto j = 1u; j < preclistli[rec].size(); j++){
			auto pp = preclistli[rec][j];
			particle[pp].pop[sett] = particle[p].pop[sett];
		}
	}
	
	if(false){
		for(auto p = 0u; p < N; p++){ 
			if(particle[p].pop[sett][0][0][0] != int(backpart[core*N+p])) emsgEC("Mpi",8);
			cout << core << ": " << " " << core*N+p << " " << particle[p].pop[sett][0][0][0] << " " << backpart[core*N+p] << " compare" << endl;
		}
		emsg("Done");
	}
	
	timer[TIME_PMCMCSWAP].stop();
}


/// Constructs a particle by gathering the pieces from different particles (from the bootstrap function)
Particle Mpi::particle_sample(const unsigned int ru, const vector < vector <unsigned int> > &backpart, vector <State> &particle, const ObservationModel &obsmodel)	
{
	timer[TIME_STATESAMPLE].start();
			
	auto N = particle.size();
	auto Ntot = backpart[0].size();
	auto nsec = obsmodel.nsection;
	
	if(false){
		for(auto sec = 1u; sec < nsec; sec++){
			cout << sec << ": "; for(auto p = 0u; p < Ntot; p++) cout << backpart[sec][p] << ",";
			cout << " Backpart" << endl;  
		}
	}

	auto p = backpart[nsec][0];
	for(int sec = nsec-1; sec >= 0; sec--){
		auto ti = obsmodel.section_ti[sec];
		auto tf = obsmodel.section_tf[sec];
		
		unsigned int co = p/N;
		if(core == 0){  
			if(co != 0){
				pack_recv(co);
				for(auto sett = ti; sett < tf; sett++){
					unpack(particle[0].pop[sett]);
					unpack(particle[0].transnum[sett]);
				}
				unpack_check();
			}
			else{
				if(p != 0){
					for(auto sett = ti; sett < tf; sett++){
						particle[0].pop[sett] = particle[p].pop[sett];
						particle[0].transnum[sett] = particle[p].transnum[sett];
					}
				}
			}
		}
		else{
			if(co == core){
				pack_initialise(0);
				for(auto sett = ti; sett < tf; sett++){
					pack(particle[p%N].pop[sett]);
					pack(particle[p%N].transnum[sett]);
				}
				pack_send(0);
			}
		}
		barrier();
		
		p = backpart[sec][p];
	}
	
//	particle_sample
	timer[TIME_STATESAMPLE].stop();
	
	return particle[0].create_particle(ru);
}


/// Gathers a double vector across all cores and returns the combined vector to core 0
vector <double> Mpi::gather(const vector <double> &vec)
{
	vector <double> vectot(vec.size()*ncore);
	
	MPI_Gather(&vec[0],vec.size(),MPI_DOUBLE,&vectot[0],vec.size(),MPI_DOUBLE,0,MPI_COMM_WORLD);
	
	return vectot;
}


/// Gathers a 2D array of doubles across all cores and returns the combined vector to core 0
vector < vector <double> > Mpi::gather(const vector < vector <double> > &array)
{
	auto Y = array.size();
	auto X = array[0].size();
	vector <double> vectrans;
	vector <double> vectot(Y*X*ncore);

	for(auto j = 0u; j < Y; j++){
		for(auto i = 0u; i < X; i++) vectrans.push_back(array[j][i]);
	}
	
	MPI_Gather(&vectrans[0],vectrans.size(),MPI_DOUBLE,&vectot[0],vectrans.size(),MPI_DOUBLE,0,MPI_COMM_WORLD);
	
	vector < vector <double> > arraytot;	

	if(core == 0){
		for(auto j = 0u; j < Y*ncore; j++){
			vector <double> vec(X);
			for(auto i = 0u; i < X; i++) vec[i] = vectot[j*X+i];
			arraytot.push_back(vec);
		}
	}
		
	return arraytot;
}


/// Gathers a long across all cores and returns the combined vector to core 0
vector <long> Mpi::gather(const long val)
{
	vector <long> valtot;
	valtot.resize(ncore);
	
	MPI_Gather(&val,1,MPI_LONG,&valtot[0],1,MPI_LONG,0,MPI_COMM_WORLD);
	
	return valtot;
}


/// Gathers a double across all cores and returns the combined vector to core 0
vector <double> Mpi::gather(const double val)
{
	vector <double> valtot;
	valtot.resize(ncore);

	MPI_Gather(&val,1,MPI_DOUBLE,&valtot[0],1,MPI_DOUBLE,0,MPI_COMM_WORLD);
	return valtot;
}


/// Scatters a double vector from core 0 across all cores 
vector <double> Mpi::scatter(const vector <double> &vectot)
{
	vector <double> vec;
	vec.resize(vectot.size()/ncore);
	
	MPI_Scatter(&vectot[0],vec.size(),MPI_DOUBLE,&vec[0],vec.size(),MPI_DOUBLE,0,MPI_COMM_WORLD);
	
	return vec;
}


/// Scatters a double vector from core 0 across all cores 
vector <unsigned int> Mpi::scatter(const vector <unsigned int> &vectot)
{
	vector <unsigned int> vec;
	vec.resize(vectot.size()/ncore);
	
	MPI_Scatter(&vectot[0],vec.size(),MPI_UNSIGNED,&vec[0],vec.size(),MPI_UNSIGNED,0,MPI_COMM_WORLD);
	
	return vec;
}


/// Sums up a long value over all cores  
long Mpi::sum(const long val)
{
	auto vec = gather(val);
	long sum = 0; for(auto val : vec) sum += val;
	return sum;	
}


/// Sums up a double value over all cores  
double Mpi::sum(const double val)
{
	auto vec = gather(val);
	auto sum = 0.0;
	if(core == 0){
		for(auto val : vec) sum += val;
	}
	return sum;	
}


/// Calculates the average of a quantities across cores
double Mpi::average(const double val) 
{
	double val_av = sum(val)/ncore;
	bcast(val_av);

	return val_av;	
}


/// Calculates the average of a vector across cores
vector <double> Mpi::average(const vector <double> &val) 
{
	vector <double> result;
	for(auto i = 0u; i < val.size(); i++){
		result.push_back(average(val[i])); 
	}
	return result;	
}


/// Gets the acceptance rate across all mpi processes
double Mpi::get_acrate(const double nac, const double ntr)
{
	return average(nac)/(average(ntr)+0.01);
}


/// Gets the ratio across all mpi processes
double Mpi::get_ratio(const double nac, const double ntr)
{
	return sum(nac)/sum(ntr);
}


/// Gets a vector of acceptance rates across all mpi processes
vector <double> Mpi::get_acrate(const vector <double> &nac, const vector <double> &ntr)
{
	vector <double> result;
	for(auto i = 0u; i < nac.size(); i++) result.push_back(get_acrate(nac[i],ntr[i]));
	return result;
}

/// Copies a variable in core 0 to all the other cores
void Mpi::bcast(double &val)
{
	MPI_Bcast(&val,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
}


/// Copies a variable in core 0 to all the other cores
void Mpi::bcast(bool &val)
{
	MPI_Bcast(&val,1,MPI_CXX_BOOL,0,MPI_COMM_WORLD);
}


/// Copies a variable in core 0 to all the other cores
void Mpi::bcast(unsigned int &val)
{
	MPI_Bcast(&val,1,MPI_UNSIGNED,0,MPI_COMM_WORLD);
}


/// Copies a vector in core 0 to all the other cores
void Mpi::bcast(vector <unsigned int> &vec)
{
	MPI_Bcast(&vec[0],vec.size(),MPI_UNSIGNED,0,MPI_COMM_WORLD);
}
	

/// Copies a vector in core 0 to all the other cores
void Mpi::bcast(vector <double> &vec)
{
	MPI_Bcast(&vec[0],vec.size(),MPI_DOUBLE,0,MPI_COMM_WORLD);
}


/// Copies a matrix in core 0 to all the other cores
void Mpi::bcast(vector < vector <double> > &M)
{
	auto n = (unsigned int)(M.size());
	bcast(n);
	for(auto j = 0u; j < n; j++) bcast(M[j]);
}


/// Copies proposals from core 0 to all the other cores
void Mpi::bcast(vector <Proposal> &vec)
{	
	unsigned int si;
	if(core == 0) si = vec.size();
	bcast(si);

	vector<unsigned int> type, num;
	if(core == 0){ for(const auto &prop : vec){ type.push_back(prop.type); num.push_back(prop.num);}}
	else{ type.resize(si); num.resize(si);}
	
	bcast(type);
	bcast(num);
	
	if(core != 0){
		vec.clear();
		for(auto i = 0u; i < type.size(); i++){ 
			Proposal p; p.type = PropType(type[i]); p.num = num[i];
			vec.push_back(p);
		}
	}
}


/*
/// This checks that all cores get to a check point with a specified number at the same time
void Mpi::check(unsigned int num)
{
	if(core == 0) cout << num << "start" << endl << flush;
	barrier();
	vector <double> vec(1);
	vec[0] = num;
	auto vectot = gather(vec);
	if(core == 0){
		for(auto i = 0u; i < vectot.size(); i++) cout << vectot[i] << ",";
		cout << " check list" << endl << flush;
		
		for(auto i = 0u; i < vectot.size(); i++){
			if( vectot[i] !=  vectot[0]) emsg("Problem");
		}
	}
	barrier();
}
*/

/// Calculates the time taken for other cores to finish what they are doing
void Mpi::barrier()
{
	timer[TIME_WAIT].start();
	MPI_Barrier(MPI_COMM_WORLD);  
	timer[TIME_WAIT].stop();
}


/// Gathers an unsigned int vector across all cores and returns the combined vector to core 0
vector <unsigned int> Mpi::gather(const vector <unsigned int> &vec)
{
	vector <unsigned int> vectot;
	vectot.resize(vec.size()*ncore);
	
	MPI_Gather(&vec[0],vec.size(),MPI_UNSIGNED,&vectot[0],vec.size(),MPI_UNSIGNED,0,MPI_COMM_WORLD);
	
	return vectot;
}


/// Combines together vector (potentially of different length) to core 0
vector <double> Mpi::combine(const vector <double> &vec)
{
	vector <double> vecout;
	
	if(core == 0){
		vecout = vec;
		for(auto co = 1u; co < ncore; co++){
			pack_recv(co);		
			vector <double> vecrec;
			unpack(vecrec);
			for(auto i = 0u; i < vecrec.size(); i++) vecout.push_back(vecrec[i]);
			unpack_check();
		}
	}
	else{
		pack_initialise(0);
		pack(vec);
		pack_send(0);
	}
	
	return vecout;
}	
	

/// Gathers together results from particles from all the cores and makes parameter and state samples on core 0
void Mpi::gather_samples(vector <ParamSample> &psamp, vector <Sample> &opsamp, const vector <Particle> &part, State &state, const string dir)
{
	barrier();
	
	if(core == 0){
		for(const auto &pa : part){
			state.initialise_from_particle(pa);
			if(details.siminf == INFERENCE) state.save(dir+"sample"+to_string(psamp.size()));
			psamp.push_back(state.create_param_sample(pa.run));
			opsamp.push_back(state.create_sample());
		}

		for(auto co = 1u; co < ncore; co++){
			unsigned int N;
			
			pack_recv(co);	
			unpack(N);
			
			for(auto i = 0u; i < N; i++){
				Particle part;
				unpack(part);
				state.initialise_from_particle(part);
				
				if(details.siminf == INFERENCE) state.save(dir+"sample"+to_string(psamp.size()));
				psamp.push_back(state.create_param_sample(part.run));
				opsamp.push_back(state.create_sample());
			}
			
			unpack_check();
		}

		if(details.siminf == INFERENCE){
			auto file = dir+"sample_information.txt";
			ofstream fout(file);
			if(!fout) emsg("Cannot open the file '"+file+"'");
			fout << "# Samples = '" << psamp.size() << "'" << endl;
			fout << "TOML file = '" << details.toml_file << "'" << endl; 
			fout << "start = '" << details.start << "'" << endl;
			fout << "end = '" << details.end << "'" << endl;
			fout << "division_per_time = '" << details.division_per_time << "'" << endl;
			fout << "timesteps = '" << details.ndivision << "'" << endl;
		}
	}
	else{
		unsigned int N = part.size();
		pack_initialise(0);
		pack(N);
		for(auto i = 0u; i < N; i++) pack(part[i]);
		pack_send(0);
	}
}


/// Gathers EF samples on different chaings 
vector < vector < vector <double> > > Mpi::gather_EF_chain_sample(const vector <Chain> &chain, const vector <Particle> &part, const unsigned int N, const unsigned int nchain, const unsigned int nrun, vector <double> &invT_total)
{
	vector < vector < vector <double> > > EF_chain_sample;
	EF_chain_sample.resize(nrun);
	for(auto run = 0u; run < nrun; run++) EF_chain_sample[run].resize(nchain);
	invT_total.resize(nchain);

	if(core == 0){
		for(auto ch = 0u; ch < N; ch++){
			auto run = part[ch].run;
			auto num = chain[ch].num;
			if(run == 0) invT_total[num] = chain[ch].invT;
			EF_chain_sample[run][num] = chain[ch].EF_samp;
		}

		for(auto co = 1u; co < ncore; co++){
			unsigned int run, num;
			double invT;
			pack_recv(co);
			for(auto ch = 0u; ch < N; ch++){
				unpack(run);
				unpack(num);
				unpack(invT);
				unpack(EF_chain_sample[run][num]);
				if(run == 0) invT_total[num] = invT;
			}
			unpack_check();
		}
	}
	else{
		pack_initialise(0);
		for(auto ch = 0u; ch < N; ch++){
			pack(part[ch].run);
			pack(chain[ch].num);
			pack(chain[ch].invT);
			pack(chain[ch].EF_samp);
		}
		pack_send(0);
	}

	return EF_chain_sample;
}


/// Gathers together results parameter samples from particles all cores onto core 0
vector <ParamSample> Mpi::gather_psamp(const vector <Particle> &part)
{
	vector <ParamSample> psamp;
	
	if(core == 0){
		for(const auto &pa : part){
			ParamSample ps; ps.run = pa.run; ps.paramval = pa.paramval; ps.EF = pa.EF;
			psamp.push_back(ps);
		}
		
		for(auto co = 1u; co < ncore; co++){
			unsigned int N;
			pack_recv(co);
			unpack(N);
			for(auto i = 0u; i < N; i++){
				ParamSample ps; 
				unpack(ps.run);
				unpack(ps.paramval);
				unpack(ps.EF);
				psamp.push_back(ps);
			}
			unpack_check();
		}
	}
	else{
		unsigned int N = part.size();
		pack_initialise(0);
		pack(N);
		for(auto i = 0u; i < N; i++){ pack(part[i].run); pack(part[i].paramval); pack(part[i].EF);}
		pack_send(0);
	}
	
	return psamp;
}


/// Gathers together results parameter samples from particles all cores onto core 0
vector <ParamSample> Mpi::gather_psamp(const vector <ParamSample> &psample)
{
	vector <ParamSample> psamp;
	
	if(core == 0){
		psamp = psample;
		
		for(auto co = 1u; co < ncore; co++){
			vector <ParamSample> psa;
			pack_recv(co);
			unpack(psa);
			for(const auto &ps : psa) psamp.push_back(ps);
			unpack_check();
		}
	}
	else{
		pack_initialise(0);
		pack(psample);
		pack_send(0);
	}
	
	return psamp;
}


/// Gathers together particles from all cores onto core 0
vector <Particle> Mpi::gather_particle(const vector <Particle> &part)
{
	vector <Particle> particle;
	
	if(core == 0){
		for(const auto &pa : part) particle.push_back(pa);
		for(auto co = 1u; co < ncore; co++){
			unsigned int N;
			pack_recv(co);
			unpack(N);
			for(auto i = 0u; i < N; i++){
				Particle pa; unpack(pa);
				particle.push_back(pa);
			}
			unpack_check();
		}
	}
	else{
		unsigned int N = part.size();
		pack_initialise(0);
		pack(N);
		for(auto i = 0u; i < N; i++) pack(part[i]);
		pack_send(0);
	}
	
	return particle;
}


/// Initialises the buffer
void Mpi::pack_initialise(const size_t size)
{
	k = 0;
	buffer.resize(size);
}

/// Checks all the values have been read from the buffer
void Mpi::unpack_check()
{
	if(k != buffer.size()) emsgEC("Mpi",9);
}


/// Returns the size of the buffer
size_t Mpi::packsize()
{
	return k;
}


/// Makes a copy of the buffer
vector <double> Mpi::copybuffer()
{
	return buffer;
}


/// Sets the buffer to a given vector
void Mpi::setbuffer(const vector <double> &vec)
{
	k = 0;
	buffer = vec;
}


/// The pointer to the buffer
double* Mpi::packbuffer()
{
	return &buffer[0];
}


/// Sends the buffer to core co
void Mpi::pack_send(const unsigned int co)
{
	unsigned int si = packsize();
	MPI_Send(&si,1,MPI_UNSIGNED,co,0,MPI_COMM_WORLD);
	MPI_Send(packbuffer(),si,MPI_DOUBLE,co,0,MPI_COMM_WORLD);
}


/// Copies the buffer to all cores
void Mpi::pack_bcast()
{
	auto si = packsize();
	MPI_Bcast(&si,1,MPI_UNSIGNED,0,MPI_COMM_WORLD);

	if(core != 0) pack_initialise(si);

	MPI_Bcast(packbuffer(),si,MPI_DOUBLE,0,MPI_COMM_WORLD);
}


/// Recieves a message from core co and places it into the buffer
void Mpi::pack_recv(const unsigned int co)
{
	unsigned int si;
	MPI_Recv(&si,1,MPI_UNSIGNED,co,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
	pack_initialise(si);
	MPI_Recv(packbuffer(),si,MPI_DOUBLE,co,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
}	

template <class T>
void Mpi::pack_item(T t)
{
	buffer.push_back(t); k++;
}

// Provide overloads for cases where default behaviour needs to be modified
void Mpi::pack_item(const string &vec)
{
	auto jmax = vec.length();
	pack_item(jmax);
	for(auto j = 0u; j < jmax; j++){
		pack_item(vec.at(j));
	}
}

void Mpi::pack(const unsigned int num)
{
	pack_item(num); 
}

void Mpi::pack(const int num)
{
	pack_item(num); 
}

void Mpi::pack(const unsigned short num)
{
	pack_item(num); 
}

void Mpi::pack(const double num)
{
	pack_item(num); 
}

void Mpi::pack(const string &vec)
{
	pack_item(vec);
}

// Use template to share implementation details between cases
template <class T>
void Mpi::pack_item(const vector<T> &vec)
{
	pack_item(vec.size());
	for (auto &item : vec) {
		pack_item(item);
	}
}

void Mpi::pack(const vector <unsigned int> &vec)
{
	pack_item(vec);
}

void Mpi::pack(const vector <unsigned short> &vec)
{
	pack_item(vec);
}

void Mpi::pack(const vector <int> &vec)
{
	pack_item(vec);
}

void Mpi::pack(const vector <double> &vec)
{
	pack_item(vec);
}

void Mpi::pack(const vector< vector <unsigned int> > &vec)
{
	pack_item(vec);
}

void Mpi::pack(const vector< vector <double> > &vec)
{
	pack_item(vec);
}

void Mpi::pack(const vector< vector <float> > &vec)
{
	pack_item(vec);
}

void Mpi::pack(const vector< vector< vector <double> > > &vec)
{
	pack_item(vec);
}

void Mpi::pack(const vector <string> &vec)
{
	pack_item(vec);
}

void Mpi::pack(const Particle &part)
{
	pack(part.paramval);
	pack(part.EF);
	pack(part.run);
		
	const auto &vec = part.transnum;
	auto imax = vec.size(); buffer.push_back(imax); k++;
	for(auto i = 0u; i < imax; i++){
		auto jmax = vec[i].size(); buffer.push_back(jmax); k++;
		for(auto j = 0u; j < jmax; j++){
			auto nmax = vec[i][j].size(); buffer.push_back(nmax); k++;
			for(auto n = 0u; n < nmax; n++){
				auto pmax = vec[i][j][n].size(); buffer.push_back(pmax); k++;
				for(auto p = 0u; p < pmax; p++){
					buffer.push_back(vec[i][j][n][p]); k++;
				}
			}
		}
	}
}

void Mpi::pack(const Observation &ob)
{           
	pack(ob.factor_spline);
	pack(ob.datatable);
	pack(ob.graph);
	pack(ob.value);
	pack(ob.sd);
	pack(ob.factor);
	pack(ob.sett_i);
	pack(ob.sett_f);
	pack(ob.area);	
	pack(ob.dp_sel);	
	pack((unsigned int) ob.obsmodel);
	pack(ob.shape);
}
		
void Mpi::pack(const Modification &cf)
{
	pack(cf.start);
	pack(cf.end);
	pack((unsigned int) cf.type);
	pack(cf.spline_name_str);
	pack(cf.trans_str);
	pack(cf.strain_str);
	pack(cf.geo_filt);
	pack(cf.democats_filt);
	pack(cf.factor);
	pack(cf.dp_sel);
	pack(cf.area);
}

void Mpi::pack(const GenerateQ &genQ)
{
	for(auto i = 0u; i < genQ.N_name.size(); i++) pack(genQ.N[i]);
	pack(genQ.M);
	pack(genQ.onlydiag);
	pack(genQ.factor);
	pack((unsigned int)(genQ.treenode.size()));
	for(auto i = 0u; i < genQ.treenode.size(); i++){
		pack(genQ.treenode[i].arearef);
		pack(genQ.treenode[i].child);
	}    
}

void Mpi::pack_item(const Area &area)
{
	pack_item(area.code);
	pack_item(area.pop_init);
	pack_item(area.pop_dp);
	pack_item(area.total_pop);
}

void Mpi::pack(const vector <Area> &vec)
{
	pack_item(vec);
}

void Mpi::pack(const vector <Proposal> &vec)
{
	pack_item(vec.size());
	for(auto si = 0u; si < vec.size(); si++){
		pack_item(vec[si].type);
		pack_item(vec[si].num);
	}
}

void Mpi::pack(const vector <GeographicMap> &geomap)
{
	unsigned int imax = geomap.size();
	pack(imax);
	for(auto i = 0u; i < imax; i++){
		pack(geomap[i].region);
		pack(geomap[i].area);
	}
}

void Mpi::pack(const vector < vector <vector <unsigned int> > > &vec)
{
	pack_item(vec);
}

void Mpi::pack(const vector < vector <vector <int> > > &vec)
{
	pack_item(vec);
}

void Mpi::pack(const Matrix &mat)
{
	pack(mat.N);
	pack(mat.ele);
}

void Mpi::pack(const SparseMatrix &mat)
{
	pack(mat.N);
	pack(mat.diag);
	pack(mat.to);
	pack(mat.val);
}

void Mpi::pack(const vector<ParamSample> &vec)
{
	unsigned int imax = vec.size();
	pack(imax);
	for(auto i = 0u; i < imax; i++){
		pack(vec[i].paramval);
		pack(vec[i].run);
		pack(vec[i].EF);
	}
}

template<class T>
void Mpi::unpack_item(T &num)
{
	num = buffer[k]; k++;
}

void Mpi::unpack_item(string &vec)
{
	unsigned int jmax;
	
	unpack_item(jmax);
	stringstream ss; for(auto j = 0u; j < jmax; j++){ ss << (char) buffer[k]; k++;}
	vec = ss.str();
}

void Mpi::unpack(unsigned int &num)
{
	unpack_item(num);
}

void Mpi::unpack(unsigned short &num)
{
	unpack_item(num);
}

void Mpi::unpack(double &num)
{
	unpack_item(num);
}

void Mpi::unpack(string &vec)
{
	unpack_item(vec);
}

template <class T>
void Mpi::unpack_item(vector<T> &vec)
{
	unsigned int size;
	unpack_item(size);
	vec.resize(size);
	for (auto &item : vec) {
		unpack_item(item);
	}
}

void Mpi::unpack(vector <unsigned int> &vec)
{
	unpack_item(vec);
}

void Mpi::unpack(vector <unsigned short> &vec)
{
	unpack_item(vec);
}

void Mpi::unpack(vector <int> &vec)
{
	unpack_item(vec);
}

void Mpi::unpack(vector <double> &vec)
{
	unpack_item(vec);
}

void Mpi::unpack(vector< vector <unsigned int> > &vec)
{
	unpack_item(vec);
}

void Mpi::unpack(vector< vector <double> > &vec)
{
	unpack_item(vec);
}

void Mpi::unpack(vector< vector <float> > &vec)
{
	unpack_item(vec);
}

void Mpi::unpack(vector< vector< vector <double> > > &vec)
{
	unpack_item(vec);
}

void Mpi::unpack(vector <string> &vec)
{
	unpack_item(vec);
}

void Mpi::unpack(Particle &part)
{
	unpack(part.paramval);
	unpack(part.EF);
	unpack(part.run);
	
	auto &vec = part.transnum;
	unsigned int imax = buffer[k]; k++; vec.resize(imax);
	for(auto i = 0u; i < imax; i++){
		unsigned int jmax = buffer[k]; k++; vec[i].resize(jmax);
		for(auto j = 0u; j < jmax; j++){
			unsigned int nmax = buffer[k]; k++; vec[i][j].resize(nmax);
			for(auto n = 0u; n < nmax; n++){
				unsigned int pmax = buffer[k]; k++; vec[i][j][n].resize(pmax);
				for(auto p = 0u; p < pmax; p++){
					vec[i][j][n][p] = buffer[k]; k++;
				}
			}
		}
	}
}

void Mpi::unpack(Observation &ob)
{
	unpack(ob.factor_spline);
	unpack(ob.datatable);
	unpack(ob.graph);
	unpack(ob.value);
	unpack(ob.sd);
	unpack(ob.factor);
	unpack(ob.sett_i);
	unpack(ob.sett_f);
	unpack(ob.area);	
	unpack(ob.dp_sel);
	unsigned int num; unpack(num); ob.obsmodel = ObsModelFunc (num);
	unpack(ob.shape);
}

void Mpi::unpack(Modification &cf)
{
	unpack(cf.start);
	unpack(cf.end);
	unsigned int num; unpack(num); cf.type = ModificationType (num);
	unpack(cf.spline_name_str);
	unpack(cf.trans_str);
	unpack(cf.strain_str);
	unpack(cf.geo_filt);
	unpack(cf.democats_filt);
	unpack(cf.factor);
	unpack(cf.dp_sel);
	unpack(cf.area);
}

void Mpi::unpack(GenerateQ &genQ)
{
	genQ.N.resize(genQ.N_name.size());
	for(auto i = 0u; i < genQ.N_name.size(); i++) unpack(genQ.N[i]);
	unpack(genQ.M);
	unpack(genQ.onlydiag);
	unpack(genQ.factor);
	unsigned int num; unpack(num); genQ.treenode.resize(num);
	for(auto i = 0u; i < genQ.treenode.size(); i++){
		unpack(genQ.treenode[i].arearef);
		unpack(genQ.treenode[i].child);
	}    
}
		
void Mpi::unpack_item(Area &area)
{
	unpack_item(area.code);
	unpack_item(area.pop_init);
	unpack_item(area.pop_dp);
	unpack_item(area.total_pop);
}

void Mpi::unpack(vector <Area> &vec)
{
	unpack_item(vec);
}

void Mpi::unpack(vector <GeographicMap> &geomap)
{
	unsigned int imax = buffer[k]; k++; geomap.resize(imax);
	for(auto i = 0u; i < geomap.size(); i++){
		unpack(geomap[i].region);
		unpack(geomap[i].area);
	}
}

void Mpi::unpack(vector < vector <vector <unsigned int> > > &vec)
{
	unpack_item(vec);
}

void Mpi::unpack(vector < vector <vector <int> > > &vec)
{
	unpack_item(vec);
}

void Mpi::unpack(Matrix &mat)
{
	unpack(mat.N);
	unpack(mat.ele);
}

void Mpi::unpack(SparseMatrix &mat)
{
	unpack(mat.N);
	unpack(mat.diag);
	unpack(mat.to);
	unpack(mat.val);
}


void Mpi::unpack(vector<ParamSample> &vec)
{
	unsigned int imax; 
	unpack(imax); vec.resize(imax);
	for(auto i = 0u; i < imax; i++){
		unpack(vec[i].paramval);
		unpack(vec[i].run);
		unpack(vec[i].EF);
	}
}