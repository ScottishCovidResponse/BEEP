/// Information and routines for transferring data using MPI

#include <sstream>

using namespace std;

#include "mpi.hh"
#include "param_prop.hh"
#include "state.hh"
#include "areatree.hh"

Mpi::Mpi()
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
void Mpi::copy_data(unsigned int &narea, vector <Area> &area, vector <GeographicMap> &region_effect, unsigned int &nobs, vector <Observation> &obs, GenerateQ &genQ, vector <CounterFact> &counterfact)
{
	if(core == 0){                              				   // Copies the above information to all the other cores
	
		pack_initialise(0);
		pack(narea);
		pack(area);
		pack(region_effect);	
		pack(nobs);	for(const auto& ob : obs)	pack(ob);
		pack(genQ);
		pack((unsigned int) counterfact.size()); 
		for(const auto& cf : counterfact) pack(cf);
	}

	pack_bcast();

	if(core != 0){
		unpack(narea);
		unpack(area);
		unpack(region_effect);
		unpack(nobs); obs.resize(nobs); for(auto& ob : obs) unpack(ob);
		unpack(genQ);
		unsigned int ncounterfact; unpack(ncounterfact); counterfact.resize(ncounterfact); 
		for(auto& cf : counterfact)	unpack(cf);
	}
}


/// Copies particles across MPI processes (ABC-MBP, ABC-MBP-GR, PAIS)
void Mpi::copy_particles(vector<Particle> &part, vector <unsigned int> &partcopy, unsigned int N, unsigned int Ntot)
{
	bcast(partcopy);
	
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
			if(partcopy[iitot] == itot && iitot/N != itot/N){
				pack_initialise(0);
				pack(part[i]);
	
				sendbuffer.push_back(copybuffer()); 
				sendbuffer_to.push_back(iitot);
			
				buffersize[i] = packsize();
				if(false) cout << "send " << itot << " -> " << iitot << " \n";
			}				
		}
	}
	
	for(auto i = 0u; i < N; i++){
		auto itot = core*N+i;
		auto iitot = partcopy[itot];
		if(iitot != UNSET){
			if(iitot/N == itot/N){
				if(i != iitot%N) part[i] = part[iitot%N]; 
				if(false) cout << "internal copy " <<iitot << " -> " << itot << "\n";
			}
			else{
				recibuffer.push_back(vector <double> ());
				recibuffer_from.push_back(iitot);
				recibuffer_to.push_back(itot);
				if(false) cout << "receive " << iitot << " -> " << itot << " \n"; 
			}
		}
	}
	
	auto buffersizetot = gather(buffersize);
	bcast(buffersizetot);
	
	auto buftot = recibuffer.size()+sendbuffer.size();
	if(buftot > 0){
		MPI_Request reqs[buftot];                 // These store information used in Isend and Irecv
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


/// Exchanges samples within a generation (used in ABC-SMC, ABC-MBP and PAIS)
void Mpi::exchange_samples(Generation &gen, unsigned int Ntot)
{
	if(core == 0){
		for(auto co = 1u; co < ncore; co++){
			pack_recv(co);		
			vector< vector <double> > paramsa;
			unpack(paramsa);
			for(auto i = 0u; i < paramsa.size(); i++) gen.param_samp.push_back(paramsa[i]);
			
			vector <double> EFsa;
			unpack(EFsa);
			for(auto i = 0u; i < EFsa.size(); i++) gen.EF_samp.push_back(EFsa[i]);
			
			vector < vector <double> > EFDT;
			unpack(EFDT);
			for(auto i = 0u; i < EFDT.size(); i++) gen.EF_datatable.push_back(EFDT[i]);
		}
		if(gen.EF_samp.size() != Ntot || gen.EF_datatable.size() != Ntot || gen.param_samp.size() != Ntot) emsgEC("Mpi",2);
	}
	else{
		pack_initialise(0);
		pack(gen.param_samp);
		pack(gen.EF_samp);
		pack(gen.EF_datatable);
		pack_send(0);
	}
	
	if(ncore > 1){
		if(core == 0){
			pack_initialise(0);
			pack(gen.param_samp);
			pack(gen.EF_samp);
		}
	
		pack_bcast();
	
		if(core != 0){ 
			unpack(gen.param_samp); 
			unpack(gen.EF_samp);
		}
	}
}


/// Gather all the particle parameter vectors to core 0
vector < vector <double> > Mpi::gather_paramval(const vector <Particle> &part)
{
	vector < vector <double> > param_tot;
	if(core == 0){
		for(const auto &pa : part) param_tot.push_back(pa.paramval);
		
		for(auto co = 1u; co < ncore; co++){
			pack_recv(co);		
			for(auto p = 0u; p < part.size(); p++){
				vector <double> vec; unpack(vec); param_tot.push_back(vec);
			}
		}
	}
	else{
		pack_initialise(0);
		for(const auto &pa : part) pack(pa.paramval);
		pack_send(0);
	}
	
	return param_tot;
}

/// At the end of an observed section, this swaps information between different particles (PMCMC)
void Mpi::end_sec_swap(vector <State> &particle, const unsigned int sett, const vector <unsigned int> &backpart, unsigned int buffersize)
{
	timer[TIME_PMCMCSWAP].start();
		
	auto N = particle.size();
	auto Ntot = backpart.size();
	
	MPI_Request reqs[Ntot];                                      // These store information used in Isend and Irecv
	MPI_Status stats[Ntot];
	
	vector < vector<double> > sendbuffer;                        // Buffers used for Isend and Irecv
	vector < vector<double> > recibuffer;
	sendbuffer.resize(N); recibuffer.resize(N);

	if(false && core == 0){
		for(auto p = 0u; p < Ntot; p++) cout << p << " " << backpart[p] << "  backpart\n";
	}
	
	auto npreclist = 0u;
	vector <unsigned int> preclist;
	vector < vector <unsigned int> > preclistli;
	
	if(false){ for(auto p = 0u; p < N; p++){ particle[p].pop[sett][0][0][0] = core*N+p;}} // For testing
	
	auto nreqs = 0u;
	auto pmin = core*N;                                          // Sets up information to be recieved
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
					nreqs++; if(nreqs > Ntot) emsgEC("Mpi",3); 
					npreclist++; if(npreclist > N) emsgEC("Mpi",4); 
				}
			}
		}
	}

	auto nsendbuf = 0u;
	for(auto p = pmin; p < pmin+N; p++){                          // Initiates information to be sent
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
				sendbuffer[nsendbuf] = copybuffer(); if(sendbuffer[nsendbuf].size() != buffersize) emsgEC("Mpi",5);
				
				for(auto j = 0u; j < ncorlist; j++){
					MPI_Isend(&sendbuffer[nsendbuf][0],buffersize,MPI_DOUBLE,corlist[j],p,MPI_COMM_WORLD,&reqs[nreqs]);
					nreqs++; if(nreqs > Ntot) emsgEC("Mpi",6); 
				}				
				nsendbuf++; if(nsendbuf > N) emsgEC("Mpi",7);
			}
		}
	}
	
	if(nreqs > 0){
		if(MPI_Waitall(nreqs,reqs,stats) != MPI_SUCCESS) emsgEC("Mpi",8);
	}
		
	for(auto rec = 0u; rec < npreclist; rec++){                           	// Unpacks the recieved information
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
			if(particle[p].pop[sett][0][0][0] != int(backpart[core*N+p])) emsgEC("Mpi",11);
			cout << core << ": " << " " << core*N+p << " " << particle[p].pop[sett][0][0][0] << " " << backpart[core*N+p] << " compare\n";
		}
		emsg("end");
	}
	
	timer[TIME_PMCMCSWAP].stop();
}


/// Constructs a particle by gathering the pieces from different particles (from the bootstrap function)
Particle Mpi::particle_sample(const vector < vector <unsigned int> > &backpart, vector <State> &particle, const ObservationModel &obsmodel)	
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
	
	return particle[0].create_particle();

	timer[TIME_STATESAMPLE].stop();
}


/// Gathers a double vector across all cores and returns the combined vector to core 0
vector <double> Mpi::gather(const vector <double> &vec)
{
	vector <double> vectot;
	vectot.resize(vec.size()*ncore);
	
	MPI_Gather(&vec[0],vec.size(),MPI_DOUBLE,&vectot[0],vec.size(),MPI_DOUBLE,0,MPI_COMM_WORLD);
	
	return vectot;
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
double Mpi::average(double val) 
{
	double val_av = sum(val)/ncore;
	bcast(val_av);

	return val_av;	
}


/// Calculates the average of a vector across cores
vector <double> Mpi::average(vector <double> val) 
{
	vector <double> result;
	for(auto i = 0u; i < val.size(); i++){
		result.push_back(average(val[i])); 
	}
	return result;	
}


/// Gets the acceptance rate across all mpi processes
double Mpi::get_acrate(unsigned int nac, unsigned int ntr)
{
	return average(nac)/(average(ntr)+0.01);
}


/// Gets a vector of acceptance rates across all mpi processes
vector <double> Mpi::get_acrate(vector <unsigned int> nac, vector <unsigned int> ntr)
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
		for(auto co = 1u; co < (unsigned int) ncore; co++){
			pack_recv(co);		
			vector <double> vecrec;
			unpack(vecrec);
			for(auto i = 0u; i < vecrec.size(); i++) vecout.push_back(vecrec[i]);
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
	if(core == 0){
		for(const auto& pa : part){
			state.initialise_from_particle(pa);
			state.save(dir+"sample"+to_string(psamp.size())+".txt");
			psamp.push_back(state.create_param_sample());
			opsamp.push_back(state.create_sample());
		}
		
		for(auto co = 1u; co < (unsigned int)ncore; co++){
			unsigned int N;
			pack_recv(co);
			unpack(N);
			for(auto i = 0u; i < N; i++){
				Particle part;
				unpack(part);
				state.initialise_from_particle(part);
				state.save(dir+"sample"+to_string(psamp.size())+".txt");
				psamp.push_back(state.create_param_sample());
				opsamp.push_back(state.create_sample());
			}
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


/// Initialises the buffer
void Mpi::pack_initialise(size_t size)
{
	k = 0;
	buffer.resize(size);
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
void Mpi::pack_send(unsigned int co)
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
void Mpi::pack_recv(unsigned int co)
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
void Mpi::pack_item(const string& vec)
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
void Mpi::pack_item(const vector<T>& vec)
{
	pack_item(vec.size());
	for (auto& item : vec) {
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
	//pack(part.ev);
	pack(part.EF);
	
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
	pack(ob.fraction_spline_ref);
	pack(ob.datatable);
	pack(ob.graph);
	pack(ob.value);
	pack(ob.sett_i);
	pack(ob.sett_f);
	pack(ob.area);	
	pack(ob.dp_sel);	
	pack(ob.logdif_offset);
	pack(ob.w);
	pack(ob.invT);
	pack((unsigned int) ob.obsmodel);
	pack(ob.shape);
}
		
void Mpi::pack(const CounterFact &cf)
{
	pack(cf.start);
	pack(cf.end);
	pack((unsigned int) cf.type);
	pack(cf.trans_str);
	pack(cf.geo);
	pack(cf.democats);
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
}

void Mpi::pack_item(const Area& area)
{
	pack_item(area.code);
	pack_item(area.region_effect);
	pack_item(area.x);
	pack_item(area.y);
	pack_item(area.pop);
	pack_item(area.covar);
	pack_item(area.ind);
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

void Mpi::pack(const vector <GeographicMap>& geomap)
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
void Mpi::unpack_item(vector<T>& vec)
{
	unsigned int size;
	unpack_item(size);
	vec.resize(size);
	for (auto& item : vec) {
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
	//unpack(part.ev);
	unpack(part.EF);

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
	unpack(ob.fraction_spline_ref);
	unpack(ob.datatable);
	unpack(ob.graph);
	unpack(ob.value);
	unpack(ob.sett_i);
	unpack(ob.sett_f);
	unpack(ob.area);	
	unpack(ob.dp_sel);
	
	unpack(ob.logdif_offset);
	unpack(ob.w);
	unpack(ob.invT);
	unsigned int num; unpack(num); ob.obsmodel = ObsModelFunc (num);
	unpack(ob.shape);
}

void Mpi::unpack(CounterFact &cf)
{
	unpack(cf.start);
	unpack(cf.end);
	unsigned int num; unpack(num); cf.type = CounterFactType (num);
	unpack(cf.trans_str);
	unpack(cf.geo);
	unpack(cf.democats);
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
}
		
void Mpi::unpack_item(Area& area)
{
	unpack_item(area.code);
	unpack_item(area.region_effect);
	unpack_item(area.x);
	unpack_item(area.y);
	unpack_item(area.pop);
	unpack_item(area.covar);
	unpack_item(area.ind);
}

void Mpi::unpack(vector <Area> &vec)
{
	unpack_item(vec);
}

void Mpi::unpack(vector <GeographicMap>& geomap)
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
