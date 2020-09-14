// These are function for packing up quantities to be sent by MPI

#include <vector>
#include <iostream>
#include <algorithm>
#include <string>
#include <sstream>  

using namespace std;

#include "consts.hh"
#include "utils.hh"
#include "data.hh"
#include "model.hh"
#include "jump.hh"
#include "mcmc.hh"
#include "timers.hh"

namespace {
	unsigned int k;

	std::vector<double> buffer;
}


/// Initialises the buffer
void pack_initialise(size_t size)
{
	k = 0;
	buffer.resize(size);
}


/// Returns the size of the buffer
size_t packsize()
{
	return k;
}


/// Makes a copy of the buffer
vector <double> copybuffer()
{
	return buffer;
}


/// Sets the buffer to a given vector
void setbuffer(const vector <double> &vec)
{
	k = 0;
	buffer = vec;
}


/// The pointer to the buffer
double *packbuffer()
{
	return &buffer[0];
}


/// Sends the buffer to core co
void pack_mpi_send(unsigned int co)
{
	unsigned int si = packsize();
	MPI_Send(&si,1,MPI_UNSIGNED,co,0,MPI_COMM_WORLD);
	MPI_Send(packbuffer(),si,MPI_DOUBLE,co,0,MPI_COMM_WORLD);
}


/// Copies the buffer to all cores
void pack_mpi_bcast()
{
	auto si = packsize();
	MPI_Bcast(&si,1,MPI_UNSIGNED,0,MPI_COMM_WORLD);

	int core;
  MPI_Comm_rank(MPI_COMM_WORLD,&core); 
	if(core != 0) pack_initialise(si);

	MPI_Bcast(packbuffer(),si,MPI_DOUBLE,0,MPI_COMM_WORLD);
}


/// Recieves a message from core co and places it into the buffer
void pack_mpi_recv(unsigned int co)
{
	unsigned int si;
	MPI_Recv(&si,1,MPI_UNSIGNED,co,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
	pack_initialise(si);
	MPI_Recv(packbuffer(),si,MPI_DOUBLE,co,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
}	


/// Gathers a double vector across all cores and returns the combined vector to core 0
vector <double> mpi_gather(const vector <double> &vec)
{
	int num;
	MPI_Comm_size(MPI_COMM_WORLD,&num);
	vector <double> vectot;
	vectot.resize(vec.size()*num);
	
	MPI_Gather(&vec[0],vec.size(),MPI_DOUBLE,&vectot[0],vec.size(),MPI_DOUBLE,0,MPI_COMM_WORLD);
	
	return vectot;
}


/// Gathers a long across all cores and returns the combined vector to core 0
vector <long> mpi_gather(const long val)
{
	int num;
	MPI_Comm_size(MPI_COMM_WORLD,&num);
	vector <long> valtot;
	valtot.resize(num);
	
	MPI_Gather(&val,1,MPI_LONG,&valtot[0],1,MPI_LONG,0,MPI_COMM_WORLD);
	
	return valtot;
}


/// Sums up a value over all cores  
long mpi_sum(const long val)
{
	auto vec = mpi_gather(val);
	long sum = 0; for(auto val : vec) sum += val;
	return sum;	
}


/// Copies a variable in core 0 to all the other cores
void mpi_bcast(double &val)
{
	MPI_Bcast(&val,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
}


/// Copies a variable in core 0 to all the other cores
void mpi_bcast(unsigned int &val)
{
	MPI_Bcast(&val,1,MPI_UNSIGNED,0,MPI_COMM_WORLD);
}


/// Copies a vector in core 0 to all the other cores
void mpi_bcast(vector <unsigned int> &vec)
{
	MPI_Bcast(&vec[0],vec.size(),MPI_UNSIGNED,0,MPI_COMM_WORLD);
}
	
	
/// Calculates the time taken for other cores to finish what they are doing
void mpi_barrier()
{
	timers.timewait -= clock();
	MPI_Barrier(MPI_COMM_WORLD);  
	timers.timewait += clock();
}


/// Gathers an unsigned int vector across all cores and returns the combined vector to core 0
vector <unsigned int> mpi_gather(const vector <unsigned int> &vec)
{
	int num;
	MPI_Comm_size(MPI_COMM_WORLD,&num);
	vector <unsigned int> vectot;
	vectot.resize(vec.size()*num);
	
	MPI_Gather(&vec[0],vec.size(),MPI_UNSIGNED,&vectot[0],vec.size(),MPI_UNSIGNED,0,MPI_COMM_WORLD);
	
	return vectot;
}

template <class T>
void pack_item(T t)
{
	buffer.push_back(t); k++;
}

// Provide overloads for cases where default behaviour needs
// to be modified
void pack_item(const string& vec)
{
	auto jmax = vec.length();
	pack_item(jmax);
	for(auto j = 0u; j < jmax; j++){
		pack_item(vec.at(j));
	}
}

void pack(const unsigned int num)
{
	pack_item(num); 
}

void pack(const int num)
{
	pack_item(num); 
}

void pack(const unsigned short num)
{
	pack_item(num); 
}

void pack(const double num)
{
	pack_item(num); 
}

void pack(const string &vec)
{
	pack_item(vec);
}

// Use template to share implementation details between cases
template <class T>
void pack_item(const vector<T>& vec)
{
	pack_item(vec.size());
	for (auto& item : vec) {
		pack_item(item);
	}
}

void pack(const vector <unsigned int> &vec)
{
	pack_item(vec);
}

void pack(const vector <unsigned short> &vec)
{
	pack_item(vec);
}

void pack(const vector <int> &vec)
{
	pack_item(vec);
}

void pack(const vector <double> &vec)
{
	pack_item(vec);
}

void pack(const vector< vector <unsigned int> > &vec)
{
	pack_item(vec);
}

void pack(const vector< vector <double> > &vec)
{
	pack_item(vec);
}

void pack(const vector< vector <float> > &vec)
{
	pack_item(vec);
}

void pack(const vector< vector< vector <double> > > &vec)
{
	pack_item(vec);
}

void pack(const vector <string> &vec)
{
	pack_item(vec);
}

void pack(const vector <Event> &vec)
{
	auto imax = vec.size(); buffer.push_back(imax); k++;
	for(auto i = 0u; i < imax; i++){
		buffer.push_back(vec[i].trans); k++;
		buffer.push_back(vec[i].ind); k++;
		buffer.push_back(vec[i].t); k++;
		buffer.push_back(vec[i].timep); k++;
	}
}

void pack(const Particle &part)
{
	pack(part.paramval);
	pack(part.ev);
	pack(part.EF);
}

void pack_item(const Area& area)
{
	pack_item(area.code);
	pack_item(area.region);
	pack_item(area.x);
	pack_item(area.y);
	pack_item(area.agepop);
	pack_item(area.pop);
	pack_item(area.covar);
	pack_item(area.ind);
}

void pack(const vector <Area> &vec)
{
	pack_item(vec);
}

void pack_item(const DataRegion& r)
{
	pack_item(r.name);
	pack_item(r.code);
}

void pack(const vector <DataRegion> &vec)
{
	pack_item(vec);
}

void pack_item(const DemographicCategory& d)
{
	pack_item(d.name);
	pack_item(d.value);
}

void pack(const vector <DemographicCategory> &vec)
{
	pack_item(vec);
}

void pack_item(const EventRef& ev)
{
	pack_item(ev.ind);
	pack_item(ev.e);
}

void pack(const Jump& jump)
{
	pack_item(jump.mbp);
	pack_item(jump.mbp_ntr);
	pack_item(jump.mbp_nac);
	pack_item(jump.stand);		
	pack_item(jump.stand_ntr);		
	pack_item(jump.stand_nac);
	pack_item(jump.naddrem);	
	pack_item(jump.standev_ntr);	
	pack_item(jump.standev_nac);	
	pack_item(jump.nvar);	
	pack_item(jump.param_not_fixed);	
}

void pack(const ChainInfo& cinfo)
{
	pack(cinfo.invT);
	pack(cinfo.ch);
	pack(cinfo.jump);
}

void pack(const vector <vector <EventRef> > &vec)
{
	pack_item(vec);
}

void pack(const vector < vector <vector <unsigned int> > > &vec)
{
	pack_item(vec);
}

template<class T>
void unpack_item(T &num)
{
	num = buffer[k]; k++;
}
void unpack_item(string &vec)
{
	unsigned int jmax;
	
	unpack_item(jmax);
	stringstream ss; for(auto j = 0u; j < jmax; j++){ ss << (char) buffer[k]; k++;}
	vec = ss.str();
}

void unpack(unsigned int &num)
{
	unpack_item(num);
}

void unpack(unsigned short &num)
{
	unpack_item(num);
}

void unpack(double &num)
{
	unpack_item(num);
}

void unpack(string &vec)
{
	unpack_item(vec);
}

template <class T>
void unpack_item(vector<T>& vec)
{
	unsigned int size;
	unpack_item(size);
	vec.resize(size);
	for (auto& item : vec) {
		unpack_item(item);
	}
}

void unpack(vector <unsigned int> &vec)
{
	unpack_item(vec);
}

void unpack(vector <unsigned short> &vec)
{
	unpack_item(vec);
}

void unpack(vector <int> &vec)
{
	unpack_item(vec);
}

void unpack(vector <double> &vec)
{
	unpack_item(vec);
}

void unpack(vector< vector <unsigned int> > &vec)
{
	unpack_item(vec);
}

void unpack(vector< vector <double> > &vec)
{
	unpack_item(vec);
}

void unpack(vector< vector <float> > &vec)
{
	unpack_item(vec);
}

void unpack(vector< vector< vector <double> > > &vec)
{
	unpack_item(vec);
}

void unpack(vector <string> &vec)
{
	unpack_item(vec);
}

void unpack(vector <Event> &vec)
{
	unsigned int imax = buffer[k]; k++; vec.resize(imax);
	for(auto i = 0u; i < imax; i++){
		vec[i].trans = buffer[k]; k++;
		vec[i].ind = buffer[k]; k++;
		vec[i].t = buffer[k]; k++;
		vec[i].timep = buffer[k]; k++;
	}
}

void unpack(Particle &part)
{
	unpack(part.paramval);
	unpack(part.ev);
	unpack(part.EF);
}

void unpack_item(Area& area)
{
	unpack_item(area.code);
	unpack_item(area.region);
	unpack_item(area.x);
	unpack_item(area.y);
	unpack_item(area.agepop);
	unpack_item(area.pop);
	unpack_item(area.covar);
	unpack_item(area.ind);
}

void unpack(vector <Area> &vec)
{
	unpack_item(vec);
}

void unpack_item(DataRegion& r)
{
	unpack_item(r.name);
	unpack_item(r.code);
}

void unpack(vector <DataRegion> &vec)
{
	unpack_item(vec);
}

void unpack_item(DemographicCategory& d)
{
	unpack_item(d.name);
	unpack_item(d.value);
}

void unpack(vector <DemographicCategory> &vec)
{
	unpack_item(vec);
}

void unpack_item(EventRef& ev)
{
	unpack_item(ev.ind);
	unpack_item(ev.e);
}

void unpack(vector <vector <EventRef> > &vec)
{
	unpack_item(vec);
}

void unpack(Jump &jump)
{
	unpack_item(jump.mbp);
	unpack_item(jump.mbp_ntr);
	unpack_item(jump.mbp_nac);
	unpack_item(jump.stand);		
	unpack_item(jump.stand_ntr);		
	unpack_item(jump.stand_nac);
	unpack_item(jump.naddrem);	
	unpack_item(jump.standev_ntr);	
	unpack_item(jump.standev_nac);	
	unpack_item(jump.nvar);
	unpack_item(jump.param_not_fixed);	
}

void unpack(ChainInfo& cinfo)
{
	unpack(cinfo.invT);
	unpack(cinfo.ch);
	unpack(cinfo.jump);
}

void unpack(vector < vector <vector <unsigned int> > > &vec)
{
	unpack_item(vec);
}
