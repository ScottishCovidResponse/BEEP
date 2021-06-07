#ifndef BEEPMBP__MPI_HH
#define BEEPMBP__MPI_HH

#include "struct.hh"

struct Mpi {
	Mpi(const Details &details);
	
public:
	unsigned int ncore;                                          // The number of cores that MPI is using
	unsigned int core;                                           // The core of the current process
	
	void copy_data(unsigned int &narea, vector <Area> &area, unsigned int &nobs, vector <Observation> &obs, GenerateQ &genQ, vector <Modification> &modification, vector <Covariate> &covar, LevelEffect &level_effect, vector <DemocatChange> &democat_change);
	void copy_particles(vector<Particle> &part, vector <unsigned int> &partcopy, const unsigned int N, const unsigned int Ntot);
	void gather_samples(vector <ParamSample> &psamp, vector <Sample> &opsamp, const vector <Particle> &part, State &state, const string dir);
	vector < vector < vector <double> > > gather_EF_chain_sample(const vector <Chain> &chain, const vector <Particle> &part, const unsigned int N, const unsigned int nchain, const unsigned int nrun, vector <double> &invT_total);
	vector <ParamSample> gather_psamp(const vector <Particle> &part);
	vector <ParamSample> gather_psamp(const vector <ParamSample> &psample);
	vector <Particle> gather_particle(const vector <Particle> &part);
	void exchange_samples(Generation &gen);
	void end_sec_swap(vector <State> &particle, const unsigned int sett, const vector <unsigned int> &backpart, const unsigned int buffersize);
	Particle particle_sample(const unsigned int ru, const vector < vector <unsigned int> > &backpart, vector <State> &particle, const ObservationModel &obsmodel);
		
	vector <double> gather(const vector <double> &vec);
	vector <unsigned int> gather(const vector <unsigned int> &vec);	
	vector <long> gather(const long val);
	vector <double> gather(const double val);
	vector <double> scatter(const vector <double> &vectot);
	vector <unsigned int> scatter(const vector <unsigned int> &vectot);
	
	void barrier();
	long sum(const long val);
	double sum(const double val);
	double average(const double val);
	vector <double> average(const vector <double> &val);
	double get_acrate(const unsigned int nac, const unsigned int ntr);
	vector <double> get_acrate(const vector <unsigned int> &nac, const vector <unsigned int> &ntr);
	double get_ratio(const double nac, const double ntr);

	void bcast(double &val);
	void bcast(bool &val);
	void bcast(unsigned int &val);
	void bcast(vector <unsigned int> &vec);
	void bcast(vector <double> &vec);
	void bcast(vector <Proposal> &vec);
	
	vector <double> combine(const vector <double> &vec);
	
private:
	vector<double> buffer;                                              // Stores packed up information to be sent between cores
	unsigned int k;                                                     // Indexes the buffer
	
	void pack_initialise(const size_t size);
	void unpack_check();
	size_t packsize();
	double *packbuffer();
	vector <double> copybuffer();
	void setbuffer(const vector <double> &vec);

	template <class T>
	void pack_item(T t);
	template <class T>
	void pack_item(const vector<T>& vec);	

	void pack_item(const string& vec);
	void pack_item(const Area& area);
		
	template<class T>
	void unpack_item(T &num);
	
	template <class T>
	void unpack_item(vector<T>& vec);

	void unpack_item(string &vec);
	void unpack_item(Area& area);

	void pack_send(const unsigned int co);
	void pack_recv(const unsigned int co);
	void pack_bcast();
	
	void pack(const int num);
	void pack(const unsigned int num);
	void pack(const unsigned short num);
	void pack(const double num);
	void pack(const string &vec);
	void pack(const vector <unsigned int> &vec);
	void pack(const vector <unsigned short> &vec);
	void pack(const vector <int> &vec);
	void pack(const vector <double> &vec);
	void pack(const vector< vector <unsigned int> > &vec);
	void pack(const vector< vector <double> > &vec);
	void pack(const vector< vector <float> > &vec);
	void pack(const vector< vector< vector <double> > > &vec);
	void pack(const vector <string> &vec);
	void pack(const Particle &part);
	void pack(const Observation &ob);
	void pack(const Modification &cf);
	void pack(const GenerateQ &genQ);
	void pack(const vector <Area> &vec);
	void pack(const vector <Proposal> &vec);
	void pack(const vector <DemographicCategory> &vec);
	void pack(const vector <GeographicMap>& geomap);
	void pack(const vector < vector <vector <unsigned int> > > &vec);
	void pack(const vector < vector <vector <int> > > &vec);
	void pack(const Matrix &mat);
	void pack(const SparseMatrix &mat);
	void pack(const vector <ParamSample> &vec);
	
	void unpack(unsigned int &num);
	void unpack(unsigned short &num);
	void unpack(double &num);
	void unpack(string &vec);
	void unpack(vector <unsigned int> &vec);
	void unpack(vector <unsigned short> &vec);
	void unpack(vector <int> &vec);
	void unpack(vector <double> &vec);
	void unpack(vector< vector <unsigned int> > &vec);
	void unpack(vector< vector <double> > &vec);
	void unpack(vector< vector <float> > &vec);
	void unpack(vector< vector< vector <double> > > &vec);
	void unpack(vector <string> &vec);
	void unpack(Particle &part);
	void unpack(Observation &ob);
	void unpack(Modification &cf);
	void unpack(GenerateQ &genQ);
	void unpack(vector <Area> &vec);
	void unpack(vector <Proposal> &vec);
	void unpack(vector <DemographicCategory> &vec);
	void unpack(vector <GeographicMap>& geomap);
	void unpack(vector < vector <vector <unsigned int> > > &vec);
	void unpack(vector < vector <vector <int> > > &vec);
	void unpack(Matrix &mat);
	void unpack(SparseMatrix &mat);
	void unpack(vector <ParamSample> &vec);
	
	const Details &details;
};
#endif
