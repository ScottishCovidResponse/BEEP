#ifndef BEEPMBP__ABC_HH
#define BEEPMBP__ABC_HH

#include "data.hh"
#include "model.hh"
#include "chain.hh"

struct ParamSample;
class Data;
class Model;
class AreaTree;
class Output;
class ObservationModel;
class Chain;
struct Sample;

class ABC
{
public:	
  ABC(const Details &details, const Data &data, const Model &model, const AreaTree &areatree, const Mpi &mpi, const Inputs &inputs, const Output &output, const ObservationModel &obsmodel);	
	void smc();
	void mbp();
	
private:
	void mcmc_updates(Generation &gen, vector <Particle> &part, Chain &chain);
	double calculate_mixing(const vector <Particle> &part, vector <unsigned int> &partcopy) const;
	void exchange_samples_mpi(Generation &gen);
	double next_generation_mpi(vector<Particle> &part, vector <unsigned int> &partcopy);
	void results_mpi(const vector <Generation> &generation, const vector <Particle> &part, Chain &chain) const;
	Sample get_sample(const Particle &part, Chain &chain) const;
	double acceptance(double rate) const;
	void calculate_particle_weight(vector <Generation> &generation, double jumpsize);
	unsigned int effective_particle_number(vector <double> w);
	
	unsigned int total_time;                 // The time the algorithm is run for
	
	unsigned int N;                          // The number of particles per core
	unsigned int Ntot;                       // The total number of particles
	
	unsigned int G;                          // The total number of generations 
	
	Chain chain;
	Jump &jump;
	
	const Details &details;
	const Data &data;
	const Model &model;
	const AreaTree &areatree;
	const Mpi &mpi;
	const Output &output;
	const ObservationModel &obsmodel;
};

#endif
