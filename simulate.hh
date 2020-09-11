#ifndef BEEPMBP__SIMULATE_HH
#define BEEPMBP__SIMULATE_HH

#include "data.hh"
#include "model.hh"

struct ParamSample;
class Data;
class Model;
class AreaTree;
class Output;
class ObservationModel;

class Simulate
{
public:	
	Simulate(const Details &details, const Data &data, const Model &model, const AreaTree &areatree, const Mpi &mpi, const Inputs &inputs, Output &output, const ObservationModel &obsmodel);	
	void run();
	void multirun();
	
private:
	void proportions(const vector< vector <Event> > &indev);
	unsigned int nsamp;                                   // The number of simulations 
	
	const Details &details;
	const Data &data;
	const Model &model;
	const AreaTree &areatree;
	const Mpi &mpi;
	const ObservationModel &obsmodel;
	Output &output;
};

#endif
