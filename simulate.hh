#ifndef BEEPMBP__SIMULATE_HH
#define BEEPMBP__SIMULATE_HH

#include "data.hh"
#include "model.hh"

class DATA;
class MODEL;
class POPTREE;
class Output;
class Obsmodel;

class Simulate
{
public:	
	Simulate(Details &details, DATA &data, MODEL &model, POPTREE &poptree, Mpi &mpi, Inputs &inputs, Output &output, Obsmodel &obsmodel, Mode mode, bool verbose);	
	
	void run();
	void multirun();
	
private:
	void proportions(const vector< vector <FEV> > &indev);

	unsigned int nsamp;                                   // The number of simulations 
	
	Details &details;
	DATA &data;
	MODEL &model;
	POPTREE &poptree;
	Mpi &mpi;
	Output &output;
	Obsmodel &obsmodel;
};

#endif
