#ifndef BEEPMBP__SIMULATE_HH
#define BEEPMBP__SIMULATE_HH

#include "data.hh"
#include "model.hh"

class DATA;
class MODEL;
class POPTREE;

class Simulate
{
public:	
	Simulate(DATA &data, MODEL &model, POPTREE &poptree, Mpi &mpi, Inputs &inputs, Mode mode, bool verbose);	
	
	void run();
	void multirun();
	
private:
	void proportions(const vector< vector <FEV> > &indev);

	unsigned int nsamp;                                   // The number of simulations 
	
	DATA &data;
	MODEL &model;
	POPTREE &poptree;
	Mpi &mpi;
};

//void simulatedata(DATA &data, MODEL &model, POPTREE &poptree, Mcmc &mcmc);

#endif
