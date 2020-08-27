#ifndef BEEPMBP__SIMULATE_HH
#define BEEPMBP__SIMULATE_HH

#include "data.hh"
#include "model.hh"

struct PARAMSAMP;
class DATA;
class MODEL;
class POPTREE;
class Output;
class Obsmodel;

class Simulate
{
public:	
	Simulate(const Details &details, DATA &data, MODEL &model, const POPTREE &poptree, const Mpi &mpi, Inputs &inputs, Output &output, Obsmodel &obsmodel);	
	void run();
	void multirun();
	
private:
	void proportions(const vector< vector <FEV> > &indev);
	unsigned int nsamp;                                   // The number of simulations 
	
	const Details &details;
	DATA &data;
	MODEL &model;
	const POPTREE &poptree;
	const Mpi &mpi;
	Output &output;
	Obsmodel &obsmodel;
};

#endif
