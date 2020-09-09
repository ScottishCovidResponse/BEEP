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
	Simulate(const Details &details, const DATA &data, const MODEL &model, const POPTREE &poptree, const Mpi &mpi, const Inputs &inputs, Output &output, const Obsmodel &obsmodel);	
	void run();
	void multirun();
	
private:
	void proportions(const vector< vector <FEV> > &indev);
	unsigned int nsamp;                                   // The number of simulations 
	
	const Details &details;
	const DATA &data;
	const MODEL &model;
	const POPTREE &poptree;
	const Mpi &mpi;
	const Obsmodel &obsmodel;
	Output &output;
};

#endif
