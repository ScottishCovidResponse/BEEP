#ifndef BEEP__IMPORTANCE_HH
#define BEEP__IMPORTANCE_HH

#include "struct.hh"
#include "mbp.hh"
#include "param_prop.hh"

class IMPORTANCE
{
public:	
  IMPORTANCE(const Details &details, const Data &data, const Model &model, Inputs &inputs, const Output &output, const ObservationModel &obsmodel, Mpi &mpi);	
	void run();
	
private:
	double invT;                             // The inverse temperature for the power posterior
	
	unsigned int N;                          // The number of particles per core
	unsigned int Ntot;                       // The total number of particles

	vector <ObsSlice> obs_slice;             // The observation at different time points
	
	State state;                             // Stores the state of the system
	
	const Details &details;
	const Data &data;
	const Model &model;
	const Output &output;
	const ObservationModel &obsmodel;
	Mpi &mpi;
};   

#endif
