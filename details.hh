#ifndef BEEPMBP__DETAILS_HH
#define BEEPMBP__DETAILS_HH

#include "inputs.hh"
#include "utils.hh"

struct Mpi
{
	Mpi();
	
	unsigned int ncore;                      // The number of cores that MPI is using
	unsigned int core;                       // The core of the current process
};

struct Details
{
	Details(Inputs &inputs);
	
	Mode mode;                               // Stores if doing simulation/inference

	
};
#endif