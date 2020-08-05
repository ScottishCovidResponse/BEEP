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
	
	unsigned int gettime(string st) const;   // Gets a time from a date
	string getdate(unsigned int t) const;    // Gets a date fri=om a time
	
	Mode mode;                               // Stores if doing simulation/inference
	string outputdir;                        // The output directory
	
	unsigned int tform;                      // The time format (e.g. times or dates)
	string tformat;                          // A description of the time format ('time' or 'date').

	unsigned int start;                      // The start time over which simulation/inference is performed
	unsigned int end;                        // The start time over which simulation/inference is performed
	unsigned int period;                     // The time over which simulation/inference is performed (e.g. in weeks)
	
};
#endif