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

struct Details                             // Provides various details used to define the simulation / inference
{
	Details(const Inputs &inputs);
	
	unsigned int gettime(string st) const;   // Gets a time from a date
	string getdate(unsigned int t) const;    // Gets a date fri=om a time
	
	Mode mode;                               // Stores if doing simulation/inference
	string outputdir;                        // The output directory
	
	unsigned int tform;                      // The time format (e.g. times or dates)
	string tformat;                          // A description of the time format ('time' or 'date').

	unsigned int start;                      // The start time over which simulation/inference is performed
	unsigned int end;                        // The start time over which simulation/inference is performed
	unsigned int period;                     // The time over which simulation/inference is performed (e.g. in weeks)
	
	unsigned int fediv;                      // # Divisions into which the global timeline is divided for events
	unsigned int fepertime;                  // # fediv per nsettime
	unsigned int settpertime;                // # nsettime per unit time

	unsigned int nsettime;                   // # Divisions into which the global timeline is divided for update of Q
	vector <double> settime;                 // The timings at which beta changes
};
#endif