#ifndef BEEPMBP__DETAILS_HH
#define BEEPMBP__DETAILS_HH

#include "inputs.hh"
#include "utils.hh"
#include "consts.hh"

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
	string output_directory;                 // The output directory
	
	TimeFormat time_format;                  // The time format (e.g. times or dates)
	string time_format_str;                  // A description of the time format ('time' or 'date').

	unsigned int start;                      // The start time over which simulation/inference is performed
	unsigned int end;                        // The start time over which simulation/inference is performed
	unsigned int period;                     // The time over which simulation/inference is performed (e.g. in weeks)
	
	unsigned int division_per_time;          // # divisions per unit time

	unsigned int ndivision;                  // # Divisions into which the global timeline is divided for update of Q
	vector <double> division_time;           // The discretised timings at which infectivity map is updated
};
#endif