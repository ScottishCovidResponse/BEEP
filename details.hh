#ifndef BEEPMBP__DETAILS_HH
#define BEEPMBP__DETAILS_HH

#include "struct.hh"
#include "inputs.hh"

struct Details                             // Provides various details used to define the simulation / inference
{
	Details(const Inputs &inputs);
	
	unsigned int gettime(string st) const;   // Gets a time from a date
	string getdate(unsigned int t) const;    // Gets a date fri=om a time
	
	Mode mode;                               // Stores the mode of operation
	SimInf siminf;                           // Stores if doing simulation/inference
	string output_directory;                 // The output directory
	
	unsigned int efoi_factor;                // Denominator when expressing external force of infection.

	TimeFormat time_format;                  // The time format (e.g. times or dates)
	string time_format_str;                  // A description of the time format ('time' or 'date').

	unsigned int start;                      // The start time over which simulation/inference is performed
	unsigned int end;                        // The start time over which simulation/inference is performed
	unsigned int period;                     // The time over which simulation/inference is performed (e.g. in weeks)
	
	vector <TimePlot> timeplot;              // Times which can be used 
	
	unsigned int division_per_time;          // # divisions per unit time
	double timestep;                         // The timestep per division

	ObsCombineType obs_combine_type;         // Determines how observations are combined 
	unsigned int pmcmc_obs_period;           // The period between particle filtering (in divisions);
	vector <bool> sec_define;                // Defines sections when performing PMCMC
	
	unsigned int ndivision;                  // # Divisions into which the global timeline is divided for update of Q
	vector <double> division_time;           // The discretised timings at which infectivity map is updated
};
#endif