#ifndef BEEPMBP__DETAILS_HH
#define BEEPMBP__DETAILS_HH

#include "struct.hh"
#include "inputs.hh"

struct Details                                                     // Provides details used in simulation / inference
{
	Details(Inputs &inputs);
	
	unsigned int gettime(const string st, const string em) const;    // Gets a time from a date
	string getdate(const unsigned int t) const;                      // Gets a date fri=om a time
	
	Mode mode;                                                       // Stores the mode of operation
	SimInf siminf;                                                   // Stores if doing simulation/inference
	string output_directory;                                         // The output directory
	
	string analysis_type;                                            // Labels "Posterior", "Simulation", "Multisim"
	
	string description;                                              // Stores a description of the analysis
	
	string toml_file;                                                // The name of the TOML file
	
	unsigned int efoi_factor;                                        // Denominator when expressing external force of infection.

	TimeFormat time_format;                                          // The time format (e.g. times or dates)
	string time_format_str;                                          // A description of the time format ('time' or 'date').

	long time_start;                                                 // The execution start time

	unsigned int start;                                              // Start time over which simulation/inference is performed
	unsigned int end;                                                // End time over which simulation/inference is performed
	unsigned int period;                                             // Period over which simulation/inference is performed
	unsigned int pred_start;                                         // The start time of prediction
	unsigned int pred_end;                                           // The end time of prediction

	unsigned int prop_max;                                           // The maximum number of proposals (ABCMBP/PAIS)

	vector <TimePlot> timeplot;                                      // Descriptive times which can be used 
	
	unsigned int division_per_time;                                  // # divisions per unit time
	double timestep;                                                 // The timestep per division

	unsigned int graph_step;                                         // The number of steps used when plotting 
		
	bool stochastic;                                                 // Determines if simulations are stochastic or not
	
	MCMCUpdate mcmc_update;                                          // Stores information about the mcmc updates
	
	bool obs_section;                                                // Set to true if observation are in sections (PMCMC)
	double trans_combine;                                            // Sets the value when combining transition observations

	unsigned int pmcmc_obs_period;                                   // The period between particle filtering (in divisions);
	vector <bool> sec_define;                                        // Defines sections when performing PMCMC
	 
	unsigned int ndivision;                                          // # Divisions into which the global timeline is divided
	vector <double> division_time;                                   // The discretised timings
};
#endif