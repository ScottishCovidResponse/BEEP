#pragma once

#include <string>
#include <vector>

using namespace std;

struct NEV {                               // Information about the immediate next events
  short type; double t;
};

struct FEV {                               // Stores information about a compartmental transition
  long trans;                              // References the transition type
	long ind;                                // The individual on which the transition happens
	double t;                                // The time of the transition
	short done;                              // Set to 1 if that transition is in the past 
};

// Stores the CPU clock times for different parts of the algorithm
struct TIMERS {
	long timetot;
	long timesim;
	long timeboot;
};
