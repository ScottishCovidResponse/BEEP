#ifndef BEEP__TIMERS_HH
#define BEEP__TIMERS_HH

#include <string>
#include <iostream>
#include <vector>

using namespace std;

#include "utils.hh"
#include "mpi.hh"

struct Timer { 
	void start();
	void stop();
	
	long val;
};

extern vector <Timer> timer;

void timersinit();
void output_timers(string file, Mpi &mpi);
#endif
