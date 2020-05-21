
#include "timers.hh"

TIMERS timers; // Stores the CPU clock times for different parts of the algorithm

void timersinit()
{
	timers.timetot = 0;
	timers.timesim = 0;
	timers.timeboot = 0;
}