// Stores the CPU clock times for different parts of the algorithm
 
#include "timers.hh"

TIMERS timers;

void timersinit()
{
	timers.timetot = 0;
	timers.timesim = 0;
	timers.timeboot = 0;
	timers.timembp = 0;
	timers.timewait = 0;
	timers.timembpQmap = 0;
	timers.timembpprop = 0;
}