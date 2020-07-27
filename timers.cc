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
	timers.timembptemp = 0;
	timers.timembptemp2 = 0;	timers.timembptemp3 = 0;	timers.timembptemp4 = 0;
	timers.timestandard = 0;	
	timers.timeparam = 0;	
	
	timers.timebetaphi = 0;
	timers.timebetaphiinit = 0;
	timers.timecovarinit = 0;
	timers.timecovar = 0;
	timers.timeaddrem = 0;	
	timers.timecompparam = 0;
	timers.timembpconRtot = 0;
	
	timers.timeswap = 0;
	timers.timeoutput = 0;
}