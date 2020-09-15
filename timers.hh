#ifndef BEEPMBP__TIMERS_HH
#define BEEPMBP__TIMERS_HH

#include <string>
#include <iostream>

using namespace std;

struct TIMERS {
	long timetot;
	long timesim;
	long timembp;
	long timewait;
	long timembpQmap;
	long timembpprop;
	long timembpinit;
	long timestandard;
	long timeparam;
	long timeswap;
	long timeabc;
	long timeabcprop;
};

extern TIMERS timers;

void timersinit();
void output_timers(string file);
#endif
