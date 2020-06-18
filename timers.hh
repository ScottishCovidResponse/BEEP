#pragma once

struct TIMERS {
	long timetot;
	long timesim;
	long timeboot;
	long timembp;
	long timewait;
	long timembpQmap;
	long timembpprop;
	long timembpinit;
};

extern TIMERS timers;

void timersinit();