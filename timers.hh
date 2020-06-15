#pragma once

struct TIMERS {
	long timetot;
	long timesim;
	long timeboot;
	long timembp;
	long timewait;
};

extern TIMERS timers;

void timersinit();