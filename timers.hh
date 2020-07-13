#ifndef BEEPMBP__TIMERS_HH
#define BEEPMBP__TIMERS_HH

struct TIMERS {
	long timetot;
	long timesim;
	long timeboot;
	long timembp;
	long timewait;
	long timembpQmap;
	long timembpprop;
	long timembpinit;
	long timembptemp;	long timembptemp2;	long timembptemp3;	long timembptemp4;
	long timestandard;
	long timeparam;
	long timebetaphiloop;
	long timecompparam;
	long timeaddrem;
};

extern TIMERS timers;

void timersinit();
#endif
