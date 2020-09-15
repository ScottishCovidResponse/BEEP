// Stores the CPU clock times for different parts of the algorithm
 
#include "timers.hh"
#include <fstream>

using namespace std;

TIMERS timers;

void timersinit()
{
	timers.timetot = 0;
	timers.timesim = 0;
	timers.timembp = 0;
	timers.timewait = 0;
	timers.timembpQmap = 0;
	timers.timembpprop = 0;
	timers.timestandard = 0;	
	timers.timeparam = 0;	
	timers.timeswap = 0;
	timers.timeabc = 0;
	timers.timeabcprop = 0;
}

void output_timers(string file)
{
	ofstream timings(file.c_str()); 

	timings << endl << "Timings for different parts of the algorithm:" << endl;
	timings << double(timers.timewait)/CLOCKS_PER_SEC << " MBP waiting time (seconds)" << endl;
	timings << double(timers.timembp)/CLOCKS_PER_SEC << " MBP time (seconds)" << endl;
	timings << double(timers.timembpinit)/CLOCKS_PER_SEC << " MBP init (seconds)" << endl;
	timings << double(timers.timembpQmap)/CLOCKS_PER_SEC << " MBP Qmap (seconds)" << endl;
	timings << double(timers.timembpprop)/CLOCKS_PER_SEC << " MBP prop (seconds)" << endl;
	timings << double(timers.timestandard)/CLOCKS_PER_SEC << " Standard (seconds)" << endl;			
	timings << double(timers.timeparam)/CLOCKS_PER_SEC << " Param (seconds)" << endl;			
	timings << double(timers.timeswap)/CLOCKS_PER_SEC << " Swap (seconds)" << endl;	
}