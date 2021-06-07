// Stores the CPU clock times for different parts of the algorithm 

#include <fstream>
#include <sstream>

using namespace std;

#include "timers.hh"
#include "mpi.hh"

vector <Timer> timer;

void Timer::start()
{
	val -= clock();
}

void Timer::stop()
{
	val += clock();
}

void timersinit()
{
	timer.resize(TIMERMAX);
	
	for(auto &ti : timer) ti.val = 0;
}


///Outputs CPU timing information to a file
void output_timers(string file, Mpi &mpi)
{
	vector <double> time_av(TIMERMAX);
	for(auto i = 0u; i < TIMERMAX; i++) time_av[i] = mpi.average(timer[i].val);
	
	if(mpi.core == 0){
		ofstream dia(file); if(!dia) emsg("Cannot open the file '"+file+"'");
	
		dia.precision(3);
			
		dia << endl << "Timings for different parts of the algorithm:" << endl;
		if(time_av[TIME_PMCMCLIKE] > 0) dia << per(time_av[TIME_PMCMCLIKE]/time_av[TIME_ALG]) << " PMCMC Likelihood" << endl;
		if(time_av[TIME_BOOTSTRAP] > 0) dia << per(time_av[TIME_BOOTSTRAP]/time_av[TIME_ALG]) << " PMCMC Bootstrap" << endl;
		if(time_av[TIME_PMCMCSWAP] > 0) dia << per(time_av[TIME_PMCMCSWAP]/time_av[TIME_ALG]) << " PMCMC Swapping state information"  << endl;
		
		if(time_av[TIME_MCMCPROP] > 0) dia << per(time_av[TIME_MCMCPROP]/time_av[TIME_ALG]) << " MCMC Proposals" << endl;
		if(time_av[TIME_TRANSNUM] > 0) dia << per(time_av[TIME_TRANSNUM]/time_av[TIME_ALG]) << " Changes to transnum " << endl;
		if(time_av[TIME_UPDATEPOP] > 0) dia << per(time_av[TIME_UPDATEPOP]/time_av[TIME_ALG]) << " Update populations" << endl;
		if(time_av[TIME_UPDATEIMAP] > 0) dia << per(time_av[TIME_UPDATEIMAP]/time_av[TIME_ALG]) << " Update infection map" << endl;
		if(time_av[TIME_SIMULATE] > 0) dia << per(time_av[TIME_SIMULATE]/time_av[TIME_ALG]) << " Simulation" << endl;
		if(time_av[TIME_SETPARAM] > 0) dia << per(time_av[TIME_SETPARAM]/time_av[TIME_ALG]) << " Setup model" << endl;
		if(time_av[TIME_INITFROMPART] > 0) dia << per(time_av[TIME_INITFROMPART]/time_av[TIME_ALG]) << " Initialise state from particle" << endl;
		if(time_av[TIME_TRANSMEAN] > 0) dia << per(time_av[TIME_TRANSMEAN]/time_av[TIME_ALG]) << " Calculating means for transitions" << endl;
	
		if(time_av[TIME_OBSPROB] > 0) dia << per(time_av[TIME_OBSPROB]/time_av[TIME_ALG]) << " Obsevation probability / Error function" << endl;
		
		if(time_av[TIME_STATESAMPLE] > 0) dia << per(time_av[TIME_STATESAMPLE]/time_av[TIME_ALG]) << " PMCMC Generating state samples" << endl;
		if(time_av[TIME_RESULTS] > 0) dia << per(time_av[TIME_RESULTS]/time_av[TIME_ALG]) << " Generating final results " << endl;
		if(time_av[TIME_WAIT] > 0) dia << per(time_av[TIME_WAIT]/time_av[TIME_ALG]) << " MPI Waiting" << endl;
		
		
		if(time_av[TIME_MVN] > 0){
			dia << endl << "Timings for different MCMC proposals:" << endl;
			if(time_av[TIME_SELF] > 0) dia << per(time_av[TIME_SELF]/time_av[TIME_ALG]) << " Self" << endl;
			if(time_av[TIME_MVN] > 0) dia << per(time_av[TIME_MVN]/time_av[TIME_ALG]) << " Univariate / MVN" << endl;
			if(time_av[TIME_SIGMA] > 0) dia << per(time_av[TIME_SIGMA]/time_av[TIME_ALG]) << " Sigma" << endl;
			if(time_av[TIME_MEANTIME] > 0) dia << per(time_av[TIME_MEANTIME]/time_av[TIME_ALG]) << " Mean Time" << endl;
			if(time_av[TIME_NEIGHBOUR] > 0) dia << per(time_av[TIME_NEIGHBOUR]/time_av[TIME_ALG]) << " Spline neighbour" << endl;
			if(time_av[TIME_JOINT] > 0) dia << per(time_av[TIME_JOINT]/time_av[TIME_ALG]) << " Spline joint" << endl;
			if(time_av[TIME_FIXEDTREE] > 0) dia << per(time_av[TIME_FIXEDTREE]/time_av[TIME_ALG]) << " Fixed tree" << endl;
			if(time_av[TIME_SLICETIME] > 0) dia << per(time_av[TIME_SLICETIME]/time_av[TIME_ALG]) << " Splice time" << endl;
		}
		 
		if(all_diagnostics == true){
			dia << endl << "Timings for different subsections:" << endl;
			if(time_av[TIME_MBP] > 0) dia << per(time_av[TIME_MBP]/time_av[TIME_ALG]) << " MBP " << endl;
			if(time_av[TIME_MBPINIT] > 0) dia << per(time_av[TIME_MBPINIT]/time_av[TIME_ALG]) << " MBP initialisation " << endl;
			if(time_av[TIME_GEN] > 0) dia << per(time_av[TIME_GEN]/time_av[TIME_ALG]) << " Gen" << endl;
		}
		
		dia << endl;
	}
}
