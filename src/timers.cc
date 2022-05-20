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
		if(time_av[TIME_LIKELIHOOD_APPROX] > 0) dia << per(time_av[TIME_LIKELIHOOD_APPROX]/time_av[TIME_ALG]) << " Likelihood" << endl;
		if(time_av[TIME_OBS_APPROX] > 0) dia << per(time_av[TIME_OBS_APPROX]/time_av[TIME_ALG]) << " Obs approx" << endl;
		
		if(time_av[TIME_CMAES] > 0) dia << per(time_av[TIME_CMAES]/time_av[TIME_ALG]) << " CMAES algorithm" << endl;
		
		if(time_av[TIME_GENERATE_SAMPLES] > 0) dia << per(time_av[TIME_GENERATE_SAMPLES]/time_av[TIME_ALG]) << " CMAES generate samples" << endl;
		
		if(time_av[TIME_TEMP1] > 0) dia << per(time_av[TIME_TEMP1]/time_av[TIME_ALG]) << " Temp1" << endl;
		if(time_av[TIME_TEMP2] > 0) dia << per(time_av[TIME_TEMP2]/time_av[TIME_ALG]) << " Temp2" << endl;
		if(time_av[TIME_TEMP3] > 0) dia << per(time_av[TIME_TEMP3]/time_av[TIME_ALG]) << " Temp3" << endl;
		if(time_av[TIME_TEMP4] > 0) dia << per(time_av[TIME_TEMP4]/time_av[TIME_ALG]) << " Temp4" << endl;
			
		if(time_av[TIME_INV_MATRIX] > 0) dia << per(time_av[TIME_INV_MATRIX]/time_av[TIME_ALG]) << " Invert matrix" << endl;
	
		if(time_av[TIME_DETERMINANT] > 0) dia << per(time_av[TIME_DETERMINANT]/time_av[TIME_ALG]) << " Finding determinant" << endl;
		
		if(time_av[TIME_LINEAR_EQ] > 0) dia << per(time_av[TIME_LINEAR_EQ]/time_av[TIME_ALG]) << " Invert matrix" << endl;
	
		if(time_av[TIME_MATRIX_MULT] > 0) dia << per(time_av[TIME_MATRIX_MULT]/time_av[TIME_ALG]) << " Matrix multiply" << endl;
	

		if(time_av[TIME_ADD_REMOVE_S] > 0) dia << per(time_av[TIME_ADD_REMOVE_S]/time_av[TIME_ALG]) << " Add remove S" << endl;
	
		if(time_av[TIME_FUTURE_OBS_APPROX] > 0) dia << per(time_av[TIME_FUTURE_OBS_APPROX]/time_av[TIME_ALG]) << " Obs approx" << endl;
		
		if(time_av[TIME_EF_CALCULATE] > 0) dia << per(time_av[TIME_EF_CALCULATE]/time_av[TIME_ALG]) << " EF calculate" << endl;
		
		if(time_av[TIME_COVAR] > 0) dia << per(time_av[TIME_COVAR]/time_av[TIME_ALG]) << " Covariance matrix update" << endl;
		if(time_av[TIME_CORRECT] > 0) dia << per(time_av[TIME_CORRECT]/time_av[TIME_ALG]) << " Correct matrix" << endl;
		if(time_av[TIME_GRAD] > 0) dia << per(time_av[TIME_GRAD]/time_av[TIME_ALG]) << " Transition gradients" << endl;
		if(time_av[TIME_POSTERIOR_SAMPLE] > 0) dia << per(time_av[TIME_POSTERIOR_SAMPLE]/time_av[TIME_ALG]) << " Generate Posterior sample" << endl;
	
		if(time_av[TIME_SCALE_COVARIANCE] > 0) dia << per(time_av[TIME_SCALE_COVARIANCE]/time_av[TIME_ALG]) << " Scale covariance" << endl;
	
		if(time_av[TIME_PMCMCLIKE] > 0) dia << per(time_av[TIME_PMCMCLIKE]/time_av[TIME_ALG]) << " PMCMC Likelihood" << endl;
		if(time_av[TIME_BOOTSTRAP] > 0) dia << per(time_av[TIME_BOOTSTRAP]/time_av[TIME_ALG]) << " PMCMC Bootstrap" << endl;
		if(time_av[TIME_PMCMCSWAP] > 0) dia << per(time_av[TIME_PMCMCSWAP]/time_av[TIME_ALG]) << " PMCMC Swapping state information"  << endl;
		
		if(time_av[TIME_MCMCPROP] > 0) dia << per(time_av[TIME_MCMCPROP]/time_av[TIME_ALG]) << " Type I MBPs" << endl;
		if(time_av[TIME_STATEPROP] > 0) dia << per(time_av[TIME_STATEPROP]/time_av[TIME_ALG]) << " Type II MBPs" << endl;
		if(time_av[TIME_PARAMPROP] > 0) dia << per(time_av[TIME_PARAMPROP]/time_av[TIME_ALG]) << " Parameter Proposals" << endl;
	
		if(time_av[TIME_UPDATE] > 0) dia << per(time_av[TIME_UPDATE]/time_av[TIME_ALG]) << " MCMC Update" << endl;
		
		if(time_av[TIME_TRANSNUM] > 0) dia << per(time_av[TIME_TRANSNUM]/time_av[TIME_ALG]) << " Changes to transnum " << endl;
		if(time_av[TIME_UPDATEPOP] > 0) dia << per(time_av[TIME_UPDATEPOP]/time_av[TIME_ALG]) << " Update populations" << endl;
		if(time_av[TIME_UPDATEIMAP] > 0) dia << per(time_av[TIME_UPDATEIMAP]/time_av[TIME_ALG]) << " Update infection map" << endl;
		if(time_av[TIME_SIMULATE] > 0) dia << per(time_av[TIME_SIMULATE]/time_av[TIME_ALG]) << " Simulation" << endl;
		if(time_av[TIME_SETPARAM] > 0) dia << per(time_av[TIME_SETPARAM]/time_av[TIME_ALG]) << " Setup model" << endl;
		if(time_av[TIME_INITFROMPART] > 0) dia << per(time_av[TIME_INITFROMPART]/time_av[TIME_ALG]) << " Initialise state from particle" << endl;
		if(time_av[TIME_TRANSMEAN] > 0) dia << per(time_av[TIME_TRANSMEAN]/time_av[TIME_ALG]) << " Calculating means for transitions" << endl;
	
		if(time_av[TIME_OBSPROB] > 0) dia << per(time_av[TIME_OBSPROB]/time_av[TIME_ALG]) << " Obsevation probability / Error function" << endl;
		
		if(time_av[TIME_STATESAMPLE] > 0) dia << per(time_av[TIME_STATESAMPLE]/time_av[TIME_ALG]) << " PMCMC Generating state samples" << endl;
		if(time_av[TIME_CREATEN] > 0) dia << per(time_av[TIME_CREATEN]/time_av[TIME_ALG]) << " Create N " << endl;
		if(time_av[TIME_BETA_FROM_R] > 0) dia << per(time_av[TIME_BETA_FROM_R]/time_av[TIME_ALG]) << " Beta from R " << endl;
		
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
			if(time_av[TIME_COVAR_AREA] > 0) dia << per(time_av[TIME_COVAR_AREA]/time_av[TIME_ALG]) << " Covar area" << endl;
			if(time_av[TIME_FIXEDTREE] > 0) dia << per(time_av[TIME_FIXEDTREE]/time_av[TIME_ALG]) << " Fixed tree" << endl;
			if(time_av[TIME_SLICETIME] > 0) dia << per(time_av[TIME_SLICETIME]/time_av[TIME_ALG]) << " Splice time" << endl;
		}
		 
		if(all_diagnostics == true || true){
			dia << endl << "Timings for different subsections:" << endl;
			if(time_av[TIME_MBP] > 0) dia << per(time_av[TIME_MBP]/time_av[TIME_ALG]) << " MBP " << endl;
			if(time_av[TIME_MBPINIT] > 0) dia << per(time_av[TIME_MBPINIT]/time_av[TIME_ALG]) << " MBP initialisation " << endl;
			if(time_av[TIME_GEN] > 0) dia << per(time_av[TIME_GEN]/time_av[TIME_ALG]) << " Gen" << endl;
		}
		
		if(time_av[TIME_CUTOFF] > 0) dia << per(time_av[TIME_CUTOFF]/time_av[TIME_ALG]) << " Cutoff" << endl;
		if(time_av[TIME_PROP] > 0) dia << per(time_av[TIME_PROP]/time_av[TIME_ALG]) << " Prop" << endl;
		
		if(time_av[TIME_MVNSETUP] > 0) dia << per(time_av[TIME_MVNSETUP]/time_av[TIME_ALG]) << " Mvn setup" << endl;
	
		if(time_av[TIME_MBPUPDATE] > 0) dia << per(time_av[TIME_MBPUPDATE]/time_av[TIME_ALG]) << " Mbp update" << endl;
		if(time_av[TIME_UPDATEPROP] > 0) dia << per(time_av[TIME_UPDATEPROP]/time_av[TIME_ALG]) << " Update prop" << endl;
	
		dia << endl;
	}
}
