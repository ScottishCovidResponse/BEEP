#ifndef BEEPMBP__MODELEVIDENCE_HH
#define BEEPMBP__MODELEVIDENCE_HH

#include <vector>

using namespace std;

#include "consts.hh"
#include "chain.hh"
#include "utils.hh"

class Chain;

class Model_Evidence                                                                    // Used to calculate the model evidence
{
	public:
		void init(const unsigned int _nchaintot, const unsigned int _quench, const double _invTmin, const double _invTmax);    
		void store(const Mpi &mpi, const vector<Chain>& chain);                             // Stores quantities used to calculate the model evidence
		double calculate() const;                                                           // Calculates the model evidence based on the store quantities
		double get_invT(unsigned int ch) const;                                             // Gives the inverse temperature on a chain
		double average_L(unsigned int ch) const;                                           // Gives the average likelihood on a chain     
		void set_invT(const unsigned int samp, vector<Chain>& chain);
		
	private:
		vector <double> invT;                                                               // The inverse temperatures of all the chains
		vector <vector <double> > L;                                                       // Likelihood samples for the chains
		
		double invTmin, invTmax;                                                            // The minimum and maximum inverse tenperatures that get run at
	
		unsigned int quench;                                                                // The number of samples over which the system is quenched
	
		unsigned int nchaintot;                                                             // The total number of chains (across all MPI processes);
};

#endif