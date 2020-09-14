#ifndef BEEPMBP__MODELEVIDENCE_HH
#define BEEPMBP__MODELEVIDENCE_HH

#include <vector>

using namespace std;

#include "consts.hh"
#include "chain.hh"
#include "utils.hh"

class Chain;

class Model_Evidence                                                       // Used to calculate the model evidence
{
	public:
		void init(const unsigned int _nchaintot, const unsigned int _quench, const double _invTmin, const double _invTmax);    
		void store(const Mpi &mpi, const vector<Chain>& chain);                    
		double calculate() const;                                               
		double get_invT(unsigned int ch) const;                                   
		double average_L(unsigned int ch) const;                                   
		void set_invT(const unsigned int samp, vector<Chain>& chain);
		
	private:
		vector <double> invT;                                                   // The inverse temperatures of all the chains
		vector <vector <double> > L;                                            // Likelihood samples from the chains
		
		double invTmin, invTmax;                                                // The minimum and maximum inverse temperatures 
	
		unsigned int quench;                                                    // Samples over which the system is quenched
	
		unsigned int nchaintot;                                                 // The number of chains (across all cores)
};

#endif