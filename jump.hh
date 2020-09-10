#ifndef BEEPMBP__JUMP_HH
#define BEEPMBP__JUMP_HH

#include <vector>

using namespace std;

#include "consts.hh"

class Jump                                   // Information about kernal for jumping in parameter space
{
public:
	vector <float> mbp;                          // The size of jumps in parameter space
	vector <unsigned int> mbp_ntr, mbp_nac;      // The number of jumps tried and accepted
	
	vector <float> stand;                        // The size of jumps in parameter space (fixed event sequence)
	vector <unsigned int> stand_ntr, stand_nac;  // The number of jumps tried and accepted

	float naddrem;                           	 // The size of adding and removing events
	unsigned int standev_ntr, standev_nac;    
	
	void init(const vector <double> &paramv);
	vector <double> mbp_prop(const vector <double> &p, unsigned int th);
	
	void setburnin(unsigned int samp, unsigned int _burnin);
	void mbp_accept(unsigned int th);
	void mbp_reject(unsigned int th);
	void stand_accept(unsigned int th);
	void stand_reject(unsigned int th);
	void standev_accept();
	void standev_reject();
	
private:
	Alter alter;
};

#endif
