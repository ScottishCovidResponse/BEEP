/// This provides all the functions relating to generating MCMC proposals in parameter space

using namespace std;

#include "jump.hh"
#include "utils.hh"

/// Quantities in jump are intialised
void Jump::init(const vector <double> &paramv)
{
	auto nparam = paramv.size();
	mbp.resize(nparam); mbp_ntr.resize(nparam); mbp_nac.resize(nparam);     
	stand.resize(nparam); stand_ntr.resize(nparam); stand_nac.resize(nparam);
	for(auto th = 0u; th < nparam; th++){
		mbp[th] = paramv[th]/2; if(mbp[th] == 0) mbp[th] = 0.1;
		mbp_ntr[th] = 0; mbp_nac[th] = 0;
		
		stand[th] = paramv[th]/10; if(stand[th] == 0) stand[th] = 0.1;
		stand_ntr[th] = 0; stand_nac[th] = 0;
	}
	
	naddrem = 20;
	standev_ntr = 0; standev_nac = 0;
}

/// Makes a proposal by chaning a given parameter from given parameter set
vector <double> Jump::mbp_prop(const vector <double> &p, unsigned int th)
{
	vector <double> paramv = p;
	paramv[th] += normal(0,mbp[th]);  
	
	return paramv;  
}

void Jump::setburnin(unsigned int samp, unsigned int burnin)
{
	if(samp < burnin){
		if(samp < 50) alter = fast;
		else alter = slow;
	}
	else alter = none;
}

/// A MBP is accepted
void Jump::mbp_accept(unsigned int th)
{
	mbp_ntr[th]++;
	mbp_nac[th]++;
	switch(alter){
		case fast: mbp[th] *= 2; break;
		case slow: mbp[th] *= 1.1; break;
		case none: break;
	}
}

/// A MBP is rejected
void Jump::mbp_reject(unsigned int th)
{
	mbp_ntr[th]++;
	switch(alter){
		case fast: mbp[th] *= 0.5*0.95; break;
		case slow: mbp[th] *= 0.95; break;
		case none: break;
	}
}

/// A standand proposal is accepted
void Jump::stand_accept(unsigned int th)
{
	stand_ntr[th]++;
	stand_nac[th]++;
	switch(alter){
		case fast: stand[th] *= 1.05; break;
		case slow: stand[th] *= 1.01; break;
		case none: break;
	}
}

/// A standard proposal is rejected
void Jump::stand_reject(unsigned int th)
{
	stand_ntr[th]++;
	switch(alter){
		case fast: stand[th] *= 0.975; break;
		case slow: stand[th] *= 0.995; break;
		case none: break;
	}
}

void Jump::standev_accept()
{
	standev_ntr++;
	standev_nac++;
	switch(alter){
		case fast: naddrem *= 1.05; break;
		case slow: naddrem *= 1.05; break;
		case none: break;
	}
}

void Jump::standev_reject()
{
	standev_ntr++;
	switch(alter){
		case fast: naddrem *= 0.95; break;
		case slow: naddrem *= 0.95;break;
		case none: break;
	}
	if(naddrem < 1) naddrem = 1;
}

	