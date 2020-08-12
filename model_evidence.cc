// Calculates the model evidence

#include <math.h>
#include <iostream>

using namespace std;

#include "model_evidence.hh"

/// Initialises quantities for calculating the model evidence
void Model_Evidence::init(const unsigned int _nchaintot, const unsigned int _quench, const double _invTmin, const double _invTmax)
{
	nchaintot = _nchaintot;
	quench = _quench;
	invTmin = _invTmin;
	invTmax = _invTmax;
	Li.resize(nchaintot); invT.resize(nchaintot);
}

void Model_Evidence::set_invT(const unsigned int samp, vector<Chain>& chain)
{
	auto pmax = 1-pow(invTmin,1.0/Tpower);
	auto pmin = 1-pow(invTmax,1.0/Tpower);
	
	for(auto p = 0u; p < nchaintot; p++){
		if(nchaintot == 1 || duplicate == 1) invT[p] = invTmax;
		else{
			auto fac=1.0;
			if(samp < quench) fac = double(samp)/quench;
		
			auto kappa = double(p)/(nchaintot-1);
			auto ppf = pmin+kappa*(pmax-pmin);
			auto ppeff = 1-fac*(1-ppf);	
			invT[p] = pow(1-ppeff,Tpower);
		}
	}
	
	for(auto& cha : chain) cha.invT = invT[cha.ch];
}

/// Stores the likelihood so the model evidence can be calculated later
void Model_Evidence::store(const Mpi &mpi, const vector<Chain>& chain) 
{
	unsigned int nchain = chain.size();
	double L[nchain], Ltot[nchaintot];
	unsigned int ch[nchain], chtot[nchaintot];
	
	for(auto p = 0u; p < nchain; p++){ L[p] = chain[p].Li; ch[p] = chain[p].ch;}
		
	MPI_Gather(L,nchain,MPI_DOUBLE,Ltot,nchain,MPI_DOUBLE,0,MPI_COMM_WORLD);
	MPI_Gather(ch,nchain,MPI_UNSIGNED,chtot,nchain,MPI_UNSIGNED,0,MPI_COMM_WORLD);

	if(mpi.core == 0){
		for(auto p = 0u; p < nchaintot; p++) Li[chtot[p]].push_back(Ltot[p]);
	}
}
	
/// Calculates the model evidence
double Model_Evidence::calculate() const
{
	unsigned int jmin, jmax;
	double max, dinvT, ME, sum;
	
	ME = 0;
	for(auto ch = 1u; ch < Li.size(); ch++){
		dinvT = invT[ch-1]-invT[ch];
		jmax = Li[ch].size(); jmin = jmax/3;
		max = -large; 
		for(auto j = jmin; j < jmax; j++){ 
			if(Li[ch][j] > max) max = Li[ch][j];
		}
		
		sum = 0; for(auto j = jmin; j < jmax; j++) sum += exp(dinvT*(Li[ch][j]-max));
		ME += dinvT*max + log(sum/(jmax-jmin));
	}
	if(ME < -large) ME = -large;
	
	return ME;
}

double Model_Evidence::get_invT(unsigned int ch) const
{
	return invT[ch];
}

double Model_Evidence::average_Li(unsigned int ch) const
{
	double av = 0; 
	auto jmax = Li[ch].size(); 
	for(auto j = 0u; j < jmax; j++) av += Li[ch][j];
	
	return av/jmax;
}
	
