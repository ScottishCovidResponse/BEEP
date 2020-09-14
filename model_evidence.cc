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
	L.resize(nchaintot); invT.resize(nchaintot);
}

void Model_Evidence::set_invT(const unsigned int samp, vector<Chain>& chain)
{
	auto pmax = 1-pow(invTmin,1.0/INVT_POWER);
	auto pmin = 1-pow(invTmax,1.0/INVT_POWER);
	
	for(auto p = 0u; p < nchaintot; p++){
		if(nchaintot == 1) invT[p] = invTmax;
		else{
			auto fac=1.0;
			if(samp < quench) fac = double(samp)/quench;
		
			auto kappa = double(p)/(nchaintot-1);
			auto ppf = pmin+kappa*(pmax-pmin);
			auto ppeff = 1-fac*(1-ppf);	
			invT[p] = pow(1-ppeff,INVT_POWER);
		}
	}
	
	for(auto& cha : chain) cha.invT = invT[cha.ch];
}

/// Stores the likelihood so the model evidence can be calculated later
void Model_Evidence::store(const Mpi &mpi, const vector<Chain> &chain) 
{
	auto nchain = chain.size();
	vector <double> Lst(nchain), Lsttot(nchaintot);
	vector <unsigned int> ch(nchain), chtot(nchaintot);
	
	for(auto p = 0u; p < nchain; p++){ Lst[p] = chain[p].initial.L; ch[p] = chain[p].ch;}
		
	MPI_Gather(&Lst[0],nchain,MPI_DOUBLE,&Lsttot[0],nchain,MPI_DOUBLE,0,MPI_COMM_WORLD);
	MPI_Gather(&ch[0],nchain,MPI_UNSIGNED,&chtot[0],nchain,MPI_UNSIGNED,0,MPI_COMM_WORLD);

	if(mpi.core == 0){
		for(auto p = 0u; p < nchaintot; p++) L[chtot[p]].push_back(Lsttot[p]);
	}
}
	
/// Calculates the model evidence
double Model_Evidence::calculate() const
{
	auto ME = 0.0;
	for(auto ch = 1u; ch < L.size(); ch++){
		auto dinvT = invT[ch-1]-invT[ch];
		auto jmax = L[ch].size();
		auto jmin = jmax/3;
		auto max = -LARGE; 
		for(auto j = jmin; j < jmax; j++){ 
			if(L[ch][j] > max) max = L[ch][j];
		}
		
		auto sum = 0.0; for(auto j = jmin; j < jmax; j++) sum += exp(dinvT*(L[ch][j]-max));
		ME += dinvT*max + log(sum/(jmax-jmin));
	}
	if(ME < -LARGE) ME = -LARGE;
	
	return ME;
}

double Model_Evidence::get_invT(unsigned int ch) const
{
	return invT[ch];
}

double Model_Evidence::average_L(unsigned int ch) const
{
	auto av = 0.0; 
	for(auto Lch : L[ch]) av += Lch;
	
	return av/L[ch].size();
}
	
