// This file contains all the functions for running a MCMC er the MBP algorithm

#include <fstream>
#include <iostream>
#include <sstream>
#include <algorithm>  
#include "stdlib.h"
#include "math.h"
#include "assert.h"

using namespace std;

#include "timers.hh"
#include "model.hh"
#include "utils.hh"
#include "chain.hh"
#include "output.hh"
#include "pack.hh"
#include "obsmodel.hh"

/// Initialises a single mcmc chain
Chain::Chain(const Details &details, const DATA &data, const MODEL &model, const POPTREE &poptree, const Obsmodel &obsmodel, unsigned int chstart) : initial(details,data,model,obsmodel), propose(details,data,model,obsmodel), comp(model.comp), lev(poptree.lev), trans(model.trans), details(details), data(data), model(model), poptree(poptree), obsmodel(obsmodel)
{
	ch = chstart;

	setuplists();
	
	indmap.resize(data.popsize);
	for(auto i = 0u; i < data.popsize; i++){
		indmap[i].resize(model.trans.size());
		for(auto tra = 0u; tra < model.trans.size(); tra++) indmap[i][tra] = 0;
	}
	
	dQmap.resize(data.narage);                                                 // Initialises vectors
	
	dQbuf.resize(data.narage);
	for(auto v = 0u; v < data.narage; v++){
		dQbuf[v].resize(data.Q.size()); for(auto q = 0u; q < data.Q.size(); q++) dQbuf[v][q] = 0;
	}
	dQbuflistv.clear(); dQbuflistq.clear();

	Rtot.resize(poptree.level); for(auto l = 0u; l < poptree.level; l++) Rtot[l].resize(lev[l].node.size()); 
	N.resize(comp.size()); 
	
	sample_from_prior();
	
	jump.init(propose.paramval);
	
	popw.resize(data.nardp);                                        // Used for event based changes
	lam.resize(data.nsettardp); lamsum.resize(data.nsettardp);
	
	if(details.mode != inf) return;
	
	initial.copy(propose);
	initial.setLPr();
}

/// Performs a "standard" set of proposals
void Chain::standard_prop(double EFcut) 
{
	initial.standard_parameter_prop(jump);
	stand_event_prop(EFcut);
}
		
/// Randomly samples from prior and generates an event sequence
void Chain::sample_from_prior()
{
	unsigned int loop, loopmax = 100;
	for(loop = 0; loop < loopmax; loop++){
		if(simulate(model.priorsamp()) == success) break;
	}
	
	if(loop == loopmax) emsg("After '"+to_string(loopmax)+"' random simulations, it was not possible to find an initial state with the number of infected individuals below the threshold 'infmax' specified in the input TOML file.");
}

/// Performs a MBP on parameter 'th'
void Chain::mbp_proposal(unsigned int th)  
{
	//timeprop -= clock();
	timers.timembpprop -= clock();

	vector <double> paramv = jump.mbp_prop(initial.paramval,th);

	double al = 0; 
	if(mbp(paramv) == success){
		propose.setLPr();
		
		al = exp(propose.Pr-initial.Pr + invT*(propose.L-initial.L));		
		if(checkon == 1) cout << al << " " << invT << " " << propose.L << " " << initial.L << " al" << endl;
	}

	if(ran() < al){
		initial.copy(propose);
		jump.mbp_accept(th);
	}
	else{
		jump.mbp_reject(th);
	}

	if(checkon == 1) initial.check();
	
	timers.timembpprop += clock();
	//timeprop += clock();	
}

/// Performs a MBP
Status Chain::mbp(const vector<double> &paramv)
{
	timers.timembpinit -= clock();
	
	if(model.inbounds(paramv) == false) return fail;
		
	if(propose.setparam(paramv) == fail) return fail;
	
	bool doev = model.dombpevents(initial.paramval,propose.paramval);

	for(auto c = 0u; c < comp.size(); c++) N[c] = 0;
	N[0] = data.popsize;
		
	propose.clear();
	
	for(auto v = 0u; v < data.narage; v++) dQmap[v] = 0;
	
	timers.timembpinit += clock();
	
	timers.timembp -= clock();
		
	unsigned int n = 0, nmax = initial.x.size();
	double t = 0;
	for(auto sett = 0u; sett < details.nsettime; sett++){
		if(details.mode == sim){
			cout  << "  Time: " << details.settime[t];
			for(auto c = 0u; c < comp.size(); c++) cout << "  " << comp[c].name << ":"	<< N[c];
			cout << endl;	
		}
		
		initial.setbetaphi(sett);	propose.setbetaphi(sett);
	
		propose.setQmapUsingdQ(sett,initial,dQmap);
	
		constructRtot(initial.Qmap[sett],propose.Qmap[sett]);
	
		double tmax = details.settime[sett+1];
		do{
			double tini;
			FEV ev;
			if(n < nmax){ ev = initial.getinfev(n); tini = ev.t;} else{ ev.ind = UNSET; tini = tmax;}
		
			double v; v = ran();
			double tinf;
			if(Rtot[0][0] <= 0) tinf = tmax; else tinf = t - log(v)/Rtot[0][0];
				
			if(tini >= tmax && tinf >= tmax){ t = tmax; break;}
			
			if(tinf < tini){  // A new infection appears
				t = tinf;
				auto c = nextinfection();
				addinfc(c,t);	
			}
			else{            // An event on initial sequence 
				t = tini;
				auto i = ev.ind;
		
				if(stat[i] == both_sus){
					auto w = data.ind[i].area*data.ndemocatpos + data.ind[i].dp;
	
					auto al = propose.lam[w]/initial.lam[w];
					if(ran() < al){                                    // Keeps the infection event
						changestat(i,not_sus,1);
						
						if(doev == true){
							model.mbpmodel(initial.indev[i],propose.indev[i],initial.paramval,propose.paramval,initial.comptrans,propose.comptrans);
							//mbpmodel(initial.indev[i],propose.indev[i]);
						}
						else propose.indev[i] = initial.indev[i];
						
						propose.addindev(i);
					}
					else changestat(i,ponly_sus,1);      // Does not keep the infection event
				}
				n++;
			}
	
			if(propose.x.size() >= model.infmax) break;
		}while(1 == 1);
		if(propose.x.size() >= model.infmax) break; 

		updatedQmap(initial.trev[sett],propose.trev[sett]);	
		if(checkon == 1) check(t,sett);
	}

	timers.timembp += clock();
		
	resetlists();

	if(propose.x.size() >= model.infmax) return fail;
	return success;
}

/// Simulates an event sequence given a parameter set (returns 1 if it is found not to be possible)
Status Chain::simulate(const vector <double>& paramv)
{
	vector <double> zero(paramv.size());
	for(auto& pval : zero) pval = 0;
	
	initial.setparam(zero);
	
	if(mbp(paramv) == fail) return fail;
	
	return success;
}

/// Constructs the tree Rtot for fast Gilelspie sampling	
void Chain::constructRtot(const vector <double> &Qmi, const vector <double> &Qmp)
{
	timers.timembpconRtot -= clock();
		
	int l = poptree.level-1;
	for(auto c = 0u; c < data.narea; c++){
		auto wmin = c*data.ndemocatpos; 
		auto wmax = wmin + data.ndemocatpos;
	
		auto faci = initial.beta*initial.areafac[c];
		auto facp = propose.beta*propose.areafac[c];
		
		double sum = 0; 
		auto dp = 0u; 
		auto v = c*data.nage; 
		for(auto w = wmin; w < wmax; w++){
			auto a = data.democatpos[dp][0];
			initial.lam[w] = initial.sus[dp]*(faci*Qmi[v+a] + initial.phi);
			propose.lam[w] = propose.sus[dp]*(facp*Qmp[v+a] + propose.phi);
			auto dlam = nindbothlist[w]*(propose.lam[w] - initial.lam[w]); if(dlam < 0) dlam = 0;
			if(std::isnan(dlam)){ emsgEC("Chain",400);}
			sum += dlam + nindponlylist[w]*propose.lam[w];
			dp++;
		}
		if(std::isnan(sum)) emsgEC("Chain",4);
		Rtot[l][c] = sum;
	}
	
	for(int l = poptree.level-2; l >= 0; l--){                                 
		auto cmax = lev[l].node.size();
		for(auto c = 0u; c < cmax; c++){
			double sum = 0; for(const auto& ch : lev[l].node[c].child) sum += Rtot[l+1][ch];
			
			Rtot[l][c] = sum;
		}
	}
	
	timers.timembpconRtot += clock();
}
	
/// Sets up lists for use with MBPs
void Chain::setuplists()
{	
	indbothlist.clear(); indponlylist.clear(); indnotlist.clear();
	indbothlist.resize(data.nardp); indponlylist.resize(data.nardp); indnotlist.resize(data.nardp);
	nindbothlist.resize(data.nardp); nindponlylist.resize(data.nardp); nindnotlist.resize(data.nardp);
	stat.resize(data.popsize); indlistref.resize(data.popsize); 

	for(auto c = 0u; c < data.narea; c++){
		for(auto dp = 0u; dp < data.ndemocatpos; dp++){
			auto w = c*data.ndemocatpos + dp;

			for(const auto& i : data.area[c].ind[dp]){
				stat[i] = both_sus;
				indlistref[i] = indbothlist[w].size();
				indbothlist[w].push_back(i);
			}
			nindbothlist[w] =  data.area[c].ind[dp].size();
			nindponlylist[w] = 0;
			nindnotlist[w] = 0;
		}
	}	
}

/// Makes a change in the status of an individual
void Chain::changestat(unsigned int i, unsigned int st, unsigned int updateR)
{
	auto c = data.ind[i].area;
	auto w = c*data.ndemocatpos + data.ind[i].dp;
	
	double dval = 0;
	if(updateR == 1){		
		auto dlam = nindbothlist[w]*(propose.lam[w] - initial.lam[w]); if(dlam < 0) dlam = 0;
		dval = -(dlam + nindponlylist[w]*propose.lam[w]);
	}
	
	int l = indlistref[i];   // Removes the einitial.xsiting entry
	int n;
	switch(stat[i]){
  case both_sus:
		if(indbothlist[w][l] != i) emsgEC("Chain",7);
		n = indbothlist[w].size();
		if(l < n-1){
			indbothlist[w][l] = indbothlist[w][n-1];
			indlistref[indbothlist[w][l]] = l;
		}
		indbothlist[w].pop_back();
		nindbothlist[w]--;
		break;
		
	case ponly_sus:
		if(indponlylist[w][l] != i) emsgEC("Chain",8);
		n = indponlylist[w].size();
		if(l < n-1){
			indponlylist[w][l] = indponlylist[w][n-1];
			indlistref[indponlylist[w][l]] = l;
		}
		indponlylist[w].pop_back();
		nindponlylist[w]--;
		break;
		
	case not_sus:
		if(indnotlist[w][l] != i) emsgEC("Chain",9);
		n = indnotlist[w].size();
		if(l < n-1){
			indnotlist[w][l] = indnotlist[w][n-1];
			indlistref[indnotlist[w][l]] = l;
		}
		indnotlist[w].pop_back();
		nindnotlist[w]--;
		break;
	
	default: emsgEC("Chain",10); break;
	}

	stat[i] = st;
	switch(stat[i]){
	case ponly_sus:
		indlistref[i] = indponlylist[w].size();
		indponlylist[w].push_back(i);
		nindponlylist[w]++;
		break;
		
	case not_sus:
		indlistref[i] = indnotlist[w].size();
		indnotlist[w].push_back(i);
		nindnotlist[w]++;
		break;
	
	case both_sus:
		indlistref[i] = indbothlist[w].size();
		indbothlist[w].push_back(i);
		nindbothlist[w]++;
		break;
		
	default: emsgEC("Chain",11); break;
	}
	
	if(updateR == 1){
		auto dlam = nindbothlist[w]*(propose.lam[w] - initial.lam[w]); if(dlam < 0) dlam = 0;
		dval += dlam + nindponlylist[w]*propose.lam[w];
		
		if(dval != 0){
			int l = poptree.level-1;
			do{
				Rtot[l][c] += dval;
				c = lev[l].node[c].parent; l--;
			}while(l >= 0);
		}
	}
}

/// Places all individuals back onto the both susceptible list 
void Chain::resetlists()
{
	for(auto w = 0u; w < data.nardp; w++){
		for(const auto& i : indponlylist[w]){
			stat[i] = both_sus;
			indlistref[i] = indbothlist[w].size();
			indbothlist[w].push_back(i);
			nindbothlist[w]++;
		}
		indponlylist[w].clear();
		nindponlylist[w] = 0;
		
		for(const auto& i : indnotlist[w]){
			stat[i] = both_sus;
			indlistref[i] = indbothlist[w].size();
			indbothlist[w].push_back(i);
			nindbothlist[w]++;
		}
		indnotlist[w].clear();
		nindnotlist[w] = 0;
	}
}

/// Updates dQmap based on events which occur in timestep sett in the initial and proposed states
void Chain::updatedQmap(const vector <EVREF> &trei, const vector <EVREF> &trep)
{
	timers.timembpQmap -= clock();
	
	auto nage = data.nage;

	for(const auto& tre : trei){
		auto i = tre.ind; 
		auto tra = initial.indev[i][tre.e].trans;
		indmap[i][tra] = 1;
	}
	
	for(const auto& tre : trep){
		auto i = tre.ind; 
		auto tra = propose.indev[i][tre.e].trans;	
		if(indmap[i][tra] == 0){
			auto v = data.ind[i].area*nage+data.democatpos[data.ind[i].dp][0];
			auto dq = trans[tra].DQ[propose.indev[i][tre.e].timep];

			if(dq != UNSET){
				for(auto loop = 0u; loop < 2; loop++){
					auto q = model.DQ[dq].q[loop];
					if(q != UNSET){
						if(dQbuf[v][q] == 0){ dQbuflistv.push_back(v); dQbuflistq.push_back(q);}
						dQbuf[v][q] += model.DQ[dq].fac[loop];
					}
				}
			}
		}
		else indmap[i][tra] = 0;
	}	
	
	for(const auto& tre : trei){
		auto i = tre.ind; 
		auto tra = initial.indev[i][tre.e].trans;
		if(indmap[i][tra] == 1){
			auto v = data.ind[i].area*nage+data.democatpos[data.ind[i].dp][0];
			auto dq = trans[tra].DQ[initial.indev[i][tre.e].timep];
			if(dq != UNSET){
				for(auto loop = 0u; loop < 2; loop++){
					auto q = model.DQ[dq].q[loop];
					if(q != UNSET){
						if(dQbuf[v][q] == 0){ dQbuflistv.push_back(v); dQbuflistq.push_back(q);}
						dQbuf[v][q] -= model.DQ[dq].fac[loop];
					}
				}
			}
			indmap[i][tra] = 0;
		}
	}
	
	if(details.mode == sim){
		for(const auto& tre : trep){
			auto tra = propose.indev[tre.ind][tre.e].trans;
			N[trans[tra].from]--;
			N[trans[tra].to]++;
		}
	}

	nage = data.nage;
	auto jmax = dQbuflistv.size();
	for(auto j = 0u; j < jmax; j++){
		auto v = dQbuflistv[j]; 
		auto q = dQbuflistq[j]; 
		auto qt = data.Q[q].Qtenref;
		
		auto fac = dQbuf[v][q];
		if(fac < -vtiny || fac > vtiny){
			auto kmax = data.genQ.Qten[qt].ntof[v];
			
			auto& cref = data.genQ.Qten[qt].tof[v];
			auto& valref = data.genQ.Qten[qt].valf[v];
			if(nage == 1){
				for(auto k = 0u; k < kmax; k++){
					dQmap[cref[k]] += fac*valref[k][0];
				}
			}
			else{
				for(auto k = 0u; k < kmax; k++){
					auto vv = cref[k]*nage;	
					for(auto a = 0u; a < nage; a++){
						dQmap[vv] += fac*valref[k][a];
						vv++;
					}
				}
			}
		}
		dQbuf[v][q] = 0;
	}
	dQbuflistv.clear(); dQbuflistq.clear(); 
	
	timers.timembpQmap += clock();
}
	
/// This samples the area in which the next infection occurs due to the MBP
unsigned int Chain::nextinfection()
{
	double sumst[4];
	
	auto l = 0u, c = 0u;                              // We start at the top level l=0 and proceed to finer and finer scales
	auto lmax = poptree.level;
	while(l < lmax-1){
		auto jmax = lev[l].node[c].child.size();
		double sum = 0;
		for(auto j = 0u; j < jmax; j++){
			sum += Rtot[l+1][lev[l].node[c].child[j]];	
			sumst[j] = sum;
		}
		
		double z = ran()*sum; auto j = 0u; while(j < jmax && z > sumst[j]) j++;
		if(j == jmax) emsgEC("Chain",12);
		
		c = lev[l].node[c].child[j];
		l++;
	};
	
	return c;
}

/// Adds an exposed individual in area
void Chain::addinfc(unsigned int c, double t)
{
	auto dpmax = data.ndemocatpos;
	
	vector <double> sumst;
	sumst.resize(dpmax);
	
	double sum = 0;                           // First selects the demographic possibility
	for(auto dp = 0u; dp < dpmax; dp++){
		auto w = c*dpmax + dp;
		double dlam = nindbothlist[w]*(propose.lam[w] - initial.lam[w]); if(dlam < 0) dlam = 0;
		sum += dlam + nindponlylist[w]*propose.lam[w];
		sumst[dp] = sum;
	}

	double z = ran()*sum; auto dp = 0u; while(dp < dpmax && z > sumst[dp]) dp++; 
	if(dp == dpmax) emsgEC("Chain",13);
	
	auto w = c*dpmax + dp;                  // Next select when individuals in the initial and proposed states are susceptible
	double dlam = nindbothlist[w]*(propose.lam[w] - initial.lam[w]); if(dlam < 0) dlam = 0;

	unsigned int i;
	if(ran() < dlam/(dlam + nindponlylist[w]*propose.lam[w])){ // Both suscetible
		auto n = indbothlist[w].size(); if(n == 0) emsgEC("Chain",14);
		i = indbothlist[w][(unsigned int)(ran()*n)];
	}
	else{                                       // Only proposed state susceptible
		auto n = indponlylist[w].size(); if(n == 0) emsgEC("Chain",15);
		i = indponlylist[w][(unsigned int)(ran()*n)];
	}
	
	changestat(i,not_sus,1);
	
	model.simmodel(propose.paramval,propose.comptrans,propose.indev[i],i,0,t);
	//propose.simmodel(i,0,t);

	propose.addindev(i);
}

/// Used for checking the code is running correctly
void Chain::check(double t, unsigned int sett) const
{
	for(auto i = 0u; i < data.popsize; i++){    // Checks stat is correct
	  auto w = data.ind[i].area*data.ndemocatpos + data.ind[i].dp;
		
		if(nindbothlist[w] != indbothlist[w].size()) emsgEC("Chain",24);
		if(nindponlylist[w] != indponlylist[w].size()) emsgEC("Chain",25);
		if(nindnotlist[w] != indnotlist[w].size()) emsgEC("Chain",26);
		
		if((initial.indev[i].size() == 0 || t < initial.indev[i][0].t) && propose.indev[i].size() == 0){
			if(stat[i] != both_sus) emsgEC("Chain",27);
			if(indbothlist[w][indlistref[i]] != i) emsgEC("Chain",28);
		}
		else{
			if((initial.indev[i].size() != 0 && t >= initial.indev[i][0].t) && propose.indev[i].size() == 0){
				if(stat[i] != ponly_sus) emsgEC("Chain",29);
				if(indponlylist[w][indlistref[i]] != i) emsgEC("Chain",30);
			}
			else{
				if(stat[i] != not_sus) emsgEC("Chain",31);
				if(indnotlist[w][indlistref[i]] != i) emsgEC("Chain",32);
			}
		}
	}
	
	auto l = poptree.level-1;
	for(auto c = 0u; c < data.narea; c++){
		auto wmin = c*data.ndemocatpos, wmax = wmin + data.ndemocatpos;
	
		double sum = 0; 
		auto dp = 0u; 
		auto v = c*data.nage; 
		for(auto w = wmin; w < wmax; w++){
			auto a = data.democatpos[dp][0];
			
			double dd;
			dd = initial.lam[w] - initial.sus[dp]*(initial.beta*initial.areafac[c]*initial.Qmap[sett][v+a] + initial.phi); if(sqrt(dd*dd) > tiny) emsgEC("Chain",33);
			dd = propose.lam[w] - propose.sus[dp]*(propose.beta*propose.areafac[c]*propose.Qmap[sett][v+a] + propose.phi); if(sqrt(dd*dd) > tiny) emsgEC("Chain",34);
	
			auto dlam = nindbothlist[w]*(propose.lam[w] - initial.lam[w]); if(dlam < 0) dlam = 0;
			sum += dlam + nindponlylist[w]*propose.lam[w];
			dp++;
		}
		auto dd = Rtot[l][c] - sum; if(sqrt(dd*dd) > tiny){ emsgEC("Chain",35);}
	}
	
	for(auto i = 0u; i < data.popsize; i++){
		for(auto tra = 0u; tra < model.trans.size(); tra++){
			if(indmap[i][tra] != 0) emsgEC("Chain",36);
		}
	}
}

/// Calculates propose.Qmap based on the initial and final sequences
void Chain::calcproposeQmap()
{	
	for(auto& dQma : dQmap) dQma = 0;
	
	for(auto sett = 0u; sett < details.nsettime; sett++){
		propose.setQmapUsingdQ(sett,initial,dQmap);
		updatedQmap(initial.trev[sett],propose.trev[sett]);	
	} 
}

/// Adds and removes infectious individuals
void Chain::stand_event_prop(double EFcut)
{	
	timers.timeaddrem -= clock();

	auto probif = 0.0, probfi = 0.0;
	
	propose.trev = initial.trev;     // Copies initial state into proposed state
	for(const auto& x : propose.x) propose.indev[x.ind].clear();
	propose.x = initial.x;
	for(const auto& x : propose.x) propose.indev[x.ind] = initial.indev[x.ind];
		 
	for(const auto& x : propose.x) changestat(x.ind,not_sus,0);
	
	if(ran() < 0.5){  // Adds individuals
		timers.timembptemp2 -= clock();
		infsampler(initial.Qmap);
		timers.timembptemp2 += clock();
		
		for(auto j = 0u; j < jump.naddrem; j++){
			if(propose.x.size() >= model.infmax){ resetlists(); return;}
			
			auto sumtot = lamsum[data.nsettardp-1]; if(sumtot == 0) emsgEC("Chain",56);
			auto z = ran()*sumtot;
			
			//k = 0; while(k < data.nsettardp && z > lamsum[k]) k += 1;
			//if(k == data.nsettardp) emsg("pr"); 
			
			auto k = 0u; auto dk = data.nsettardp/10; if(dk == 0) dk = 1;
			do{
				while(k < data.nsettardp && z > lamsum[k]) k += dk;
				if(dk == 1) break;
				if(k >= dk) k -= dk;
				dk /= 10; if(dk == 0) dk = 1;
			}while(1 == 1);
			if(k >= data.nsettardp) emsgEC("Chain",57);
		
			if(k > 0){
				if(k >= lamsum.size()) emsgEC("Chain",58);
				if(!(z < lamsum[k] && z > lamsum[k-1])) emsgEC("Chain",59);
			}
			
			probif += log(lam[k]/sumtot);

			auto sett = k/data.nardp, w = k%data.nardp;
			
			if(nindbothlist[w] == 0){ resetlists(); return;}
			auto i = indbothlist[w][int(ran()*nindbothlist[w])];
			probif += log(1.0/nindbothlist[w]);
		
			changestat(i,not_sus,0);
			
			auto dt = details.settime[sett+1]-details.settime[sett];
			auto t = details.settime[sett] + ran()*dt;
			probif += log(1.0/dt);
			
			model.simmodel(initial.paramval,propose.comptrans,propose.indev[i],i,0,t);
			//initial.simmodel(i,0,t);
		
			propose.addindev(i);
			
			probfi += log(1.0/propose.x.size());
		}
		
		timers.timembptemp3 -= clock();
		propose.sortx();
		calcproposeQmap();
		timers.timembptemp3 += clock();
	}
	else{    // Removes individuals
		vector <int> kst;
		for(auto j = 0u; j < jump.naddrem; j++){
			if(propose.x.size() == 0){ resetlists(); return;}
			
			auto l = int(ran()*propose.x.size());
			probif += log(1.0/propose.x.size());
			auto i = propose.x[l].ind;
			auto sett = (unsigned int)(details.nsettime*initial.indev[i][propose.x[l].e].t/details.period); 

			auto c = data.ind[i].area;
			auto w = c*data.ndemocatpos + data.ind[i].dp;
	
			auto dt = details.settime[sett+1]-details.settime[sett];

			probfi += log(1.0/dt);
			kst.push_back(sett*data.nardp + w);
			
			propose.indev[i].clear();
			
			propose.x[l] = propose.x[propose.x.size()-1];
			propose.x.pop_back();
			
			changestat(i,both_sus,0);
			
			probfi += log(1.0/nindbothlist[w]);
		}
		
		for(auto& trev : propose.trev){  // Removes events in propose.trev
			auto j = 0u;
			auto jmax = trev.size();
			while(j < jmax){
				if(propose.indev[trev[j].ind].size() == 0){
					jmax--;
					trev[j] = trev[jmax];
					trev.pop_back();
				}
				else j++;
			}
		}
		
		timers.timembptemp3 -= clock();
		propose.sortx();
		calcproposeQmap();
		timers.timembptemp3 += clock();
		
		timers.timembptemp2 -= clock();
		infsampler(propose.Qmap);
		timers.timembptemp2 += clock();
		
		auto sumtot = lamsum[data.nsettardp-1]; 
		for(const auto& ks : kst) probfi += log(lam[ks]/sumtot);
	}
	
	timers.timembptemp4 -= clock();
	propose.setlikelihood();
	
	double al;
	if(details.mode == abcmbp){
		propose.EF = obsmodel.Lobs(propose.trev,propose.indev);
		if(propose.EF > EFcut) al = 0;
		else al = exp(propose.Lev-initial.Lev + probfi - probif);
		if(checkon == 1) cout << al << " " << initial.EF << " " << propose.EF << " " << initial.Lev << " " << propose.Lev <<  " " << EFcut << "al" << endl;		
	}
	else{
		propose.L = obsmodel.Lobs(propose.trev,propose.indev);
		al = exp(invT*(propose.L-initial.L) + propose.Lev-initial.Lev + probfi - probif);
		if(checkon == 1) cout << al << " " << initial.L << " " << propose.L << " " << initial.Lev << " " << propose.Lev << "al" << endl;		
	}
	timers.timembptemp4 += clock();
	
	if(ran() < al){
		initial.Lev = propose.Lev;
		if(details.mode == abcmbp) initial.EF = propose.EF;
		else initial.L = propose.L;
		
		initial.trev = propose.trev;
		initial.Qmap = propose.Qmap;
		
		for(const auto& x : initial.x) initial.indev[x.ind].clear();
		for(const auto& x : propose.x) initial.indev[x.ind] = propose.indev[x.ind];
		
		initial.x = propose.x;
		jump.standev_accept();
	}
	else{
		jump.standev_reject();
	}
	
	resetlists();

	if(checkon == 1) initial.check();
	timers.timeaddrem += clock();
}

/// Generates a sampler for adding infected individuals into the system
void Chain::infsampler(const vector< vector<double> > &Qmap)
{
	auto sum = 0.0;
	for(auto sett = 0u; sett < details.nsettime; sett++){
		auto phi = initial.disc_spline[model.phispline_ref][sett]; 
		auto beta = initial.disc_spline[model.betaspline_ref][sett];
	
		for(auto c = 0u; c < data.narea; c++){
			auto fac = beta*initial.areafac[c];
			for(auto dp = 0u; dp < data.ndemocatpos; dp++){
				auto w = c*data.ndemocatpos + dp;
				auto tot = sett*data.nardp + w;
				auto v = c*data.nage + data.democatpos[dp][0];
				
				auto val = nindbothlist[w]*initial.sus[dp]*(fac*Qmap[sett][v] + phi);
				sum += val;
				
				lam[tot] = val;				
				lamsum[tot] = sum;
			}
		}
	}
}

/// Compresses the events to take up as little memory as possible (used for abcmbp)
vector <FEV> Chain::event_compress(const vector < vector <FEV> > &indev) const
{
	vector <FEV> store;
	for(const auto& inde : indev){
		for(const auto& ev : inde) store.push_back(ev);
	}
	
	return store;
}

/// Generates a particle (used for abcmbp)
void Chain::generate_particle(Particle &part) const
{
	part.EF = initial.EF;
	part.paramval = initial.paramval;
	part.ev = event_compress(initial.indev);
}

Status Chain::abcmbp_proposal(const vector <double> &paramv, double EFcut)  
{
	auto al = 1.0;
	for(auto th = 0u; th < model.param.size(); th++){
		if(paramv[th] < model.param[th].min || paramv[th] > model.param[th].max) al = 0;
		if(paramv[th] < model.param[th].min || paramv[th] > model.param[th].max) al = 0;
	}

	if(al == 1){

		if(mbp(propose.paramval) == fail) al = 0;
		else{
			propose.EF = obsmodel.Lobs(propose.trev,propose.indev);
			if(propose.EF >= EFcut) al = 0;
			else{
				propose.Pr = model.prior(propose.paramval);
				al = exp(propose.Pr-initial.Pr);
				//cout << al << " " << propose.Pr << " " << initial.Pr <<  "al\n";
			}
		}
	}
	
	if(ran() < al){
		initial.copy(propose);
		return success;
	}

	return fail;
}

