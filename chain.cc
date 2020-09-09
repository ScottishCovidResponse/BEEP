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
//#include "state.hh"

struct EVREFT {                
	unsigned int ind;                   
	unsigned int e;	              
	double t;	 
};

static bool compEVREFT(EVREFT lhs, EVREFT rhs)
{
	return lhs.t < rhs.t;
};

struct LCONT {                
	unsigned int w;       
	unsigned int num;       	
	double betafac;	              
	double phifac;	 
};

/// Initialises a single mcmc chain
Chain::Chain(const Details &details, const DATA &data, MODEL &model, const POPTREE &poptree, Obsmodel &obsmodel, unsigned int chstart) : initial(model), propose(model), comp(model.comp), lev(poptree.lev), trans(model.trans), details(details), data(data), model(model), poptree(poptree), obsmodel(obsmodel)
{
	ch = chstart;

	initial.disc_spline.resize(model.spline.size());
	propose.disc_spline.resize(model.spline.size());
	
	//initial.init();
	//propose.init();
		
	initial.x.clear();
	initial.trev.clear(); initial.trev.resize(details.nsettime); 
	initial.indev.clear(); initial.indev.resize(data.popsize);
	propose.indev.clear(); propose.indev.resize(data.popsize);
	
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
	
	initial.Qmap.resize(details.nsettime); propose.Qmap.resize(details.nsettime); 
	for(auto sett = 0u; sett < details.nsettime; sett++){
		initial.Qmap[sett].resize(data.narage); for(auto v = 0u; v < data.narage; v++) initial.Qmap[sett][v] = 0;
		propose.Qmap[sett].resize(data.narage); 
	}

	initial.lam.resize(data.nardp); propose.lam.resize(data.nardp);
	Rtot.resize(poptree.level); for(auto l = 0u; l < poptree.level; l++) Rtot[l].resize(lev[l].node.size()); 
	N.resize(comp.size()); 
	
	sample_from_prior();

	proposal_init();

	popw.resize(data.nardp);                                        // Used for event based changes
	lam.resize(data.nsettardp); lamsum.resize(data.nsettardp);
	
	if(details.mode != inf) return;
	
	initial.trev = propose.trev;
	initial.Qmap = propose.Qmap;	
	initial.indev = propose.indev;
	initial.x = propose.x;
	
	initial.L = obsmodel.Lobs(initial.trev,initial.indev);
	initial.Pr = model.prior(paramval);

	initial.disc_spline = propose.disc_spline;
	
	setinitialQmap(1);
}

/// Randomly samples from prior and generates an event sequence
void Chain::sample_from_prior()
{
	unsigned int loop, loopmax = 100;
	for(loop = 0; loop < loopmax; loop++){
		do{	paramval = model.priorsamp();}while(model.setup(paramval) == 1);  
			
		if(simulate(paramval) == 0) break;
	}
	
	if(loop == loopmax) emsg("After '"+to_string(loopmax)+"' random simulations, it was not possible to find an initial state with the number of infected individuals below the threshold 'infmax' specified in the input TOML file.");
}

/// Simulates an event sequence given a parameter set
unsigned int Chain::simulate(const vector <double>& paramv)
{
	if(model.setup(paramv) == 1) return 1;
	
	paramval = paramv;
	vector <double> paramvalinit = paramval;
	for(auto& pval : paramvalinit) pval = 0;
	
	model.setup(paramvalinit);                                       // To generate initial state mbp is used to simulate
	model.copyi(paramvalinit);
	for(auto sp = 0u; sp < model.spline.size(); sp++) initial.disc_spline[sp] = model.create_disc_spline(sp,paramvalinit);
	initial.sus = model.create_sus(paramvalinit);    
	
	model.setup(paramval);
	model.copyp(paramval);
	for(auto sp = 0u; sp < model.spline.size(); sp++) propose.disc_spline[sp] = model.create_disc_spline(sp,paramval);
	propose.sus = model.create_sus(paramval);  
	
	if(mbp() == 1) return 1;

	return 0;
}

void Chain::proposal_init()
{
	auto nparam = model.param.size();
	paramjump.resize(nparam); ntr.resize(nparam); nac.resize(nparam);         // Initialises proposal and diagnostic information
	paramjumpstand.resize(nparam); ntrstand.resize(nparam); nacstand.resize(nparam);
	for(auto th = 0u; th < nparam; th++){
		paramjump[th] = paramval[th]/2; if(paramjump[th] == 0) paramjump[th] = 0.1;
		ntr[th] = 0; nac[th] = 0;
		
		paramjumpstand[th] = paramval[th]/10; if(paramjumpstand[th] == 0) paramjumpstand[th] = 0.1;
		ntrstand[th] = 0; nacstand[th] = 0;
	}
	
	sigmajump = 0.01;
		
	timeprop = 0;
	
	numaddrem = 20;
	ntr_addrem = 0; nac_addrem = 0;
}

/// Performs a MBP
unsigned int Chain::mbp()
{
	timers.timembpinit -= clock();

	unsigned int doev = model.dombpevents();

	for(auto c = 0u; c < comp.size(); c++) N[c] = 0;
	N[0] = data.popsize;
		
	for(const auto& i : propose.x) propose.indev[i.ind].clear();
	//propose.indev.clear(); propose.indev.resize(data.popsize);	
	
	propose.x.clear();
	propose.trev.clear(); propose.trev.resize(details.nsettime);
	
	for(auto v = 0u; v < data.narage; v++) dQmap[v] = 0;
	
	timers.timembpinit += clock();
	
	timers.timembp -= clock();
		
	unsigned int n = 0;
	double t = 0;
	for(auto sett = 0u; sett < details.nsettime; sett++){
		if(details.mode == sim){
			cout  << "  Time: " << details.settime[t];
			for(auto c = 0u; c < comp.size(); c++) cout << "  " << comp[c].name << ":"	<< N[c];
			cout << endl;	
		}
		
		initial.phi = initial.disc_spline[model.phispline_ref][sett]; propose.phi = propose.disc_spline[model.phispline_ref][sett];	
		initial.beta = initial.disc_spline[model.betaspline_ref][sett]; propose.beta = propose.disc_spline[model.betaspline_ref][sett];	
	
		for(auto v = 0u; v < data.narage; v++){
			double val = initial.Qmap[sett][v] + dQmap[v];
			if(val < -tiny){ cout << val << "val\n"; emsgEC("Chain",1);}
			if(val < 0) val = 0;	
			propose.Qmap[sett][v] = val;
		}

		constructRtot(initial.Qmap[sett],propose.Qmap[sett]);
	
		double tmax = details.settime[sett+1];
		do{
			double tini;
			FEV ev;
			if(n < initial.x.size()){ ev = initial.indev[initial.x[n].ind][initial.x[n].e]; tini = ev.t;} else{ ev.ind = UNSET; tini = tmax;}
		
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
						
						if(doev == 1) model.mbpmodel(initial.indev[i],propose.indev[i]);
						else propose.indev[i] = initial.indev[i];
						
						addindev(i,propose.indev[i],propose.x,propose.trev);
					}
					else changestat(i,ponly_sus,1);      // Does not keep the infection event
				}
				n++;
			}
	
			if(propose.x.size() >= model.infmax) break;
		}while(1 == 1);
		if(propose.x.size() >= model.infmax) break; 

		updatedQmap(initial.trev[sett],propose.trev[sett]);	
		if(checkon == 1) check(1,t,sett);
	}

	timers.timembp += clock();
		
	resetlists();

	if(propose.x.size() >= model.infmax) return 1;
	return 0;
}

/// Adds an individual event sequence
void Chain::addindev(unsigned int i, vector <FEV> &indev, vector <EVREF> &x, vector <vector <EVREF> > &trev)
{
	unsigned int e, emax, se;
	EVREF evref;
	
	emax = indev.size();
	if(emax == 0) return;
	
	evref.ind = i; evref.e = 0;
	x.push_back(evref);
	for(e = 0; e < indev.size(); e++){
		evref.e = e;
		se = (unsigned int)(details.nsettime*indev[e].t/details.period); 
		if(se < details.nsettime) trev[se].push_back(evref);
	}
}
		
/// Based on the the event sequence in initial.x, this sets initial.Qmap
void Chain::setinitialQmap(unsigned int check)
{
	for(auto& dQma : dQmap) dQma = 0;

	auto nage = data.nage;
	for(auto sett = 0u; sett < details.nsettime; sett++){
		for(auto v = 0u; v < data.narage; v++){
			auto val = dQmap[v];
			if(check == 1){
				if(val < -tiny) emsgEC("Chain",2);
				if(val < initial.Qmap[sett][v]-tiny || val > initial.Qmap[sett][v]+tiny) emsgEC("Chain",3);
			}
			if(val < 0){ val = 0; dQmap[v] = 0;}	
			
			initial.Qmap[sett][v] = val;
		}
		
		for(const auto& trev : initial.trev[sett]){
			auto i = trev.ind;
			FEV fev = initial.indev[i][trev.e];

			auto v = data.ind[i].area*data.nage+data.democatpos[data.ind[i].dp][0];
			auto dq = trans[fev.trans].DQ[fev.timep];
			if(dq != UNSET){
				for(auto loop = 0u; loop < 2; loop++){
					auto q = model.DQ[dq].q[loop];
					if(q != UNSET){
						auto fac = model.DQ[dq].fac[loop];
						
						auto qt = data.Q[q].Qtenref;
						auto kmax = data.genQ.Qten[qt].ntof[v];
						auto& cref = data.genQ.Qten[qt].tof[v];
						auto& valref = data.genQ.Qten[qt].valf[v];
						if(nage == 1){
							for(auto k = 0u; k < kmax; k++){
								dQmap[cref[k]*nage] += fac*valref[k][0];
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
				}
			}
		}
	}
}

/// Constructs the tree Rtot for fast Gilelspie sampling	
void Chain::constructRtot(vector <double> &Qmi, vector <double> &Qmp)
{
	timers.timembpconRtot -= clock();
		
	int l = poptree.level-1;
	for(auto c = 0u; c < data.narea; c++){
		auto wmin = c*data.ndemocatpos; 
		auto wmax = wmin + data.ndemocatpos;
	
		auto faci = initial.beta*model.areafaci[c];
		auto facp = propose.beta*model.areafacp[c];
		
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
	
/// Performs a MBP on parameter 'th'
void Chain::proposal(unsigned int th, unsigned int samp, unsigned int burnin)  
{
	timeprop -= clock();
	timers.timembpprop -= clock();
	
	model.setup(paramval);
	model.copyi(paramval);
	for(auto sp = 0u; sp < model.spline.size(); sp++) initial.disc_spline[sp] = model.create_disc_spline(sp,paramval);
	initial.sus = model.create_sus(paramval);
			
	auto valst = paramval[th];

	paramval[th] += normal(0,paramjump[th]);               // Makes a change to a parameter

	double al; 
	propose.L = initial.L;
	propose.Pr = initial.Pr;
	if(paramval[th] < model.param[th].min || paramval[th] > model.param[th].max) al = 0;
	else{
		if(model.setup(paramval) == 1) al = 0;
		else{
			model.copyp(paramval);
			for(auto sp = 0u; sp < model.spline.size(); sp++) propose.disc_spline[sp] = model.create_disc_spline(sp,paramval);
			propose.sus = model.create_sus(paramval);
				
			if(mbp() == 1) al = 0;
			else{
				propose.L = obsmodel.Lobs(propose.trev,propose.indev);
				propose.Pr = model.prior(paramval);
				
				al = exp(propose.Pr-initial.Pr + invT*(propose.L-initial.L));		
				if(checkon == 1) cout << al << " " << invT << " " << propose.L << " " << initial.L << " al" << endl;
			}
		}
	}
	
	ntr[th]++;
	if(ran() < al){
		initial.L = propose.L;
		initial.Pr = propose.Pr;
		initial.trev = propose.trev;
		initial.Qmap = propose.Qmap;
		
		initial.disc_spline = propose.disc_spline;
		initial.sus = propose.sus;
		
		for(const auto& x : initial.x) initial.indev[x.ind].clear();
		for(const auto& x : propose.x) initial.indev[x.ind] = propose.indev[x.ind];
		//initial.indev = propose.indev;
		
		initial.x = propose.x;
		nac[th]++;
		if(samp < burnin){ if(samp < 50) paramjump[th] *= 2; else paramjump[th] *= 1.1;}
	}
	else{
		paramval[th] = valst;
		if(samp < burnin){ if(samp < 50) paramjump[th] *= 0.5; paramjump[th] *= 0.95;}
	}

	if(checkon == 1){
		model.setup(paramval);
		double dd;
		dd = initial.L - obsmodel.Lobs(initial.trev,initial.indev); if(sqrt(dd*dd) > tiny) emsgEC("Chain",5);
		dd = initial.Pr - model.prior(paramval); if(sqrt(dd*dd) > tiny) emsgEC("Chain",6);
	}
	
	timers.timembpprop += clock();
	timeprop += clock();
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
void Chain::updatedQmap(vector <EVREF> &trei, vector <EVREF> &trep)
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

/// Adds an epropose.xosed indivdual in area
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
	
	model.simmodel(paramval,propose.indev[i],i,0,t);

	addindev(i,propose.indev[i],propose.x,propose.trev);
}

/// Used for checking the code is running correctly
void Chain::check(unsigned int /* num */, double t, unsigned int sett)
{
	for(auto j = 0u; j < propose.x.size(); j++){ // Checks order
		auto i = propose.x[j].ind, e = propose.x[j].e;
		if(i >= propose.indev.size()) emsgEC("Chain",16);
		if(e >= propose.indev[i].size()) emsgEC("Chain",17);
		if(j < propose.x.size()-1){
			if(propose.indev[i][e].t > propose.indev[propose.x[j+1].ind][propose.x[j+1].e].t) emsgEC("Chain",18);
		}
	}

	for(const auto& indev : propose.indev){
		auto emax = indev.size();
		if(emax > 0){
			auto c = 0u; 
			double tt = 0;
			for(const auto& ev : indev){
				auto ttt = ev.t; if(ttt <= tt){ model.oe("here",indev); emsgEC("Chain",19);}
				auto tra = ev.trans;
				if(trans[tra].from != c) emsgEC("Chain",20);
				c = trans[tra].to; tt = ttt;
			}
			if(comp[c].trans.size() != 0) emsgEC("Chain",21);
			
			for(auto timep = 0u; timep < model.ntimeperiod; timep++){
				tt = model.timeperiod[timep].tend;
				if(tt > indev[0].t && tt < indev[emax-1].t){
					unsigned int e;
					for(e = 0u; e < emax; e++) if(indev[e].t == tt) break;
					if(timep <  model.ntimeperiod-1){
						if(e == emax) emsgEC("Chain",22);
					}
					else{
						if(e != emax) emsgEC("Chain",23);
					}
				}
			}
		}
	}
	
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
			dd = initial.lam[w] - initial.sus[dp]*(initial.beta*model.areafaci[c]*initial.Qmap[sett][v+a] + initial.phi); if(sqrt(dd*dd) > tiny) emsgEC("Chain",33);
			dd = propose.lam[w] - propose.sus[dp]*(propose.beta*model.areafacp[c]*propose.Qmap[sett][v+a] + propose.phi); if(sqrt(dd*dd) > tiny) emsgEC("Chain",34);
	
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

/// Checks that quantities used when adding and removing events are correctly updated
void Chain::check_addrem()
{
	for(auto j = 0u; j < initial.x.size(); j++){
		auto i = initial.x[j].ind, e = initial.x[j].e;
		if(i >= initial.indev.size()) emsgEC("Chain",37);
		if(e >= initial.indev[i].size()) emsgEC("Chain",38);
		if(j < initial.x.size()-1){
			if(initial.indev[i][e].t > initial.indev[initial.x[j+1].ind][initial.x[j+1].e].t) emsgEC("Chain",39);
		}
	}

	vector < vector <int> > done;
	done.resize(initial.indev.size());
	auto num = 0u;
	for(auto i = 0u; i < initial.indev.size(); i++){
		if(initial.indev[i].size() > 0){
			num++;
			done[i].resize(initial.indev[i].size());
			for(auto e = 0u; e < initial.indev[i].size(); e++) done[i][e] = 0;
		}
	}
	if(num != initial.x.size()) emsgEC("Chain",40);

	for(auto sett = 0u; sett < details.nsettime; sett++){
		for(auto j = 0u; j < initial.trev[sett].size(); j++){
			auto i = initial.trev[sett][j].ind, e = initial.trev[sett][j].e;
			if(e >= initial.indev[i].size()) emsgEC("Chain",41);
			
			auto se = (unsigned int)(details.nsettime*initial.indev[i][e].t/details.period); 
			if(se != sett) emsgEC("Chain",42);
			if(done[i][e] != 0) emsgEC("Chain",43);
			done[i][e] = 1;
		}
	}
	
	for(auto i = 0u; i < initial.indev.size(); i++){
		for(auto e = 0u; e < initial.indev[i].size(); e++){
			if(initial.indev[i][e].t < details.period){
				if(done[i][e] != 1) emsgEC("Chain",44);
			}
			else{
				if(done[i][e] != 0) emsgEC("Chain",45);
			}
		}
	}
	
	setinitialQmap(1);
}

/// Calculates the likelihood in the initial state
double Chain::likelihood(vector < vector<double> > &Qmap, vector <EVREF> &x, vector <vector<FEV> > &indev, 	vector < vector <double> > &disc_spline, vector <double> &sus)
{	
	model.setup(paramval);
		
	for(auto c = 0u; c < data.narea; c++){
		for(auto dp = 0u; dp < data.ndemocatpos; dp++){
			auto w = c*data.ndemocatpos + dp;
			popw[w] = data.area[c].ind[dp].size();
		}
	}		
			
	auto L = 0.0;
	auto t = 0.0; 
	auto n = 0u;
	for(auto sett = 0u; sett < details.nsettime; sett++){
		auto phi = disc_spline[model.phispline_ref][sett]; 
		auto beta = disc_spline[model.betaspline_ref][sett];
	
		auto tmax = details.settime[sett+1];
		
		for(auto c = 0u; c < data.narea; c++){
			auto fac = beta*model.areafac[c];
			for(auto dp = 0u; dp < data.ndemocatpos; dp++){
				auto w = c*data.ndemocatpos + dp;
				auto v = c*data.nage + data.democatpos[dp][0];
				initial.lam[w] = sus[dp]*(fac*Qmap[sett][v] + phi);		
				if(initial.lam[w] < 0) emsgEC("Chain",46);
				
				L -= initial.lam[w]*popw[w]*(tmax-t);
			}
		}
		
		while(n < x.size()){
			auto i = x[n].ind;
			auto ev = indev[i][x[n].e];
			auto tt = ev.t;
			if(tt >= tmax) break;
	
			t = tt;
			
			auto c = data.ind[i].area;
			auto w = c*data.ndemocatpos + data.ind[i].dp;
			L += log(initial.lam[w]);
			if(std::isnan(L)) emsgEC("Chain",47);
			popw[w]--;
			n++;
			
			L += initial.lam[w]*(tmax-t);
		}
		t = tmax;
	} 
	
	return L;
}


/// Calculates propose.Qmap based on the initial and final sequences
void Chain::calcproposeQmap()
{	
	for(auto& dQma : dQmap) dQma = 0;
	
	for(auto sett = 0u; sett < details.nsettime; sett++){
		for(auto v = 0u; v < data.narage; v++){
			auto val = initial.Qmap[sett][v] + dQmap[v];
			if(val < -tiny){ emsgEC("Chain",48);}
			if(val < 0) val = 0;	
			propose.Qmap[sett][v] = val;
		}
		updatedQmap(initial.trev[sett],propose.trev[sett]);	
	} 
}

/// This incorporates standard proposals which adds and removes events as well as changes parameters
void Chain::standard_prop(unsigned int samp, unsigned int burnin, double EFcut)
{
	timers.timestandard -= clock();

	model.setup(paramval);
	
	if(checkon == 1){ double dd = initial.Pr - model.prior(paramval); if(sqrt(dd*dd) > tiny){ emsgEC("Chainbegin",51);}}
	
	timers.timembptemp -= clock();
	initial.Lev = likelihood(initial.Qmap,initial.x,initial.indev,initial.disc_spline,initial.sus);
	//initial.Lev = initial.likelihood();
	timers.timembptemp += clock();
	
	
	timers.timeparam -= clock();
	betaphi_prop(samp,burnin);
	area_prop(samp,burnin);
	model.compparam_prop(samp,burnin,initial.x,initial.indev,paramval,paramjumpstand,ntrstand,nacstand,initial.Pr);
	if(model.regioneffect == 1) fixarea_prop(samp,burnin);
	timers.timeparam += clock();
		
	if(checkon == 1){ double dd = likelihood(initial.Qmap,initial.x,initial.indev,initial.disc_spline,initial.sus) - initial.Lev; 
	 //double dd = initial.likelihood() - initial.Lev; 
	if(dd*dd > tiny) emsgEC("Chain",49);}

	timers.timeaddrem -= clock();
	addrem_prop(samp,burnin,EFcut);
	timers.timeaddrem += clock();
	
	if(checkon == 1){
		double dd = likelihood(initial.Qmap,initial.x,initial.indev,initial.disc_spline,initial.sus) - initial.Lev; 
//	double dd = initial.likelihood() - initial.Lev; 
if(dd*dd > tiny) emsgEC("Chain",50);}

	if(checkon == 1){ double dd = initial.Pr - model.prior(paramval); if(sqrt(dd*dd) > tiny){ emsgEC("Chain",51);}}

	timers.timestandard += clock();
}

/// Makes proposal to beta and phi
void Chain::betaphi_prop(unsigned int samp, unsigned int burnin)
{	
	timers.timebetaphiinit -= clock();
	
	vector <unsigned int> map;
	map.resize(data.nardp); for(auto& ma : map) ma = 0;
	
	for(auto c = 0u; c < data.narea; c++){
		for(auto dp = 0u; dp < data.ndemocatpos; dp++){
			auto w = c*data.ndemocatpos + dp;
			popw[w] = data.area[c].ind[dp].size();
		}
	}		
			
	model.setup(paramval);
	
	vector <double> betafac, phifac;
	betafac.resize(details.nsettime); phifac.resize(details.nsettime);
	
	auto t = 0.0; 
	auto n = 0u;
	
	vector<	vector <LCONT> > lc;
	for(auto sett = 0u; sett < details.nsettime; sett++){
		auto tmax = details.settime[sett+1];
		vector <LCONT> lcontlist;
		
		auto betasum = 0.0, phisum = 0.0;
		for(auto c = 0u; c < data.narea; c++){
			auto fac = model.areafac[c];
			for(auto dp = 0u; dp < data.ndemocatpos; dp++){
				auto w = c*data.ndemocatpos + dp;
				auto v = c*data.nage + data.democatpos[dp][0];	
				betasum -= fac*initial.sus[dp]*initial.Qmap[sett][v]*popw[w]*(tmax-t);
				phisum -= initial.sus[dp]*popw[w]*(tmax-t);
			}
		}
		
		while(n < initial.x.size()){
			auto i = initial.x[n].ind;
			FEV ev = initial.indev[i][initial.x[n].e];
			auto tt = ev.t;
			if(tt >= tmax) break;
	
			t = tt;
			
			auto c = data.ind[i].area;
			auto fac = model.areafac[c];
			
			auto dp = data.ind[i].dp;
			auto w = c*data.ndemocatpos + dp;
			auto v = c*data.nage + data.democatpos[dp][0]; 
			
			if(map[w] == 0){
				LCONT lcont;
				lcont.w = w; lcont.betafac = fac*initial.sus[dp]*initial.Qmap[sett][v]; lcont.phifac = initial.sus[dp]; lcont.num = UNSET;
				lcontlist.push_back(lcont);
			}
			map[w]++;
		
			betasum += fac*initial.sus[dp]*initial.Qmap[sett][v]*(tmax-t);
			phisum += initial.sus[dp]*(tmax-t);
			popw[w]--;
			n++;
		}
		
		betafac[sett] = betasum; phifac[sett] = phisum;
					
		for(auto& lcont : lcontlist){
			auto w = lcont.w;
			lcont.num = map[w];
			map[w] = 0;
		}
					
		lc.push_back(lcontlist);
		
		t = tmax;
	} 
	
	vector <unsigned int> parampos;
	for(const auto& spl : model.spline[model.betaspline_ref]) parampos.push_back(spl.param);
	for(const auto& spl : model.spline[model.phispline_ref]) parampos.push_back(spl.param);
		
	timers.timebetaphiinit += clock();
		
	unsigned int loopmax = 12/parampos.size(); if(loopmax == 0) loopmax = 1;
	
	timers.timebetaphi -= clock();
	for(auto loop = 0u; loop < loopmax; loop++){
		for(auto th : parampos){
			if(model.param[th].min != model.param[th].max){
				auto valst = paramval[th];	
				paramval[th] += normal(0,paramjumpstand[th]);               // Makes a change to a parameter

				propose.Pr = initial.Pr; 
				propose.Lev = initial.Lev;

				double al;
				if(paramval[th] < model.param[th].min || paramval[th] > model.param[th].max) al = 0;
				else{
					model.setup(paramval);
					for(auto sp = 0u; sp < model.spline.size(); sp++) propose.disc_spline[sp] = model.create_disc_spline(sp,paramval);

					propose.Lev = 0; 
					for(auto sett = 0u; sett < details.nsettime; sett++){
						auto phi = propose.disc_spline[model.phispline_ref][sett]; 
						auto beta = propose.disc_spline[model.betaspline_ref][sett];
	
						
						propose.Lev += betafac[sett]*beta + phifac[sett]*phi;
						
						for(const auto& l : lc[sett]) propose.Lev += l.num*log(l.betafac*beta + l.phifac*phi);
						if(std::isnan(propose.Lev)) emsgEC("Chain",52);
					}
					
					if(smooth_spline == 1) propose.Pr = model.prior(paramval);
					al = exp(propose.Pr - initial.Pr + propose.Lev-initial.Lev);
				}
			
				ntrstand[th]++;
				if(ran() < al){
					initial.Pr = propose.Pr;
					initial.Lev = propose.Lev;
					initial.disc_spline = propose.disc_spline;
					nacstand[th]++;
					if(samp < burnin){ if(samp < 50) paramjumpstand[th] *= 1.05; else paramjumpstand[th] *= 1.01;}
				}
				else{
					paramval[th] = valst;
					if(samp < burnin){ if(samp < 50) paramjumpstand[th] *= 0.975; else  paramjumpstand[th] *= 0.995;}
				}
			}
		}
	}
	timers.timebetaphi += clock();
}

/// Makes proposal to change factors affecting transmission rates in areas
void Chain::area_prop(unsigned int samp, unsigned int burnin)
{	
	timers.timecovarinit -= clock();
	
	model.setup(paramval);

	vector <double> lamareafac, lamphifac;
	lamareafac.resize(data.nardp);
	lamphifac.resize(data.nardp);
	
	vector < vector <double> >	mult, add;
	mult.resize(data.narea);
	add.resize(data.narea);
	
	vector <double> areasum;
	areasum.resize(data.narea);
	
	for(auto c = 0u; c < data.narea; c++){
		areasum[c] = 0;
		for(auto dp = 0u; dp < data.ndemocatpos; dp++){
			auto w = c*data.ndemocatpos + dp;
			popw[w] = data.area[c].ind[dp].size();
		}
	}		
	
	auto L0 = 0.0, t = 0.0;
	auto n = 0u;
	for(auto sett = 0u; sett < details.nsettime; sett++){
		auto phi = initial.disc_spline[model.phispline_ref][sett]; 
		auto beta = initial.disc_spline[model.betaspline_ref][sett];
	
		auto tmax = details.settime[sett+1];
		
		for(auto c = 0u; c < data.narea; c++){
			for(auto dp = 0u; dp < data.ndemocatpos; dp++){
				auto w = c*data.ndemocatpos + dp;
				auto v = c*data.nage + data.democatpos[dp][0];
				lamareafac[w] = initial.sus[dp]*beta*initial.Qmap[sett][v];
				lamphifac[w] = initial.sus[dp]*phi;
				
				areasum[c] -= lamareafac[w]*popw[w]*(tmax-t);
				L0 -= lamphifac[w] *popw[w]*(tmax-t);
			}
		}
		
		while(n < initial.x.size()){
			auto i = initial.x[n].ind;
			FEV ev = initial.indev[i][initial.x[n].e];
			auto tt = ev.t;
			if(tt >= tmax) break;
	
			t = tt;
			
			auto c = data.ind[i].area;
			auto w = c*data.ndemocatpos + data.ind[i].dp;
			
			mult[c].push_back(lamareafac[w]);
			add[c].push_back(lamphifac[w]);
		
			popw[w]--;
			n++;
			
 			areasum[c] += lamareafac[w]*(tmax-t);
			L0 += lamphifac[w]*(tmax-t);
		}
		
		t = tmax;
	} 
	
	timers.timecovarinit += clock();
		
	timers.timecovar -= clock();
	
	auto num = model.covar_param.size(); if(model.regioneffect == 1) num += data.nregion;
	unsigned int loopmax = 12/num; if(loopmax == 0) loopmax = 1;
	
	for(auto loop = 0u; loop < loopmax; loop++){
		for(auto th : model.covar_param){ 
			if(model.param[th].min != model.param[th].max) area_prop2(samp,burnin,th,L0,areasum,mult,add);
		}
		
		if(model.regioneffect == 1){
			for(auto th : model.regioneff_param){ 
				if(model.param[th].min != model.param[th].max) area_prop2(samp,burnin,th,L0,areasum,mult,add);
			}
		
			auto th = model.sigma_param;
			if(model.param[th].min != model.param[th].max) area_prop2(samp,burnin,th,L0,areasum,mult,add);
		}
	}
	timers.timecovar += clock();
}

void Chain::area_prop2(unsigned int samp, unsigned int burnin, unsigned int th, double L0, vector <double> &areasum, vector < vector <double> >&mult, vector < vector <double> > &add)
{
	auto valst = paramval[th];	
	paramval[th] += normal(0,paramjumpstand[th]);               // Makes a change to a parameter

	propose.Lev=initial.Lev;
	propose.Pr = initial.Pr;
	double al;
	if(paramval[th] < model.param[th].min || paramval[th] > model.param[th].max) al = 0;
	else{
		model.setup(paramval);
		propose.sus = model.create_sus(paramval);

		propose.Lev = L0;
		for(auto c = 0u; c < data.narea; c++){
			auto fac = model.areafac[c];
			propose.Lev += areasum[c]*fac;
			auto kmax = mult[c].size();
			for(auto k = 0u; k < kmax; k++) propose.Lev += log(mult[c][k]*fac + add[c][k]);
		}
		if(std::isnan(propose.Lev)) emsgEC("Chain",53);
	
	  propose.Pr = model.prior(paramval);
		al = exp(propose.Pr-initial.Pr + propose.Lev-initial.Lev);
	}

	ntrstand[th]++;
	if(ran() < al){
		initial.Lev = propose.Lev;
		initial.Pr = propose.Pr;
		initial.sus = propose.sus;
		nacstand[th]++;
		if(samp < burnin){ if(samp < 50) paramjumpstand[th] *= 1.05; else paramjumpstand[th] *= 1.01;}
	}
	else{
		paramval[th] = valst;
		if(samp < burnin){ if(samp < 50) paramjumpstand[th] *= 0.975; else paramjumpstand[th] *= 0.995;}
	}
}

/// Makes fast proposals whilst fixing area factor 
void Chain::fixarea_prop(unsigned int samp, unsigned int burnin)
{
	const auto loopmax=10;
	for(auto loop = 0u; loop < loopmax; loop++){	
		auto th = model.sigma_param;
		auto valst = paramval[th];
		
		paramval[th] += normal(0,sigmajump);
		
		propose.Pr = initial.Pr;
		double al;
		if(paramval[th] < model.param[th].min || paramval[th] > model.param[th].max) al = 0;
		else{
			propose.Pr = model.prior(paramval);
			al = exp(propose.Pr-initial.Pr);
		}

		if(ran() < al){
			initial.Pr = propose.Pr;
			if(samp < burnin){ if(samp < 50) sigmajump *= 1.05; else sigmajump *= 1.01;}
		}
		else{
			paramval[th] = valst;
			if(samp < burnin){ if(samp < 50) sigmajump *= 0.975; else sigmajump *= 0.995;}
		}
	}
}

/// Time orders x
void Chain::sortx(vector <EVREF> &x, vector <vector <FEV> > &indev)
{
	vector <EVREFT> xt;
	for(const auto& xx : x){
		EVREFT evreft;
		evreft.ind = xx.ind; evreft.e = xx.e;
		if(indev[xx.ind].size() == 0) emsgEC("Chain",54);
		
		evreft.t = indev[xx.ind][xx.e].t;
		xt.push_back(evreft);	
	}
	sort(xt.begin(),xt.end(),compEVREFT);

	for(auto i = 0u; i < x.size(); i++){
		x[i].ind = xt[i].ind;	x[i].e = xt[i].e;
	}
}

/// Adds and removes infectious individuals
void Chain::addrem_prop(unsigned int samp, unsigned int burnin, double EFcut)
{
	model.setup(paramval);
		
	if(checkon == 1){
		auto dd = likelihood(initial.Qmap,initial.x,initial.indev,initial.disc_spline,initial.sus) - initial.Lev; 
//auto dd = initial.likelihood() - initial.Lev; 
if(dd*dd > tiny) emsgEC("Chain",55);}

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
		
		for(auto j = 0u; j < numaddrem; j++){
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
			
			model.simmodel(paramval,propose.indev[i],i,0,t);
			addindev(i,propose.indev[i],propose.x,propose.trev);
			
			probfi += log(1.0/propose.x.size());
		}
		
		timers.timembptemp3 -= clock();
		sortx(propose.x,propose.indev);
		calcproposeQmap();
		timers.timembptemp3 += clock();
	}
	else{    // Removes individuals
		vector <int> kst;
		for(auto j = 0u; j < numaddrem; j++){
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
		sortx(propose.x,propose.indev);
		calcproposeQmap();
		timers.timembptemp3 += clock();
		
		timers.timembptemp2 -= clock();
		infsampler(propose.Qmap);
		timers.timembptemp2 += clock();
		
		auto sumtot = lamsum[data.nsettardp-1]; 
		for(const auto& ks : kst) probfi += log(lam[ks]/sumtot);
	}
	
	timers.timembptemp4 -= clock();
	propose.Lev = likelihood(propose.Qmap,propose.x,propose.indev,propose.disc_spline,propose.sus);
	//propose.Lev = propose.likelihood();
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
	
	ntr_addrem++;
	if(ran() < al){
		initial.Lev = propose.Lev;
		if(details.mode == abcmbp) initial.EF = propose.EF;
		else initial.L = propose.L;
		
		initial.trev = propose.trev;
		initial.Qmap = propose.Qmap;
		
		for(const auto& x : initial.x) initial.indev[x.ind].clear();
		for(const auto& x : propose.x) initial.indev[x.ind] = propose.indev[x.ind];
		//initial.indev = propose.indev;
		
		initial.x = propose.x;
		nac_addrem++;
		if(samp < burnin) numaddrem *= 1.05;
	}
	else{
		if(samp < burnin){ numaddrem *= 0.95; if(numaddrem < 1) numaddrem = 1;}
	}
	
	resetlists();

	if(checkon == 1) check_addrem();
}

/// Generates a sampler for adding infected individuals into the system
void Chain::infsampler(vector< vector<double> > &Qmap)
{
	auto sum = 0.0;
	for(auto sett = 0u; sett < details.nsettime; sett++){
		auto phi = initial.disc_spline[model.phispline_ref][sett]; 
		auto beta = initial.disc_spline[model.betaspline_ref][sett];
	
		for(auto c = 0u; c < data.narea; c++){
			auto fac = beta*model.areafac[c];
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

/// Initialises the chain based on a particle (used for abcmbp)
void Chain::initialise_from_particle(const Particle &part)
{
	paramval = part.paramval;
	model.setup(paramval);
	
	for(const auto& x : initial.x) initial.indev[x.ind].clear();   // Removes the einitial.xsiting initial sequence 
	initial.x.clear();
	
	initial.trev.clear(); initial.trev.resize(details.nsettime); 
	
	vector <int> indlist;
	for(const auto& ev : part.ev){
		int i = ev.ind;
		if(initial.indev[i].size() == 0) indlist.push_back(i);
		initial.indev[i].push_back(ev);
	}	
	
	for(auto i : indlist) addindev(i,initial.indev[i],initial.x,initial.trev);
	
	sortx(initial.x,initial.indev);
	
	setinitialQmap(0);

	initial.EF = part.EF;
	initial.Pr = model.prior(paramval);
	
	if(initial.EF != obsmodel.Lobs(initial.trev,initial.indev)) emsg("Observation does not agree");
}

/// Generates a particle (used for abcmbp)
void Chain::generate_particle(Particle &part) const
{
	part.EF = initial.EF;
	part.paramval = paramval;
	part.ev = event_compress(initial.indev);
}

int Chain::abcmbp_proposal(const vector <double> &param_propose, double EFcut)  
{
	model.setup(paramval);
	model.copyi(paramval);
	for(auto sp = 0u; sp < model.spline.size(); sp++) initial.disc_spline[sp] = model.create_disc_spline(sp,paramval);
	initial.sus = model.create_sus(paramval);
		
	vector <double> valst = paramval;
	paramval = param_propose;
	
	auto al = 1.0;
	for(auto th = 0u; th < model.param.size(); th++){
		if(paramval[th] < model.param[th].min || paramval[th] > model.param[th].max) al = 0;
		if(paramval[th] < model.param[th].min || paramval[th] > model.param[th].max) al = 0;
	}

	propose.EF = initial.EF;
	propose.Pr = initial.Pr;
	if(al == 1){
		if(model.setup(paramval) == 1) al = 0;
		else{
			model.copyp(paramval);
			for(auto sp = 0u; sp < model.spline.size(); sp++) propose.disc_spline[sp] = model.create_disc_spline(sp,paramval);
			propose.sus = model.create_sus(paramval);
				
			if(mbp() == 1) al = 0;
			else{
				propose.EF = obsmodel.Lobs(propose.trev,propose.indev);
				if(propose.EF >= EFcut) al = 0;
				else{
					propose.Pr = model.prior(paramval);
					al = exp(propose.Pr-initial.Pr);
					//cout << al << " " << propose.Pr << " " << initial.Pr <<  "al\n";
				}
			}
		}
	}
	
	if(ran() < al){
		initial.EF = propose.EF;
		initial.Pr = propose.Pr;
		initial.trev = propose.trev;
		initial.Qmap = propose.Qmap;
		initial.sus = propose.sus;
		
		for(const auto& x : initial.x) initial.indev[x.ind].clear();
		for(const auto& x : propose.x) initial.indev[x.ind] = propose.indev[x.ind];
		//initial.indev = propose.indev;
		
		initial.x = propose.x;
		return 1;
	}
	else{
		paramval = valst;
	}
	return 0;
}

