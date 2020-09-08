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
Chain::Chain(const Details &details, const DATA &data, MODEL &model, const POPTREE &poptree, Obsmodel &obsmodel, unsigned int chstart) : comp(model.comp), lev(poptree.lev), trans(model.trans), details(details), data(data), model(model), poptree(poptree), obsmodel(obsmodel)
{
	ch = chstart;

	xi.clear();
	trevi.clear(); trevi.resize(details.nsettime); 
	indevi.clear(); indevi.resize(data.popsize);
	indevp.clear(); indevp.resize(data.popsize);
	
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
	
	Qmapi.resize(details.nsettime); Qmapp.resize(details.nsettime); 
	for(auto sett = 0u; sett < details.nsettime; sett++){
		Qmapi[sett].resize(data.narage); for(auto v = 0u; v < data.narage; v++) Qmapi[sett][v] = 0;
		Qmapp[sett].resize(data.narage); 
	}

	lami.resize(data.nardp); lamp.resize(data.nardp);
	Rtot.resize(poptree.level); for(auto l = 0u; l < poptree.level; l++) Rtot[l].resize(lev[l].node.size()); 
	N.resize(comp.size()); 
	
	sample_from_prior();

	proposal_init();

	popw.resize(data.nardp);                                        // Used for event based changes
	lam.resize(data.nsettardp); lamsum.resize(data.nsettardp);
	
	if(details.mode != inf) return;
	
	trevi = trevp;
	Qmapi = Qmapp;	
	indevi = indevp;
	xi = xp;
	
	Li = obsmodel.Lobs(trevi,indevi);
	Pri = model.prior(paramval);

	setQmapi(1);
}

/// Randomly samples from prior and generates an event sequence
void Chain::sample_from_prior()
{
	const auto loopmax = 100;
	auto loop = 0;
	do{
		//if(details.mode == MODE_INF) cout << ch << "Initialisation try: " << loop << endl;
		do{	paramval = model.priorsamp(); }while(model.setup(paramval) == 1);             // Randomly samples parameters from the prior	

		vector <double> paramvalinit = paramval;

		 // Sets the initial state to zero force of infection
		for(auto &spl : model.betaspline) paramvalinit[spl.param] = 0; 
		for(auto &spl : model.phispline) paramvalinit[spl.param] = 0;
			
		model.setup(paramvalinit);                                       // To generate initial state mbp is used to simulate
		model.copyi(paramvalinit);
		model.setup(paramval);
		model.copyp(paramval);

		if(mbp() == 0) break;

		loop++;
	}while(loop < loopmax);                                           // Checks not too many infected individuals (based on prior)
	if(loop == loopmax) emsg("After '"+to_string(loopmax)+"' random simulations, it was not possible to find an initial state with the number of infected individuals below the threshold 'infmax' specified in the input TOML file.");
	
}

/// Simulates an event sequence given a parameter set
unsigned int Chain::simulate(const vector <double>& paramv)
{
	if(model.setup(paramv) == 1) return 1;
	
	auto nparam = model.param.size();                   
	paramval.resize(nparam); for(auto th = 0u; th < nparam; th++) paramval[th] = paramv[th];
	vector <double> paramvalinit = paramval;

	// Sets the initial state to zero force of infection
	for(auto &spl : model.betaspline) paramvalinit[spl.param] = 0; 
	for(auto &spl : model.phispline) paramvalinit[spl.param] = 0;
		
	model.setup(paramvalinit);                                       // To generate initial state mbp is used to simulate
	model.copyi(paramvalinit);
	model.setup(paramval);
	model.copyp(paramval);

	if(mbp() == 1) return 1;

	return 0;
}

void Chain::proposal_init()
{
	auto nparam = model.param.size();
	paramjump.resize(nparam); ntr.resize(nparam); nac.resize(nparam);         // Initialises proposal and diagnostic information
	paramjumpxi.resize(nparam); ntrxi.resize(nparam); nacxi.resize(nparam);
	for(auto th = 0u; th < nparam; th++){
		paramjump[th] = paramval[th]/2; if(paramjump[th] == 0) paramjump[th] = 0.1;
		ntr[th] = 0; nac[th] = 0;
		
		paramjumpxi[th] = paramval[th]/10; if(paramjumpxi[th] == 0) paramjumpxi[th] = 0.1;
		ntrxi[th] = 0; nacxi[th] = 0;
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
		
	for(auto &i : xp) indevp[i.ind].clear();
	//indevp.clear(); indevp.resize(data.popsize);	
	
	xp.clear();
	trevp.clear(); trevp.resize(details.nsettime);
	
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
		
		phii = model.phii[sett]; phip = model.phip[sett];	
		betai = model.betai[sett]; betap = model.betap[sett];

		for(auto v = 0u; v < data.narage; v++){
			double val = Qmapi[sett][v] + dQmap[v];
			if(val < -tiny){ cout << val << "val\n"; emsgEC("Chain",1);}
			if(val < 0) val = 0;	
			Qmapp[sett][v] = val;
		}

		constructRtot(Qmapi[sett],Qmapp[sett]);
	
		double tmax = details.settime[sett+1];
		do{
			double txi;
			FEV ev;
			if(n < xi.size()){ ev = indevi[xi[n].ind][xi[n].e]; txi = ev.t;} else{ ev.ind = UNSET; txi = tmax;}
		
			double v; v = ran();
			double tinf;
			if(Rtot[0][0] <= 0) tinf = tmax; else tinf = t - log(v)/Rtot[0][0];
				
			if(txi >= tmax && tinf >= tmax){ t = tmax; break;}
			
			if(tinf < txi){  // A new infection appears
				t = tinf;
				auto c = nextinfection();
				addinfc(c,t);	
			}
			else{            // An event on initial sequence 
				t = txi;
				auto i = ev.ind;
		
				if(stat[i] == both_sus){
					auto w = data.ind[i].area*data.ndemocatpos + data.ind[i].dp;
	
					auto al = lamp[w]/lami[w];
					if(ran() < al){                                    // Keeps the infection event
						changestat(i,not_sus,1);
						
						if(doev == 1) model.mbpmodel(indevi[i],indevp[i]);
						else indevp[i] = indevi[i];
						
						addindev(i,indevp[i],xp,trevp);
					}
					else changestat(i,ponly_sus,1);      // Does not keep the infection event
				}
				n++;
			}
	
			if(xp.size() >= model.infmax) break;
		}while(1 == 1);
		if(xp.size() >= model.infmax) break; 

		updatedQmap(trevi[sett],trevp[sett]);	
		if(checkon == 1) check(1,t,sett);
	}

	timers.timembp += clock();
		
	resetlists();

	if(xp.size() >= model.infmax) return 1;
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
		
/// Based on the the event sequence in xi, this sets Qmapi
void Chain::setQmapi(unsigned int check)
{
	for(auto &dQma : dQmap) dQma = 0;

	auto nage = data.nage;
	for(auto sett = 0u; sett < details.nsettime; sett++){
		for(auto v = 0u; v < data.narage; v++){
			auto val = dQmap[v];
			if(check == 1){
				if(val < -tiny) emsgEC("Chain",2);
				if(val < Qmapi[sett][v]-tiny || val > Qmapi[sett][v]+tiny) emsgEC("Chain",3);
			}
			if(val < 0){ val = 0; dQmap[v] = 0;}	
			
			Qmapi[sett][v] = val;
		}
		
		for(auto &trev : trevi[sett]){
			auto i = trev.ind;
			FEV fev = indevi[i][trev.e];

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
	//unsigned int c, cmax, wmin, wmax, dp, v, a, w, j, jmax;
	//int l;
	//double sum, dlam, faci, facp;
	
	timers.timembpconRtot -= clock();
		
	int l = poptree.level-1;
	for(auto c = 0u; c < data.narea; c++){
		auto wmin = c*data.ndemocatpos; 
		auto wmax = wmin + data.ndemocatpos;
	
		auto faci = betai*model.areafaci[c];
		auto facp = betap*model.areafacp[c];
		
		double sum = 0; 
		auto dp = 0u; 
		auto v = c*data.nage; 
		for(auto w = wmin; w < wmax; w++){
			auto a = data.democatpos[dp][0];
			lami[w] = model.susi[dp]*(faci*Qmi[v+a] + phii);
			lamp[w] = model.susp[dp]*(facp*Qmp[v+a] + phip);
			auto dlam = nindbothlist[w]*(lamp[w] - lami[w]); if(dlam < 0) dlam = 0;
			if(std::isnan(dlam)){ emsgEC("Chain",400);}
			sum += dlam + nindponlylist[w]*lamp[w];
			dp++;
		}
		if(std::isnan(sum)) emsgEC("Chain",4);
		Rtot[l][c] = sum;
	}
	
	for(int l = poptree.level-2; l >= 0; l--){                                 
		auto cmax = lev[l].node.size();
		for(auto c = 0u; c < cmax; c++){
			double sum = 0; for(auto &ch : lev[l].node[c].child) sum += Rtot[l+1][ch];
			
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
			
	auto valst = paramval[th];

	paramval[th] += normal(0,paramjump[th]);               // Makes a change to a parameter

	double al, Lp=Li, Prp=Pri;
	if(paramval[th] < model.param[th].min || paramval[th] > model.param[th].max) al = 0;
	else{
		if(model.setup(paramval) == 1) al = 0;
		else{
			model.copyp(paramval);
			if(mbp() == 1) al = 0;
			else{
				Lp = obsmodel.Lobs(trevp,indevp);
				Prp = model.prior(paramval);
				
				al = exp(Prp-Pri + invT*(Lp-Li));		
				if(checkon == 1) cout << al << " " << invT << " " << Lp << " " << Li << " al" << endl;
			}
		}
	}
	
	ntr[th]++;
	if(ran() < al){
		Li = Lp;
		Pri = Prp;
		trevi = trevp;
		Qmapi = Qmapp;
		
		for(auto &x : xi) indevi[x.ind].clear();
		for(auto &x : xp) indevi[x.ind] = indevp[x.ind];
		//indevi = indevp;
		
		xi = xp;
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
		dd = Li - obsmodel.Lobs(trevi,indevi); if(sqrt(dd*dd) > tiny) emsgEC("Chain",5);
		dd = Pri - model.prior(paramval); if(sqrt(dd*dd) > tiny) emsgEC("Chain",6);
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

			for(auto &i : data.area[c].ind[dp]){
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
		auto dlam = nindbothlist[w]*(lamp[w] - lami[w]); if(dlam < 0) dlam = 0;
		dval = -(dlam + nindponlylist[w]*lamp[w]);
	}
	
	int l = indlistref[i];   // Removes the exisiting entry
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
		auto dlam = nindbothlist[w]*(lamp[w] - lami[w]); if(dlam < 0) dlam = 0;
		dval += dlam + nindponlylist[w]*lamp[w];
		
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
		for(auto &i : indponlylist[w]){
			stat[i] = both_sus;
			indlistref[i] = indbothlist[w].size();
			indbothlist[w].push_back(i);
			nindbothlist[w]++;
		}
		indponlylist[w].clear();
		nindponlylist[w] = 0;
		
		for(auto &i : indnotlist[w]){
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

	for(auto &tre : trei){
		auto i = tre.ind; 
		auto tra = indevi[i][tre.e].trans;
		indmap[i][tra] = 1;
	}
	
	for(auto &tre : trep){
		auto i = tre.ind; 
		auto tra = indevp[i][tre.e].trans;	
		if(indmap[i][tra] == 0){
			auto v = data.ind[i].area*nage+data.democatpos[data.ind[i].dp][0];
			auto dq = trans[tra].DQ[indevp[i][tre.e].timep];

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
	
	for(auto &tre : trei){
		auto i = tre.ind; 
		auto tra = indevi[i][tre.e].trans;
		if(indmap[i][tra] == 1){
			auto v = data.ind[i].area*nage+data.democatpos[data.ind[i].dp][0];
			auto dq = trans[tra].DQ[indevi[i][tre.e].timep];
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
		for(auto &tre : trep){
			auto tra = indevp[tre.ind][tre.e].trans;
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

/// Adds an exposed indivdual in area
void Chain::addinfc(unsigned int c, double t)
{
	auto dpmax = data.ndemocatpos;
	
	vector <double> sumst;
	sumst.resize(dpmax);
	
	double sum = 0;                           // First selects the demographic possibility
	for(auto dp = 0u; dp < dpmax; dp++){
		auto w = c*dpmax + dp;
		double dlam = nindbothlist[w]*(lamp[w] - lami[w]); if(dlam < 0) dlam = 0;
		sum += dlam + nindponlylist[w]*lamp[w];
		sumst[dp] = sum;
	}

	double z = ran()*sum; auto dp = 0u; while(dp < dpmax && z > sumst[dp]) dp++; 
	if(dp == dpmax) emsgEC("Chain",13);
	
	auto w = c*dpmax + dp;                  // Next select when individuals in the initial and proposed states are susceptible
	double dlam = nindbothlist[w]*(lamp[w] - lami[w]); if(dlam < 0) dlam = 0;

	unsigned int i;
	if(ran() < dlam/(dlam + nindponlylist[w]*lamp[w])){ // Both suscetible
		auto n = indbothlist[w].size(); if(n == 0) emsgEC("Chain",14);
		i = indbothlist[w][(unsigned int)(ran()*n)];
	}
	else{                                       // Only proposed state susceptible
		auto n = indponlylist[w].size(); if(n == 0) emsgEC("Chain",15);
		i = indponlylist[w][(unsigned int)(ran()*n)];
	}
	
	changestat(i,not_sus,1);
	
	model.simmodel(paramval,indevp[i],i,0,t);

	addindev(i,indevp[i],xp,trevp);
}

/// Used for checking the code is running correctly
void Chain::check(unsigned int /* num */, double t, unsigned int sett)
{
	for(auto j = 0u; j < xp.size(); j++){ // Checks order
		auto i = xp[j].ind, e = xp[j].e;
		if(i >= indevp.size()) emsgEC("Chain",16);
		if(e >= indevp[i].size()) emsgEC("Chain",17);
		if(j < xp.size()-1){
			if(indevp[i][e].t > indevp[xp[j+1].ind][xp[j+1].e].t) emsgEC("Chain",18);
		}
	}

	for(auto &indev : indevp){
		auto emax = indev.size();
		if(emax > 0){
			auto c = 0u; 
			double tt = 0;
			for(auto &ev : indev){
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
		
		if((indevi[i].size() == 0 || t < indevi[i][0].t) && indevp[i].size() == 0){
			if(stat[i] != both_sus) emsgEC("Chain",27);
			if(indbothlist[w][indlistref[i]] != i) emsgEC("Chain",28);
		}
		else{
			if((indevi[i].size() != 0 && t >= indevi[i][0].t) && indevp[i].size() == 0){
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
			dd = lami[w] - model.susi[dp]*(betai*model.areafaci[c]*Qmapi[sett][v+a] + phii); if(sqrt(dd*dd) > tiny) emsgEC("Chain",33);
			dd = lamp[w] - model.susp[dp]*(betap*model.areafacp[c]*Qmapp[sett][v+a] + phip); if(sqrt(dd*dd) > tiny) emsgEC("Chain",34);
	
			auto dlam = nindbothlist[w]*(lamp[w] - lami[w]); if(dlam < 0) dlam = 0;
			sum += dlam + nindponlylist[w]*lamp[w];
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
	for(auto j = 0u; j < xi.size(); j++){
		auto i = xi[j].ind, e = xi[j].e;
		if(i >= indevi.size()) emsgEC("Chain",37);
		if(e >= indevi[i].size()) emsgEC("Chain",38);
		if(j < xi.size()-1){
			if(indevi[i][e].t > indevi[xi[j+1].ind][xi[j+1].e].t) emsgEC("Chain",39);
		}
	}

	vector < vector <int> > done;
	done.resize(indevi.size());
	auto num = 0u;
	for(auto i = 0u; i < indevi.size(); i++){
		if(indevi[i].size() > 0){
			num++;
			done[i].resize(indevi[i].size());
			for(auto e = 0u; e < indevi[i].size(); e++) done[i][e] = 0;
		}
	}
	if(num != xi.size()) emsgEC("Chain",40);

	for(auto sett = 0u; sett < details.nsettime; sett++){
		for(auto j = 0u; j < trevi[sett].size(); j++){
			auto i = trevi[sett][j].ind, e = trevi[sett][j].e;
			if(e >= indevi[i].size()) emsgEC("Chain",41);
			
			auto se = (unsigned int)(details.nsettime*indevi[i][e].t/details.period); 
			if(se != sett) emsgEC("Chain",42);
			if(done[i][e] != 0) emsgEC("Chain",43);
			done[i][e] = 1;
		}
	}
	
	for(auto i = 0u; i < indevi.size(); i++){
		for(auto e = 0u; e < indevi[i].size(); e++){
			if(indevi[i][e].t < details.period){
				if(done[i][e] != 1) emsgEC("Chain",44);
			}
			else{
				if(done[i][e] != 0) emsgEC("Chain",45);
			}
		}
	}
	
	setQmapi(1);
}

/// Calculates the likelihood in the initial state
double Chain::likelihood(vector < vector<double> > &Qmap, vector <EVREF> &x, vector <vector<FEV> > &indev)
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
		auto beta = model.beta[sett], phi = model.phi[sett];
		auto tmax = details.settime[sett+1];
		
		for(auto c = 0u; c < data.narea; c++){
			auto fac = beta*model.areafac[c];
			for(auto dp = 0u; dp < data.ndemocatpos; dp++){
				auto w = c*data.ndemocatpos + dp;
				auto v = c*data.nage + data.democatpos[dp][0];
				lami[w] = model.sus[dp]*(fac*Qmap[sett][v] + phi);		
				if(lami[w] < 0) emsgEC("Chain",46);
				
				L -= lami[w]*popw[w]*(tmax-t);
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
			L += log(lami[w]);
			if(std::isnan(L)) emsgEC("Chain",47);
			popw[w]--;
			n++;
			
			L += lami[w]*(tmax-t);
		}
		t = tmax;
	} 
	
	return L;
}

/// Calculates Qmapp based on the initial and final sequences
void Chain::calcQmapp()
{	
	for(auto &dQma : dQmap) dQma = 0;
	
	for(auto sett = 0u; sett < details.nsettime; sett++){
		for(auto v = 0u; v < data.narage; v++){
			auto val = Qmapi[sett][v] + dQmap[v];
			if(val < -tiny){ emsgEC("Chain",48);}
			if(val < 0) val = 0;	
			Qmapp[sett][v] = val;
		}
		updatedQmap(trevi[sett],trevp[sett]);	
	} 
}

/// This incorporates standard proposals which adds and removes events as well as changes parameters
void Chain::standard_prop(unsigned int samp, unsigned int burnin, double EFcut)
{
	timers.timestandard -= clock();

	model.setup(paramval);
	
	if(checkon == 1){ double dd = Pri - model.prior(paramval); if(sqrt(dd*dd) > tiny){ emsgEC("Chainbegin",51);}}
	
	timers.timembptemp -= clock();
	Levi = likelihood(Qmapi,xi,indevi);
	timers.timembptemp += clock();
	
	
	timers.timeparam -= clock();
	betaphi_prop(samp,burnin);
	area_prop(samp,burnin);
	model.compparam_prop(samp,burnin,xi,indevi,paramval,paramjumpxi,ntrxi,nacxi,Pri);
	if(model.regioneffect == 1) fixarea_prop(samp,burnin);
	timers.timeparam += clock();
		
	if(checkon == 1){ double dd = likelihood(Qmapi,xi,indevi) - Levi; if(dd*dd > tiny) emsgEC("Chain",49);}

	timers.timeaddrem -= clock();
	addrem_prop(samp,burnin,EFcut);
	timers.timeaddrem += clock();
	
	if(checkon == 1){ double dd = likelihood(Qmapi,xi,indevi) - Levi; if(dd*dd > tiny) emsgEC("Chain",50);}

	if(checkon == 1){ double dd = Pri - model.prior(paramval); if(sqrt(dd*dd) > tiny){ emsgEC("Chain",51);}}

	timers.timestandard += clock();
}

/// Makes proposal to beta and phi
void Chain::betaphi_prop(unsigned int samp, unsigned int burnin)
{	
	timers.timebetaphiinit -= clock();
	
	vector <unsigned int> map;
	map.resize(data.nardp); for(auto &ma : map) ma = 0;
	
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
				betasum -= fac*model.sus[dp]*Qmapi[sett][v]*popw[w]*(tmax-t);
				phisum -= model.sus[dp]*popw[w]*(tmax-t);
			}
		}
		
		while(n < xi.size()){
			auto i = xi[n].ind;
			FEV ev = indevi[i][xi[n].e];
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
				lcont.w = w; lcont.betafac = fac*model.sus[dp]*Qmapi[sett][v]; lcont.phifac = model.sus[dp]; lcont.num = UNSET;
				lcontlist.push_back(lcont);
			}
			map[w]++;
		
			betasum += fac*model.sus[dp]*Qmapi[sett][v]*(tmax-t);
			phisum += model.sus[dp]*(tmax-t);
			popw[w]--;
			n++;
		}
		
		betafac[sett] = betasum; phifac[sett] = phisum;
					
		for(auto &lcont : lcontlist){
			auto w = lcont.w;
			lcont.num = map[w];
			map[w] = 0;
		}
					
		lc.push_back(lcontlist);
		
		t = tmax;
	} 
	
	vector <unsigned int> parampos;
	for(auto &spl : model.betaspline) parampos.push_back(spl.param);
	for(auto &spl : model.phispline) parampos.push_back(spl.param);
		
	timers.timebetaphiinit += clock();
		
	unsigned int loopmax = 12/parampos.size(); if(loopmax == 0) loopmax = 1;
	
	timers.timebetaphi -= clock();
	for(auto loop = 0u; loop < loopmax; loop++){
		for(auto &th : parampos){
			if(model.param[th].min != model.param[th].max){
				auto valst = paramval[th];	
				paramval[th] += normal(0,paramjumpxi[th]);               // Makes a change to a parameter

				auto Prp = Pri, Levp = Levi;
				double al;
				if(paramval[th] < model.param[th].min || paramval[th] > model.param[th].max) al = 0;
				else{
					model.setup(paramval);

					Levp = 0; 
					for(auto sett = 0u; sett < details.nsettime; sett++){
						auto beta = model.beta[sett], phi = model.phi[sett];
						Levp += betafac[sett]*beta + phifac[sett]*phi;
						
						for(auto &l : lc[sett]) Levp += l.num*log(l.betafac*beta + l.phifac*phi);
						if(std::isnan(Levp)) emsgEC("Chain",52);
					}
					
					if(smooth_spline == 1) Prp = model.prior(paramval);
					al = exp(Prp - Pri + Levp-Levi);
				}
			
				ntrxi[th]++;
				if(ran() < al){
					Pri = Prp;
					Levi = Levp;
					nacxi[th]++;
					if(samp < burnin){ if(samp < 50) paramjumpxi[th] *= 1.05; else paramjumpxi[th] *= 1.01;}
				}
				else{
					paramval[th] = valst;
					if(samp < burnin){ if(samp < 50) paramjumpxi[th] *= 0.975; else  paramjumpxi[th] *= 0.995;}
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
		auto beta = model.beta[sett], phi = model.phi[sett];
		auto tmax = details.settime[sett+1];
		
		for(auto c = 0u; c < data.narea; c++){
			for(auto dp = 0u; dp < data.ndemocatpos; dp++){
				auto w = c*data.ndemocatpos + dp;
				auto v = c*data.nage + data.democatpos[dp][0];
				lamareafac[w] = model.sus[dp]*beta*Qmapi[sett][v];
				lamphifac[w] = model.sus[dp]*phi;
				
				areasum[c] -= lamareafac[w]*popw[w]*(tmax-t);
				L0 -= lamphifac[w] *popw[w]*(tmax-t);
			}
		}
		
		while(n < xi.size()){
			auto i = xi[n].ind;
			FEV ev = indevi[i][xi[n].e];
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
		for(auto &th : model.covar_param){ 
			if(model.param[th].min != model.param[th].max) area_prop2(samp,burnin,th,L0,areasum,mult,add);
		}
		
		if(model.regioneffect == 1){
			for(auto &th : model.regioneff_param){ 
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
	paramval[th] += normal(0,paramjumpxi[th]);               // Makes a change to a parameter

	auto Levp=Levi, Prp=Pri;
	double al;
	if(paramval[th] < model.param[th].min || paramval[th] > model.param[th].max) al = 0;
	else{
		model.setup(paramval);

		Levp = L0;
		for(auto c = 0u; c < data.narea; c++){
			auto fac = model.areafac[c];
			Levp += areasum[c]*fac;
			auto kmax = mult[c].size();
			for(auto k = 0u; k < kmax; k++) Levp += log(mult[c][k]*fac + add[c][k]);
		}
		if(std::isnan(Levp)) emsgEC("Chain",53);
	
	  Prp = model.prior(paramval);
		al = exp(Prp-Pri + Levp-Levi);
	}

	ntrxi[th]++;
	if(ran() < al){
		Levi = Levp;
		Pri = Prp;
		nacxi[th]++;
		if(samp < burnin){ if(samp < 50) paramjumpxi[th] *= 1.05; else paramjumpxi[th] *= 1.01;}
	}
	else{
		paramval[th] = valst;
		if(samp < burnin){ if(samp < 50) paramjumpxi[th] *= 0.975; else paramjumpxi[th] *= 0.995;}
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
		
		auto Prp=Pri;
		double al;
		if(paramval[th] < model.param[th].min || paramval[th] > model.param[th].max) al = 0;
		else{
			auto sd = paramval[model.sigma_param];
			Prp = 0; for(auto &th : model.regioneff_param) Prp += normalprob(paramval[th],0,sd*sd);
			al = exp(Prp-Pri);
		}

		if(ran() < al){
			Pri = Prp;
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
	for(auto &xx : x){
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
		
	if(checkon == 1){ auto dd = likelihood(Qmapi,xi,indevi) - Levi; if(dd*dd > tiny) emsgEC("Chain",55);}

	auto probif = 0.0, probfi = 0.0;
	
	trevp = trevi;     // Copies initial state into proposed state
	for(auto &x : xp) indevp[x.ind].clear();
	xp = xi;
	for(auto &x : xp) indevp[x.ind] = indevi[x.ind];
		 
	for(auto &x : xp) changestat(x.ind,not_sus,0);
	
	if(ran() < 0.5){  // Adds individuals
		timers.timembptemp2 -= clock();
		infsampler(Qmapi);
		timers.timembptemp2 += clock();
		
		for(auto j = 0u; j < numaddrem; j++){
			if(xp.size() >= model.infmax){ resetlists(); return;}
			
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
			
			model.simmodel(paramval,indevp[i],i,0,t);
			addindev(i,indevp[i],xp,trevp);
			
			probfi += log(1.0/xp.size());
		}
		
		timers.timembptemp3 -= clock();
		sortx(xp,indevp);
		calcQmapp();
		timers.timembptemp3 += clock();
	}
	else{    // Removes individuals
		vector <int> kst;
		for(auto j = 0u; j < numaddrem; j++){
			if(xp.size() == 0){ resetlists(); return;}
			
			auto l = int(ran()*xp.size());
			probif += log(1.0/xp.size());
			auto i = xp[l].ind;
			auto sett = (unsigned int)(details.nsettime*indevi[i][xp[l].e].t/details.period); 

			auto c = data.ind[i].area;
			auto w = c*data.ndemocatpos + data.ind[i].dp;
	
			auto dt = details.settime[sett+1]-details.settime[sett];

			probfi += log(1.0/dt);
			kst.push_back(sett*data.nardp + w);
			
			indevp[i].clear();
			
			xp[l] = xp[xp.size()-1];
			xp.pop_back();
			
			changestat(i,both_sus,0);
			
			probfi += log(1.0/nindbothlist[w]);
		}
		
		for(auto &trev : trevp){  // Removes events in trevp
			auto j = 0u;
			auto jmax = trev.size();
			while(j < jmax){
				if(indevp[trev[j].ind].size() == 0){
					jmax--;
					trev[j] = trev[jmax];
					trev.pop_back();
				}
				else j++;
			}
		}
		
		timers.timembptemp3 -= clock();
		sortx(xp,indevp);
		calcQmapp();
		timers.timembptemp3 += clock();
		
		timers.timembptemp2 -= clock();
		infsampler(Qmapp);
		timers.timembptemp2 += clock();
		
		auto sumtot = lamsum[data.nsettardp-1]; 
		for(auto &ks : kst) probfi += log(lam[ks]/sumtot);
	}
	
	timers.timembptemp4 -= clock();
	auto Levp = likelihood(Qmapp,xp,indevp);
	
	double al, Lp=0, EFp=0;
	if(details.mode == abcmbp){
		EFp = obsmodel.Lobs(trevp,indevp);
		if(EFp > EFcut) al = 0;
		else al = exp(Levp-Levi + probfi - probif);
		if(checkon == 1) cout << al << " " << EF << " " << EFp << " " << Levi << " " << Levp <<  " " << EFcut << "al" << endl;		
	}
	else{
		Lp = obsmodel.Lobs(trevp,indevp);
		al = exp(invT*(Lp-Li) + Levp-Levi + probfi - probif);
		if(checkon == 1) cout << al << " " << Li << " " << Lp << " " << Levi << " " << Levp << "al" << endl;		
	}
	timers.timembptemp4 += clock();
	
	ntr_addrem++;
	if(ran() < al){
		Levi = Levp;
		if(details.mode == abcmbp) EF = EFp;
		else Li = Lp;
		
		trevi = trevp;
		Qmapi = Qmapp;
		
		for(auto &x : xi) indevi[x.ind].clear();
		for(auto &x : xp) indevi[x.ind] = indevp[x.ind];
		//indevi = indevp;
		
		xi = xp;
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
		auto beta = model.beta[sett], phi = model.phi[sett];
	
		for(auto c = 0u; c < data.narea; c++){
			auto fac = beta*model.areafac[c];
			for(auto dp = 0u; dp < data.ndemocatpos; dp++){
				auto w = c*data.ndemocatpos + dp;
				auto tot = sett*data.nardp + w;
				auto v = c*data.nage + data.democatpos[dp][0];
				
				auto val = nindbothlist[w]*model.sus[dp]*(fac*Qmap[sett][v] + phi);
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
	for(auto &inde : indev){
		for(auto &ev : inde) store.push_back(ev);
	}
	
	return store;
}

/// Initialises the chain based on a particle (used for abcmbp)
void Chain::initialise_from_particle(const Particle &part)
{
	paramval = part.paramval;
	model.setup(paramval);
	
	for(auto &x : xi) indevi[x.ind].clear();   // Removes the exisiting initial sequence 
	xi.clear();
	
	trevi.clear(); trevi.resize(details.nsettime); 
	
	vector <int> indlist;
	for(auto &ev : part.ev){
		int i = ev.ind;
		if(indevi[i].size() == 0) indlist.push_back(i);
		indevi[i].push_back(ev);
	}	
	
	for(auto i : indlist) addindev(i,indevi[i],xi,trevi);
	
	sortx(xi,indevi);
	
	setQmapi(0);

	EF = part.EF;
	Pri = model.prior(paramval);
	
	if(EF != obsmodel.Lobs(trevi,indevi)) emsg("Observation does not agree");
}

/// Generates a particle (used for abcmbp)
void Chain::generate_particle(Particle &part) const
{
	part.EF = EF;
	part.paramval = paramval;
	part.ev = event_compress(indevi);
}

int Chain::abcmbp_proposal(const vector <double> &param_propose, double EFcut)  
{
	model.setup(paramval);
	model.copyi(paramval);
			
	vector <double> valst = paramval;
	paramval = param_propose;
	
	auto al = 1.0;
	for(auto th = 0u; th < model.param.size(); th++){
		if(paramval[th] < model.param[th].min || paramval[th] > model.param[th].max) al = 0;
	}

	double EFp = EF, Prp = Pri;
	if(al == 1){
		if(model.setup(paramval) == 1) al = 0;
		else{
			model.copyp(paramval);
			if(mbp() == 1) al = 0;
			else{
				EFp = obsmodel.Lobs(trevp,indevp);
				if(EFp >= EFcut) al = 0;
				else{
					Prp = model.prior(paramval);
					al = exp(Prp-Pri);
					//cout << al << " " << Prp << " " << Pri <<  "al\n";
				}
			}
		}
	}
	
	if(ran() < al){
		EF = EFp;
		Pri = Prp;
		trevi = trevp;
		Qmapi = Qmapp;
		
		for(auto &x : xi) indevi[x.ind].clear();
		for(auto &x : xp) indevi[x.ind] = indevp[x.ind];
		//indevi = indevp;
		
		xi = xp;
		return 1;
	}
	else{
		paramval = valst;
	}
	return 0;
}

