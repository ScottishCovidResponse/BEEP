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
Chain::Chain(const Details &details, const DATA &data, MODEL &model, const POPTREE &poptree, Obsmodel &obsmodel, unsigned int chstart) : comp(model.comp), lev(poptree.lev), trans(model.trans), details(details), data(data), model(model), poptree(poptree), obsmodel(obsmodel)
{
	unsigned int v, q, sett, i, tra;
	int l;

	ch = chstart;

	xi.clear();
	trevi.clear(); trevi.resize(details.nsettime); 
	indevi.clear(); indevi.resize(data.popsize);
	indevp.clear(); indevp.resize(data.popsize);
	
	setuplists();
	
	indmap.resize(data.popsize);
	for(i = 0; i < data.popsize; i++){
		indmap[i].resize(model.trans.size());
		for(tra = 0; tra < model.trans.size(); tra++) indmap[i][tra] = 0;
	}
	
	dQmap.resize(data.narage);                                                 // Initialises vectors
	
	dQbuf.resize(data.narage);
	for(v = 0; v < data.narage; v++){
		dQbuf[v].resize(data.Q.size()); for(q = 0; q < data.Q.size(); q++) dQbuf[v][q] = 0;
	}
	dQbuflistv.clear(); dQbuflistq.clear();
	
	Qmapi.resize(details.nsettime); Qmapp.resize(details.nsettime); 
	for(sett = 0; sett < details.nsettime; sett++){
		Qmapi[sett].resize(data.narage); for(v = 0; v < data.narage; v++) Qmapi[sett][v] = 0;
		Qmapp[sett].resize(data.narage); 
	}

	lami.resize(data.nardp); lamp.resize(data.nardp);
	Rtot.resize(poptree.level); for(l = 0; l < (int)poptree.level; l++) Rtot[l].resize(lev[l].node.size()); 
	N.resize(comp.size()); 
	
	sample_from_prior();

	if(details.mode != inf) return;
	
	trevi = trevp;
	Qmapi = Qmapp;	
	indevi = indevp;
	xi = xp;
	
	Li = obsmodel.Lobs(trevi,indevi);
	Pri = model.prior();

	setQmapi(1);

	proposal_init();

	popw.resize(data.nardp);                                        // Used for event based changes
	lam.resize(data.nsettardp); lamsum.resize(data.nsettardp);
}

/// Randomly samples from prior and generates an event sequence
void Chain::sample_from_prior()
{
	const auto loopmax = 100;
	auto loop = 0;
	do{
		//if(details.mode == MODE_INF) cout << ch << "Initialisation try: " << loop << endl;
		do{	model.priorsamp(); }while(model.setup(model.paramval) == 1);             // Randomly samples parameters from the prior	

		auto nparam = model.param.size();                   
		paramval.resize(nparam); for(auto th = 0u; th < nparam; th++) paramval[th] = model.paramval[th];
		vector <double> paramvalinit = paramval;

		 // Sets the initial state to zero force of infection
		for(auto j = 0u; j < model.betaspline.size(); j++) paramvalinit[model.betaspline[j].param] = 0; 
		for(auto j = 0u; j < model.phispline.size(); j++) paramvalinit[model.phispline[j].param] = 0;
			
		model.setup(paramvalinit);                                       // To generate initial state mbp is used to simulate
		model.copyi();
		model.setup(paramval);
		model.copyp();

		if(mbp() == 0) break;

		loop++;
	}while(loop < loopmax);                          // Checks not too many infected (based on prior)
	if(loop == loopmax) emsg("After '"+to_string(loopmax)+"' random simulations, it was not possible to find an initial state with the number of infected individuals below the threshold 'infmax' specified in the input TOML file.");
	
}

/// Simulates an event sequence given a 
unsigned int Chain::simulate(const vector <double>& paramv)
{
	if(model.setup(paramv) == 1) return 1;
	
	auto nparam = model.param.size();                   
	paramval.resize(nparam); for(auto th = 0u; th < nparam; th++) paramval[th] = model.paramval[th];
	vector <double> paramvalinit = paramval;

	 // Sets the initial state to zero force of infection
	for(auto j = 0u; j < model.betaspline.size(); j++) paramvalinit[model.betaspline[j].param] = 0; 
	for(auto j = 0u; j < model.phispline.size(); j++) paramvalinit[model.phispline[j].param] = 0;
		
	model.setup(paramvalinit);                                       // To generate initial state mbp is used to simulate
	model.copyi();
	model.setup(paramval);
	model.copyp();

	if(mbp() == 1) return 1;

	return 0;
}

void Chain::proposal_init()
{
	auto nparam = model.param.size();
	paramjump.resize(nparam); ntr.resize(nparam); nac.resize(nparam);         // Initialises proposal and diagnostic information
	paramjumpxi.resize(nparam); ntrxi.resize(nparam); nacxi.resize(nparam);
	for(auto th = 0u; th < nparam; th++){
		paramval[th] = model.paramval[th];
		paramjump[th] = paramval[th]/2; if(paramjump[th] == 0) paramjump[th] = 0.1;
		ntr[th] = 0; nac[th] = 0;
		
		paramjumpxi[th] = paramval[th]/10; if(paramjumpxi[th] == 0) paramjumpxi[th] = 0.1;
		ntrxi[th] = 0; nacxi[th] = 0;
	}
	
	logbetajump = 0.01;
	sigmajump = 0.01;
		
	timeprop = 0;
	
	numaddrem = 20;
	ntr_addrem = 0; nac_addrem = 0;
}

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

/// Performs a MBP
unsigned int Chain::mbp()
{
	unsigned int j, jmax, c, v, sett, n, i, w, doev;
	double t, tmax, val, txi, tinf, al;
	FEV ev;

	timers.timembpinit -= clock();

	doev = model.dombpevents();

	for(c = 0; c < comp.size(); c++) N[c] = 0;
	N[0] = data.popsize;
		
	jmax = xp.size(); for(j = 0; j < jmax; j++) indevp[xp[j].ind].clear();
	//indevp.clear(); indevp.resize(data.popsize);	
	
	xp.clear();
	trevp.clear(); trevp.resize(details.nsettime);
	
	for(v = 0; v < data.narage; v++) dQmap[v] = 0;
	
	timers.timembpinit += clock();
	
	timers.timembp -= clock();
		
	t = 0; n = 0;
	for(sett = 0; sett < details.nsettime; sett++){
		if(details.mode == sim){
			cout  << "  Time: " << details.settime[t];
			for(c = 0; c < comp.size(); c++) cout << "  " << comp[c].name << ":"	<< N[c];
			cout << endl;	
		}
		
		phii = model.phii[sett]; phip = model.phip[sett];	
		betai = model.betai[sett]; betap = model.betap[sett];

		for(v = 0; v < data.narage; v++){
			val = Qmapi[sett][v] + dQmap[v];
			if(val < -tiny){ cout << val << "val\n"; emsgEC("Chain",1);}
			if(val < 0) val = 0;	
			Qmapp[sett][v] = val;
		}

		constructRtot(Qmapi[sett],Qmapp[sett]);
	
		tmax = details.settime[sett+1];
		do{
			if(n < xi.size()){ ev = indevi[xi[n].ind][xi[n].e]; txi = ev.t;} else{ ev.ind = UNSET; txi = tmax;}
			
			double v; v = ran();
			if(Rtot[0][0] <= 0) tinf = tmax; else tinf = t - log(v)/Rtot[0][0];
				
			if(txi >= tmax && tinf >= tmax){ t = tmax; break;}
			
			if(tinf < txi){  // A new infection appears
				t = tinf;
				c = nextinfection();
				addinfc(c,t);	
			}
			else{            // An event on initial sequence 
				t = txi;
				i = ev.ind;
			
				if(stat[i] == both_sus){
					w = data.ind[i].area*data.ndemocatpos + data.ind[i].dp;
		
					al = lamp[w]/lami[w];
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
	unsigned int v, dq, q, j, jmax, k, kmax, i, sett, a, nage, vv, loop, qt;
	double val, fac;
	FEV fev;

	for(v = 0; v < data.narage; v++) dQmap[v] = 0;

	nage = data.nage;
	for(sett = 0; sett < details.nsettime; sett++){
		for(v = 0; v < data.narage; v++){
			val = dQmap[v];
			if(check == 1){
				if(val < -tiny) emsgEC("Chain",2);
				if(val < Qmapi[sett][v]-tiny || val > Qmapi[sett][v]+tiny) emsgEC("Chain",3);
			}
			if(val < 0){ val = 0; dQmap[v] = 0;}	
			
			Qmapi[sett][v] = val;
		}
		
		jmax = trevi[sett].size();
		for(j = 0; j < jmax; j++){
			i = trevi[sett][j].ind;
			fev = indevi[i][trevi[sett][j].e];

			v = data.ind[i].area*data.nage+data.democatpos[data.ind[i].dp][0];
			dq = trans[fev.trans].DQ[fev.timep];
			if(dq != UNSET){
				for(loop = 0; loop < 2; loop++){
					q = model.DQ[dq].q[loop];
					if(q != UNSET){
						fac = model.DQ[dq].fac[loop];
						
						qt = data.Q[q].Qtenref;
						kmax = data.genQ.Qten[qt].ntof[v];
						auto& cref = data.genQ.Qten[qt].tof[v];
						auto& valref = data.genQ.Qten[qt].valf[v];
						if(nage == 1){
							for(k = 0; k < kmax; k++){
								dQmap[cref[k]*nage] += fac*valref[k][0];
							}
						}
						else{
							for(k = 0; k < kmax; k++){
								vv = cref[k]*nage;	
								for(a = 0; a < nage; a++){
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
	unsigned int c, cmax, wmin, wmax, dp, v, a, w, j, jmax;
	int l;
	double sum, dlam, faci, facp;
	
	timers.timembpconRtot -= clock();
		
	l = poptree.level-1;
	for(c = 0; c < data.narea; c++){
		wmin = c*data.ndemocatpos; wmax = wmin + data.ndemocatpos;
	
		faci = betai*model.areafaci[c];
		facp = betap*model.areafacp[c];
		
		sum = 0; dp = 0; v = c*data.nage; 
		for(w = wmin; w < wmax; w++){
			a = data.democatpos[dp][0];
			lami[w] = model.susi[dp]*(faci*Qmi[v+a] + phii);
			lamp[w] = model.susp[dp]*(facp*Qmp[v+a] + phip);
			dlam = nindbothlist[w]*(lamp[w] - lami[w]); if(dlam < 0) dlam = 0;
			
				if(std::isnan(dlam)){ cout <<   model.susi[dp] << " " << faci << " " << Qmi[v+a] << " " << phii << " g " <<  model.susp[dp] << " " << facp << " " << Qmp[v+a] << " " << phip << " g1\n"; emsgEC("Chain",400);}// zz
				
					if(std::isnan(lamp[w])){ cout <<   model.susi[dp] << " " << faci << " " << Qmi[v+a] << " " << phii << " g " <<  model.susp[dp] << " " << facp << " " << Qmp[v+a] << " " << phip << " g2\n"; emsgEC("Chain",401);}
			sum += dlam + nindponlylist[w]*lamp[w];
			dp++;
		}
		if(std::isnan(sum)) emsgEC("Chain",4);
		Rtot[l][c] = sum;
	}
	
	for(l = poptree.level-2; l >= 0; l--){                                 
		cmax = lev[l].node.size();
		for(c = 0; c < cmax; c++){
			jmax = lev[l].node[c].child.size();
			sum = 0; for(j = 0; j < jmax; j++) sum += Rtot[l+1][lev[l].node[c].child[j]];
			
			Rtot[l][c] = sum;
		}
	}
	
	timers.timembpconRtot += clock();
}
	
/// Performs a MBP on parameter 'th'
void Chain::proposal(unsigned int th, unsigned int samp, unsigned int burnin)  
{
	unsigned int j, jmax;
	double al, valst, Lp=0, Prp=0, dd;

	timeprop -= clock();
	timers.timembpprop -= clock();
	
	model.setup(paramval);
	model.copyi();
			
	valst = paramval[th];

	paramval[th] += normal(0,paramjump[th]);               // Makes a change to a parameter

	if(paramval[th] < model.param[th].min || paramval[th] > model.param[th].max) al = 0;
	else{
		if(model.setup(paramval) == 1) al = 0;
		else{
			model.copyp();
			if(mbp() == 1) al = 0;
			else{
				Lp = obsmodel.Lobs(trevp,indevp);
				Prp = model.prior();
				
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
		
		jmax = xi.size(); for(j = 0; j < jmax; j++) indevi[xi[j].ind].clear();
		jmax = xp.size(); for(j = 0; j < jmax; j++) indevi[xp[j].ind] = indevp[xp[j].ind];
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
		dd = Li - obsmodel.Lobs(trevi,indevi); if(sqrt(dd*dd) > tiny) emsgEC("Chain",5);
		dd = Pri - model.prior(); if(sqrt(dd*dd) > tiny) emsgEC("Chain",6);
	}
	
	timers.timembpprop += clock();
	timeprop += clock();
}

/// Sets up lists for use with MBPs
void Chain::setuplists()
{	
	unsigned int c, dp, j, jmax, w, i;

	indbothlist.clear(); indponlylist.clear(); indnotlist.clear();
	indbothlist.resize(data.nardp); indponlylist.resize(data.nardp); indnotlist.resize(data.nardp);
	nindbothlist.resize(data.nardp); nindponlylist.resize(data.nardp); nindnotlist.resize(data.nardp);
	stat.resize(data.popsize); indlistref.resize(data.popsize); 

	for(c = 0; c < data.narea; c++){
		for(dp = 0; dp < data.ndemocatpos; dp++){
			w = c*data.ndemocatpos + dp;

			jmax = data.area[c].ind[dp].size();
			for(j = 0; j < jmax; j++){
				i = data.area[c].ind[dp][j];
				stat[i] = both_sus;
				indlistref[i] = indbothlist[w].size();
				indbothlist[w].push_back(i);
			}
			nindbothlist[w] = jmax;
			nindponlylist[w] = 0;
			nindnotlist[w] = 0;
		}
	}	
}

/// Makes a change in the status of an individual
void Chain::changestat(unsigned int i, unsigned int st, unsigned int updateR)
{
	unsigned int c, w;
	int l, n;
	double dlam, dval=0;

	c = data.ind[i].area;
	w = c*data.ndemocatpos + data.ind[i].dp;
		
	if(updateR == 1){		
		dlam = nindbothlist[w]*(lamp[w] - lami[w]); if(dlam < 0) dlam = 0;
		dval = -(dlam + nindponlylist[w]*lamp[w]);
	}
	
	l = indlistref[i];   // Removes the exisiting entry
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
		dlam = nindbothlist[w]*(lamp[w] - lami[w]); if(dlam < 0) dlam = 0;
		dval += dlam + nindponlylist[w]*lamp[w];
		
		if(dval != 0){
			l = poptree.level-1;
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
	unsigned int w, j, jmax, i;
	
	for(w = 0; w < data.nardp; w++){
		jmax = indponlylist[w].size();
		for(j = 0; j < jmax; j++){
			i = indponlylist[w][j];
			stat[i] = both_sus;
			indlistref[i] = indbothlist[w].size();
			indbothlist[w].push_back(i);
			nindbothlist[w]++;
		}
		indponlylist[w].clear();
		nindponlylist[w] = 0;
		
		jmax = indnotlist[w].size();
		for(j = 0; j < jmax; j++){
			i = indnotlist[w][j];
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
	unsigned int j, jmax, k, kmax, i, tra, v, dq, q, vv, a, nage, loop, qt;
	double fac;

	timers.timembpQmap -= clock();
	
	nage = data.nage;

	jmax = trei.size();
	for(j = 0; j < jmax; j++){
		i = trei[j].ind; 
		tra = indevi[i][trei[j].e].trans;
		indmap[i][tra] = 1;
	}
	
	jmax = trep.size(); 
	for(j = 0; j < jmax; j++){
		i = trep[j].ind; 
		tra = indevp[i][trep[j].e].trans;	
		if(indmap[i][tra] == 0){
			v = data.ind[i].area*nage+data.democatpos[data.ind[i].dp][0];
			dq = trans[tra].DQ[indevp[i][trep[j].e].timep];

			if(dq != UNSET){
				for(loop = 0; loop < 2; loop++){
					q = model.DQ[dq].q[loop];
					if(q != UNSET){
						if(dQbuf[v][q] == 0){ dQbuflistv.push_back(v); dQbuflistq.push_back(q);}
						dQbuf[v][q] += model.DQ[dq].fac[loop];
					}
				}
			}
		}
		else indmap[i][tra] = 0;
	}	
	
	jmax = trei.size(); 
	for(j = 0; j < jmax; j++){
		i = trei[j].ind; 
		tra = indevi[i][trei[j].e].trans;
		if(indmap[i][tra] == 1){
			
			v = data.ind[i].area*nage+data.democatpos[data.ind[i].dp][0];
			dq = trans[tra].DQ[indevi[i][trei[j].e].timep];
			if(dq != UNSET){
				for(loop = 0; loop < 2; loop++){
					q = model.DQ[dq].q[loop];
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
		jmax = trep.size(); 
		for(j = 0; j < jmax; j++){
			tra = indevp[trep[j].ind][trep[j].e].trans;
		
			N[trans[tra].from]--;
			N[trans[tra].to]++;
		}
	}

	nage = data.nage;
	jmax = dQbuflistv.size();
	for(j = 0; j < jmax; j++){
		v = dQbuflistv[j]; q = dQbuflistq[j]; qt = data.Q[q].Qtenref;
		fac = dQbuf[v][q];
		if(fac < -vtiny || fac > vtiny){
			kmax = data.genQ.Qten[qt].ntof[v];
			
			auto& cref = data.genQ.Qten[qt].tof[v];
			auto& valref = data.genQ.Qten[qt].valf[v];
			if(nage == 1){
				for(k = 0; k < kmax; k++){
					dQmap[cref[k]] += fac*valref[k][0];
				}
			}
			else{
				for(k = 0; k < kmax; k++){
					vv = cref[k]*nage;	
					for(a = 0; a < nage; a++){
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
	unsigned int l, lmax, c, j, jmax;
	double z, sum, sumst[4];
	
	l = 0; c = 0;                              // We start at the top level l=0 and proceed to finer and finer scales
	lmax = poptree.level;
	while(l < lmax-1){
		jmax = lev[l].node[c].child.size();
		sum = 0;
		for(j = 0; j < jmax; j++){
			sum += Rtot[l+1][lev[l].node[c].child[j]];	
			sumst[j] = sum;
		}
		
		z = ran()*sum; j = 0; while(j < jmax && z > sumst[j]) j++;
		if(j == jmax) emsgEC("Chain",12);
		
		c = lev[l].node[c].child[j];
		l++;
	};
	
	return c;
}

/// Adds an exposed indivdual in area
void Chain::addinfc(unsigned int c, double t)
{
	unsigned int dp, dpmax, w, i, n;
	double sum, dlam, z;
	vector <double> sumst;
	
	dpmax = data.ndemocatpos;
	sumst.resize(dpmax);
	
	sum = 0;                           // First selects the demographic possibility
	for(dp = 0; dp < dpmax; dp++){
		w = c*dpmax + dp;
		dlam = nindbothlist[w]*(lamp[w] - lami[w]); if(dlam < 0) dlam = 0;
		sum += dlam + nindponlylist[w]*lamp[w];
		sumst[dp] = sum;
	}

	z = ran()*sum; dp = 0; while(dp < dpmax && z > sumst[dp]) dp++; 
	if(dp == dpmax) emsgEC("Chain",13);
	
	w = c*dpmax + dp;                  // Next select when individuals in the initial and proposed states are susceptible
	dlam = nindbothlist[w]*(lamp[w] - lami[w]); if(dlam < 0) dlam = 0;
	if(ran() < dlam/(dlam + nindponlylist[w]*lamp[w])){ // Both suscetible
		n = indbothlist[w].size(); if(n == 0) emsgEC("Chain",14);
		i = indbothlist[w][(unsigned int)(ran()*n)];
	}
	else{                                       // Only proposed state susceptible
		n = indponlylist[w].size(); if(n == 0) emsgEC("Chain",15);
		i = indponlylist[w][(unsigned int)(ran()*n)];
	}
	
	changestat(i,not_sus,1);
	
	model.simmodel(indevp[i],i,0,t);

	addindev(i,indevp[i],xp,trevp);
}

/// Used for checking the code is running correctly
void Chain::check(unsigned int /* num */, double t, unsigned int sett)
{
	unsigned int c, l, i, w, dp, a, wmin, wmax, v, j, e, emax, tra, timep;
	double dd, dlam, sum, tt, ttt;
	vector <double> lai,lap;

	for(j = 0; j < xp.size(); j++){ // Checks order
		i = xp[j].ind; e = xp[j].e;
		if(i >= indevp.size()) emsgEC("Chain",16);
		if(e >= indevp[i].size()) emsgEC("Chain",17);
		if(j < xp.size()-1){
			if(indevp[i][e].t > indevp[xp[j+1].ind][xp[j+1].e].t) emsgEC("Chain",18);
		}
	}

	for(i = 0; i < data.popsize; i++){
		emax = indevp[i].size();
		if(emax > 0){
			c = 0; tt = 0;
			for(e = 0; e < emax; e++){
				ttt = indevp[i][e].t; if(ttt <= tt){ model.oe("here",indevp[i]); emsgEC("Chain",19);}
				tra = indevp[i][e].trans;
				if(trans[tra].from != c) emsgEC("Chain",20);
				c = trans[tra].to; tt = ttt;
			}
			if(comp[c].trans.size() != 0) emsgEC("Chain",21);
			
			for(timep = 0; timep < model.ntimeperiod; timep++){
				tt = model.timeperiod[timep].tend;
				if(tt > indevp[i][0].t && tt < indevp[i][emax-1].t){
					for(e = 0; e < emax; e++) if(indevp[i][e].t == tt) break;
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
	
	for(i = 0; i < data.popsize; i++){    // Checks stat is correct
	  w = data.ind[i].area*data.ndemocatpos + data.ind[i].dp;
		
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
	
	l = poptree.level-1;
	for(c = 0; c < data.narea; c++){
		wmin = c*data.ndemocatpos; wmax = wmin + data.ndemocatpos;
	
		sum = 0; dp = 0; v = c*data.nage; 
		for(w = wmin; w < wmax; w++){
			a = data.democatpos[dp][0];
			dd = lami[w] - model.susi[dp]*(betai*model.areafaci[c]*Qmapi[sett][v+a] + phii); if(sqrt(dd*dd) > tiny) emsgEC("Chain",33);
			dd = lamp[w] - model.susp[dp]*(betap*model.areafacp[c]*Qmapp[sett][v+a] + phip); if(sqrt(dd*dd) > tiny) emsgEC("Chain",34);
	
			dlam = nindbothlist[w]*(lamp[w] - lami[w]); if(dlam < 0) dlam = 0;
			sum += dlam + nindponlylist[w]*lamp[w];
			dp++;
		}
		dd = Rtot[l][c] - sum; if(sqrt(dd*dd) > tiny){ emsgEC("Chain",35);}
	}
	
	for(i = 0; i < data.popsize; i++){
		for(tra = 0; tra < model.trans.size(); tra++){
			if(indmap[i][tra] != 0) emsgEC("Chain",36);
		}
	}
}

/// Checks that quantities used when adding and removing events are correctly updated
void Chain::check_addrem()
{
	unsigned int j, i, e, num, sett, se;
	vector < vector <int> > done;
	
	for(j = 0; j < xi.size(); j++){
		i = xi[j].ind; e = xi[j].e;
		if(i >= indevi.size()) emsgEC("Chain",37);
		if(e >= indevi[i].size()) emsgEC("Chain",38);
		if(j < xi.size()-1){
			if(indevi[i][e].t > indevi[xi[j+1].ind][xi[j+1].e].t) emsgEC("Chain",39);
		}
	}

	done.resize(indevi.size());
	num = 0;
	for(i = 0; i < indevi.size(); i++){
		if(indevi[i].size() > 0){
			num++;
			done[i].resize(indevi[i].size());
			for(e = 0; e < indevi[i].size(); e++) done[i][e] = 0;
		}
	}
	if(num != xi.size()) emsgEC("Chain",40);

	for(sett = 0; sett < details.nsettime; sett++){
		for(j = 0; j < trevi[sett].size(); j++){
			i = trevi[sett][j].ind; e = trevi[sett][j].e;
			if(e >= indevi[i].size()) emsgEC("Chain",41);
			se = (unsigned int)(details.nsettime*indevi[i][e].t/details.period); 
			if(se != sett) emsgEC("Chain",42);
			if(done[i][e] != 0) emsgEC("Chain",43);
			done[i][e] = 1;
		}
	}
	
	for(i = 0; i < indevi.size(); i++){
		for(e = 0; e < indevi[i].size(); e++){
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
	unsigned int c, dp, w, v, i, n, sett;
	double L, t, tt, tmax, beta, phi, fac;
	FEV ev;
	
	model.setup(paramval);
		
	for(c = 0; c < data.narea; c++){
		for(dp = 0; dp < data.ndemocatpos; dp++){
			w = c*data.ndemocatpos + dp;
			popw[w] = data.area[c].ind[dp].size();
		}
	}		
			
	L = 0;
	
	t = 0; n = 0;
	for(sett = 0; sett < details.nsettime; sett++){
		beta = model.beta[sett]; phi = model.phi[sett];
		tmax = details.settime[sett+1];
		
		for(c = 0; c < data.narea; c++){
			fac = beta*model.areafac[c];
			for(dp = 0; dp < data.ndemocatpos; dp++){
				w = c*data.ndemocatpos + dp;
				v = c*data.nage + data.democatpos[dp][0];
				lami[w] = model.sus[dp]*(fac*Qmap[sett][v] + phi);		
				if(lami[w] < 0) emsgEC("Chain",46);
				
				L -= lami[w]*popw[w]*(tmax-t);
			}
		}
		
		while(n < x.size()){
			i = x[n].ind;
			ev = indev[i][x[n].e];
			tt = ev.t;
			if(tt >= tmax) break;
	
			t = tt;
			
			c = data.ind[i].area;
			w = c*data.ndemocatpos + data.ind[i].dp;
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
	unsigned int v, sett;
	double val;
	
	for(v = 0; v < data.narage; v++) dQmap[v] = 0;
	
	for(sett = 0; sett < details.nsettime; sett++){
		for(v = 0; v < data.narage; v++){
			val = Qmapi[sett][v] + dQmap[v];
			if(val < -tiny){ emsgEC("Chain",48);}
			if(val < 0) val = 0;	
			Qmapp[sett][v] = val;
		}
		updatedQmap(trevi[sett],trevp[sett]);	
	} 
}

/// This incorporates standard proposals which adds and removes events as well as changes parameters
void Chain::standard_prop(unsigned int samp, unsigned int burnin)
{
	unsigned int loop, loopmax = 1;
	
	timers.timestandard -= clock();
	
	model.setup(paramval);
	
	timers.timembptemp -= clock();
	Levi = likelihood(Qmapi,xi,indevi);
	timers.timembptemp += clock();
	
	for(loop = 0; loop < loopmax; loop++){
		timers.timeparam -= clock();
		betaphi_prop(samp,burnin);
		area_prop(samp,burnin);
		model.compparam_prop(samp,burnin,xi,indevi,paramval,paramjumpxi,ntrxi,nacxi,Pri);
		if(model.regioneffect == 1) fixarea_prop(samp,burnin);
		timers.timeparam += clock();
			
		if(checkon == 1){ double dd = likelihood(Qmapi,xi,indevi) - Levi; if(dd*dd > tiny) emsgEC("Chain",49);}

		timers.timeaddrem -= clock();
		if(loop%2 == 0) addrem_prop(samp,burnin);
		timers.timeaddrem += clock();
		
		if(checkon == 1){ double dd = likelihood(Qmapi,xi,indevi) - Levi; if(dd*dd > tiny) emsgEC("Chain",50);}
	}
	
	if(checkon == 1){ double dd = Pri - model.prior(); if(sqrt(dd*dd) > tiny){ emsgEC("Chain",51);}}

	timers.timestandard += clock();
}

/// Makes proposal to beta and phi
void Chain::betaphi_prop(unsigned int samp, unsigned int burnin)
{	
	unsigned int c, dp, w, v, i, j, jmax, n, sett, loop, loopmax, th, pos;
	double t, tt, tmax, betasum, phisum, al, Levp, valst, fac;
	vector <unsigned int> map;
	FEV ev;
	vector <double> betafac, phifac;
	LCONT lcont;
	vector <LCONT> lcontlist;
	vector<	vector <LCONT> > lc;
	vector <unsigned int> parampos;
			
	timers.timebetaphiinit -= clock();
	
	map.resize(data.nardp); for(w = 0; w < data.nardp; w++) map[w] = 0;
	
	for(c = 0; c < data.narea; c++){
		for(dp = 0; dp < data.ndemocatpos; dp++){
			w = c*data.ndemocatpos + dp;
			popw[w] = data.area[c].ind[dp].size();
		}
	}		
			
	model.setup(paramval);
	betafac.resize(details.nsettime); phifac.resize(details.nsettime);
	
	t = 0; n = 0;
	for(sett = 0; sett < details.nsettime; sett++){
		//beta = model.beta[sett]; phi = model.phi[sett];
		tmax = details.settime[sett+1];

		lcontlist.clear();
		
		betasum = 0; phisum = 0;
		for(c = 0; c < data.narea; c++){
			fac = model.areafac[c];
			for(dp = 0; dp < data.ndemocatpos; dp++){
				w = c*data.ndemocatpos + dp;
				v = c*data.nage + data.democatpos[dp][0];	
				betasum -= fac*model.sus[dp]*Qmapi[sett][v]*popw[w]*(tmax-t);
				phisum -= model.sus[dp]*popw[w]*(tmax-t);
			}
		}
		
		while(n < xi.size()){
			i = xi[n].ind;
			ev = indevi[i][xi[n].e];
			tt = ev.t;
			if(tt >= tmax) break;
	
			t = tt;
			
			c = data.ind[i].area;
			fac = model.areafac[c];
			
			dp = data.ind[i].dp;
			w = c*data.ndemocatpos + dp;
			v = c*data.nage + data.democatpos[dp][0]; 
			
			if(map[w] == 0){
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
					
		for(j = 0; j < lcontlist.size(); j++){
			w = lcontlist[j].w;
			lcontlist[j].num = map[w];
			map[w] = 0;
		}
					
		lc.push_back(lcontlist);
		
		t = tmax;
	} 
	
	for(j = 0; j < model.betaspline.size(); j++) parampos.push_back(model.betaspline[j].param);
	for(j = 0; j < model.phispline.size(); j++) parampos.push_back(model.phispline[j].param);
		
	timers.timebetaphiinit += clock();
		
	loopmax = 12/parampos.size(); if(loopmax == 0) loopmax = 1;
	
	timers.timebetaphi -= clock();
	for(loop = 0; loop < loopmax; loop++){
		for(pos = 0; pos < parampos.size(); pos++){
			th = parampos[pos];

			if(model.param[th].min != model.param[th].max){
				valst = paramval[th];	
				paramval[th] += normal(0,paramjumpxi[th]);               // Makes a change to a parameter

				if(paramval[th] < model.param[th].min || paramval[th] > model.param[th].max){ al = 0; Levp = -large;}
				else{
					model.setup(paramval);

					Levp = 0; 
					for(sett = 0; sett < details.nsettime; sett++){
						double beta = model.beta[sett], phi = model.phi[sett];
						Levp += betafac[sett]*beta + phifac[sett]*phi;
						
						jmax = lc[sett].size();
				
						for(j = 0; j < jmax; j++){
							Levp += lc[sett][j].num*log(lc[sett][j].betafac*beta + lc[sett][j].phifac*phi);
						}
						if(std::isnan(Levp)) emsgEC("Chain",52);
					}
					al = exp(Levp-Levi);
				}
			
				ntrxi[th]++;
				if(ran() < al){
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
	unsigned int c, dp, w, v, i, n, sett, j, th, loop, loopmax, num;
	double L0, t, tt, tmax, beta, phi;
	FEV ev;
	
	vector <double> areasum, lamareafac, lamphifac;
	vector < vector <double> >	mult, add;

	timers.timecovarinit -= clock();
	
	model.setup(paramval);

	lamareafac.resize(data.nardp);
	lamphifac.resize(data.nardp);
	mult.resize(data.narea);
	add.resize(data.narea);
	
	areasum.resize(data.narea);
	for(c = 0; c < data.narea; c++){
		areasum[c] = 0;
		for(dp = 0; dp < data.ndemocatpos; dp++){
			w = c*data.ndemocatpos + dp;
			popw[w] = data.area[c].ind[dp].size();
		}
	}		
	
	L0 = 0;
	t = 0; n = 0;
	for(sett = 0; sett < details.nsettime; sett++){
		beta = model.beta[sett]; phi = model.phi[sett];
		tmax = details.settime[sett+1];
		
		for(c = 0; c < data.narea; c++){
			for(dp = 0; dp < data.ndemocatpos; dp++){
				w = c*data.ndemocatpos + dp;
				v = c*data.nage + data.democatpos[dp][0];
				lamareafac[w] = model.sus[dp]*beta*Qmapi[sett][v];
				lamphifac[w] = model.sus[dp]*phi;
				
				areasum[c] -= lamareafac[w]*popw[w]*(tmax-t);
				L0 -= lamphifac[w] *popw[w]*(tmax-t);
			}
		}
		
		while(n < xi.size()){
			i = xi[n].ind;
			ev = indevi[i][xi[n].e];
			tt = ev.t;
			if(tt >= tmax) break;
	
			t = tt;
			
			c = data.ind[i].area;
			w = c*data.ndemocatpos + data.ind[i].dp;
			
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
	
	num = model.covar_param.size(); if(model.regioneffect == 1) num += data.nregion;
	loopmax = 12/num; if(loopmax == 0) loopmax = 1;
	
	for(loop = 0; loop < loopmax; loop++){
		th = model.logbeta_param;
		if(model.param[th].min != model.param[th].max) area_prop2(samp,burnin,th,L0,areasum,mult,add);
		
		for(j = 0; j < model.covar_param.size(); j++){ 
			th = model.covar_param[j];
			if(model.param[th].min != model.param[th].max) area_prop2(samp,burnin,th,L0,areasum,mult,add);
		}
		
		if(model.regioneffect == 1){
			for(j = 0; j < data.nregion; j++){ 
				th = model.regioneff_param[j];
				if(model.param[th].min != model.param[th].max) area_prop2(samp,burnin,th,L0,areasum,mult,add);
			}
		
			th = model.sigma_param;
			if(model.param[th].min != model.param[th].max) area_prop2(samp,burnin,th,L0,areasum,mult,add);
		}
	}
	timers.timecovar += clock();
}

void Chain::area_prop2(unsigned int samp, unsigned int burnin, unsigned int th, double L0, vector <double> &areasum, vector < vector <double> >&mult, vector < vector <double> > &add)
{
	unsigned int c, k, kmax;
	double valst, al, Levp, fac, Prp;
	
	valst = paramval[th];	
	paramval[th] += normal(0,paramjumpxi[th]);               // Makes a change to a parameter

	if(paramval[th] < model.param[th].min || paramval[th] > model.param[th].max){ al = 0; Levp = -large; Prp = -large;}
	else{
		model.setup(paramval);

		Levp = L0;
		for(c = 0; c < data.narea; c++){
			fac = model.areafac[c];
			Levp += areasum[c]*fac;
			kmax = mult[c].size();
			for(k = 0; k < kmax; k++) Levp += log(mult[c][k]*fac + add[c][k]);
		}
		if(std::isnan(Levp)) emsgEC("Chain",53);
	
	  Prp = model.prior();
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
	unsigned int th, r, loop, loopmax=10;
	double al, Prp, sd, dlogbeta, valst;
	
	for(loop = 0; loop < loopmax; loop++){
		th = model.logbeta_param;
		dlogbeta = normal(0,logbetajump);
		paramval[th] += dlogbeta;
		for(r = 0; r < data.nregion; r++) paramval[model.regioneff_param[r]] -= dlogbeta;
	
		if(paramval[th] < model.param[th].min || paramval[th] > model.param[th].max){ al = 0; Prp = Pri;}
		else{
			sd = paramval[model.sigma_param]; Prp = 0; for(r = 0; r < data.nregion; r++) Prp += normalprob(paramval[model.regioneff_param[r]],0,sd*sd);
			al = exp(Prp-Pri);
		}

		if(ran() < al){
			Pri = Prp;
			if(samp < burnin){ if(samp < 50) logbetajump *= 1.05; else logbetajump *= 1.01;}
		}
		else{
			paramval[th] -= dlogbeta;
			for(r = 0; r < data.nregion; r++) paramval[model.regioneff_param[r]] += dlogbeta;
			if(samp < burnin){ if(samp < 50) logbetajump *= 0.975; else logbetajump *= 0.995;}
		}
		
		th = model.sigma_param;
		valst = paramval[th];
		paramval[th] += normal(0,sigmajump);
		if(paramval[th] < model.param[th].min || paramval[th] > model.param[th].max){ al = 0; Prp = Pri;}
		else{
			sd = paramval[model.sigma_param]; Prp = 0; for(r = 0; r < data.nregion; r++) Prp += normalprob(paramval[model.regioneff_param[r]],0,sd*sd);
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
	unsigned int i;
	EVREFT evreft;         
	vector <EVREFT> xt;

	for(i = 0; i < x.size(); i++){
		evreft.ind = x[i].ind; evreft.e = x[i].e;
		if(indev[x[i].ind].size() == 0) emsgEC("Chain",54);
		
		evreft.t = indev[x[i].ind][x[i].e].t;
		xt.push_back(evreft);	
	}
	sort(xt.begin(),xt.end(),compEVREFT);

	for(i = 0; i < x.size(); i++){
		x[i].ind = xt[i].ind;	x[i].e = xt[i].e;
	}
}

/// Adds and removes infectious individuals
void Chain::addrem_prop(unsigned int samp, unsigned int burnin)
{
	unsigned int j, jmax, k, dk, sett, w, c, i, l;
	double z, probif, probfi, Levp, Lp, t, dt, al, dd, sumtot;
	vector <int> kst;
	
	model.setup(paramval);
		
	if(checkon == 1){ dd = likelihood(Qmapi,xi,indevi) - Levi; if(dd*dd > tiny) emsgEC("Chain",55);}

	probif = 0; probfi = 0;
	
	trevp = trevi;     // Copies initial state into proposed state
	jmax = xp.size(); for(j = 0; j < jmax; j++) indevp[xp[j].ind].clear();
	xp = xi;
	jmax = xp.size(); for(j = 0; j < jmax; j++) indevp[xp[j].ind] = indevi[xp[j].ind];
		 
	for(j = 0; j < xp.size(); j++) changestat(xp[j].ind,not_sus,0);
	
	if(ran() < 0.5){  // Adds individuals
		timers.timembptemp2 -= clock();
		infsampler(Qmapi);
		timers.timembptemp2 += clock();
		
		for(j = 0; j < numaddrem; j++){
			if(xp.size() >= model.infmax){ resetlists(); return;}
			
			sumtot = lamsum[data.nsettardp-1]; if(sumtot == 0) emsgEC("Chain",56);
			z = ran()*sumtot;
			
			k = 0; dk = data.nsettardp/10; if(dk == 0) dk = 1;
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

			sett = k/data.nardp; w = k%data.nardp;
			
			if(nindbothlist[w] == 0){ resetlists(); return;}
			i = indbothlist[w][int(ran()*nindbothlist[w])];
			probif += log(1.0/nindbothlist[w]);
		
			changestat(i,not_sus,0);
			
			dt = details.settime[sett+1]-details.settime[sett];
			t = details.settime[sett] + ran()*dt;
			probif += log(1.0/dt);
			
			model.simmodel(indevp[i],i,0,t);
			addindev(i,indevp[i],xp,trevp);
			
			probfi += log(1.0/xp.size());
		}
		
		timers.timembptemp3 -= clock();
		sortx(xp,indevp);
		calcQmapp();
		timers.timembptemp3 += clock();
	}
	else{    // Removes individuals
		for(j = 0; j < numaddrem; j++){
			if(xp.size() == 0){ resetlists(); return;}
			
			l = int(ran()*xp.size());
			probif += log(1.0/xp.size());
			i = xp[l].ind;
			sett = (unsigned int)(details.nsettime*indevi[i][xp[l].e].t/details.period); 

			c = data.ind[i].area;
			w = c*data.ndemocatpos + data.ind[i].dp;
	
			dt = details.settime[sett+1]-details.settime[sett];

			probfi += log(1.0/dt);
			kst.push_back(sett*data.nardp + w);
			
			indevp[i].clear();
			
			xp[l] = xp[xp.size()-1];
			xp.pop_back();
			
			changestat(i,both_sus,0);
			
			probfi += log(1.0/nindbothlist[w]);
		}
		
		for(sett = 0; sett < details.nsettime; sett++){  // Removes events in trevp
			j = 0; jmax = trevp[sett].size();
			while(j < jmax){
				if(indevp[trevp[sett][j].ind].size() == 0){
					jmax--;
					trevp[sett][j] = trevp[sett][jmax];
					trevp[sett].pop_back();
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
		
		for(j = 0; j < numaddrem; j++) probfi += log(lam[kst[j]]/lamsum[data.nsettardp-1]);
	}
	
	timers.timembptemp4 -= clock();
	Levp = likelihood(Qmapp,xp,indevp);
	
	Lp = obsmodel.Lobs(trevp,indevp);
	timers.timembptemp4 += clock();
		
	al = exp(invT*(Lp-Li) + Levp-Levi + probfi - probif);
	if(checkon == 1) cout << al << " " << Li << " " << Lp << " " << Levi << " " << Levp << "al" << endl;		
	 
	ntr_addrem++;
	if(ran() < al){
		Levi = Levp;
		Li = Lp;
		trevi = trevp;
		Qmapi = Qmapp;
		
		jmax = xi.size(); for(j = 0; j < jmax; j++) indevi[xi[j].ind].clear();
		jmax = xp.size(); for(j = 0; j < jmax; j++) indevi[xp[j].ind] = indevp[xp[j].ind];
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
	unsigned int sett, tot, c, dp, v, w;
	double val, sum, beta, phi, fac;
		
	sum = 0;
	for(sett = 0; sett < details.nsettime; sett++){
		beta = model.beta[sett]; phi = model.phi[sett];
	
		for(c = 0; c < data.narea; c++){
			fac = beta*model.areafac[c];
			for(dp = 0; dp < data.ndemocatpos; dp++){
				w = c*data.ndemocatpos + dp;
				tot = sett*data.nardp + w;
				v = c*data.nage + data.democatpos[dp][0];
				
				val = nindbothlist[w]*model.sus[dp]*(fac*Qmap[sett][v] + phi);
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
	
	for(auto i = 0u; i < indev.size(); i++){
		for(auto e = 0u; e < indev[i].size(); e++){
			store.push_back(indev[i][e]);
		}
	}
	
	return store;
}

/// Initialises the chain based on a particle (used for abcmbp)
void Chain::initialise_from_particle(const Particle &part)

{
	auto jmax = xi.size(); for(auto j = 0u; j < jmax; j++) indevi[xi[j].ind].clear();   // Removes the exisiting initial sequence 
	xi.clear();
	
	trevi.clear(); trevi.resize(details.nsettime); 
	
	vector <int> indlist;
	for(auto e = 0u; e < part.ev.size(); e++){
		int i = part.ev[e].ind;
		if(indevi[i].size() == 0) indlist.push_back(i);
		indevi[i].push_back(part.ev[e]);
	}	
	
	for(int i : indlist) addindev(i,indevi[i],xi,trevi);
	
	sortx(xi,indevi);
	
	setQmapi(0);

	EF = part.EF;
	Pri = model.prior();
	
	if(EF != obsmodel.Lobs(trevi,indevi)) emsg("Observation does not agree");
	
	paramval = part.paramval;
}

/// Generates a particle (used for abcmbp)
void Chain::generate_particle(Particle &part) const
{
	part.EF = EF;
	part.paramval = paramval;
	part.ev = event_compress(indevi);
}

int Chain::abcmbp_proposal(const vector <double> param_propose, double EFcut)  
{
	model.setup(paramval);
	model.copyi();
			
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
			model.copyp();
			if(mbp() == 1) al = 0;
			else{
				EFp = obsmodel.Lobs(trevp,indevp);
				if(EFp >= EFcut) al = 0;
				else{
					Prp = model.prior();
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
		
		auto jmax = xi.size(); for(auto j = 0u; j < jmax; j++) indevi[xi[j].ind].clear();
		jmax = xp.size(); for(auto j = 0u; j < jmax; j++) indevi[xp[j].ind] = indevp[xp[j].ind];
		//indevi = indevp;
		
		xi = xp;
		return 1;
	}
	else{
		paramval = valst;
	}
	return 0;
}

