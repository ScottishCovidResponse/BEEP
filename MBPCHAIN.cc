// This file contains all the functions for running a chain under the MBP algorithm

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
#include "MBPCHAIN.hh"
#include "output.hh"
#include "pack.hh"
#include "obsmodel.hh"

MBPCHAIN::MBPCHAIN(DATA &data, MODEL &model, POPTREE &poptree, double invTstart, unsigned int chstart) : data(data), model(model), poptree(poptree), trans(model.trans), comp(model.comp), lev(poptree.lev)
{
	unsigned int th, nparam, v, q, j, sett, i, tra, loop, loopmax=1000;
	int l;
	vector <double> paramvalinit;

	invT = invTstart;	invTtrue = invTstart;
	ch = chstart;

	xi.clear();
	trevi.clear(); trevi.resize(data.nsettime); 
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
	
	Qmapi.resize(data.nsettime); Qmapp.resize(data.nsettime); 
	for(sett = 0; sett < data.nsettime; sett++){
		Qmapi[sett].resize(data.narage); for(v = 0; v < data.narage; v++) Qmapi[sett][v] = 0;
		Qmapp[sett].resize(data.narage); 
	}

	lami.resize(data.nardp); lamp.resize(data.nardp);
	Rtot.resize(poptree.level); for(l = 0; l < (int)poptree.level; l++) Rtot[l].resize(lev[l].node.size()); 
	N.resize(comp.size()); 
	
	loop = 0;
	do{
		//if(data.mode == MODE_INF) cout << "Initialisation try: " << loop << endl;
		do{	model.priorsamp(); }while(model.setup(model.paramval) == 1);             // Randomly samples parameters from the prior	

		nparam = model.param.size();                   
		paramval.resize(nparam); for(th = 0; th < nparam; th++) paramval[th] = model.paramval[th];
		paramvalinit = paramval;
		
		 // Sets the initial state to zero force of infection
		for(j = 0; j < model.betaspline.size(); j++) paramvalinit[model.betaspline[j].param] = 0; 
		for(j = 0; j < model.phispline.size(); j++) paramvalinit[model.phispline[j].param] = 0;
			
		model.setup(paramvalinit);                                       // To generate initial state mbp is used to simulate
		model.copyi();
		model.setup(paramval);
		model.copyp();
	
		if(mbp() == 0) break;

		loop++;
	}while(loop < loopmax);                          // Checks not too many infected (based on prior)

	if(loop == loopmax) emsg("Cannot find initial state with number of events under infmax");
	
	trevi = trevp;
	Qmapi = Qmapp;	
	indevi = indevp;
	xi = xp;
	
	if(data.mode != MODE_INF) return;
	
	Li = Lobs(data,model,poptree,trevi,indevi);
	//emsg("P");
	Pri = model.prior();

	setQmapi(1);

	paramjump.resize(nparam); ntr.resize(nparam); nac.resize(nparam);         // Initialises proposal and diagnostic information
	paramjumpxi.resize(nparam); ntrxi.resize(nparam); nacxi.resize(nparam);
	for(th = 0; th < nparam; th++){
		paramval[th] = model.paramval[th];
		paramjump[th] = paramval[th]/2; if(paramjump[th] == 0) paramjump[th] = 0.1;
		ntr[th] = 0; nac[th] = 0;
		
		paramjumpxi[th] = paramval[th]/10; if(paramjumpxi[th] == 0) paramjumpxi[th] = 0.1;
		ntrxi[th] = 0; nacxi[th] = 0;
	}
	timeprop = 0;
	
	numaddrem = 20;
	ntr_addrem = 0; nac_addrem = 0;
	
	popw.resize(data.nardp);                                        // Used for event based changes
	lam.resize(data.nsettardp); lamsum.resize(data.nsettardp);
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
unsigned int MBPCHAIN::mbp()
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
	trevp.clear(); trevp.resize(data.nsettime);
	
	for(v = 0; v < data.narage; v++) dQmap[v] = 0;
	
	timers.timembpinit += clock();
	
	timers.timembp -= clock();
		
	t = 0; n = 0;
	for(sett = 0; sett < data.nsettime; sett++){
		if(data.mode == MODE_SIM){
			cout  << "  Time: " << data.settime[t];
			for(c = 0; c < comp.size(); c++) cout << "  " << comp[c].name << ":"	<< N[c];
			cout << endl;	
		}
		
		phii = model.phii[sett]; phip = model.phip[sett];	
		betai = model.betai[sett]; betap = model.betap[sett];

		for(v = 0; v < data.narage; v++){
			val = Qmapi[sett][v] + dQmap[v];
			if(val < -tiny) emsg("MBPchain: EC31");
			if(val < 0) val = 0;	
			Qmapp[sett][v] = val;
		}

		constructRtot(Qmapi[sett],Qmapp[sett]);
	
		tmax = data.settime[sett+1];
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
			
				if(stat[i] == BOTH){
					w = data.ind[i].area*data.ndemocatpos + data.ind[i].dp;
		
					al = lamp[w]/lami[w];
					if(ran() < al){                                    // Keeps the infection event
						changestat(i,NOT,1);
						
						if(doev == 1) model.mbpmodel(indevi[i],indevp[i]);
						else indevp[i] = indevi[i];
						
						addindev(i,indevp[i],xp,trevp);
					}
					else changestat(i,PONLY,1);      // Does not keep the infection event
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
void MBPCHAIN::addindev(unsigned int i, vector <FEV> &indev, vector <EVREF> &x, vector <vector <EVREF> > &trev)
{
	unsigned int e, emax, se;
	EVREF evref;
	
	emax = indev.size();
	if(emax == 0) return;
	
	evref.ind = i; evref.e = 0;
	x.push_back(evref);
	for(e = 0; e < indev.size(); e++){
		evref.e = e;
		se = (unsigned int)(data.nsettime*indev[e].t/data.period); 
		if(se < data.nsettime) trev[se].push_back(evref);
	}
}
		
/// Based on the the event sequence in xi, this sets Qmapi
void MBPCHAIN::setQmapi(unsigned int check)
{
	unsigned int v, dq, q, j, jmax, k, kmax, i, sett, a, nage, vv, loop, qt;
	double val, fac;
	unsigned short *cref;
	float **valref;	
	FEV fev;

	for(v = 0; v < data.narage; v++) dQmap[v] = 0;

	nage = data.nage;
	for(sett = 0; sett < data.nsettime; sett++){
		for(v = 0; v < data.narage; v++){
			val = dQmap[v];
			if(check == 1){
				if(val < -tiny) emsg("MBPchain: EC17a");
				if(val < Qmapi[sett][v]-tiny || val > Qmapi[sett][v]+tiny) emsg("MBPchain: EC17b");
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
						cref = data.genQ.Qten[qt].tof[v];
						valref = data.genQ.Qten[qt].valf[v];
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
void MBPCHAIN::constructRtot(vector <double> &Qmi, vector <double> &Qmp)
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
			sum += dlam + nindponlylist[w]*lamp[w];
			dp++;
		}
		if(std::isnan(sum)) emsg("nan here");
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
	
/// Performs an MBP on parameter 'th'
void MBPCHAIN::proposal(unsigned int th, unsigned int samp, unsigned int burnin)  
{
	unsigned int j, jmax;
	double al, valst, Lp=0, Prp=0, dd;

	timeprop -= clock();
	timers.timembpprop -= clock();
	
	model.setup(paramval);
	model.copyi();
			
	valst = paramval[th];

	//cout << model.param[th].name << "up\n";
	paramval[th] += normal(0,paramjump[th]);               // Makes a change to a parameter

	if(paramval[th] < model.param[th].min || paramval[th] > model.param[th].max) al = 0;
	else{
		if(model.setup(paramval) == 1) al = 0;
		else{
			model.copyp();
			if(mbp() == 1) al = 0;
			else{
				Lp = Lobs(data,model,poptree,trevp,indevp);
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
		dd = Li - Lobs(data,model,poptree,trevi,indevi); if(sqrt(dd*dd) > tiny) emsg("MBPchain: EC34");
		dd = Pri - model.prior(); if(sqrt(dd*dd) > tiny) emsg("MBPchain: EC35");
	}
	
	timers.timembpprop += clock();
	timeprop += clock();
}

/// Sets up lists for use with MBPs
void MBPCHAIN::setuplists()
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
				stat[i] = BOTH;
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
void MBPCHAIN::changestat(unsigned int i, unsigned int st, unsigned int updateR)
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
  case BOTH:
		if(indbothlist[w][l] != i) emsg("MBPchain: EC43b");
		n = indbothlist[w].size();
		if(l < n-1){
			indbothlist[w][l] = indbothlist[w][n-1];
			indlistref[indbothlist[w][l]] = l;
		}
		indbothlist[w].pop_back();
		nindbothlist[w]--;
		break;
		
	case PONLY:
		if(indponlylist[w][l] != i) emsg("MBPchain: EC44");
		n = indponlylist[w].size();
		if(l < n-1){
			indponlylist[w][l] = indponlylist[w][n-1];
			indlistref[indponlylist[w][l]] = l;
		}
		indponlylist[w].pop_back();
		nindponlylist[w]--;
		break;
		
	case NOT:
		if(indnotlist[w][l] != i) emsg("MBPchain: EC44");
		n = indnotlist[w].size();
		if(l < n-1){
			indnotlist[w][l] = indnotlist[w][n-1];
			indlistref[indnotlist[w][l]] = l;
		}
		indnotlist[w].pop_back();
		nindnotlist[w]--;
		break;
	
	default: emsg("MBPChain EC46"); break;
	}

	stat[i] = st;
	switch(stat[i]){
	case PONLY:
		indlistref[i] = indponlylist[w].size();
		indponlylist[w].push_back(i);
		nindponlylist[w]++;
		break;
		
	case NOT:
		indlistref[i] = indnotlist[w].size();
		indnotlist[w].push_back(i);
		nindnotlist[w]++;
		break;
	
	case BOTH:
		indlistref[i] = indbothlist[w].size();
		indbothlist[w].push_back(i);
		nindbothlist[w]++;
		break;
		
	default: emsg("MBPChain EC47"); break;
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
void MBPCHAIN::resetlists()
{
	unsigned int w, j, jmax, i;
	
	for(w = 0; w < data.nardp; w++){
		jmax = indponlylist[w].size();
		for(j = 0; j < jmax; j++){
			i = indponlylist[w][j];
			stat[i] = BOTH;
			indlistref[i] = indbothlist[w].size();
			indbothlist[w].push_back(i);
			nindbothlist[w]++;
		}
		indponlylist[w].clear();
		nindponlylist[w] = 0;
		
		jmax = indnotlist[w].size();
		for(j = 0; j < jmax; j++){
			i = indnotlist[w][j];
			stat[i] = BOTH;
			indlistref[i] = indbothlist[w].size();
			indbothlist[w].push_back(i);
			nindbothlist[w]++;
		}
		indnotlist[w].clear();
		nindnotlist[w] = 0;
	}
}

/// Updates dQmap based on events which occur in timestep sett in the initial and proposed states
void MBPCHAIN::updatedQmap(vector <EVREF> &trei, vector <EVREF> &trep)
{
	unsigned int j, jmax, k, kmax, i, tra, v, dq, q, vv, a, nage, loop, qt;
	double fac;
	unsigned short *cref;
	float **valref;	

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
			v = data.ind[i].area*data.nage+data.democatpos[data.ind[i].dp][0];
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
			
			v = data.ind[i].area*data.nage+data.democatpos[data.ind[i].dp][0];
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
	
	if(data.mode == MODE_SIM){
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
			
			cref = data.genQ.Qten[qt].tof[v];
			valref = data.genQ.Qten[qt].valf[v];
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
unsigned int MBPCHAIN::nextinfection()
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
		if(j == jmax) emsg("MBPchain: EC15");
		
		c = lev[l].node[c].child[j];
		l++;
	};
	
	return c;
}

/// Adds an exposed indivdual in area
void MBPCHAIN::addinfc(unsigned int c, double t)
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
	if(dp == dpmax) emsg("MBPchain: EC56");
	
	w = c*dpmax + dp;                  // Next select when individuals in the initial and proposed states are susceptible
	dlam = nindbothlist[w]*(lamp[w] - lami[w]); if(dlam < 0) dlam = 0;
	if(ran() < dlam/(dlam + nindponlylist[w]*lamp[w])){ // Both suscetible
		n = indbothlist[w].size(); if(n == 0) emsg("MBPchain: EC43");
		i = indbothlist[w][(unsigned int)(ran()*n)];
	}
	else{                                       // Only proposed state susceptible
		n = indponlylist[w].size(); if(n == 0) emsg("MBPchain: EC44b");
		i = indponlylist[w][(unsigned int)(ran()*n)];
	}
	
	changestat(i,NOT,1);
	
	model.simmodel(indevp[i],i,0,t);

	addindev(i,indevp[i],xp,trevp);
}

/// Used for checking the code is running correctly
void MBPCHAIN::check(unsigned int /* num */, double t, unsigned int sett)
{
	unsigned int c, l, i, w, dp, a, wmin, wmax, v, j, e, emax, tra, timep;
	double dd, dlam, sum, tt, ttt;
	vector <double> lai,lap;

	for(j = 0; j < xp.size(); j++){ // Checks order
		i = xp[j].ind; e = xp[j].e;
		// Unsigned quantities always >= 0
		if( /* i < 0 || */ i >= indevp.size()) emsg("MBPchain: EC57");
		if( /* e < 0 || */ e >= indevp[i].size()) emsg("MBPchain: EC58");
		if(j < xp.size()-1){
			if(indevp[i][e].t > indevp[xp[j+1].ind][xp[j+1].e].t) emsg("MBPchain: EC59");
		}
	}

	for(i = 0; i < data.popsize; i++){
		emax = indevp[i].size();
		if(emax > 0){
			c = 0; tt = 0;
			for(e = 0; e < emax; e++){
				ttt = indevp[i][e].t; if(ttt <= tt){ model.oe("here",indevp[i]); emsg("MBPchain: EC70");}
				tra = indevp[i][e].trans;
				if(trans[tra].from != c) emsg("MBPchain: EC71");
				c = trans[tra].to; tt = ttt;
			}
			if(comp[c].trans.size() != 0) emsg("MBPchain: EC72");
			
			for(timep = 0; timep < model.ntimeperiod; timep++){
				tt = model.timeperiod[timep].tend;
				if(tt > indevp[i][0].t && tt < indevp[i][emax-1].t){
					for(e = 0; e < emax; e++) if(indevp[i][e].t == tt) break;
					if(timep <  model.ntimeperiod-1){
						if(e == emax) emsg("MBPchain: EC73");
					}
					else{
						if(e != emax) emsg("MBPchain: EC73b");
					}
				}
			}
		}
	}
	
	for(i = 0; i < data.popsize; i++){    // Checks stat is correct
	  w = data.ind[i].area*data.ndemocatpos + data.ind[i].dp;
		
		if(nindbothlist[w] != indbothlist[w].size()) emsg("MBPchain: EC22");
		if(nindponlylist[w] != indponlylist[w].size()) emsg("MBPchain: EC22");
		if(nindnotlist[w] != indnotlist[w].size()) emsg("MBPchain: EC22");
		
		if((indevi[i].size() == 0 || t < indevi[i][0].t) && indevp[i].size() == 0){
			if(stat[i] != BOTH) emsg("MBPchain: EC22");
			if(indbothlist[w][indlistref[i]] != i) emsg("MBPchain: EC23");
		}
		else{
			if((indevi[i].size() != 0 && t >= indevi[i][0].t) && indevp[i].size() == 0){
				if(stat[i] != PONLY) emsg("MBPchain: EC24");
				if(indponlylist[w][indlistref[i]] != i) emsg("MBPchain: EC25");
			}
			else{
				if(stat[i] != NOT) emsg("MBPchain: EC26");
				if(indnotlist[w][indlistref[i]] != i) emsg("MBPchain: EC27");
			}
		}
	}
	
	l = poptree.level-1;
	for(c = 0; c < data.narea; c++){
		wmin = c*data.ndemocatpos; wmax = wmin + data.ndemocatpos;
	
		sum = 0; dp = 0; v = c*data.nage; 
		for(w = wmin; w < wmax; w++){
			a = data.democatpos[dp][0];
			dd = lami[w] - model.susi[dp]*(betai*model.areafaci[c]*Qmapi[sett][v+a] + phii); if(sqrt(dd*dd) > tiny) emsg("MBPchain: EC67");
			dd = lamp[w] - model.susp[dp]*(betap*model.areafacp[c]*Qmapp[sett][v+a] + phip); if(sqrt(dd*dd) > tiny) emsg("MBPchain: EC68");
	
			dlam = nindbothlist[w]*(lamp[w] - lami[w]); if(dlam < 0) dlam = 0;
			sum += dlam + nindponlylist[w]*lamp[w];
			dp++;
		}
		dd = Rtot[l][c] - sum; if(sqrt(dd*dd) > tiny){ emsg("MBPchain: EC69");}
	}
	
	for(i = 0; i < data.popsize; i++){
		for(tra = 0; tra < model.trans.size(); tra++){
			if(indmap[i][tra] != 0) emsg("MBPchain: EC65");
		}
	}
}

/// Checks that quantities used when adding and removing events are correctly updated
void MBPCHAIN::check_addrem()
{
	unsigned int j, i, e, num, sett, se;
	vector < vector <int> > done;
	
	for(j = 0; j < xi.size(); j++){
		i = xi[j].ind; e = xi[j].e;
		if(i >= indevi.size()) emsg("MBPchain: EC57");
		if(e >= indevi[i].size()) emsg("MBPchain: EC58");
		if(j < xi.size()-1){
			if(indevi[i][e].t > indevi[xi[j+1].ind][xi[j+1].e].t) emsg("MBPchain: EC59");
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
	if(num != xi.size()) emsg("MBPchain: EC60");

	for(sett = 0; sett < data.nsettime; sett++){
		for(j = 0; j < trevi[sett].size(); j++){
			i = trevi[sett][j].ind; e = trevi[sett][j].e;
			if(e >= indevi[i].size()) emsg("MBPchain: EC61");
			se = (unsigned int)(data.nsettime*indevi[i][e].t/data.period); 
			if(se != sett) emsg("MBPchain: EC62");
			if(done[i][e] != 0) emsg("MBPchain: EC62");
			done[i][e] = 1;
		}
	}
	
	for(i = 0; i < indevi.size(); i++){
		for(e = 0; e < indevi[i].size(); e++){
			if(indevi[i][e].t < data.period){
				if(done[i][e] != 1) emsg("MBPchain: EC63");
			}
			else{
				if(done[i][e] != 0) emsg("MBPchain: EC64");
			}
		}
	}
	
	setQmapi(1);
}

/// Calculates the likelihood in the initial state
double MBPCHAIN::likelihood(vector < vector<double> > &Qmap, vector <EVREF> &x, vector <vector<FEV> > &indev)
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
	for(sett = 0; sett < data.nsettime; sett++){
		beta = model.beta[sett]; phi = model.phi[sett];
		tmax = data.settime[sett+1];
		
		for(c = 0; c < data.narea; c++){
			fac = beta*model.areafac[c];
			for(dp = 0; dp < data.ndemocatpos; dp++){
				w = c*data.ndemocatpos + dp;
				v = c*data.nage + data.democatpos[dp][0];
				lami[w] = model.sus[dp]*(fac*Qmap[sett][v] + phi);		
				if(lami[w] < 0) emsg("n");
				
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
			if(std::isnan(L)) emsg("MBPchain: EC87");
			popw[w]--;
			n++;
			
			L += lami[w]*(tmax-t);
		}
		t = tmax;
	} 
	
	return L;
}

/// Calculates Qmapp based on the initial and final sequences
void MBPCHAIN::calcQmapp()
{	
	unsigned int v, sett;
	double val;
	
	for(v = 0; v < data.narage; v++) dQmap[v] = 0;
	
	for(sett = 0; sett < data.nsettime; sett++){
		for(v = 0; v < data.narage; v++){
			val = Qmapi[sett][v] + dQmap[v];
			if(val < -tiny){ emsg("MBPchain: EC31b");}
			if(val < 0) val = 0;	
			Qmapp[sett][v] = val;
		}
		updatedQmap(trevi[sett],trevp[sett]);	
	} 
}

/// This incorporates standard proposals which adds and removes events as well as changes parameters
void MBPCHAIN::standard_prop(unsigned int samp, unsigned int burnin)
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
		covar_prop(samp,burnin);
		model.compparam_prop(samp,burnin,xi,indevi,paramval,paramjumpxi,ntrxi,nacxi,Pri);
		timers.timeparam += clock();
			
		if(checkon == 1){ double dd = likelihood(Qmapi,xi,indevi) - Levi; if(dd*dd > tiny) emsg("MBPchain: EC24b");}

		timers.timeaddrem -= clock();
		if(loop%2 == 0) addrem_prop(samp,burnin);
		timers.timeaddrem += clock();
		
		if(checkon == 1){ double dd = likelihood(Qmapi,xi,indevi) - Levi; if(dd*dd > tiny) emsg("MBPchain: EC24c");}
	}
	
	if(checkon == 1){ double dd = Pri - model.prior(); if(sqrt(dd*dd) > tiny){ emsg("MBPchain: EC24d");}}
		
	timers.timestandard += clock();
}

/// Makes proposal to beta and phi
void MBPCHAIN::betaphi_prop(unsigned int samp, unsigned int burnin)
{	
	unsigned int c, dp, w, v, i, j, jmax, n, sett, loop, loopmax, th, pos;
	double t, tt, tmax, betasum, phisum, beta, phi, al, Levp, valst, fac;
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
	betafac.resize(data.nsettime); phifac.resize(data.nsettime);
	
	t = 0; n = 0;
	for(sett = 0; sett < data.nsettime; sett++){
		beta = model.beta[sett]; phi = model.phi[sett];
		tmax = data.settime[sett+1];

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
				lcont.w = w; lcont.betafac = fac*model.sus[dp]*Qmapi[sett][v]; lcont.phifac = model.sus[dp];
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
					for(sett = 0; sett < data.nsettime; sett++){
						beta = model.beta[sett]; phi = model.phi[sett];
						Levp += betafac[sett]*beta + phifac[sett]*phi;
						
						jmax = lc[sett].size();
				
						for(j = 0; j < jmax; j++){
							Levp += lc[sett][j].num*log(lc[sett][j].betafac*beta + lc[sett][j].phifac*phi);
						}
						if(std::isnan(Levp)) emsg("MBPchain: EC77b");
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

/// Makes proposal to change covariates
void MBPCHAIN::covar_prop(unsigned int samp, unsigned int burnin)
{	
	unsigned int c, dp, w, v, i, n, sett, k, kmax, j, th, loop, loopmax;
	double L0, t, tt, tmax, beta, phi, fac, valst, al, Levp;
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
	for(sett = 0; sett < data.nsettime; sett++){
		beta = model.beta[sett]; phi = model.phi[sett];
		tmax = data.settime[sett+1];
		
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
	loopmax = 12/model.covar_param.size(); if(loopmax == 0) loopmax = 1;
	for(loop = 0; loop < loopmax; loop++){
		for(j = 0; j < model.covar_param.size(); j++){ 
			th = model.covar_param[j];
		
			if(model.param[th].min != model.param[th].max){
				valst = paramval[th];	
				paramval[th] += normal(0,paramjumpxi[th]);               // Makes a change to a parameter

				if(paramval[th] < model.param[th].min || paramval[th] > model.param[th].max){ al = 0; Levp = -large;}
				else{
					model.setup(paramval);

					Levp = L0;
					for(c = 0; c < data.narea; c++){
						fac = model.areafac[c];
						Levp += areasum[c]*fac;
						kmax = mult[c].size();
						for(k = 0; k < kmax; k++) Levp += log(mult[c][k]*fac + add[c][k]);
					}
					if(std::isnan(Levp)) emsg("MBPchain: EC77b");
	
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
	timers.timecovar += clock();
}

/// Time orders x
void MBPCHAIN::sortx(vector <EVREF> &x, vector <vector <FEV> > &indev)
{
	unsigned int i;
	EVREFT evreft;         
	vector <EVREFT> xt;

	for(i = 0; i < x.size(); i++){
		evreft.ind = x[i].ind; evreft.e = x[i].e;
		if(indev[x[i].ind].size() == 0) emsg("MBPchain: EC87");
		
		evreft.t = indev[x[i].ind][x[i].e].t;
		xt.push_back(evreft);	
	}
	sort(xt.begin(),xt.end(),compEVREFT);

	for(i = 0; i < x.size(); i++){
		x[i].ind = xt[i].ind;	x[i].e = xt[i].e;
	}
}

/// Adds and removes infectious individuals
void MBPCHAIN::addrem_prop(unsigned int samp, unsigned int burnin)
{
	unsigned int j, jmax, k, dk, sett, w, c, i, l;
	double z, probif, probfi, Levp, Lp, t, dt, al, dd, sumtot;
	vector <int> kst;
	
	model.setup(paramval);
		
	if(checkon == 1){ dd = likelihood(Qmapi,xi,indevi) - Levi; if(dd*dd > tiny) emsg("MBPchain: EC24");}

	probif = 0; probfi = 0;
	
	trevp = trevi;     // Copies initial state into proposed state
	jmax = xp.size(); for(j = 0; j < jmax; j++) indevp[xp[j].ind].clear();
	xp = xi;
	jmax = xp.size(); for(j = 0; j < jmax; j++) indevp[xp[j].ind] = indevi[xp[j].ind];
		 
	for(j = 0; j < xp.size(); j++) changestat(xp[j].ind,NOT,0);
	
	if(ran() < 0.5){  // Adds individuals
		timers.timembptemp2 -= clock();
		infsampler(Qmapi);
		timers.timembptemp2 += clock();
		
		for(j = 0; j < numaddrem; j++){
			if(xp.size() >= model.infmax){ resetlists(); return;}
			
			sumtot = lamsum[data.nsettardp-1]; if(sumtot == 0) emsg("MBPchain: EC32a");
			z = ran()*sumtot;
			
			k = 0; dk = data.nsettardp/10; if(dk == 0) dk = 1;
			do{
				while(k < data.nsettardp && z > lamsum[k]) k += dk;
				if(dk == 1) break;
				if(k >= dk) k -= dk;
				dk /= 10; if(dk == 0) dk = 1;
			}while(1 == 1);
			if(k >= data.nsettardp) emsg("MBPchain: EC32b");
		
			if(k > 0){
				if(k >= lamsum.size()) emsg("MBPchain: EC32c");
				if(!(z < lamsum[k] && z > lamsum[k-1])) emsg("MBPchain: EC33");
			}
			
			probif += log(lam[k]/sumtot);

			sett = k/data.nardp; w = k%data.nardp;
			
			if(nindbothlist[w] == 0){ resetlists(); return;}
			i = indbothlist[w][int(ran()*nindbothlist[w])];
			probif += log(1.0/nindbothlist[w]);
		
			changestat(i,NOT,0);
			
			dt = data.settime[sett+1]-data.settime[sett];
			t = data.settime[sett] + ran()*dt;
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
			sett = (unsigned int)(data.nsettime*indevi[i][xp[l].e].t/data.period); 

			c = data.ind[i].area;
			w = c*data.ndemocatpos + data.ind[i].dp;
	
			dt = data.settime[sett+1]-data.settime[sett];

			probfi += log(1.0/dt);
			kst.push_back(sett*data.nardp + w);
			
			indevp[i].clear();
			
			xp[l] = xp[xp.size()-1];
			xp.pop_back();
			
			changestat(i,BOTH,0);
			
			probfi += log(1.0/nindbothlist[w]);
		}
		
		for(sett = 0; sett < data.nsettime; sett++){  // Removes events in trevp
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
	
	Lp = Lobs(data,model,poptree,trevp,indevp);
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
void MBPCHAIN::infsampler(vector< vector<double> > &Qmap)
{
	unsigned int sett, tot, c, dp, v, w;
	double val, sum, beta, phi, fac;
		
	sum = 0;
	for(sett = 0; sett < data.nsettime; sett++){
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
