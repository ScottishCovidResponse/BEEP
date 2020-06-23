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
#include "PART.hh"
#include "MBPCHAIN.hh"
#include "output.hh"
#include "pack.hh"
#include "obsmodel.hh"

MBPCHAIN::MBPCHAIN(DATA &data, MODEL &model, POPTREE &poptree) : data(data), model(model), poptree(poptree), trans(model.trans), comp(model.comp), lev(poptree.lev)
{
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

/// Initialises an MCMC chain
void MBPCHAIN::init(DATA &data, MODEL &model, POPTREE &poptree, double invTstart, vector < vector <FEV> > &indev, unsigned int chstart)
{
	unsigned int th, nparam, i, v, q, d, j, sett;
	int l;
	EVREF evref;
	
	invT = invTstart;
	ch = chstart;

	nparam = model.param.size();
	paramval.resize(nparam); paramjump.resize(nparam); ntr.resize(nparam); nac.resize(nparam);
	for(th = 0; th < nparam; th++){
		paramval[th] = model.paramval[th];
		paramjump[th] = paramval[th]/10;
		ntr[th] = 0; nac[th] = 0;
	}

	indevi = indev;
	indevp.clear(); indevp.resize(data.popsize);	

	trevi.resize(data.nsettime); 
	for(i = 0; i < data.popsize; i++) addindev(i,indevi[i],xi,trevi);

	EVREFT evreft;
	vector <EVREFT> xit;
	for(i = 0; i < xi.size(); i++){
		evreft.ind = xi[i].ind;	evreft.e = xi[i].e;	evreft.t = indevi[xi[i].ind][xi[i].e].t;
		xit.push_back(evreft);	
	}
	sort(xit.begin(),xit.end(),compEVREFT);

	for(i = 0; i < xi.size(); i++){
		xi[i].ind = xit[i].ind;	xi[i].e = xit[i].e;
	}
	
	Li = Lobs_mbp(data,model,poptree,trevi,indevi);

	indinfi.clear(); for(i = 0; i < indevi.size(); i++){ if(indevi[i].size() > 0) indinfi.push_back(i);}
	
	timeprop = 0;
		
	setuplists();
	
	dQmap.resize(data.narage);                                                 // Initialises Qmapi and event buffer
	
	dQbuf.resize(data.narage);
	for(v = 0; v < data.narage; v++){
		dQbuf[v].resize(model.DQnum); for(q = 0; q < model.DQnum; q++) dQbuf[v][q] = 0;
	}
	dQbuflistv.clear(); dQbuflistq.clear();
	
	Qmapi.resize(data.nsettime); Qmapp.resize(data.nsettime); 
	for(sett = 0; sett < data.nsettime; sett++){
		Qmapi[sett].resize(data.narage); Qmapp[sett].resize(data.narage); 
	}

	setQmapi();
	
	lami.resize(data.nardp); lamp.resize(data.nardp);
	Rtot.resize(poptree.level); for(l = 0; l < poptree.level; l++) Rtot[l].resize(lev[l].node.size()); 
	
	popw.resize(data.nardp);
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
		se = (unsigned int)(data.nsettime*indev[e].t/data.period + tiny); 
		if(se < data.nsettime) trev[se].push_back(evref);
	}
}
		
/// Based on the the event sequence in xi, this sets Qmapi
void MBPCHAIN::setQmapi()
{
	unsigned int v, q, j, jmax, k, kmax, i, d, sett, a, nage, vv;
	double val;
	FEV fev;

	for(v = 0; v < data.narage; v++) dQmap[v] = 0;

	nage = data.nage;
	for(sett = 0; sett < data.nsettime; sett++){
		for(v = 0; v < data.narage; v++){
			val = dQmap[v];
			if(val < -tiny) emsg("MBPchain: EC17");
			if(val < 0){ val = 0; dQmap[v] = 0;}	
			Qmapi[sett][v] = val;
		}
		
		jmax = trevi[sett].size();
		for(j = 0; j < jmax; j++){
			i = trevi[sett][j].ind;
			fev = indevi[i][trevi[sett][j].e];

			v = data.ind[i].area*data.nage+data.democatpos[data.ind[i].dp][0];
			q = trans[fev.trans].DQ[fev.timep];
			if(q != UNSET){
				kmax = model.DQto[q][v].size();
				for(k = 0; k < kmax; k++){
					vv = model.DQto[q][v][k]*nage;	
					for(a = 0; a < nage; a++){
						dQmap[vv] += model.DQval[q][v][k][a];
						vv++;
					}
				}
			}
		}
	}
}

/// Constructs the tree Rtot for fast Gilelspie sampling	
void MBPCHAIN::constructRtot(unsigned int sett)
{
	unsigned int c, cmax, wmin, wmax, dp, v, a, w, j, jmax;
	int l;
	double sum, dlam;
	
	l = poptree.level-1;
	for(c = 0; c < data.narea; c++){
		wmin = c*data.ndemocatpos; wmax = wmin + data.ndemocatpos;
	
		sum = 0; dp = 0; v = c*data.nage; 
		for(w = wmin; w < wmax; w++){
			a = data.democatpos[dp][0];
			lami[w] = model.susi[dp]*(betai*Qmapi[sett][v+a] + phii);
			lamp[w] = model.susp[dp]*(betap*Qmapp[sett][v+a] + phip);

			dlam = nindbothlist[w]*(lamp[w] - lami[w]); if(dlam < 0) dlam = 0;
			sum += dlam + nindponlylist[w]*lamp[w];
			dp++;
		}
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
}
	
/// Performs an MBP on parameter 'th'
void MBPCHAIN::proposal(DATA &data, MODEL &model, POPTREE &poptree, unsigned int th, unsigned int samp, unsigned int burnin)  
{
	unsigned int j, jmax;
	double al, valst, Lp=0;
	
	timeprop -= clock();
	timers.timembpprop -= clock();
	
	model.paramval = paramval;
	model.betaspline(data);
	model.setsus(data);
	model.parami = paramval; model.betai = model.beta; model.susi = model.sus;
			
	valst = paramval[th];

	paramval[th] += normal(0,paramjump[th]);               // Makes a change to a parameter

	if(paramval[th] < model.param[th].min || paramval[th] > model.param[th].max) al = 0;
	else{
		model.paramval = paramval;
		model.betaspline(data);
		model.setsus(data);
		model.paramp = paramval; model.betap = model.beta; model.susp = model.sus;
		
		if(model.settransprob() == 0) al = 0;
		else{
			if(mbp() == 1) al = 0;
			else{
				Lp = Lobs_mbp(data,model,poptree,trevp,indevp);
			
				al = exp(invT*(Lp-Li));	//cout << al << " al\n";/
			}
		}
	}
	
	ntr[th]++;
	if(ran() < al){
		Li = Lp;
		xi = xp;
		trevi = trevp;
		
		Qmapi = Qmapp;
		
		jmax = indinfi.size(); for(j = 0; j < jmax; j++) indevi[indinfi[j]].clear();
		jmax = indinfp.size(); for(j = 0; j < jmax; j++) indevi[indinfp[j]] = indevp[indinfp[j]];
		//indevi = indevp;
		indinfi = indinfp;
		nac[th]++;
		if(samp < burnin) paramjump[th] *= 1.1;
	}
	else{
		paramval[th] = valst;
		if(samp < burnin) paramjump[th] *= 0.95;
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
void MBPCHAIN::changestat(unsigned int i, unsigned int st)
{
	unsigned int c, w, n;
	int l;
	double dlam, dval;
	
	c = data.ind[i].area;
	w = c*data.ndemocatpos + data.ind[i].dp;;
		
	dlam = nindbothlist[w]*(lamp[w] - lami[w]); if(dlam < 0) dlam = 0;
	dval = -(dlam + nindponlylist[w]*lamp[w]);
		
	l = indlistref[i];   // Removes the exisiting entry
	switch(stat[i]){
  case BOTH:
		if(indbothlist[w][l] != i) emsg("MBPchain: EC43");
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
	
	default: emsg("MBPChain EC47"); break;
	}
	
	dlam = nindbothlist[w]*(lamp[w] - lami[w]); if(dlam < 0) dlam = 0;
	dval += dlam + nindponlylist[w]*lamp[w];
	
	if(dval != 0){
		l = poptree.level-1;
		if(dval != 0){
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

/// Performs a MBP
unsigned int MBPCHAIN::mbp()
{
	unsigned int j, jmax, c, l, v, sett, n, i, w;
	double t, tmax, val, txi, tinf, al;
	FEV ev;

	timers.timembpinit -= clock();
		
	N.resize(comp.size()); for(c = 0; c < comp.size(); c++) N[c] = 0;
 	N[0] = data.popsize;
		
	jmax = indinfp.size(); for(j = 0; j < jmax; j++) indevp[indinfp[j]].clear();
	indinfp.clear();
	//indevp.clear(); indevp.resize(data.popsize);	
	
	xp.clear();
	trevp.clear(); trevp.resize(data.nsettime);
	
	for(v = 0; v < data.narage; v++) dQmap[v] = 0;
	
	timers.timembpinit += clock();
	
	timers.timembp -= clock();
		
	t = 0; n = 0;
	for(sett = 0; sett < data.nsettime; sett++){
		phii = model.parami[model.phiparam]; phip = model.paramp[model.phiparam];	
		betai = model.betai[sett]; betap = model.betap[sett];

		for(v = 0; v < data.narage; v++){
			val = Qmapi[sett][v] + dQmap[v];
			if(val < -tiny) emsg("MBPchain: EC31");
			if(val < 0) val = 0;	
			Qmapp[sett][v] = val;
		}
			
		constructRtot(sett);
		
		tmax = data.settime[sett+1];
		do{
			if(n < xi.size()){ ev = indevi[xi[n].ind][xi[n].e]; txi = ev.t;} else txi = tmax;
			
			if(Rtot[0][0] <= 0) tinf = tmax; else tinf = t - log(ran())/Rtot[0][0];
				
			if(txi >= tmax && tinf >= tmax) break;
			
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
						changestat(i,NOT);
						
						indevp[i] = indevi[i];
						//model.mbpmodel(indevi[i],indevp[i]);
						addindev(i,indevp[i],xp,trevp);

						indinfp.push_back(i);
					}
					else changestat(i,PONLY);      // Does not keep the infection event
				}
				n++;
			}
			
			if(indinfp.size() >= INFMAX) return 1; 
		}while(1 == 1);
		
		updatedQmap(sett);
		
		if(checkon == 1) check(1,t,sett);
	}

	timers.timembp += clock();
		
	resetlists();
	
	return 0;
}

/// Updates dQmap based on events which occur in timestep sett in the initial and proposed states
void MBPCHAIN::updatedQmap(unsigned int sett)
{
	unsigned int j, jmax, k, kmax, i, v, q, vv, a, nage;
	int num;
	FEV fev;
	TRANS tr;
	
	nage = data.nage;
	
	jmax = trevi[sett].size();
	for(j = 0; j < jmax; j++){
		i = trevi[sett][j].ind;
		fev = indevi[i][trevi[sett][j].e];

		v = data.ind[i].area*data.nage+data.democatpos[data.ind[i].dp][0];
		q = trans[fev.trans].DQ[fev.timep];
		if(q != UNSET){
			if(dQbuf[v][q] == 0){ dQbuflistv.push_back(v); dQbuflistq.push_back(q);}
			dQbuf[v][q]--;
		}
	}
	
	jmax = trevp[sett].size();
	for(j = 0; j < jmax; j++){
		i = trevp[sett][j].ind;
		fev = indevp[i][trevp[sett][j].e];

		tr = trans[fev.trans];
		N[tr.from]--; if(N[tr.from] < 0) emsg("MBPCHAIN: EC6"); 
		N[tr.to]++;
	
		v = data.ind[i].area*data.nage+data.democatpos[data.ind[i].dp][0];
		q = trans[fev.trans].DQ[fev.timep];
		if(q != UNSET){
			if(dQbuf[v][q] == 0){ dQbuflistv.push_back(v); dQbuflistq.push_back(q);}
			dQbuf[v][q]++;
		}
	}
	
	timers.timembpQmap -= clock();
	
	nage = data.nage;
	jmax = dQbuflistv.size();
	for(j = 0; j < jmax; j++){
		v = dQbuflistv[j]; q = dQbuflistq[j]; 
		num = dQbuf[v][q];
		if(num != 0){
			kmax = model.DQto[q][v].size();
			for(k = 0; k < kmax; k++){
				vv = model.DQto[q][v][k]*nage;	
				for(a = 0; a < nage; a++){
					dQmap[vv] += num*model.DQval[q][v][k][a];
					vv++;
				}
			}
			dQbuf[v][q] = 0;
		}
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
		if(j == jmax) emsg("Part: EC15");
		
		c = lev[l].node[c].child[j];
		l++;
	};
	
	return c;
}

/// Adds an exposed indivdual in area
void MBPCHAIN::addinfc(unsigned int c, double t)
{
	unsigned int dp, dpmax, w, i, j, jmax, n;
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
	
	changestat(i,NOT);
	
	model.simmodel(indevp[i],i,0,t);
	addindev(i,indevp[i],xp,trevp);

	indinfp.push_back(i);
}

/// Used for checking the code is running correctly
void MBPCHAIN::check(unsigned int num, double t, unsigned int sett)
{
	unsigned int c, d, l, i, j, k, w, dp, a, wmin, wmax, v;
	double dd, dlam, sui, sup, sum;
	vector <double> lai,lap;
	
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
			dd = lami[w] - model.susi[dp]*(betai*Qmapi[sett][v+a] + phii); if(sqrt(dd*dd) > tiny) emsg("MBPchain: EC67");
			dd = lamp[w] - model.susp[dp]*(betap*Qmapp[sett][v+a] + phip); if(sqrt(dd*dd) > tiny) emsg("MBPchain: EC68");
	
			dlam = nindbothlist[w]*(lamp[w] - lami[w]); if(dlam < 0) dlam = 0;
			sum += dlam + nindponlylist[w]*lamp[w];
			dp++;
		}
		dd = Rtot[l][c] - sum; if(sqrt(dd*dd) > tiny){ emsg("MBPchain: EC69");}
	}
}

/// This generates parameter proposals based on fixed event sequence
void MBPCHAIN::param_prop()
{
/*

	vector <vector <double> > Linf_beta, Linf_phi;
	popw.resize(data.nardp);
	Linf_beta.resize(data.nsettime); Linf_phi.resize(data.nsettime);
	*/
	
	model.paramval = paramval;
	model.betaspline(data);
	model.setsus(data);
	model.parami = paramval; model.betai = model.beta; model.susi = model.sus;
		
	cout << "Lik\n";
	cout << likelihood() << " " << " L\n";
	emsg("cl");
}

/// Calculates the likelihood
double MBPCHAIN::likelihood()
{	
/*
	unsigned int c, dp, w, v, d, i, j;
	double L, t, tt;
			cout << "st\n";
	for(c = 0; c < data.narea; c++){
		for(dp = 0; dp < data.ndemocatpos; dp++){
			w = c*data.ndemocatpos + dp;
			popw[w] = data.area[c].ind[dp].size();
		}
	}		
			
	phii = model.paramval[model.phiparam];
		
	L = 0;
	
	for(sett = 0; sett < data.nsettime; sett++){
		
	}
	
	t = 0;
	for(sett = 0; sett < data.nsettime; sett++){
		cout << sett << "s\n";
		betai = model.betai[sett];
	  
		for(c = 0; c < data.narea; c++){
			for(dp = 0; dp < data.ndemocatpos; dp++){
				w = c*data.ndemocatpos + dp;
				v = c*data.nage + data.democatpos[dp][0];
				lami[w] = model.susi[dp]*(betai*Qmapi[sett][v] + phii);
				if(lami[w] < 0) emsg("n");
			}
		}
		
		for(d = sett*data.fepertime; d < (sett+1)*data.fepertime; d++){
			for(j = 0; j < xi[d].size(); j++){
				if(xi[d][j].trans == 0){
					tt = xi[d][j].t;
					for(w = 0; w < data.nardp; w++) L -= lami[w]*popw[w]*(tt-t);
					t = tt;
					
					i = xi[d][j].ind;
					c = data.ind[i].area;
					dp = data.ind[i].dp;
					w = c*data.ndemocatpos + dp;
					v = c*data.nage + data.democatpos[dp][0]; 
					cout << L << " " << lami[w] << " kk\n";
					L += log(lami[w]);
					if(isnan(L)) emsg("P");
					popw[w]--;
				}
			}
		}
		
		tt = data.settime[sett+1]; cout << tt << " tt\n";
		for(w = 0; w < data.nardp; w++) L -= lami[w]*popw[w]*(tt-t);
		t = tt;
	} 
	
	return L;
	*/
}

/*
/// Calculates the likelihood
void MBPCHAIN::param_prop()
{
	L = 0;
		numinf[sett] = 0;
	
	inf_beta.clear(); inf_phi.clear();

	for(d = sett*data.fepertime; d < (sett+1)*data.fepertime; d++){
			for(j = 0; j < xi[d].size(); j++){
				if(xi[d][j].trans == 0){
					tt = xi[d][j].t;
					
					for(w = 0; w < data.nardp; w++){
						L -= lami[w]*popw[w]*(tt-t);
					}
					
					for(c = 0; c < data.narea; c++){
						for(dp = 0; dp < data.ndemocatpos; dp++){
							w = c*data.ndemocatpos+dp;
							v = c*data.nage + data.democatpos[dp][0]; 
							beta_dt[w] += popw[w]*model.susi[dp]*Qmapi[sett][v]*(tt-t);
							phi_dt[w] += popw[w]*model.susi[dp]*(tt-t);
						}
					}
					
					i = xi[d][j].ind;
					c = data.ind[i].area;
					dp = data.ind[i].dp;
					w = c*data.ndemocatpos + dp;
					v = c*data.nage + data.democatpos[dp][0]; 
					L += log(lami[w]);
					Linf_beta.push_back(model.susi[dp]*Qmapi[sett][v]);
					Linf_phi.push_back(model.susi[dp]);
				}
			}
		}
	} 
			for(d = sett*data.fepertime; d < (sett+1)*data.fepertime; d++){
			for(j = 0; j < xi[d].size(); j++){
				if(xi[d][j].trans == 0){
					tt = xi[d][j].t;
					
					for(w = 0; w < data.nardp; w++) L -= lami[w]*popw[w]*(tt-t);
					
					i = xi[d][j].ind;
					c = data.ind[i].area;
					dp = data.ind[i].dp;
					w = c*data.ndemocatpos + dp;
					v = c*data.nage + data.democatpos[dp][0]; 
					L += log(lami[w]);
					Linf_beta.push_back(model.susi[dp]*Qmapi[sett][v]);
					Linf_phi.push_back(model.susi[dp]);
				}
			}
		}
	} 
}
*/