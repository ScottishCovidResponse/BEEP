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

MBPCHAIN::MBPCHAIN(DATA &data, MODEL &model, POPTREE &poptree) : data(data), model(model), comp(model.comp), trans(model.trans), poptree(poptree), lev(poptree.lev)
{
}

/// Initialises an MCMC chain
void MBPCHAIN::init(DATA &data, MODEL &model, POPTREE &poptree, double invTstart, vector < vector <FEV> > &xistart, vector < vector <FEV> > &indev, int chstart)
{
	int th, nparam, d, j, i;

	invT = invTstart;
	ch = chstart;
	fediv = data.fediv;
		
	nparam = model.param.size();
	paramval.resize(nparam); paramjump.resize(nparam); ntr.resize(nparam); nac.resize(nparam);
	for(th = 0; th < nparam; th++){
		paramval[th] = model.paramval[th];
		paramjump[th] = paramval[th]/10;
		ntr[th] = 0; nac[th] = 0;
	}

	xi = xistart; 
	Li = Lobstot(data,model,poptree,xi,1);
	
	indevi = indev;
	
	ninftot = 0; for(i = 0; i < indevi.size(); i++){ if(indevi[i].size() > 0) ninftot++;}
	
	timeprop = 0;
		
	setuplists();
}

/// Performs an MBP on parameter 'th'
void MBPCHAIN::proposal(DATA &data, MODEL &model, POPTREE &poptree, int th, int samp, int burnin)  
{
	double al, valst, Lp;
	
	timeprop -= clock();
	timers.timembpprop -= clock();
	
	model.paramval = paramval;
	model.betaspline(data.period);
	model.setsus(data);
	model.parami = paramval; model.betai = model.beta; model.susi = model.sus;
			
	valst = paramval[th];

	paramval[th] += normal(0,paramjump[th]);               // Makes a change to a parameter

	if(paramval[th] < model.param[th].min || paramval[th] > model.param[th].max) al = 0;
	else{
		model.paramval = paramval;
		model.betaspline(data.period);
		model.setsus(data);
		model.paramp = paramval; model.betap = model.beta; model.susp = model.sus;
		
		if(model.settransprob() == 0) al = 0;
		else{
			timers.timembpinit -= clock();
			mbpinit();
			timers.timembpinit += clock();
			
			timers.timembp -= clock();
			mbp();
			timers.timembp += clock();
				
			if(ninftotprop >= INFMAX) al = 0;
			else{
				Lp = Lobstot(data,model,poptree,xp,1);
			
				al = exp(invT*(Lp-Li)); cout << al << " al\n";
			}
		}
	}
	
	ntr[th]++;
	if(ran() < al){
		Li = Lp;
		xi = xp;
		indevi = indevp;
		ninftot = ninftotprop;
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
	int c, dp, j, jmax, w, i;

	indbothlist.resize(data.nardp); indponlylist.resize(data.nardp); indnotlist.resize(data.nardp);
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
		}
	}	
}

/// Makes a change in the status of an individual
void MBPCHAIN::changestat(int i, int st)
{
	int w, l, dp, n;
	double sui, sup;
	
	dp = data.ind[i].dp;
	w = data.ind[i].area*data.ndemocatpos + dp;
	sui = model.susi[dp];
	sup = model.susp[dp];
		
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
		susbothi[w] -= sui;
		susbothp[w] -= sup;
		break;
		
	case PONLY:
		if(indponlylist[w][l] != i) emsg("MBPchain: EC44");
		n = indponlylist[w].size();
		if(l < n-1){
			indponlylist[w][l] = indponlylist[w][n-1];
			indlistref[indponlylist[w][l]] = l;
		}
		indponlylist[w].pop_back();
		susponly[w] -= sup;
		break;
		
	default: emsg("MBPChain EC46"); break;
	}

	stat[i] = st;
	switch(stat[i]){
	case PONLY:
		indlistref[i] = indponlylist[w].size();
		indponlylist[w].push_back(i);
		susponly[w] += sup;
		break;
		
	case NOT:
		indlistref[i] = indnotlist[w].size();
		indnotlist[w].push_back(i);
		break;
	
	default: emsg("MBPChain EC47"); break;
	}
}

/// Places all individuals back onto the both susceptible list 
void MBPCHAIN::resetlists()
{
	int w, j, jmax, i;
	
	for(w = 0; w < data.nardp; w++){
		jmax = indponlylist[w].size();
		for(j = 0; j < jmax; j++){
			i = indponlylist[w][j];
			stat[i] = BOTH;
			indlistref[i] = indbothlist[w].size();
			indbothlist[w].push_back(i);
		}
		indponlylist[w].clear();
		
		jmax = indnotlist[w].size();
		for(j = 0; j < jmax; j++){
			i = indnotlist[w][j];
			stat[i] = BOTH;
			indlistref[i] = indbothlist[w].size();
			indbothlist[w].push_back(i);
		}
		indnotlist[w].clear();
	}
}

/// Initialises a MBP
void MBPCHAIN::mbpinit()
{
	int c, cmax, l, i, v, w, dp;
		
	fediv = data.fediv;
	
	xp.clear(); xp.resize(fediv);
	ninftotprop = 0;
		
	N.resize(comp.size()); for(c = 0; c < comp.size(); c++) N[c] = 0;
 	N[0] = data.popsize;
	
	sett = 0;
		
	phii = model.parami[model.phiparam]; phip = model.paramp[model.phiparam];	
	
	Qmapi.resize(data.narage); Qmapp.resize(data.narage);
	for(v = 0; v < data.narage; v++){ Qmapi[v] = 0; Qmapp[v] = 0;}
	
	sussame = 1;
	
	susbothi.resize(data.nardp); susbothp.resize(data.nardp); susponly.resize(data.nardp);
	lami.resize(data.nardp); lamp.resize(data.nardp);
	for(w = 0; w < data.nardp; w++){
		c = w/data.ndemocatpos; dp = w%data.ndemocatpos;
		
		if(data.area[c].pop[dp] != indbothlist[w].size()) emsg("MBPchain: EC21");
		susbothi[w] = model.susi[dp]*data.area[c].pop[dp];
		susbothp[w] = model.susp[dp]*data.area[c].pop[dp];
		if(susbothi[w] != susbothp[w]) sussame = 0;
		susponly[w] = 0;
		lami[w] = susbothi[w]*phii; lamp[w] = susbothp[w]*phip; 
	}
	
	Rtot.resize(poptree.level);
	for(l = 0; l < poptree.level; l++){
		cmax = lev[l].node.size();
		Rtot[l].resize(cmax); 
		for(c = 0; c < cmax; c++) Rtot[l][c] = 0; 
	}
	
	indevp.clear(); indevp.resize(data.popsize);
	
	xitdnext = 0; xitdfnext = 0;
	while(xitdnext < fediv && xi[xitdnext].size() == 0) xitdnext++;	
	
	xptdnext = fediv;
}
	
/// Performs a MBP
void MBPCHAIN::mbp()
{
	int td, j, jsel, c, cmax, l;
	double t, tmax = data.period;
	NEV n;  
	vector <NEV> nev;

	cmax = data.narea;
	if(sett == nsettime) emsg("MBP: EC1");

	betai = model.betai[sett]; betap = model.betap[sett];

	t = 0; 
	do{
		nev.clear();                           // We decide what event is next
		n.t = model.settime[sett];
		n.type = SET_EV;
		nev.push_back(n); 
		 
		if(xitdnext < fediv) n.t = xi[xitdnext][xitdfnext].t; else n.t = tmax;
		n.type = XIFEV_EV;
		nev.push_back(n);

		if(xptdnext < fediv) n.t = xp[xptdnext][xptdfnext].t; else n.t = tmax;
		n.type = XPFEV_EV;
		nev.push_back(n);
	
		if(Rtot[0][0] < tiny) n.t = tmax; else n.t = t - log(ran())/Rtot[0][0];
		n.type = INF_EV;
		nev.push_back(n);
		
		t = tmax; for(j = 0; j < nev.size(); j++){ if(nev[j].t < t){ t = nev[j].t; jsel = j;}}
		if(t == tmax) break;
	
		switch(nev[jsel].type){
		case SET_EV:                 // These are "settime" events which allow the value of beta to change in time
			sett++; if(sett >= nsettime) emsg("MBP: EC2");
			betai = model.betai[sett]; betap = model.betap[sett];
	
			for(l = 0; l < poptree.level; l++) lev[l].donelist.clear();
			for(c = 0; c < cmax; c++) updateRtot(c,0);
			uptree();
			break;
		
		case INF_EV:                 // These are new infection events introduced by the MBP
			ninftotprop++;	
			c = nextinfection();
			addinfc(c,t);	
			break;
			
		case XIFEV_EV:               // Events on xi
			xidofe();
			break;
			
		case XPFEV_EV:               // Events on xp
			xpdofe(1);
			break;
		
		default: emsg("MBP: EC4"); break;
		}
		
		if(ninftotprop >= INFMAX) return;  // This limits the total number of infections within the system
	}while(t < tmax);

	if(checkon == 1) check(1,t);
	
	resetlists();
}

/// This samples the area in which the next infection occurs due to the MBP
int MBPCHAIN::nextinfection()
{
	int l, lmax, c, j, jmax;
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

// Updates Rtot so samping can be performed for generating new infection events
void MBPCHAIN::updateRtot(int c, int single) 
{
	int l, w, wmin, wmax, dp, v, a, cc;
	double dlam, sum, dval;
	
	wmin = c*data.ndemocatpos; wmax = wmin + data.ndemocatpos;
	
	sum = 0; dp = 0; v = c*data.nage; 
	for(w = wmin; w < wmax; w++){
		a = data.democatpos[dp][0];
		lami[w] = betai*Qmapi[v+a] + phii;
		lamp[w] = betap*Qmapp[v+a] + phip;
	
		dlam = susbothp[w]*lamp[w] - susbothi[w]*lami[w]; if(dlam < 0) dlam = 0;
		sum += dlam + susponly[w]*lamp[w];
		dp++;
	}
	
	l = poptree.level-1;
	dval = sum-Rtot[l][c];
	if(single == 1){
		do{
			Rtot[l][c] += dval;
			c = lev[l].node[c].parent; l--;
		}while(l >= 0);
	}
	else{
		Rtot[l][c] = sum;
		cc = lev[l].node[c].parent;
		if(lev[l-1].add[cc] == 0){ lev[l-1].donelist.push_back(cc);}		
		lev[l-1].add[cc] += dval;
	}
}

/// Makes compartmental transitions in the proposed event sequence
void MBPCHAIN::xpdofe(int update)
{
	int tra;
	TRANS tr;

	tra = xp[xptdnext][xptdfnext].trans;
	tr = trans[tra];
	
	if(update == 1) Qmapupdate(xp[xptdnext][xptdfnext].ind,tra,xp[xptdnext][xptdfnext].t,0,1);
		
	N[tr.from]--; if(N[tr.from] < 0) emsg("MBPCHAIN: EC6"); 
	N[tr.to]++;
	
	xptdfnext++;
	if(xptdfnext == xp[xptdnext].size()){
		xptdnext++; xptdfnext = 0; 
		while(xptdnext < fediv && xp[xptdnext].size() == 0) xptdnext++;
	}
}

/// Goes thorough events on the exisitng event sequence xi and decide to keep them or not
void MBPCHAIN::xidofe()
{
	int i, c, tra, j, jmax, dp, w;
	double al, t;
	TRANS tr;

	i = xi[xitdnext][xitdfnext].ind;
	c = data.ind[i].area;

	tra = xi[xitdnext][xitdfnext].trans;
	tr = trans[tra];
	
	t = xi[xitdnext][xitdfnext].t;                         	// This part decides to copy over event on xi into xp or not
	if(tr.from == 0){
		if(indevp[i].size() == 0){
			dp = data.ind[i].dp;
			w = c*data.ndemocatpos + dp;
		
			al = (model.susp[dp]*lamp[w])/(model.susi[dp]*lami[w]);
			if(ran() < al){                                    // Keeps the infection event
				changestat(i,NOT);
				ninftotprop++;
				
				indevp[i] = indevi[i];
				//model.mbpmodel(indevi[i],indevp[i]);
				jmax = indevp[i].size(); for(j = 0; j < jmax; j++) addxp(indevp[i][j],data.period,t);
			}
			else{ 
				changestat(i,PONLY);
			}                               // Does not keep the infection event
			updateRtot(c,1); 
		}	
	}
	else{
		if(xptdnext < fediv && xp[xptdnext][xptdfnext].t == t){
			xpdofe(0);
			Qmapupdate(i,tra,t,1,1);
		}
		else{
			Qmapupdate(i,tra,t,1,0);
		}
	}
	
	xitdfnext++;
	if(xitdfnext == xi[xitdnext].size()){
		xitdnext++; xitdfnext = 0; 
		while(xitdnext < fediv && xi[xitdnext].size() == 0) xitdnext++;
	}
}

/// Updates the force of infection generated by the mixing matrix
void MBPCHAIN::Qmapupdate(int i, int tra, double t, int upMIi, int upMIp)
{
	int c, cc, ccc, l, k, kmax, j, jmax, q, timep, dp, v, vv, a, nage, w, wmin, wmax, per, nper, ndemocatpos;
	TRANS tr;
	double val, val2, dval, dlam, sum;
		
	nper = data.ndemocatposperage;
	nage = data.nage;
	ndemocatpos = data.ndemocatpos;
	
	tr = trans[tra];
	
	timep = 0; while(timep < data.ntimeperiod && t > data.timeperiod[timep]) timep++;
	if(timep >= tr.DQ.size()) emsg("MBPchain: EC66");
	
	q = tr.DQ[timep]; if(q == -1) return;
	
	c = data.ind[i].area;
	dp = data.ind[i].dp;
	v = c*nage+data.democatpos[dp][0];
	
	timers.timembpQmap -= clock();
	
	for(l = 0; l < poptree.level; l++) lev[l].donelist.clear();

	l = poptree.level-1;                                                             // Makes change to Qmap
	kmax = model.nDQ[q][v];
	for(k = 0; k < kmax; k++){
		cc = model.DQto[q][v][k];
		vv = cc*nage; w = cc*ndemocatpos; 
		
		if(upMIi == 1 && upMIp == 1 && sussame == 1 && betai == betap){                // This is used to speed up things
			dval = 0;
			for(a = 0; a < nage; a++){
				val = model.DQval[q][v][k][a];
				Qmapi[vv] += val;	Qmapp[vv] += val;	vv++;
				val2 = betai*val;
				for(per = 0; per < nper; per++){ lami[w] += val2; lamp[w] += val2; dval += susponly[w]*val2; w++;}
			}
		
			Rtot[l][cc] += dval;
		}
		else{
			if(upMIi == 1){ 
				if(upMIp == 1){
					for(a = 0; a < nage; a++){
						val = model.DQval[q][v][k][a];
						Qmapi[vv] += val;	Qmapp[vv] += val;	vv++;
						for(per = 0; per < nper; per++){ lami[w] += betai*val; lamp[w] += betap*val; w++;}
					}
				}
				else{
					for(a = 0; a < nage; a++){
						val = model.DQval[q][v][k][a];
						Qmapi[vv] += val; vv++;
						for(per = 0; per < nper; per++){ lami[w] += betai*val; w++;}
					}
				}
			}
			else{
				if(upMIp == 1){
					for(a = 0; a < nage; a++){
						val = model.DQval[q][v][k][a];
						Qmapp[vv] += val; vv++;
						for(per = 0; per < nper; per++){ lamp[w] += betap*val; w++;}
					}
				}
			}

			wmin = cc*ndemocatpos; wmax = wmin + ndemocatpos;
		
			sum = 0;
			for(w = wmin; w < wmax; w++){
				dlam = susbothp[w]*lamp[w] - susbothi[w]*lami[w]; if(dlam < 0) dlam = 0;
				sum += dlam + susponly[w]*lamp[w];
			}
			
			dval = sum-Rtot[l][cc];
			Rtot[l][cc] = sum;
		}
	
		ccc = lev[l].node[cc].parent;
		if(lev[l-1].add[ccc] == 0){ lev[l-1].donelist.push_back(ccc);}		
		lev[l-1].add[ccc] += dval;
	}
	
	uptree();
		
	timers.timembpQmap += clock();
}

/// Propages changes in Rtot up the tree
void MBPCHAIN::uptree()
{
	int l, j, jmax, c, cc;
	double dval;
	
	for(l = poptree.level-2; l >= 0; l--){                                 
		jmax = lev[l].donelist.size();
		for(j = 0; j < jmax; j++){
			c = lev[l].donelist[j];
			
			dval = lev[l].add[c];
			Rtot[l][c] += dval;
			lev[l].add[c] = 0;
			 
			if(l > 0){
				cc = lev[l].node[c].parent;
				if(lev[l-1].add[cc] == 0){ lev[l-1].donelist.push_back(cc);}		
				lev[l-1].add[cc] += dval;
			}
		}
	}
}

/// Adds an exposed indivdual in area
void MBPCHAIN::addinfc(int c, double t)
{
	int dp, dpmax, w, i, j, jmax, n;
	double sum, dlam, z;
	vector <double> sumst;
	
	dpmax = data.ndemocatpos;
	sumst.resize(dpmax);
	
	sum = 0;                           // First selects the demographic possibility
	for(dp = 0; dp < dpmax; dp++){
		w = c*dpmax + dp;
		dlam = susbothp[w]*lamp[w] - susbothi[w]*lami[w]; if(dlam < 0) dlam = 0;
		sum += dlam + susponly[w]*lamp[w];
		sumst[dp] = sum;
	}
	
	z = ran()*sum; dp = 0; while(dp < dpmax && z > sumst[dp]) dp++; 
	if(dp == dpmax) emsg("MBPchain: EC56");
	
	w = c*dpmax + dp;                  // Next select when individuals in the initial and proposed states are susceptible
	dlam = susbothp[w]*lamp[w] - susbothi[w]*lami[w]; if(dlam < 0) dlam = 0;
	if(ran() < dlam/(dlam +  susponly[w]*lamp[w])){ // Both suscetible
		n = indbothlist[w].size(); if(n == 0) emsg("MBPchain: EC43");
		i = indbothlist[w][int(ran()*n)];
	}
	else{                                       // Only proposed state susceptible
		n = indponlylist[w].size(); if(n == 0) emsg("MBPchain: EC44");
		i = indponlylist[w][int(ran()*n)];
	}
	
	changestat(i,NOT);
	
	updateRtot(c,1); 
	
	model.simmodel(indevp[i],i,0,t);
	jmax = indevp[i].size(); for(j = 0; j < jmax; j++) addxp(indevp[i][j],data.period,t);
}

/// Used for checking the code is running correctly
void MBPCHAIN::check(int num, double t)
{
	int c, cmax, l, j, jmax, i, w, dp, a;
	double d, val, dlam, sui, sup, sum, dd;
	vector <double> lai,lap;

	for(i = 0; i < data.popsize; i++){    // checks stat is correct
	  w = data.ind[i].area*data.ndemocatpos + data.ind[i].dp;
		
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
	
	for(c = 0; c < data.narea; c++){
		for(w = 0; w < data.nardp; w++){	
			dp = w%data.ndemocatpos;
			sui = model.susi[dp];
			sup = model.susp[dp];
	
			dd = susbothi[w] - sui*indbothlist[w].size(); if(d*d > tiny) emsg("MBPchain: EC28");
			dd = susbothp[w] - sup*indbothlist[w].size(); if(d*d > tiny) emsg("MBPchain: EC29");
			dd = susponly[w] - sup*indponlylist[w].size(); if(d*d > tiny) emsg("MBPchain: EC30");
		}
	}
		
	lai.resize(data.nardp); lap.resize(data.nardp);
		
	l = poptree.level-1;
	for(c = 0; c < data.narea; c++){
		sum = 0; 
		for(dp = 0; dp < data.ndemocatpos; dp++){
			a = data.democatpos[dp][0];
			w = c*data.ndemocatpos + dp;
			
			lai[w] = betai*Qmapi[c*data.nage+a] + phii;
			lap[w] = betap*Qmapp[c*data.nage+a] + phip;
	
			dd = lai[w] - lami[w]; if(d*d > tiny) emsg("MBPchain: EC31");
			dd = lap[w] - lamp[w]; if(d*d > tiny) emsg("MBPchain: EC32");
	
			dlam = susbothp[w]*lamp[w] - susbothi[w]*lami[w]; if(dlam < 0) dlam = 0;
			sum += dlam + susponly[w]*lamp[w];
		}
		dd = sum - Rtot[l][c]; if(d*d > tiny) emsg("MBPchain: EC33");
	}
	
	for(l = 0; l < poptree.level; l++){
		for(c = 0; c < lev[l].node.size(); c++){
			if( lev[l].add[c] != 0) emsg("Pe2");
		}
	}
}

/// Adds a future event to the timeline
void MBPCHAIN::addxp(FEV fe, double period, double tnow)
{
	int d, j, jmax;
	double t;
	
	t = fe.t; if(t < tnow) emsg("MBPCHAIN: EC10");
	if(t >= period) return;
	
	d = int((t/period)*xp.size());
	j = 0; jmax = xp[d].size();
	if(t != tnow){ while(j < jmax && t >= xp[d][j].t) j++;}
	else{ while(j < jmax && t > xp[d][j].t) j++;}
	
	if(j == jmax) xp[d].push_back(fe);
	else xp[d].insert(xp[d].begin()+j,fe);
	
	if(t != tnow){
		if(d == xptdnext){ if(j < xptdfnext) xptdfnext = j;}
		if(d < xptdnext){ xptdnext = d; xptdfnext = j;}
	}
	else{
		TRANS tr = trans[fe.trans];
		N[tr.from]--; if(N[tr.from] < 0) emsg("Part: EC12"); 
		N[tr.to]++;
		
		if(d == xptdnext) xptdfnext++;
	}
}
