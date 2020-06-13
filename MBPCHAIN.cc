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
	timeprop = 0;
		
	nparam = model.param.size();
	paramval.resize(nparam); paramjump.resize(nparam); ntr.resize(nparam); nac.resize(nparam);
	for(th = 0; th < nparam; th++){
		paramval[th] = model.paramval[th];
		paramjump[th] =paramval[th]/10;
		ntr[th] = 0; nac[th] = 0;
	}

	xi = xistart; 
	Li = Lobstot(data,model,poptree,xi,1);
	
	indevi = indev;
	
	ninftot = 0; for(i = 0; i < indevi .size(); i++){ if(indevi[i].size() > 0) ninftot++;}
}

/// Performs an MBP on parameter "th"
void MBPCHAIN::proposal(DATA &data, MODEL &model, POPTREE &poptree, int th, int samp, int burnin)  
{
	double al, valst, Lp;
	
	timeprop -= clock();
	
	model.paramval = paramval;
	model.betaspline(data.period);
	//if(model.param[p].suschange == 1) poptree.setsus(model);
	//if(model.param[p].infchange == 1) poptree.setinf(model);
	model.parami = paramval; model.betai = model.beta;
			
	valst = paramval[th];
					
	paramval[th] += normal(0,paramjump[th]);               // Makes a change to a parameter
				
	if(paramval[th] < model.param[th].min || paramval[th] > model.param[th].max) al = 0;
	else{
		model.paramval = paramval;
		model.betaspline(data.period);
		//if(model.param[p].suschange == 1) poptree.setsus(model);
		//if(model.param[p].infchange == 1) poptree.setinf(model);
		model.paramp = paramval; model.betap = model.beta;
		
		if(model.settransprob() == 0) al = 0;
		else{
			timers.timembp -= clock();
			mbpinit();
			mbp();
			timers.timembp += clock();
				
			if(ninftotprop >= INFMAX) al = 0;
			else{
				Lp = Lobstot(data,model,poptree,xp,1);
			
				al = exp(invT*(Lp-Li)); //cout << al << " al\n";
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
	
	timeprop += clock();
}
	
/// Initialises a MBP
void MBPCHAIN::mbpinit()
{
	int c, cmax, l, i;
		
	fediv = data.fediv;
	
	xp.clear(); xp.resize(fediv);
	ninftotprop = 0;
		
	N.resize(comp.size()); for(c = 0; c < comp.size(); c++) N[c] = 0;
 	N[0] = data.popsize;
	
	sett = 0;
		
	phii = model.parami[model.phiparam]; phip = model.paramp[model.phiparam];	
	
	MIi.resize(poptree.Cfine); MIp.resize(poptree.Cfine);
	susboth.resize(poptree.Cfine); susp.resize(poptree.Cfine);
	lami.resize(poptree.Cfine); lamp.resize(poptree.Cfine);
	
	l = poptree.level-1;
	for(c = 0; c < poptree.Cfine; c++){
		MIi[c] = 0; MIp[c] = 0;
		susboth[c] = lev[l].node[c].sussum; susp[c] = 0;
		lami[c] = susboth[c]*phii; lamp[c] = susboth[c]*phip; 
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
	
// Performs a MBP
void MBPCHAIN::mbp()
{
	int td, j, jsel, c, cmax;
	double t, tmax = data.period;
	
	NEV n;  
	vector <NEV> nev;

	cmax = poptree.Cfine;
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
			for(c = 0; c < cmax; c++) updateRtot(c); 
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

	if(checkon == 1) check(0);
}

/// This samples the node on the fine scale in which the next infection occurs due to the MBP
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
void MBPCHAIN::updateRtot(int c) 
{
	int l;
	double dlam, val, dval;
	
	lami[c] = betai*MIi[c] + phii;
	lamp[c] = betap*MIp[c] + phip;
			
	dlam = (lamp[c]-lami[c])*susboth[c]; if(dlam < 0) dlam = 0;
	val = dlam + lamp[c]*susp[c];
	
	l = poptree.level-1;
	dval = val-Rtot[l][c];
	do{
		Rtot[l][c] += dval;
		c = lev[l].node[c].parent; l--;
	}while(l >= 0);
}

/// Makes compartmental transitions in the proposed event sequence
void MBPCHAIN::xpdofe(int update)
{
	int tra;
	TRANS tr;

	tra = xp[xptdnext][xptdfnext].trans;
	tr = trans[tra];
	
	if(update == 1) MIupdate(xp[xptdnext][xptdfnext].ind,tra,0,1);
		
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
	int i, c, tra, j, jmax;
	double sus, al, t;
	TRANS tr;

	i = xi[xitdnext][xitdfnext].ind;
	c = poptree.ind[i].noderef;

	tra = xi[xitdnext][xitdfnext].trans;
	tr = trans[tra];
	
	t = xi[xitdnext][xitdfnext].t;                         	// This part decides to copy over event on xi into xp or not
	if(tr.from == 0){
		sus = poptree.ind[i].sus;
	
		if(indevp[i].size() == 0){
			susboth[c] -= sus;
		
			al = lamp[c]/lami[c];
			if(ran() < al){                                    // Keeps the infection event
				ninftotprop++;
				
				model.mbpmodel(indevi[i],indevp[i]);
				//indevp[i] = indevi[i];
				jmax = indevp[i].size(); for(j = 0; j < jmax; j++) addxp(indevp[i][j],data.period,t);
			}
			else susp[c] += sus;                               // Does not keep the infection event
			updateRtot(c); 
		}	
	}
	else{
		if(xptdnext < fediv && xp[xptdnext][xptdfnext].t == t){
			xpdofe(0);
			MIupdate(i,tra,1,1);
		}
		else{
			MIupdate(i,tra,1,0);
		}
	}
	
	xitdfnext++;
	if(xitdfnext == xi[xitdnext].size()){
		xitdnext++; xitdfnext = 0; 
		while(xitdnext < fediv && xi[xitdnext].size() == 0) xitdnext++;
	}
}

/// Updates the force of infection generated by the mixing matrix
void MBPCHAIN::MIupdate(int i, int tra, int upMIi, int upMIp)
{
	int c, cc, l, k, kmax, j, jmax;
	TRANS tr;
	double fac, val, dval, dlam;
		
	tr = trans[tra];
	fac = poptree.ind[i].inf*(comp[tr.to].infectivity - comp[tr.from].infectivity);
	if(fac == 0) return;

	c = poptree.ind[i].noderef;

	l = poptree.level-1;
	kmax = poptree.nMfineval[c];
	for(k = 0; k < kmax; k++){
		cc = poptree.Mfinenoderef[c][k];
		val = fac*poptree.Mfineval[c][k];
		
		if(upMIi == 1) MIi[cc] += val;	
		if(upMIp == 1) MIp[cc] += val;
		
		//updateRtot(cc); 
		
		lami[cc] = betai*MIi[cc] + phii;         // This is a faster method than using updateRtot on every node
		lamp[cc] = betap*MIp[cc] + phip;
				
		dlam = (lamp[cc]-lami[cc])*susboth[cc]; if(dlam < 0) dlam = 0;
		val = dlam + lamp[cc]*susp[cc];
		dval = val-Rtot[l][cc];
		Rtot[l][cc] += dval;
		lev[l-1].add[lev[l].node[cc].parent] += dval; 
	}
	
	for(l = poptree.level-2; l >= 0; l--){    // Adds the contribution up the tree
		kmax = poptree.Mfineadd[c][l].size();
		for(k = 0; k < kmax; k++){
			cc = poptree.Mfineadd[c][l][k];
			dval = lev[l].add[cc];
			Rtot[l][cc] += dval;
			lev[l].add[cc] = 0;
			if(l > 0) lev[l-1].add[lev[l].node[cc].parent] += dval; 
		}
	}
}

/// Adds an exposed indivdual on node c on the finest scale (i.e. level-1)
void MBPCHAIN::addinfc(int c, double t)
{
	int l, i, j, jmax, cc, k, kmax;
	double dR, sum, sus, z, dlam;
	vector <double> sumst;
	
	jmax = poptree.subpop[c].size(); 
	sumst.resize(jmax);
	sum = 0; 
	for(j = 0; j < jmax; j++){
		i = poptree.subpop[c][j];
		if(indevp[i].size() == 0){
			if(indevi[i].size() == 0 || t < indevi[i][0].t){
				dlam = (lamp[c]-lami[c])*poptree.ind[i].sus; if(dlam < 0) dlam = 0;
				sum += dlam;
			}
			else{
				sum += lamp[c]*poptree.ind[i].sus;
			}
		}
		sumst[j] = sum;
	}
	
	z = ran()*sum;       
	j = 0; while(j < jmax && z > sumst[j]) j++; if(j == jmax) emsg("MBP: EC11");
	i = poptree.subpop[c][j];

	sus = poptree.ind[i].sus;
	if(indevi[i].size() == 0 || t < indevi[i][0].t) susboth[c] -= sus;
	else susp[c] -= sus;
	updateRtot(c); 
	
	model.simmodel(indevp[i],i,0,t);
			
	jmax = indevp[i].size(); for(j = 0; j < jmax; j++) addxp(indevp[i][j],data.period,t);
}

/// Used for checking the code is running correctly
void MBPCHAIN::check(int num)
{
	int c, cmax, l, j, jmax, i;
	double d, val, dlam, sumboth, sump;

	cmax = poptree.Cfine;      

	for(c = 0; c < cmax; c++){
		sumboth = 0; sump = 0;
		jmax = poptree.subpop[c].size(); 
		for(j = 0; j < jmax; j++){
			i = poptree.subpop[c][j];
			if(indevi[i].size() == 0 && indevp[i].size() == 0) sumboth += poptree.ind[i].sus; 
			if(indevi[i].size() != 0 && indevp[i].size() == 0) sump += poptree.ind[i].sus; 
		}
		d = sumboth - susboth[c]; if(d*d > tiny) emsg("both prob");
		d = sump - susp[c]; if(d*d > tiny) emsg("sump prob");
	}
		
	l = poptree.level-1;
	for(c = 0; c < cmax; c++){
		lami[c] = betai*MIi[c] + phii;
		lamp[c] = betap*MIp[c] + phip;
			
		dlam = (lamp[c]-lami[c])*susboth[c];
		if(dlam < 0) dlam = 0;
		val = dlam + lamp[c]*susp[c];
		d = val-Rtot[l][c]; if(d*d > tiny) emsg("MBPCHAIN: EC3");
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
