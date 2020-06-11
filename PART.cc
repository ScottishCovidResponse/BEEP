// This contains all the code related to particles in PMCMC

#include <iostream>
#include <algorithm>

using namespace std;

#include "math.h"

#include "timers.hh"
#include "utils.hh"
#include "PART.hh"
#include "output.hh"
#include "consts.hh"
#include "pack.hh"

PART::PART(DATA &data, MODEL &model, POPTREE &poptree) : data(data), model(model), comp(model.comp), trans(model.trans), poptree(poptree), lev(poptree.lev)
{
}

/// Initialises a particle
void PART::partinit(int p)
{
	int c, cmax, cc, k, kmax, h, i, imax, j, jmax, l, loop;
	
	fediv = data.fediv;
		
	pst = p;
	N.resize(comp.size()); for(c = 0; c < comp.size(); c++) N[c] = 0;
	N[0] = data.popsize;
 
	ffine.clear(); 
	if(checkon == 1){ ffine.resize(poptree.Cfine); for(c = 0; c < poptree.Cfine; c++) ffine[c] = 0;}

	indinf.clear(); indinf.resize(poptree.Cfine);
	
	fev.clear(); fev.resize(fediv);

	Rtot.resize(poptree.level); addlater.resize(poptree.level); sussum.resize(poptree.level);
	for(l = 0; l < poptree.level; l++){
		cmax = lev[l].node.size();
		Rtot[l].resize(cmax); addlater[l].resize(cmax); sussum[l].resize(cmax);
		for(c = 0; c < cmax; c++){ 
			Rtot[l][c] = 0; addlater[l][c] = 0;
			sussum[l][c] = lev[l].node[c].sussum;
			//sussumheff[l][c] = 0;
		}
	}
	
	sett = 0;
	
	tdnext = fediv;
}

/// Adds an exposed indivdual on node c on the finest scale (i.e. level-1)
void PART::addinfc(int c, double t)
{
	int l, i, j, jmax, cc, k, kmax;
	double dR, sum, sus, z;
	vector <double> sumst;
	
	jmax = poptree.subpop[c].size(); 
	sumst.resize(jmax);
	sum = 0; 
	for(j = 0; j < jmax; j++){
		sum += poptree.ind[poptree.subpop[c][j]].sus;
		sumst[j] = sum;
	}
	
	kmax = indinf[c].size();
	do{
	  z = ran()*sum;                                           // Samples in proportion to individual susceptibility
		j = 0; while(j < jmax && z > sumst[j]) j++; if(j == jmax) emsg("Part: EC1");
		i = poptree.subpop[c][j];
		
		for(k = 0; k < kmax; k++) if(indinf[c][k] == i) break;   // Checks selected individual is not infected
	}while(k < kmax);
	indinf[c].push_back(i);
	
	l = poptree.level-1; cc = c;
	sus = poptree.ind[i].sus;
	dR = -sus*Rtot[l][cc]/sussum[l][cc];
	do{
		Rtot[l][cc] += dR;
		sussum[l][cc] -= sus;
		cc = lev[l].node[cc].parent; l--;
	}while(l >= 0);
	
	model.simmodel(fev,tdnext,tdfnext,i,t,data.period,N);
}

/// Performs the modified Gillespie algorithm between times ti and tf 
void PART::gillespie(double ti, double tf, int outp)
{
	int td, j, c;
	double t, tpl;
	NEV n;
	vector <NEV> nev;
	
	if(sett == nsettime) emsg("Part: EC4");
	
	t = ti; tpl = t;
	do{
		nev.clear();                     // First we decide what event is next
		n.t = model.settime[sett];
		n.type = SET_EV;
		nev.push_back(n); 
		 
		if(tdnext < fediv) n.t = fev[tdnext][tdfnext].t; else n.t = tf;
		n.type = FEV_EV;
		nev.push_back(n);
	
		if(Rtot[0][0] < tiny) n.t = tf; else n.t = t - log(ran())/(model.beta[sett]*Rtot[0][0]);
		n.type = INF_EV;
		nev.push_back(n);
		
		n.t = t - log(ran())/(sussum[0][0]*model.paramval[model.phiparam]);
		n.type = EXT_EV;
		nev.push_back(n);
		
		sort(nev.begin(),nev.end(),compNEV);
		
		if(outp == 1){
			while(t > tpl){ 
				cout  << "Time: " << tpl;
				for(c =0; c < comp.size(); c++) cout << "  " << comp[c].name << ":"	<< N[c];
				cout << endl;
				tpl++;
			}
		}
	
		t = nev[0].t; if(t >= tf) break;

		switch(nev[0].type){
		case SET_EV:                 // These are "settime" events which allow the value of beta to change in time
			sett++; if(sett >= nsettime) emsg("Part: EC5");
			break;
		
		case INF_EV:                 // These are infection events within the system
 		case EXT_EV:                 // These are external infection events

	  	c = nextinfection(nev[0].type);
			addinfc(c,t);	
			break;
			
		case FEV_EV:                 // These correspond to other compartmental transitions (e.g. E->A, E->I etc...)
			dofe();
			break;
		
		default: emsg("Part: EC6"); break;
		}
		
		if(checkon == 1) check(0);
	}while(t < tf);
}

/// Used to check that various quantities are being correctly updated
void PART::check(int num)
{
	int l, c, cmax, cc;
	double dd;
	for(l = 0; l < poptree.level; l++){
		cmax = lev[l].add.size();
		for(c = 0; c < cmax; c++){ if(lev[l].add[c] != 0) emsg("Part: EC7");}
	}	
	
	double sum=0;
	for(cc = 0; cc < poptree.Cfine; cc++) sum += ffine[cc]*sussum[poptree.level-1][cc];    
	dd = Rtot[0][0] - sum;
	if(dd < -tiny || dd > tiny) emsg("Part: EC8");
}

/// Makes changes corresponding to a compartmental transition in one of the individuals
void PART::dofe()
{
	int i, c, cmax, cc, ccc, j, jmax, h, ii, k, kmax, l, ll;
	double fac, val, num, dd, ffnew;
	TRANS tr;

	int **&nMval(poptree.nMval);
	int ***&Mnoderef(poptree.Mnoderef);
	float ***&Mval(poptree.Mval);
	
	if(checkon == 1){
		if(trans[fev[tdnext][tdfnext].trans].type == INFECTION) emsg("Part: EC8a");
		
		for(l = 0; l < poptree.level; l++){
			cmax = lev[l].add.size();
			for(c = 0; c < cmax; c++){ if(lev[l].add[c] != 0) emsg("Part: EC9");}
		}	
		
		double sum=0;
		for(cc = 0; cc < poptree.Cfine; cc++) sum += ffine[cc]*sussum[poptree.level-1][cc];    
		dd = Rtot[0][0] - sum;
		if(dd < -tiny || dd > tiny)	emsg("Part: EC10");
	}
	
	if(fev[tdnext][tdfnext].done != 0) emsg("Part: EC11");
	
	i = fev[tdnext][tdfnext].ind; 
	fev[tdnext][tdfnext].done = 1;
	c = poptree.ind[i].noderef;

	tr = trans[fev[tdnext][tdfnext].trans];
	N[tr.from]--; if(N[tr.from] < 0){ cout << tr.from << " " <<  N[tr.from] << " fr\n"; emsg("Part: EC12"); }
	N[tr.to]++;
	
	fac = poptree.ind[i].inf*(comp[tr.to].infectivity - comp[tr.from].infectivity);

	tdfnext++;
	if(tdfnext == fev[tdnext].size()){
		tdnext++; tdfnext = 0; 
		while(tdnext < fediv && fev[tdnext].size() == 0) tdnext++;
	}
		
	if(fac == 0) return;

	/*
	h = poptree.ind[i].houseref;                               // Updates household effect
	jmax = poptree.house[h].ind.size();
	for(j = 0; j < jmax; j++){
		ii = poptree.house[h].ind[j];
		if(ii != i){
			kmax = indinf[c].size(); k = 0; while(k < kmax && indinf[c][k] != ii) k++;
			if(k == kmax){
				
			}
		}
	}
	*/
	
	if(checkon == 1){                           // These are checks to see if the algorithm is working properly
		for(l = poptree.level-1; l >= 0; l--){
			kmax = nMval[c][l];
			for(k = 0; k < kmax; k++){
				cc = Mnoderef[c][l][k];
				val = fac*Mval[c][l][k];
			
				num = val*sussum[l][cc];
				jmax = lev[l].node[cc].fine.size();
				for(j = 0; j < jmax; j++){
					ccc = lev[l].node[cc].fine[j];
					ffnew = ffine[ccc]+val; if(ffnew < 0){ if(ffnew < -tiny) emsg("Part: EC13"); ffnew = 0;}
					ffine[ccc] = ffnew;
				}
			}
		}
	}
		
	for(l = poptree.level-1; l >= 0; l--){              // This updates the different levels of Rtot (used for sampling later) 
		kmax = nMval[c][l];
		for(k = 0; k < kmax; k++){
			cc = Mnoderef[c][l][k];
			val = fac*Mval[c][l][k]*sussum[l][cc];
			lev[l].add[cc] = val;
			Rtot[l][cc] += val;
		}
		
		kmax = poptree.naddnoderef[c][l];
		for(k = 0; k < kmax; k++){
			cc = poptree.addnoderef[c][l][k];

			jmax = lev[l].node[cc].child.size();
			val = 0; for(j = 0; j < jmax; j++) val += lev[l+1].add[lev[l].node[cc].child[j]];
			lev[l].add[cc] = val;
			Rtot[l][cc] += val;
		}

		if(l < poptree.level-1){
			kmax = nMval[c][l];
			for(k = 0; k < kmax; k++){
				cc = Mnoderef[c][l][k];
				val = fac*Mval[c][l][k];
				addlater[l][cc] += val;
			}
		}
	}
	
	for(l = poptree.level-1; l >= 0; l--){
		kmax = nMval[c][l]; for(k = 0; k < kmax; k++) lev[l].add[Mnoderef[c][l][k]] = 0;
		kmax = poptree.naddnoderef[c][l]; for(k = 0; k < kmax; k++) lev[l].add[poptree.addnoderef[c][l][k]] = 0;
	}
		
	if(checkon == 1) check(1);
}

/// This samples the node on the fine scale in which the next infection occurs
int PART::nextinfection(int type)
{
	int l, lmax, c, cc, j, jmax;
	double z, sum, sumst[4], val, dd, Rnew;
	
	l = 0; c = 0;                              // We start at the top level l=0 and proceed to fine and finer scales
	lmax = poptree.level;
	while(l < lmax-1){
		val = addlater[l][c]; addlater[l][c] = 0;
	
		jmax = lev[l].node[c].child.size();
		sum = 0;
		for(j = 0; j < jmax; j++){
			cc = lev[l].node[c].child[j];
			
			if(val != 0){
				Rnew = Rtot[l+1][cc]+val*sussum[l+1][cc];
				if(Rnew < 0){ if(Rnew < -tiny) emsg("Part: EC14"); Rnew = 0;}
				Rtot[l+1][cc] = Rnew;
				if(l < lmax-2) addlater[l+1][cc] += val;
			}	
			if(type == INF_EV) sum += Rtot[l+1][cc];
			else sum += sussum[l+1][cc];
				
			sumst[j] = sum;
		}
		
		z = ran()*sum; j = 0; while(j < jmax && z > sumst[j]) j++;
		if(j == jmax) emsg("Part: EC15");
		
		c = lev[l].node[c].child[j];
		l++;
	};
	
	if(checkon == 1){
		dd = ffine[c]*sussum[l][c] - Rtot[l][c];
		if(dd < -tiny || dd > tiny) emsg("Part: EC16"); 
	}
	
	return c;
}

/// Packs up all the particle information (from time fedivmin until the end) to the be sent by MPI
void PART::partpack(int fedivmin)
{
	int l, c, cmax, cc, j, jmax;
	double val;

	for(l = 0; l < poptree.level; l++){        // Process all the add later requests
		cmax = lev[l].node.size();
		for(c = 0; c < cmax; c++){
			val = lev[l].add[c]; lev[l].add[c] = 0;
			Rtot[l][c] += val*sussum[l][c];
			if(l < poptree.level-1){
				val += addlater[l][c]; addlater[l][c] = 0;
				for(j = 0; j < 4; j++) lev[l+1].add[lev[l].node[c].child[j]] += val;
			}
		}			
	}

	l = poptree.level-1;
	
	packinit();
	if(checkon == 1) pack(ffine);
	pack(indinf);
	pack(Rtot[l]);
	pack(fev,fedivmin,fediv);
	pack(N);
	pack(tdnext); 
	pack(tdfnext);
	pack(sett);
}

/// Unpacks particle 
void PART::partunpack(int fedivmin)
{
	int k, kmax, j, l, c, cmax, cc;
	double val, val2;
	
	l = poptree.level-1;
	
	packinit();
	if(checkon == 1) unpack(ffine);
	unpack(indinf);
	unpack(Rtot[l]);
	unpack(fev,fedivmin,fediv);
	unpack(N);
	unpack(tdnext);
	unpack(tdfnext);
	unpack(sett);

	cmax = poptree.Cfine;
	for(c = 0; c < cmax; c++){
		val = lev[l].node[c].sussum;
		kmax = indinf[c].size();
		for(k = 0; k < kmax; k++) val -= poptree.ind[indinf[c][k]].sus;
		sussum[l][c] = val;
	}
		
	for(l = poptree.level-2; l >= 0; l--){
		cmax = lev[l].node.size();
		for(c = 0; c < cmax; c++){
			val = 0; val2 = 0;
			for(j = 0; j < 4; j++){
				cc = lev[l].node[c].child[j];
				val += sussum[l+1][cc];
				val2 += Rtot[l+1][cc];
			}
			sussum[l][c] = val;
			Rtot[l][c] = val2;
			addlater[l][c] = 0;
		}
	}

	if(checkon == 1){
		for(l = 0; l < poptree.level; l++){
			cmax = lev[l].node.size();
			for(c = 0; c < cmax; c++){
				if(addlater[l][c] != 0) emsg("Part: EC30");
				if(lev[l].add[c] != 0) emsg("Part: EC31");
			}
		}
	}
}

/// Copies in all the information from another particle
void PART::copy(const PART &other, int fedivmin)
{
	int d;
	
	if(checkon == 1) ffine = other.ffine;
	indinf = other.indinf;
	Rtot = other.Rtot; 
	sussum = other.sussum;
	addlater = other.addlater;
	for(d = fedivmin; d < fediv; d++) fev[d] = other.fev[d];
	N = other.N;
	tdnext = other.tdnext;
	tdfnext = other.tdfnext;
	sett = other.sett;
}
