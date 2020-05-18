
#include <iostream>
#include <algorithm>

using namespace std;

#include "math.h"

#include "types.hh"
#include "var.hh"
#include "functions.hh"
#include "PART.hh"


// Initialises a particle
void PART::partinit(long p)
{
	long c, cmax, cc, k, kmax, h, i, imax, j, jmax, l, popu, loop;
	
	pst = p;
	N.resize(comp.size()); for(c = 0; c < comp.size(); c++) N[c] = 0;
 
	ffine.clear(); 
	ffine.resize(Cfine); for(c = 0; c < Cfine; c++) ffine[c] = 0;

	indinf.clear(); indinf.resize(Cfine);
	
	fev.clear(); fev.resize(fediv);

	Rtot.resize(level); addlater.resize(level); pop.resize(level);
	for(l = 0; l < level; l++){
		cmax = lev[l].node.size();
		Rtot[l].resize(cmax); addlater[l].resize(cmax); pop[l].resize(cmax);
		for(c = 0; c < cmax; c++){ Rtot[l][c] = 0; addlater[l][c] = 0; pop[l][c] = lev[l].node[c].popu;}
	}
	
	sett = 0;
	
	tdnext = fediv;
	// For simplicity we assume three randomly distributed initally exposed individuals
	// This will be changed in the proper analysis
	for(loop = 0; loop < 3; loop++){
		do{ c = long(ran()*Cfine);}while(long(subpop[c].size()) - long(indinf[c].size()) == 0);
		addinfc(c,0);
	}
}

// Copies in all the information from another particle
void PART::copy(long pfrom)
{
	short c;
	
	ffine = part[pfrom]->ffine;
	indinf = part[pfrom]->indinf;
	Rtot = part[pfrom]->Rtot; 
	pop = part[pfrom]->pop;
	addlater = part[pfrom]->addlater;
	fev = part[pfrom]->fev;
	for(c = 0; c < comp.size(); c++) N[c] = part[pfrom]->N[c];
	tdnext = part[pfrom]->tdnext;
	tdfnext = part[pfrom]->tdfnext;
	sett = part[pfrom]->sett;
}

// Returns the number of transitions for individuals going from compartment "from" to compartment "to" 
// in different regions over the time range ti - tf
vector <long> PART::getnumtrans(string from, string to, short ti, short tf)
{
	long d, k, r, tra;
	FEV fe;
	vector <long> num;
	
	tra = 0; while(tra < trans.size() && !(comp[trans[tra].from].name == from && comp[trans[tra].to].name == to)) tra++;
	if(tra == trans.size()) emsg("Finescale: Cannot find transition");
	
	for(r = 0; r < nregion; r++) num.push_back(0);

	for(d = long(fediv*double(ti)/tmax); d <= long(fediv*double(tf)/tmax); d++){
		if(d < fediv){
			for(k = 0; k < fev[d].size(); k++){
				fe = fev[d][k];
				if(fe.t > tf) break;
				if(fe.t > ti && fe.trans == tra) num[ind[fe.ind].region]++;
			}
		}
	}
	
	return num;
}

// Adds an exposed indivdual on node c on the finest scale (i.e. level-1)
void PART::addinfc(long c, double t)
{
	long l, i, cc, k, kmax;
	double dR, sum;
	
	kmax = indinf[c].size();
	do{
		l = long(ran()*long(subpop[c].size()));
		i = subpop[c][l];
		for(k = 0; k < kmax; k++) if(indinf[c][k] == i) break;
	}while(k < kmax);
	indinf[c].push_back(i);
	
	l = level-1; cc = c;
	dR = -Rtot[l][c]/pop[l][cc];
	do{
		Rtot[l][cc] += dR;
		pop[l][cc]--;
		cc = lev[l].node[cc].parent; l--;
	}while(l >= 0);
	simmodel(i,0,t);
}



// Once an individual goes into the exposed class, this function simulates all the subsequent future events
void PART::simmodel(long i, short enter, double t)
{
	short c, k, kmax, tra;
	double tnext, tnextmin, mean, sd;
	TRANS tr;
	
	N[enter]++;

	double dt;
	
	c = enter;
	do{
		kmax = comp[c].trans.size();
		if(kmax == 0) break;
		
		tnextmin = large;
		for(k = 0; k < kmax; k++){
			tr = trans[comp[c].trans[k]];
			switch(tr.type){
				case EXP_DIST:
					tnext = t - log(ran())/param[tr.param1].val;
					break;
				
				case GAMMA_DIST:
					mean = param[tr.param1].val; sd = param[tr.param2].val;
					dt = gammasamp(mean*mean/(sd*sd),mean/(sd*sd));
					tnext = t + dt; 
					break;
			}
			
			if(tnext < tnextmin){ tnextmin = tnext; tra = comp[c].trans[k];}
		}
		if(tnextmin == large) emsg("Model: EC2");
		
		addfev(tnextmin,tra,i);

		c = trans[tra].to; t = tnextmin;
	}while(1 == 1);
}

// Adds a future event to the timeline
void PART::addfev(double t, long trans, long i)
{
	long d, j, jmax;
	
	if(t >= tmax) return;
	
	FEV fe; fe.t = t; fe.trans = trans; fe.ind = i; fe.done = 0;
	
	d = long((t/tmax)*fediv);
	j = 0; jmax = fev[d].size();
	while(j < jmax && t > fev[d][j].t) j++;
	if(j == jmax) fev[d].push_back(fe);
	else fev[d].insert(fev[d].begin()+j,fe);
	
	if(d == tdnext){ if(j < tdfnext) tdfnext = j;}
	if(d < tdnext){ tdnext = d; tdfnext = j;}
}

// Performs the modified Gillespie algorithm between times ti and tf 
void PART::gillespie(double ti, double tf, short siminf)
{
	long td, j, c, NIfine[Cfine];
	double t, tpl;
	NEV n;
	vector <NEV> nev;
	
	if(sett == nsettime) emsg("Simulate: EC1");
	
	t = ti; tpl = t;
	do{
		nev.clear();                     // First we decide what event is next
		n.t = settime[sett];
		n.type = SET_EV;
		nev.push_back(n); 
		 
		if(tdnext < fediv) n.t = fev[tdnext][tdfnext].t; else n.t = tf;
		n.type = FEV_EV;
		nev.push_back(n);
	
		if(Rtot[0][0] < tiny) n.t = tf; else n.t = t - log(ran())/(beta[sett]*Rtot[0][0]);
		n.type = INF_EV;
		nev.push_back(n);
		
		sort(nev.begin(),nev.end(),compNEV);
		
		if(siminf == 1){
			while(t > tpl){ 
				cout  << "Time: " << tpl;
				for(c =0; c < comp.size(); c++) cout << "  " << comp[c].name << ":"	<< N[c];
				cout << "\n";
				tpl++;
			}
		}
	
		t = nev[0].t; if(t >= tf) break;
	
		switch(nev[0].type){
		case SET_EV:                 // These are "settime" events which allow the value of beta to change in time
			sett++; if(sett >= nsettime) emsg("Simulate: EC1a");
			break;
		
		case INF_EV:                 // These are infection events
			c = nextinfection();
			addinfc(c,t);	
			break;
			
		case FEV_EV:                 // These correspond to other compartmental transitions (e.g. E->A, E->I etc...)
			dofe();
			break;
			
		default: emsg("Simulate: EC2"); break;
		}
	}while(t < tf);
}

// Makes changes corresponding to a compartmental transition in one of the individuals
void PART::dofe()
{
	long i, c, cmax, cc, ccc, j, jmax, k, kmax, l, ll;
	double fac, val, num, dd, ffnew;
	TRANS tr;

	i = fev[tdnext][tdfnext].ind; if(fev[tdnext][tdfnext].done != 0) emsg("Simulate: EC3");
	fev[tdnext][tdfnext].done = 1;
	c = ind[i].noderef;

	tr = trans[fev[tdnext][tdfnext].trans];
	N[tr.from]--; if(N[tr.from] < 0) emsg("Simulate: EC4"); 
	N[tr.to]++;
	
	fac = comp[tr.to].infectivity - comp[tr.from].infectivity;

	tdfnext++;
	if(tdfnext == fev[tdnext].size()){
		tdnext++; tdfnext = 0; 
		while(tdnext < fediv && fev[tdnext].size() == 0) tdnext++;
	}
		
	if(fac == 0) return;
	
	if(checkon == 1){                           // These are checks to see if the algorithm is working properly
		for(l = level-1; l >= 0; l--){
			kmax = nMval[c][l];
			for(k = 0; k < kmax; k++){
				cc = Mnoderef[c][l][k];
				val = fac*Mval[c][l][k];
			
				num = val*pop[l][cc];
				jmax = lev[l].node[cc].fine.size();
				for(j = 0; j < jmax; j++){
					ccc = lev[l].node[cc].fine[j];
					ffnew = ffine[ccc]+val; if(ffnew < 0){ if(ffnew < -tiny) emsg("Simulate: EC5"); ffnew = 0;}
					ffine[ccc] = ffnew;
				}
			}
		}
	}
		
	for(l = level-1; l >= 0; l--){              // This updates the different levels of Rtot (used for sampling later) 
		kmax = nMval[c][l];
		for(k = 0; k < kmax; k++){
			cc = Mnoderef[c][l][k];
			val = fac*Mval[c][l][k]*pop[l][cc];
			lev[l].add[cc] = val;
			Rtot[l][cc] += val;
		}
		
		kmax = naddnoderef[c][l];
		for(k = 0; k < kmax; k++){
			cc = addnoderef[c][l][k];
			jmax = lev[l].node[cc].child.size();
			val = 0; for(j = 0; j < jmax; j++) val += lev[l+1].add[lev[l].node[cc].child[j]];
			lev[l].add[cc] = val;
			Rtot[l][cc] += val;
		}

		if(l < level-1){
			kmax = nMval[c][l];
			for(k = 0; k < kmax; k++){
				cc = Mnoderef[c][l][k];
				val = fac*Mval[c][l][k];
				addlater[l][cc] += val;
			}
		}
	}
	
	for(l = level-1; l >= 0; l--){
		kmax = nMval[c][l]; for(k = 0; k < kmax; k++) lev[l].add[Mnoderef[c][l][k]] = 0;
		kmax = naddnoderef[c][l]; for(k = 0; k < kmax; k++) lev[l].add[addnoderef[c][l][k]] = 0;
	}
		
	if(checkon == 1){
		for(l = 0; l < level; l++){
			cmax = lev[l].add.size();
			for(c = 0; c < cmax; c++){ if(lev[l].add[c] != 0) emsg("Simulate: EC6");}
		}	
		
		double sum=0;
		for(cc = 0; cc < Cfine; cc++) sum += ffine[cc]*pop[level-1][cc];    
		dd = Rtot[0][0] - sum;
		if(dd < -tiny || dd > tiny)	emsg("Simulate: EC7");
	}
}

// This samples the node on the fine scale in which the next infection occurs
long PART::nextinfection()
{
	long l, c, cc, j, jmax;
	double z, sum, sumst[4], val, dd, Rnew;
	
	l = 0; c = 0;                              // We start at the top level l=0 and proceed to fine and finer scales
	while(l < level-1){
		val = addlater[l][c]; addlater[l][c] = 0;
	
		jmax = lev[l].node[c].child.size();
		sum = 0;
		for(j = 0; j < jmax; j++){
			cc = lev[l].node[c].child[j];
			
			if(val != 0){
				Rnew = Rtot[l+1][cc]+val*pop[l+1][cc]; if(Rnew < 0){ if(Rnew < -tiny) emsg("Simulate: EC8"); Rnew = 0;}
				Rtot[l+1][cc] = Rnew;
				addlater[l+1][cc] += val;
			}	
			sum += Rtot[l+1][cc];
		
			sumst[j] = sum;
		}
		
		z = ran()*sum; j = 0; while(j < jmax && z > sumst[j]) j++;
		if(j == jmax) emsg("Simulate: EC9");
		
		c = lev[l].node[c].child[j];
		l++;
	};
	
	if(checkon == 1){
		dd = ffine[c]*pop[l][c] - Rtot[l][c];
		if(dd < -tiny || dd > tiny) emsg("Simulate: EC10"); 
	}
	
	return c;
}

// Measures how well the particle agrees with the observations within a given time period
// (which in this case is weekly hospitalised case data)
void PART::Lobs(short ti, short tf)
{
	short tt, r;
	double mean, var;
	vector <long> num;

	Li = 0;	
	for(tt = ti; tt < tf; tt += 7){
		num = getnumtrans("I","H",tt,tt+7);
		for(r = 0; r < nregion; r++){
			mean = ncase[r][tt/7];
			var = mean; if(var < 5) var = 5;
			Li += invT*(-0.5*log(2*3.141592654*var) - (mean-num[r])*(mean-num[r])/(2*var));
		}
	}
}
