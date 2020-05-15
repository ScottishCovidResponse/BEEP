
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
