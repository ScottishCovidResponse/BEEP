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
	int c, cmax, j, jmax, l, dp, a;
	double sum, val;

	fediv = data.fediv;
		
	pst = p;
	N.resize(comp.size()); for(c = 0; c < comp.size(); c++) N[c] = 0;
	N[0] = data.popsize;
 
	indinf.clear(); indinf.resize(data.narea);
	
	fev.clear(); fev.resize(fediv);

	Rtot.resize(poptree.level); sussum.resize(poptree.level);
	for(l = 0; l < poptree.level; l++){
		cmax = lev[l].node.size();
		Rtot[l].resize(cmax); sussum[l].resize(cmax);
		for(c = 0; c < cmax; c++) Rtot[l][c] = 0;
	}
	
	l = poptree.level-1;
	susage.resize(data.narea);
	Qmap.resize(data.narea);
	for(c = 0; c < data.narea; c++){ 
		susage[c].resize(data.nage); for(a = 0; a < data.nage; a++) susage[c][a] = 0;
		Qmap[c].resize(data.nage); for(a = 0; a < data.nage; a++) Qmap[c][a] = 0;
			
		sum = 0; 
		for(dp = 0; dp < data.ndemocatpos; dp++){
			val = model.sus[dp]*data.area[c].pop[dp];
			susage[c][data.democatpos[dp][0]] += val;
			sum += val;
		}
		sussum[l][c] = sum;
	}
	
	for(l = poptree.level-2; l >= 0; l--){                                                // Propages sussum up the tree
		cmax = poptree.lev[l].node.size();
		for(c = 0; c < cmax; c++){
			jmax =  poptree.lev[l].node[c].child.size(); 
			sum = 0; for(j = 0; j < jmax; j++) sum += sussum[l+1][poptree.lev[l].node[c].child[j]];
			sussum[l][c] = sum;
		}
	}
	sussumst = sussum;
	susagest = susage;
	
	sett = 0;

	tdnext = fediv;
	
	indev.resize(data.popsize);
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
 		case EXT_EV:                 // These are external infection event
	  	c = nextinfection(nev[0].type);
			addinfc(c,t);	
			break;
			
		case FEV_EV:                 // These correspond to other compartmental transitions (e.g. E->A, E->I etc...)
			dofe();
			break;
		
		default: emsg("Part: EC6"); break;
		}		 
	}while(t < tf);
	
	if(checkon == 1) check(0,t);
}

/// Adds an exposed indivdual in area c
void PART::addinfc(int c, double t)
{
	int l, i, dp, j, jmax, k, kmax, a;
	double dR, sum, sus, z;
	vector <double> sumst;
	vector <FEV> evlist;
	
	sum = 0; sumst.resize(data.ndemocatpos);
	for(dp = 0; dp < data.ndemocatpos; dp++){ 
		sum += data.area[c].pop[dp]*model.sus[dp];
		sumst[dp] = sum;
	}
	
	kmax = indinf[c].size();
	do{
	  z = ran()*sum;                                           // Samples in proportion to individual susceptibility
		dp = 0; while(dp < data.ndemocatpos && z > sumst[dp]) dp++; if(dp == data.ndemocatpos) emsg("Part: EC1");
		i = data.area[c].ind[dp][int(ran()*data.area[c].pop[dp])];
		
		for(k = 0; k < kmax; k++) if(indinf[c][k] == i) break;   // Checks selected individual is not infected
	}while(k < kmax);
	indinf[c].push_back(i);
	
	sus = model.sus[dp];
	
	a = data.democatpos[dp][0];
	susage[c][a] -= sus;
	dR = -sus*Qmap[c][a];
	
	l = poptree.level-1; 
	do{
		Rtot[l][c] += dR;
		sussum[l][c] -= sus;
		c = lev[l].node[c].parent; l--;
	}while(l >= 0);
	
	model.simmodel(indev[i],i,0,t);
	jmax = indev[i].size(); for(j = 0; j < jmax; j++) addfev(indev[i][j],data.period,t);
}

/// Used to check that various quantities are being correctly updated
void PART::check(int num, double t)
{
	int l, c, cmax, cc, k, j, dp, i, a, aa, v, timep, q;
	double dd, sum, sum2, val, inf;
	vector <double> susag;
	vector <vector <double> > Qma;
	
	susag.resize(data.nage);
	
	for(c = 0; c < data.narea; c++){
		for(a = 0; a < data.nage; a++) susag[a] = 0;
		 
		sum = 0; 
		for(dp = 0; dp < data.ndemocatpos; dp++){
			val = model.sus[dp]*data.area[c].pop[dp];
			susag[data.democatpos[dp][0]] += val;
			sum += val;
		}
		
		for(j = 0; j < indinf[c].size(); j++){
			i = indinf[c][j];
			dp = data.ind[i].dp;
			
			sum -= model.sus[dp];
			susag[data.democatpos[dp][0]] -= model.sus[dp];
		}		
			
		dd = sum - sussum[poptree.level-1][c]; if(dd*dd > tiny) emsg("Part: EC20");
		
		for(dp = 0; dp < data.nage; dp++){
			dd = susag[dp] - susage[c][dp]; 
			if(dd*dd > tiny) emsg("Part: EC21");
		}
	}
		
	Qma.resize(data.narea);
	for(c = 0; c < data.narea; c++){ Qma[c].resize(data.nage); for(a = 0; a < data.nage; a++) Qma[c][a] = 0;}
	
	timep = 0; while(timep < data.ntimeperiod && t > data.timeperiod[timep]) timep++;
	
	for(c = 0; c < data.narea; c++){
		for(j = 0; j < indinf[c].size(); j++){
			i = indinf[c][j];
			k = 0; cc = 0; while(k < indev[i].size() && t >= indev[i][k].t){ cc = trans[indev[i][k].trans].to; k++;}
			
			q = 0; while(q < data.Qnum && !(data.Qcomp[q] == comp[cc].name && data.Qtimeperiod[q] == timep)) q++;
			if(q < data.Qnum){ 			
				inf = comp[cc].infectivity;
				
				dp = data.ind[i].dp;
				a = data.democatpos[dp][0];
				v = c*data.nage + a;
				for(k = 0; k < data.nQ[q][v]; k++){
					cc = data.Qto[q][v][k];
					for(aa = 0; aa < data.nage; aa++){
						Qma[cc][aa] += inf*data.Qval[q][v][k][aa];
					}
				}
			}
		}
	}
		
	for(c = 0; c < data.narea; c++){
		for(a = 0; a < data.nage; a++){
			dd = Qma[c][a] - Qmap[c][a]; if(dd*dd > tiny){ cout << Qma[c][a] << " " << Qmap[c][a] << "\n";  emsg("Part: EC22");}
		}
	}
	
	sum = 0; 
	for(c = 0; c < data.narea; c++){
		sum2 = 0; for(a = 0; a < data.nage; a++) sum2 += Qma[c][a]*susage[c][a];	
		dd = sum2 - Rtot[poptree.level-1][c]; if(dd*dd > tiny) emsg("Part: EC22b");
		sum += sum2;
	}
	
	for(l = 0; l < poptree.level; l++){
		cmax = lev[l].add.size();
		sum2 = 0; for(c = 0; c < cmax; c++) sum2 += Rtot[l][c];
		dd = sum - sum2; if(dd*dd > tiny) emsg("Part: EC23");
	}
	
	for(l = 0; l < poptree.level; l++){
		cmax = lev[l].add.size();
		for(c = 0; c < cmax; c++){ if(lev[l].add[c] != 0) emsg("Part: EC7");}
	}	
}

/// Makes changes corresponding to a compartmental transition in one of the individuals
void PART::dofe()
{
	int i, c, cc, ccc, dp, v, a, l, k, kmax, j, jmax, q, timep;
	double sum, val, t;
	TRANS tr;
		 
	if(fev[tdnext][tdfnext].done != 0) emsg("Part: EC11");
	
	i = fev[tdnext][tdfnext].ind; 
	t = fev[tdnext][tdfnext].t; 
	fev[tdnext][tdfnext].done = 1;
	c = data.ind[i].area;
		 
	tr = trans[fev[tdnext][tdfnext].trans];
	N[tr.from]--; if(N[tr.from] < 0){ cout << tr.from << " " <<  N[tr.from] << " fr\n"; emsg("Part: EC12"); }
	N[tr.to]++;

	tdfnext++;
	if(tdfnext == fev[tdnext].size()){
		tdnext++; tdfnext = 0; 
		while(tdnext < fediv && fev[tdnext].size() == 0) tdnext++;
	}
		
	timep = 0; while(timep < data.ntimeperiod && t > data.timeperiod[timep]) timep++;
	if(timep >= tr.DQ.size()) emsg("Part: EC66");
	q = tr.DQ[timep]; if(q == -1) return;
	
	dp = data.ind[i].dp;
	v = c*data.nage+data.democatpos[dp][0];
	
	for(l = 0; l < poptree.level; l++) lev[l].donelist.clear();
	
	l = poptree.level-1;                                                             // Makes change to Rtot
	kmax = model.nDQ[q][v];
	for(k = 0; k < kmax; k++){
		cc = model.DQto[q][v][k];
		
		sum = 0; 
		for(a = 0; a < data.nage; a++){
			val = model.DQval[q][v][k][a];
			Qmap[cc][a] += val;
			sum += val*susage[cc][a];
		}
		Rtot[l][cc] += sum;
		
		ccc = lev[l].node[cc].parent;
		if(lev[l-1].add[ccc] == 0){ lev[l-1].donelist.push_back(ccc);}		
		lev[l-1].add[ccc] += sum;
	}
	
	for(l = poptree.level-2; l >= 0; l--){                                        // Propages change up the tree
		jmax = lev[l].donelist.size();
		for(j = 0; j < jmax; j++){
			c = lev[l].donelist[j];
			
			sum = lev[l].add[c];
			Rtot[l][c] += sum;
			lev[l].add[c] = 0;
			 
			if(l > 0){
				cc = lev[l].node[c].parent;
				if(lev[l-1].add[cc] == 0){ lev[l-1].donelist.push_back(cc);}		
				lev[l-1].add[cc] += sum;
			}
		}
	}
}

/// This samples the node on the fine scale in which the next infection occurs
int PART::nextinfection(int type)
{
	int l, lmax, c, cc, j, jmax;
	double z, sum, sumst[4], dd, Rnew;
	
	l = 0; c = 0;                              // We start at the top level l=0 and proceed to fine and finer scales
	lmax = poptree.level;
	while(l < lmax-1){
		jmax = lev[l].node[c].child.size();
		sum = 0;
		for(j = 0; j < jmax; j++){
			cc = lev[l].node[c].child[j];
	
			if(type == INF_EV) sum += Rtot[l+1][cc];
			else sum += sussum[l+1][cc];
				
			sumst[j] = sum;
		}
		
		z = ran()*sum; j = 0; while(j < jmax && z > sumst[j]) j++;
		if(j == jmax) emsg("Part: EC15");
		
		c = lev[l].node[c].child[j];
		l++;
	};
	
	return c;
}

/// Packs up all the particle information (from time fedivmin until the end) to the be sent by MPI
void PART::partpack(int fedivmin)
{
	int l, c, cmax, cc, j, jmax;

	l = poptree.level-1;
	
	packinit();
	pack(indinf);
	pack(Qmap);
	pack(fev,fedivmin,fediv);
	pack(N);
	pack(tdnext); 
	pack(tdfnext);
}

/// Unpacks particle 
void PART::partunpack(int fedivmin)
{
	int k, kmax, j, jmax, l, c, cmax, cc, dp, a;
	double val, val2, sus, sum;
	
	l = poptree.level-1;
	
	packinit();
	unpack(indinf);
	unpack(Qmap);
	unpack(fev,fedivmin,fediv);
	unpack(N);
	unpack(tdnext); 
	unpack(tdfnext);

	sussum = sussumst;
	susage = susagest;
	
	cmax = data.narea;
	for(c = 0; c < cmax; c++){
		val = sussumst[l][c]; 
		kmax = indinf[c].size();
		for(k = 0; k < kmax; k++){
			dp = data.ind[indinf[c][k]].dp;
			sus = model.sus[dp];
			sussum[l][c] -= sus;
			susage[c][data.democatpos[dp][0]] -= sus;
		}
		
		sum = 0; for(a = 0; a < data.nage; a++) sum += susage[c][a]*Qmap[c][a];
		Rtot[l][c] = sum;
	}
		
	for(l = poptree.level-2; l >= 0; l--){
		cmax = lev[l].node.size();
		for(c = 0; c < cmax; c++){
			val = 0; val2 = 0;
			jmax = lev[l].node[c].child.size();
			for(j = 0; j < jmax; j++){
				cc = lev[l].node[c].child[j];
				val += sussum[l+1][cc];
				val2 += Rtot[l+1][cc];
			}
			sussum[l][c] = val;
			Rtot[l][c] = val2;
		}
	}
}

/// Copies in all the information from another particle
void PART::copy(const PART &other, int fedivmin)
{
	int d;
	
	indinf = other.indinf;
	Rtot = other.Rtot; 
	Qmap = other.Qmap; 
	sussum = other.sussum;
	susage = other.susage;
	for(d = fedivmin; d < fediv; d++) fev[d] = other.fev[d];
	N = other.N;
	tdnext = other.tdnext;
	tdfnext = other.tdfnext;
}

/// Adds a future event to the timeline
void PART::addfev(FEV fe, double period, double tnow)
{
	int d, j, jmax;
	double t;
	
	t = fe.t; if(t < tnow) emsg("MBPCHAIN: EC10");
	if(t >= period) return;
	
	d = int((t/period)*fev.size());
	j = 0; jmax = fev[d].size();
	if(t != tnow){ while(j < jmax && t >= fev[d][j].t) j++;}
	else{ while(j < jmax && t > fev[d][j].t) j++;}
	
	if(j == jmax) fev[d].push_back(fe);
	else fev[d].insert(fev[d].begin()+j,fe);
	
	if(t != tnow){
		if(d == tdnext){ if(j < tdfnext) tdfnext = j;}
		if(d < tdnext){ tdnext = d; tdfnext = j;}
	}
	else{
		TRANS tr = trans[fe.trans];
		N[tr.from]--; if(N[tr.from] < 0) emsg("Part: EC12"); 
		N[tr.to]++;
		
		if(d == tdnext) tdfnext++;
	}
}
