
#include <iostream>
#include <algorithm>

using namespace std;

#include "math.h"

#include "timers.hh"
#include "utils.hh"
#include "PART.hh"

struct NEV {                               // Information about the immediate next events
  short type; double t;
};

struct FEV {                               // Stores information about a compartmental transition
  long trans;                              // References the transition type
	long ind;                                // The individual on which the transition happens
	double t;                                // The time of the transition
	short done;                              // Set to 1 if that transition is in the past 
};

static bool compNEV(NEV lhs, NEV rhs)
{
	return lhs.t < rhs.t;
};

PART::PART(MODEL &model, POPTREE &poptree) : model(model), comp(model.comp), trans(model.trans), poptree(poptree), lev(poptree.lev)
{
}

/// Initialises a particle
void PART::partinit(long p)
{
	long c, cmax, cc, k, kmax, h, i, imax, j, jmax, l, loop;
	
	pst = p;
	N.resize(comp.size()); for(c = 0; c < comp.size(); c++) N[c] = 0;
 
	ffine.clear(); 
	ffine.resize(poptree.Cfine); for(c = 0; c < poptree.Cfine; c++) ffine[c] = 0;

	indinf.clear(); indinf.resize(poptree.Cfine);
	
	fev.clear(); fev.resize(fediv);

	Rtot.resize(poptree.level); addlater.resize(poptree.level); sussum.resize(poptree.level);
	for(l = 0; l < poptree.level; l++){
		cmax = lev[l].node.size();
		Rtot[l].resize(cmax); addlater[l].resize(cmax); sussum[l].resize(cmax);
		for(c = 0; c < cmax; c++){ 
			Rtot[l][c] = 0; addlater[l][c] = 0;
			sussum[l][c] = lev[l].node[c].sussum;
		}
	}
	
	sett = 0;
	
	tdnext = fediv;
}

/// Copies in all the information from another particle
void PART::copy(const PART &other)
{
	short c;
	
	ffine = other.ffine;
	indinf = other.indinf;
	Rtot = other.Rtot; 
	sussum = other.sussum;
	addlater = other.addlater;
	fev = other.fev;
	for(c = 0; c < comp.size(); c++) N[c] = other.N[c];
	tdnext = other.tdnext;
	tdfnext = other.tdfnext;
	sett = other.sett;
}

/// Returns the number of transitions for individuals going from compartment "from" to compartment "to" 
/// in different regions over the time range ti - tf
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
				if(fe.t > ti && fe.trans == tra) num[poptree.ind[fe.ind].region]++;
			}
		}
	}
	
	return num;
}

/// Adds an exposed indivdual on node c on the finest scale (i.e. level-1)
void PART::addinfc(long c, double t)
{
	long l, i, j, jmax, cc, k, kmax;
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
		j = 0; while(j < jmax && z > sumst[j]) j++; if(j == jmax) emsg("Part: EC9");
		i = poptree.subpop[c][j];
		
		for(k = 0; k < kmax; k++) if(indinf[c][k] == i) break;   // Checks selected individual is not infected
	}while(k < kmax);
	indinf[c].push_back(i);
	
	kmax = indinf[c].size();
	do{
		l = long(ran()*long(poptree.subpop[c].size()));
		i = poptree.subpop[c][l];
		for(k = 0; k < kmax; k++) if(indinf[c][k] == i) break;
	}while(k < kmax);
	indinf[c].push_back(i);
	
	l = poptree.level-1; cc = c;
	sus = poptree.ind[i].sus;
	dR = -sus*Rtot[l][c]/sussum[l][cc];
	do{
		Rtot[l][cc] += dR;
		sussum[l][cc] -= sus;
		cc = lev[l].node[cc].parent; l--;
	}while(l >= 0);
	simmodel(i,0,t);
}

/// Draws a sample from the gamma distribution x^(a-1)*exp(-b*x)
static double gammasamp(double a, double b)
{
  if(a < 0 || b < 0) emsg("Model: EC1");

  if(a < 1){
    double u = ran();
    return gammasamp(1.0 + a, b) * pow (u, 1.0 / a);
  }
  else{
    double x, v, u;
    double d = a - 1.0 / 3.0;
    double c = (1.0 / 3.0) / sqrt (d);
 
    while(1 == 1){
      do{
        x = sqrt(-2*log(ran()))*cos(2*M_PI*ran());
        v = 1.0 + c * x;
      }while (v < 0);

      v = v*v*v;
      u = ran();

      if (u < 1 - 0.0331*x*x*x*x) break;

      if (log(u) < 0.5*x*x + d*(1 - v + log(v))) break;
    }

    return d*v/b;
  }
}

/// Once an individual goes into the exposed class, this function simulates all the subsequent future events
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
					tnext = t - log(ran())/model.param[tr.param1].val;
					break;
				
				case GAMMA_DIST:
					mean = model.param[tr.param1].val; sd = model.param[tr.param2].val;
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

/// Adds a future event to the timeline
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

/// Performs the modified Gillespie algorithm between times ti and tf 
void PART::gillespie(double ti, double tf, short siminf)
{
	long td, j, c, NIfine[poptree.Cfine];
	double t, tpl;
	NEV n;
	vector <NEV> nev;
	
	if(sett == nsettime) emsg("Simulate: EC1");
	
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
		
		n.t = t - log(ran())/(sussum[0][0]*model.param[model.phiparam].val);
		n.type = EXT_EV;
		nev.push_back(n);
		
		sort(nev.begin(),nev.end(),compNEV);
		
		if(siminf == 1){
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
			sett++; if(sett >= nsettime) emsg("Simulate: EC1a");
			break;
		
		case INF_EV:                 // These are infection events
			c = nextinfection();
			addinfc(c,t);	
			break;
			
		case FEV_EV:                 // These correspond to other compartmental transitions (e.g. E->A, E->I etc...)
			dofe();
			break;
			
		case EXT_EV:
			c = externalinfection();
			addinfc(c,t);	
			break;
			
		default: emsg("Simulate: EC2"); break;
		}
	}while(t < tf);
}

/// Makes changes corresponding to a compartmental transition in one of the individuals
void PART::dofe()
{
	long i, c, cmax, cc, ccc, j, jmax, k, kmax, l, ll;
	double fac, val, num, dd, ffnew;
	TRANS tr;

	long **&nMval(poptree.nMval);
	long ***&Mnoderef(poptree.Mnoderef);
	float ***&Mval(poptree.Mval);
	
	i = fev[tdnext][tdfnext].ind; if(fev[tdnext][tdfnext].done != 0) emsg("Simulate: EC3");
	fev[tdnext][tdfnext].done = 1;
	c = poptree.ind[i].noderef;

	tr = trans[fev[tdnext][tdfnext].trans];
	N[tr.from]--; if(N[tr.from] < 0) emsg("Simulate: EC4"); 
	N[tr.to]++;
	
	fac = poptree.ind[i].inf*(comp[tr.to].infectivity - comp[tr.from].infectivity);

	tdfnext++;
	if(tdfnext == fev[tdnext].size()){
		tdnext++; tdfnext = 0; 
		while(tdnext < fediv && fev[tdnext].size() == 0) tdnext++;
	}
		
	if(fac == 0) return;

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
					ffnew = ffine[ccc]+val; if(ffnew < 0){ if(ffnew < -tiny) emsg("Simulate: EC5"); ffnew = 0;}
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
		
	if(checkon == 1){
		for(l = 0; l < poptree.level; l++){
			cmax = lev[l].add.size();
			for(c = 0; c < cmax; c++){ if(lev[l].add[c] != 0) emsg("Simulate: EC6");}
		}	
		
		double sum=0;
		for(cc = 0; cc < poptree.Cfine; cc++) sum += ffine[cc]*sussum[poptree.level-1][cc];    
		dd = Rtot[0][0] - sum;
		if(dd < -tiny || dd > tiny)	emsg("Simulate: EC7");
	}
}

long PART::externalinfection()       // An infection is caused by something external to the system (i.e. from abroad)
{
	long l, c, cc, j, jmax;
	double z, sum, sumst[4];
	
	l = 0; c = 0;                              // We start at the top level l=0 and proceed to fine and finer scales
	while(l < poptree.level-1){
		jmax = lev[l].node[c].child.size();
		sum = 0;
		for(j = 0; j < jmax; j++){
			cc = lev[l].node[c].child[j];
			sum += sussum[l+1][cc];
			sumst[j] = sum;
		}
		
		z = ran()*sum; j = 0; while(j < jmax && z > sumst[j]) j++;
		if(j == jmax) emsg("Simulate: EC9");
		
		c = lev[l].node[c].child[j];
		l++;
	};
	
	return c;
}
	
/// This samples the node on the fine scale in which the next infection occurs
long PART::nextinfection()
{
	long l, c, cc, j, jmax;
	double z, sum, sumst[4], val, dd, Rnew;
	
	l = 0; c = 0;                              // We start at the top level l=0 and proceed to fine and finer scales
	while(l < poptree.level-1){
		val = addlater[l][c]; addlater[l][c] = 0;
	
		jmax = lev[l].node[c].child.size();
		sum = 0;
		for(j = 0; j < jmax; j++){
			cc = lev[l].node[c].child[j];
			
			if(val != 0){
				Rnew = Rtot[l+1][cc]+val*sussum[l+1][cc]; if(Rnew < 0){ if(Rnew < -tiny) emsg("Simulate: EC8"); Rnew = 0;}
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
		dd = ffine[c]*sussum[l][c] - Rtot[l][c];
		if(dd < -tiny || dd > tiny) emsg("Simulate: EC10"); 
	}
	
	return c;
}

/// Measures how well the particle agrees with the observations within a given time period
/// (which in this case is weekly hospitalised case data)
void PART::Lobs(short ti, short tf, long ncase[nregion][tmax/7+1])
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
