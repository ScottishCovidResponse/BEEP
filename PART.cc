
#include <iostream>
#include <algorithm>

using namespace std;

#include "math.h"

#include "timers.hh"
#include "utils.hh"
#include "PART.hh"
#include "output.hh"

struct NEV {                               // Information about the immediate next events
  short type; double t;
};

static bool compNEV(NEV lhs, NEV rhs)
{
	return lhs.t < rhs.t;
};

PART::PART(DATA &data, MODEL &model, POPTREE &poptree) : data(data), model(model), comp(model.comp), trans(model.trans), poptree(poptree), lev(poptree.lev)
{
}

/// Initialises a particle
void PART::partinit(long p)
{
	long c, cmax, cc, k, kmax, h, i, imax, j, jmax, l, loop;
	
	pst = p;
	N.resize(comp.size()); for(c = 0; c < comp.size(); c++) N[c] = 0;
 
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
	
	simmodel(i,0,t);
	
}

/// Draws a sample from the gamma distribution x^(a-1)*exp(-b*x)
static double gammasamp(double a, double b)
{
  if(a < 0 || b < 0) emsg("Part: EC2");

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
		if(tnextmin == large) emsg("Part: EC3");
		
		addfev(tnextmin,tra,i);

		c = trans[tra].to; t = tnextmin;
	}while(1 == 1);
}

/// Adds a future event to the timeline
void PART::addfev(double t, long trans, long i)
{
	long d, j, jmax;
	
	if(t >= data.tmax) return;
	
	FEV fe; fe.t = t; fe.trans = trans; fe.ind = i; fe.done = 0;
	
	d = long((t/data.tmax)*fediv);
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

void PART::check(short num)
{
	long l, c, cmax, cc;
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
	long i, c, cmax, cc, ccc, j, jmax, h, ii, k, kmax, l, ll;
	double fac, val, num, dd, ffnew;
	TRANS tr;

	long **&nMval(poptree.nMval);
	long ***&Mnoderef(poptree.Mnoderef);
	float ***&Mval(poptree.Mval);
	
	if(checkon == 1){
		for(l = 0; l < poptree.level; l++){
			cmax = lev[l].add.size();
			for(c = 0; c < cmax; c++){ if(lev[l].add[c] != 0) emsg("Part: EC9");}
		}	
		
		double sum=0;
		for(cc = 0; cc < poptree.Cfine; cc++) sum += ffine[cc]*sussum[poptree.level-1][cc];    
		dd = Rtot[0][0] - sum;
		if(dd < -tiny || dd > tiny)	emsg("Part: EC10");
	}
	
	
	i = fev[tdnext][tdfnext].ind; if(fev[tdnext][tdfnext].done != 0) emsg("Part: EC11");
	fev[tdnext][tdfnext].done = 1;
	c = poptree.ind[i].noderef;

	tr = trans[fev[tdnext][tdfnext].trans];
	N[tr.from]--; if(N[tr.from] < 0) emsg("Part: EC12"); 
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
long PART::nextinfection(short type)
{
	long l, lmax, c, cc, j, jmax;
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

/// Measures how well the particle agrees with the observations within a given time period
/// (which in this case is weekly hospitalised case data)
void PART::Lobs(short w, double invT)
{
	short r;
	double mean, var;
	vector <long> num;
	
	Li = 0;	
	num = getnumtrans(data,model,poptree,fev,"I","H",w*timestep,(w+1)*timestep);
	for(r = 0; r < data.nregion; r++){
		mean = data.ncase[r][w];
		var = mean; if(var < 5) var = 5;
		//Li += -0.5*log(2*3.141592654*var) - (mean-num[r])*(mean-num[r])/(2*var);
		Li +=  - invT*(mean-num[r])*(mean-num[r])/(2*var);
	}
}

/// Copies in all the information from another particle

void PART::pack(vector <double> &pac, long num)
{
	pac.push_back(num);
}

void PART::pack(vector <double> &pac, vector <long> &vec)
{
	long imax, i; imax = vec.size(); pac.push_back(imax); for(i = 0; i < imax; i++) pac.push_back(vec[i]);
}

void PART::pack(vector <double> &pac, vector <double> &vec)
{
	long imax, i; imax = vec.size(); pac.push_back(imax); for(i = 0; i < imax; i++) pac.push_back(vec[i]);
}

void PART::pack(vector <double> &pac, vector< vector <long> > &vec)
{
	long imax, i, jmax, j;
	imax = vec.size(); pac.push_back(imax); 
	for(i = 0; i < imax; i++){
		jmax = vec[i].size(); pac.push_back(jmax); for(j = 0; j < jmax; j++) pac.push_back(vec[i][j]);
	}
}

void PART::pack(vector <double> &pac, vector< vector <double> > &vec)
{
	long imax, i, jmax, j;
	imax = vec.size(); pac.push_back(imax); 
	for(i = 0; i < imax; i++){
		jmax = vec[i].size(); pac.push_back(jmax); for(j = 0; j < jmax; j++) pac.push_back(vec[i][j]);
	}
}

void PART::pack(vector <double> &pac, vector< vector <FEV> > &vec, short fedivmin, short fedivmax)
{
	long imax, i, jmax, j;
	imax = vec.size(); pac.push_back(imax); 
	for(i = fedivmin; i < fedivmax; i++){
		jmax = vec[i].size(); pac.push_back(jmax); 
		for(j = 0; j < jmax; j++){
			pac.push_back(vec[i][j].trans);
			pac.push_back(vec[i][j].ind);
			pac.push_back(vec[i][j].t);
			pac.push_back(vec[i][j].done);
		}
	}
}

void PART::unpack(long &k, double *buffer, long &num)
{
	num = buffer[k]; k++;
}

void PART::unpack(long &k, double *buffer, vector <long> &vec)
{
	long imax, i; imax = buffer[k]; k++; vec.resize(imax); for(i = 0; i < imax; i++){ vec[i] = buffer[k]; k++;}
}

void PART::unpack(long &k, double *buffer, vector <double> &vec)
{
	long imax, i; imax = buffer[k]; k++; vec.resize(imax); for(i = 0; i < imax; i++){ vec[i] = buffer[k]; k++;}
}

void PART::unpack(long &k, double *buffer, vector< vector <long> > &vec)
{
	long imax, i, jmax, j;
	imax = buffer[k]; k++; vec.resize(imax);
	for(i = 0; i < imax; i++){
		jmax = buffer[k]; k++; vec[i].resize(jmax); for(j = 0; j < jmax; j++){ vec[i][j] = buffer[k]; k++;}
	}
}

void PART::unpack(long &k, double *buffer, vector< vector <double> > &vec)
{
	long imax, i, jmax, j;
	imax = buffer[k]; k++; vec.resize(imax);
	for(i = 0; i < imax; i++){
		jmax = buffer[k]; k++; vec[i].resize(jmax); for(j = 0; j < jmax; j++){ vec[i][j] = buffer[k]; k++;}
	}
}

void PART::unpack(long &k, double *buffer, vector< vector <FEV> > &vec, short fedivmin, short fedivmax)
{
	long imax, i, jmax, j;
	imax = buffer[k]; k++; vec.resize(imax);
	for(i = fedivmin; i < fedivmax; i++){
		jmax = buffer[k]; k++; vec[i].resize(jmax);
		for(j = 0; j < jmax; j++){ 
			vec[i][j].trans = buffer[k]; k++;
			vec[i][j].ind = buffer[k]; k++;
			vec[i][j].t = buffer[k]; k++;
			vec[i][j].done = buffer[k]; k++;
		}
	}
}

long PART::fevpack(double *buffer, short fedivmin, short fedivmax)
{
	long j, jmax;
	vector <double> pac; 
	
	pack(pac,fev,fedivmin,fedivmax);
	
	jmax = pac.size();
	for(j = 0; j < jmax; j++) buffer[j] = pac[j];
	
	return jmax;
}

long PART::partpack(double *buffer, short fedivmin)
{
	long l, c, cmax, cc, j, jmax;
	double val;
	vector <double> pac; 
	
	for(l = 0; l < poptree.level; l++){        // Process all the add later requests
		cmax = lev[l].node.size();
		for(c = 0; c < cmax; c++){
			val = lev[l].add[c]; lev[l].add[c] = 0;
			Rtot[l][c] += val*sussum[l][c];
			if(l <  poptree.level-1){
				val += addlater[l][c]; addlater[l][c] = 0;
				for(j = 0; j < 4; j++) lev[l+1].add[ lev[l].node[c].child[j]] += val;
			}
		}			
	}
	for(l = 0; l < poptree.level; l++){
		cmax = lev[l].node.size();
		for(c = 0; c < cmax; c++){
			if(addlater[l][c] != 0){ cout << l << " " << addlater[l][c] << "p\n";  emsg("Not ze");}
			if(lev[l].add[c] != 0) emsg("PP");
		}
	}
	
	l = poptree.level-1;
	
	if(checkon == 1) pack(pac,ffine);// cout << pac.size() <<  "f1\n";
	pack(pac,indinf); //cout << pac.size() <<  "f2\n";
	pack(pac,Rtot[l]); //cout << pac.size() <<  "f3\n";
	//pack(pac,Rtot); 
	//pack(pac,sussum); //cout << pac.size() <<  "f4\n";
	//pack(pac,addlater); //cout << pac.size() <<  "f5\n";
	pack(pac,fev,fedivmin,fediv); // cout << pac.size() <<  "f6\n";
	pack(pac,N);
	//cout << pac.size() <<  "f7\n";
	pack(pac,tdnext); //cout << pac.size() <<  "f8\n";
	pack(pac,tdfnext); //cout << pac.size() <<  "f9\n";
	pack(pac,sett); //cout << pac.size() <<  "f10\n";

	jmax = pac.size();
	for(j = 0; j < jmax; j++) buffer[j] = pac[j];

	return jmax;
}

void PART::partunpack(double *buffer, int max, short fedivmin)
{
	long k=0, kmax, j, l, c, cmax, cc;
	double val, val2;
	
	l = poptree.level-1;
	
	if(checkon == 1) unpack(k,buffer,ffine);
	unpack(k,buffer,indinf);
	unpack(k,buffer,Rtot[l]);
	unpack(k,buffer,fev,fedivmin,fediv);
	unpack(k,buffer,N);
	unpack(k,buffer,tdnext);
	unpack(k,buffer,tdfnext);
	unpack(k,buffer,sett);
	if(k != max) emsg("PART: EC17");
		
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
				if(addlater[l][c] != 0) emsg("Not ze");
				if(lev[l].add[c] != 0) emsg("PP");
			}
		}
	}
}

void PART::copy(const PART &other, short fedivmin)
{
	short d;
	
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
