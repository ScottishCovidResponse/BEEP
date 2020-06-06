// The PMBP algorithm

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
#include "MBPPART.hh"
#include "output.hh"
#include "pack.hh"

MBPPART::MBPPART(DATA &data, MODEL &model, POPTREE &poptree) : data(data), model(model), comp(model.comp), trans(model.trans), poptree(poptree), lev(poptree.lev)
{
}

/// Initialises a particle for a MBP
void MBPPART::mbpinit(long p, vector < vector <FEV> > &xi)
{
	long c, cmax, l, i;
		
	fev.clear(); fev.resize(fediv);
		
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
	
	statilist1.clear(); statplist1.clear(); statplist2.clear();
	stati.resize(data.popsize); statp.resize(data.popsize);
	for(i = 0; i < data.popsize; i++){ stati[i] = 0;	statp[i] = 0;}
	
	tdnext = fediv;
	
	xitdnext = 0; xitdfnext = 0;
	while(xitdnext < fediv && xi[xitdnext].size() == 0) xitdnext++;	
}
	
void MBPPART::mbp(double ti, double tf, vector < vector <FEV> > &xi)
{
	long td, j, c, cmax;
	double t;
	
	NEV n;
	vector <NEV> nev;

	cmax = poptree.Cfine;
	if(sett == nsettime) emsg("MBP: EC1");

	betai = model.betai[sett]; betap = model.betap[sett];

	t = ti; 
	do{
		//cout << t << " t" << endl;
		
		nev.clear();                     // First we decide what event is next
		n.t = model.settime[sett];
		n.type = SET_EV;
		nev.push_back(n); 
		 
		if(tdnext < fediv) n.t = fev[tdnext][tdfnext].t; else n.t = tf;
		n.type = FEV_EV;
		nev.push_back(n);
	
		if(xitdnext < fediv) n.t = xi[xitdnext][xitdfnext].t; else n.t = tf;
		n.type = XIFEV_EV;
		nev.push_back(n);
		
		if(Rtot[0][0] < tiny) n.t = tf; else n.t = t - log(ran())/Rtot[0][0];
		n.type = INF_EV;
		nev.push_back(n);
		
		sort(nev.begin(),nev.end(),compNEV);
		
		t = nev[0].t; if(t >= tf) break;
	
		switch(nev[0].type){
		case SET_EV:                 // These are "settime" events which allow the value of beta to change in time
			sett++; if(sett >= nsettime) emsg("MBP: EC2");
			betai = model.betai[sett]; betap = model.betap[sett];
			for(c = 0; c < cmax; c++) updateRtot(c); 
			break;
		
		case INF_EV:                 // These are infection events within the system		
			c = mbpnextinfection();
			mbpaddinfc(c,t);	
			break;
			
		case FEV_EV:                 // These correspond to other compartmental transitions (e.g. E->A, E->I etc...)
			mbpdofe(1);
			break;
			
		case XIFEV_EV:
			mbpxidofe(xi);
			break;
		
		default: emsg("MBP: EC4"); break;
		}
	}while(t < tf);

	if(checkon == 1) check(0);
}

/// Recreates xi for the reverese transition 
void MBPPART::recreate(double ti, double tf, vector < vector <FEV> > &xi, vector < vector <FEV> > &xp)
{
	long td, j, jmax, d, c, cmax, xptdnext, xptdfnext, tra, i, fedivmin, fedivmax, doxi, doxp;
	double timeset, timexife, timexpfe, sus;
	TRANS tr;
	FEV fe;
	
	fedivmin = (fediv*ti)/data.tmax; fedivmax = (fediv*tf)/data.tmax;

	cmax = poptree.Cfine;	
	if(sett == nsettime) emsg("MBP: EC1");

	betai = model.betai[sett]; betap = model.betap[sett];

	xptdnext = fedivmin; xptdfnext = 0;
	while(xptdnext < fediv && xp[xptdnext].size() == 0) xptdnext++;	
	
	do{
		timeset = model.settime[sett];	 
		if(xptdnext < fediv) timexpfe = xp[xptdnext][xptdfnext].t; else timexpfe = tf;
		if(xitdnext < fediv) timexife = xi[xitdnext][xitdfnext].t; else timexife = tf;
		
		if(timeset >= tf && timexpfe >= tf && timexife >= tf) break;
		
		if(timeset < timexpfe && timeset < timexife){
			sett++; if(sett >= nsettime) emsg("MBP: EC2");
			betai = model.betai[sett]; betap = model.betap[sett];
		}
		else{
			doxi = 0; doxp = 0;
			if(timexpfe < timexife) doxp = 1;
			else{
				if(timexife < timexpfe) doxi = 1;
				else{ doxi = 1; doxp = 1; if(timexife != timexpfe) emsg("MBPART: EC9"); }
			}

			if(doxp == 1){   // event on xp
				tra = xp[xptdnext][xptdfnext].trans;
				tr = trans[tra];
				i = xp[xptdnext][xptdfnext].ind;
				c = poptree.ind[i].noderef;

				N[tr.from]--; if(N[tr.from] < 0) emsg("Part: EC12"); 
				N[tr.to]++;
		
				if(tr.type == INFECTION){
					sus = poptree.ind[i].sus;
					if(doxi == 0){
						if(stati[i] == 0) susboth[c] -= sus;
						else susp[c] -= sus;
						statp[i] = 1; statplist1.push_back(i);	
					}
					else{
						susboth[c] -= sus;
						statp[i] = 2; statplist2.push_back(i);
					}
				}
				else{
					if(doxi == 0) MIupdate(i,tra,0,1);
				}
				
				xptdfnext++;
				if(xptdfnext == xp[xptdnext].size()){
					xptdnext++; xptdfnext = 0; 
					while(xptdnext < fediv && xp[xptdnext].size() == 0) xptdnext++;
				}
			}
			
			if(doxi == 1){
				if(doxp == 0){
					tra = xi[xitdnext][xitdfnext].trans;
					tr = trans[tra];
					i = xi[xitdnext][xitdfnext].ind;
					c = poptree.ind[i].noderef;
				}
				else{
					if(tra != xi[xitdnext][xitdfnext].trans) emsg("MBPART: EC10");
					if(i != xi[xitdnext][xitdfnext].ind) emsg("MBPART: EC11");
				}
					
				if(tr.type == INFECTION){
					stati[i] = 1; statilist1.push_back(i);	
						
					if(doxp == 0 && statp[i] == 0){
						sus = poptree.ind[i].sus;
						susboth[c] -= sus; susp[c] += sus;}
					}
				}
				else{
					if(doxp == 0) MIupdate(i,tra,1,0);
					else MIupdate(i,tra,1,1);
				}
				
				xitdfnext++;
				if(xitdfnext == xi[xitdnext].size()){
					xitdnext++; xitdfnext = 0; 
					while(xitdnext < fediv && xi[xitdnext].size() == 0) xitdnext++;
				}
			}
		}
	}while(1 == 1);
	
	for(d = fedivmin; d < fedivmax; d++) fev[d] = xp[d];   // Recreates fev
	for(d = fedivmax; d < fediv; d++){
		fev[d].clear();
		jmax = xp[d].size();
		for(j = 0; j < jmax; j++){
			fe = xp[d][j];
			if(statp[fe.ind] == 1) fev[d].push_back(fe);
		}
	}
	
	tdnext = fedivmax; tdfnext = 0;
	while(tdnext < fediv && fev[tdnext].size() == 0) tdnext++;	
	
	for(c = 0; c < cmax; c++) updateRtot(c); 
		
	if(checkon == 1) check(1);
}

/// This samples the node on the fine scale in which the next infection occurs due to the MBP
long MBPPART::mbpnextinfection()
{
	long l, lmax, c, j, jmax;
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

// Updates Rtot so samping can be performed
void MBPPART::updateRtot(long c) 
{
	short l;
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

void MBPPART::mbpdofe(short update)
{
	long tra;
	TRANS tr;

	tra = fev[tdnext][tdfnext].trans;
	tr = trans[tra];
	
	if(update == 1) MIupdate(fev[tdnext][tdfnext].ind,tra,0,1);
		
	N[tr.from]--; if(N[tr.from] < 0) emsg("MBPPART: EC6"); 
	N[tr.to]++;
	
	tdfnext++;
	if(tdfnext == fev[tdnext].size()){
		tdnext++; tdfnext = 0; 
		while(tdnext < fediv && fev[tdnext].size() == 0) tdnext++;
	}
}

void MBPPART::MIupdate(long i, long tra, short upMIi, short upMIp)
{
	long c, cc, l, k, kmax, j, jmax;
	long **&nMval(poptree.nMval);
	long ***&Mnoderef(poptree.Mnoderef);
	float ***&Mval(poptree.Mval);
	TRANS tr;
	double fac, val;
		
	tr = trans[tra];
	fac = poptree.ind[i].inf*(comp[tr.to].infectivity - comp[tr.from].infectivity);
	if(fac == 0) return;

	c = poptree.ind[i].noderef;

	l = poptree.level-1;
	kmax = nMval[c][l];
	for(k = 0; k < kmax; k++){
		cc = Mnoderef[c][l][k];
		val = fac*Mval[c][l][k];
		
		if(upMIi == 1) MIi[cc] += val;	
		if(upMIp == 1) MIp[cc] += val;
		updateRtot(cc); 
	}
}

void MBPPART::mbpxidofe(vector < vector <FEV> > &xi)
{
	long i, c, tra;
	double sus, al, t;
	TRANS tr;

	i = xi[xitdnext][xitdfnext].ind;
	c = poptree.ind[i].noderef;

	tra = xi[xitdnext][xitdfnext].trans;
	tr = trans[tra];
	
	t = xi[xitdnext][xitdfnext].t;                         	// This part decides to copy over event on xi into xp or not
	if(tr.type == INFECTION){
		sus = poptree.ind[i].sus;
		stati[i] = 1; statilist1.push_back(1);
	
		if(statp[i] == 0){
			susboth[c] -= sus;
		
			al = lamp[c]/lami[c]; 
			if(ran() < al){                                    // Keeps the infection event
				statp[i] = 2; statplist2.push_back(i);
			
				addfev(t,tra,i,1);
			}
			else susp[c] += sus;                               // Does not keep the infection event
			updateRtot(c); 
		}	
	}
	else{
		if(statp[i] == 2){                                   // Keeps the non infection event if infection happened at the same time
			addfev(t,tra,i,0);
			mbpdofe(0);
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

/// Adds an exposed indivdual on node c on the finest scale (i.e. level-1)
void MBPPART::mbpaddinfc(long c, double t)
{
	long l, i, j, jmax, cc, k, kmax;
	double dR, sum, sus, z, dlam;
	vector <double> sumst;
	
	jmax = poptree.subpop[c].size(); 
	sumst.resize(jmax);
	sum = 0; 
	for(j = 0; j < jmax; j++){
		i = poptree.subpop[c][j];
		if(statp[i] == 0){
			if(stati[i] == 0){
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
	if(stati[i] == 0) susboth[c] -= sus;
	else susp[c] -= sus;
	updateRtot(c); 
	
	statp[i] = 1; statplist1.push_back(i);
	
	simmodel(i,t);
}


// COPIED FROM PART

/// This function simulates events as in individual passes through the compartmental model
void MBPPART::simmodel(long i, double t)
{
	short c, k, kmax, tra;
	double tnext, tnextmin, mean, sd, dt;
	TRANS tr;
	
	switch(modelsel){
	case MOD_IRISH:
		if(ran() < model.param[model.afracparam].val) tra = 0;
		else tra = 1;
		break;
		
	case MOD_OLD:
		tra = 0;
		break;
	}
	addfev(t,tra,i,1);
	c = trans[tra].to;
	
	do{
		kmax = comp[c].trans.size();
		if(kmax == 0) break;
		
		tnextmin = data.tmax;
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
					
				case LOGNORM_DIST:
					mean = model.param[tr.param1].val; sd = model.param[tr.param2].val;
					dt = exp(normal(mean,sd));
					tnext = t + dt; 
					break;
					
				default: emsg("Part: EC2b"); break;
			}
			
			if(tnext < tnextmin){ tnextmin = tnext; tra = comp[c].trans[k];}
		}
		
		if(tnextmin == data.tmax) break; 
		
		addfev(tnextmin,tra,i,0);

		c = trans[tra].to; t = tnextmin;
	}while(1 == 1);
}

/// Adds a future event to the timeline
void MBPPART::addfev(double t, long tra, long i, short done)
{
	long d, j, jmax;
	
	if(t >= data.tmax) return;
	
	FEV fe; fe.t = t; fe.trans = tra; fe.ind = i; //fe.done = done;
	
	d = long((t/data.tmax)*fediv);
	j = 0; jmax = fev[d].size();
	if(done == 0){ while(j < jmax && t >= fev[d][j].t) j++;}
	else{ while(j < jmax && t >= fev[d][j].t) j++;}
	
	if(j == jmax) fev[d].push_back(fe);
	else fev[d].insert(fev[d].begin()+j,fe);
	
	if(done == 0){
		if(d == tdnext){ if(j < tdfnext) tdfnext = j;}
		if(d < tdnext){ tdnext = d; tdfnext = j;}
	}
	else{
		TRANS tr = trans[tra];
		N[tr.from]--; if(N[tr.from] < 0) emsg("Part: EC12"); 
		N[tr.to]++;
		
		if(d == tdnext) tdfnext++;
	}
}

/// Measures how well the particle agrees with the observations for a given week w
/// (which in this case is weekly hospitalised case data)
void MBPPART::Lobs(short w, vector < vector <FEV> > &fev)
{
	short r;
	double mean, var;
	vector <long> num;
	
	Li = 0;	
	num = getnumtrans(data,model,poptree,fev,"I","H",w*timestep,(w+1)*timestep);
	for(r = 0; r < data.nregion; r++){
		mean = data.ncase[r][w];
		var = mean; if(var < 5) var = 5;
		//var *= 4;

		//Li += -0*0.5*log(2*3.141592654*var) - (mean-num[r])*(mean-num[r])/(2*var);
		Li += -(mean-num[r])*(mean-num[r])/(2*var);
	}
}

void MBPPART::check(short num)
{
	long c, cmax, l, j, jmax, i;
	double d, val, dlam, sumboth, sump;

	cmax = poptree.Cfine;      

	for(c = 0; c < cmax; c++){
		sumboth = 0; sump = 0;
		jmax = poptree.subpop[c].size(); 
		for(j = 0; j < jmax; j++){
			i = poptree.subpop[c][j];
			if(stati[i] == 0 && statp[i] == 0) sumboth += poptree.ind[i].sus; 
			if(stati[i] != 0 && statp[i] == 0) sump += poptree.ind[i].sus; 
		}
		d = sumboth - susboth[c]; if(d*d > tiny) emsg("both prob");
		d = sump - susp[c]; if(d*d > tiny) emsg("sump prob");
	}
		
	l = poptree.level-1;    // checks that Rtot is correct
	for(c = 0; c < cmax; c++){
		lami[c] = betai*MIi[c] + phii;
		lamp[c] = betap*MIp[c] + phip;
			
		dlam = (lamp[c]-lami[c])*susboth[c];
		if(dlam < 0) dlam = 0;
		val = dlam + lamp[c]*susp[c];
		d = val-Rtot[l][c]; if(d*d > tiny) emsg("Mbppart: EC3");
	}
}

/// Copies in all the information from another particle
void MBPPART::copy(const MBPPART &other, short fedivmin)
{
	short d;
	 
	susboth = other.susboth;
	susp = other.susp;
	lami = other.lami;
	lamp = other.lamp;
	
	Rtot = other.Rtot; 
	
	stati = other.stati;
	statp = other.statp;
	//
	
	MIi = other.MIi;
	MIp = other.MIp;
	statilist1 = other.statilist1; statplist1 = other.statplist1; statplist2 = other.statplist2; 
	
	for(d = fedivmin; d < fediv; d++) fev[d] = other.fev[d];
	
	N = other.N;
	tdnext = other.tdnext;
	tdfnext = other.tdfnext;
	xitdnext = other.xitdnext;
	xitdfnext = other.xitdfnext;
	sett = other.sett;
}

/// Packs up all the particle information (from time fedivmin until the end) to the be sent by MPI
void MBPPART::partpack(short fedivmin)
{
	long l, c, cmax, cc, j, jmax;
	double val;

	/*
	packinit();
	pack(susboth);
	pack(susp);
	pack(lami);
	pack(lamp);
	
	pack(Rtot); 
	
	pack(stati);
	pack(statp);

	pack(MIi);
	pack(MIp);
	pack(statilist1); pack(statplist1); pack(statplist2); 
	
	pack(fev,fedivmin,fediv);
	
	pack(N);
	pack(tdnext);
	pack(tdfnext);
	pack(xitdnext);
	pack(xitdfnext);
	pack(sett);
	*/
}

/// Unpacks particle 
void MBPPART::partunpack(short fedivmin)
{
	/*
	packinit();
	unpack(susboth);
	unpack(susp);
	unpack(lami);
	unpack(lamp);
	
	unpack(Rtot); 
	
	unpack(stati);
	unpack(statp);

	unpack(MIi);
	unpack(MIp);
	unpack(statilist1); unpack(statplist1); unpack(statplist2); 
	
	unpack(fev,fedivmin,fediv);
	
	unpack(N);
	unpack(tdnext);
	unpack(tdfnext);
	unpack(xitdnext);
	unpack(xitdfnext);
	unpack(sett);
	*/
}
