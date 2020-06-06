// Performs particle MBPs

#include <fstream>
#include <iostream>
#include <sstream>
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
#include "consts.hh"

MPI_Request mbpreqs[partmax];               // These are information used Isend and Irecv
MPI_Status mbpstats[partmax];

double mbpsendbuffer[SENDRECMAX][BUFMAX];
double mbprecibuffer[SENDRECMAX][BUFMAX];

static double mbpbootstrap(short core, short ncore, short npart, short fedivmin, short w, short *backpart, short dir);
static double pmbpsample(MODEL &model, short core, short ncore, short npart, double tmax, vector < vector <FEV> > &xi, vector < vector <FEV> > &xp, double invT);
static void mbpgeneventsample(short core, short ncore, short npart, short nweek, short *backpart, vector < vector <FEV> > &xp);
static void checkrecreate(short nweek, vector < vector <FEV> > &xi);

MBPPART* mbppart[partmax]; 

void PMBP(DATA &data, MODEL &model, POPTREE &poptree, long nsamp, short core, short ncore, short npart)
{
	long p, pp, d, j, nparam, st, w, burnin = nsamp/3, samp, accept, nweek = data.tmax/timestep, r;
	double valst, al, L;
	vector <SAMPLE> opsamp;            // Stores sample for generating satatistics later
	vector <long> num;
	PART *part;
	double invT =  1;                  // The inverse temperature (used to relax the observation model)

	vector <PARAM> &param(model.param);
		
	vector < vector <FEV> > xi, xp;
	
	srand(core);
		
	part = new PART(data,model,poptree);
	part->partinit(0);
	part->gillespie(0,data.tmax, 1 /* simulating */);

	xi = part->fev;     // start

	/*
	outputsimulateddata(data,model,poptree,part->fev);
	for(w = 0; w < data.nweek; w++){  // Fixes data to simulation
		num = getnumtrans(data,model,poptree,xi,"I","H",w*timestep,(w+1)*timestep);
		for(r = 0; r < data.nregion; r++) data.ncase[r][w] = num[r];
	}
	*/
	
	if(core == 0) outputinit(model,ncore*npart);
	
	for(p = 0; p < npart; p++){ mbppart[p] = new MBPPART(data,model,poptree);}
	 
	nparam = model.param.size();
	model.parami.resize(nparam); model.paramp.resize(nparam);
	
	for(samp = 0; samp < nsamp; samp++){
		if(core == 0 && samp%1 == 0) cout << " Sample: " << samp << endl;
	
		for(p = 0; p < 8; p++){
		  //for(p = 0; p < long(param.size()); p++){
			if(param[p].min != param[p].max){
				for(pp = 0; pp < nparam; pp++) model.parami[pp] = model.param[pp].val;  // Sets the initial parameter set
				for(st = 0; st < nsettime; st++) model.betai[st] = model.beta[st];
	
				valst = param[p].val;
				
				if(core == 0) param[p].val += normal(0,param[p].jump);                  // Makes a change to a parameter
				MPI_Bcast(&param[p].val,1,MPI_DOUBLE,0,MPI_COMM_WORLD);

				if(param[p].betachange == 1) model.betaspline(data.tmax);
				if(param[p].suschange == 1) poptree.setsus(model);
				if(param[p].infchange == 1) poptree.setinf(model);
				
				if(param[p].val < param[p].min || param[p].val > param[p].max) al = 0;
				else{
					for(pp = 0; pp < nparam; pp++) model.paramp[pp] = model.param[pp].val;// Sets the proposed parameters
					for(st = 0; st < nsettime; st++) model.betap[st] = model.beta[st];
				
					al = pmbpsample(model,core,ncore,npart,data.tmax,xi,xp,invT);         // Performs a pmbp proposal
				}
				if(core == 0){ if(ran() < al) accept = 1; else accept = 0;}
				
				MPI_Bcast(&accept,1,MPI_LONG,0,MPI_COMM_WORLD);
					
				param[p].ntr++;
				if(accept == 1){
					if(core == 0){
						param[p].nac++;
						xi = xp;
						if(samp < burnin) param[p].jump *= 1.1;
					}
				}
				else{
					param[p].val = valst;
					if(param[p].betachange == 1) model.betaspline(data.tmax);
					if(param[p].suschange == 1) poptree.setsus(model);
					if(param[p].infchange == 1) poptree.setinf(model);
				
					if(core == 0){
						if(samp < burnin) param[p].jump *= 0.95;
					}
				}
			}
		}
	
		if(core == 0){
			L = 0; for(w = 0; w < nweek; w++){ mbppart[0]->Lobs(w,xi); L += mbppart[0]->Li;}
			opsamp.push_back(outputsamp(1,samp,L,data,model,poptree,xi));
		}			
	}
	
	if(core == 0) outputresults(data,model,opsamp,0,ncore*npart);
}

/// This modifies the event sequence using the change change in parameters
/// Particles which account for the observation probability 
static double pmbpsample(MODEL &model, short core, short ncore, short npart, double tmax, vector < vector <FEV> > &xi, vector < vector <FEV> > &xp, double invT)
{
	short p, w, nweek = tmax/timestep, backpart[ncore*npart*nweek], st;
	double al, Liav, Lpav, tmp;
	vector <double> temp;
	
	if(checkon == 1 && npart > 1) checkrecreate(nweek,xi);
	
	for(p = 0; p < npart; p++) mbppart[p]->mbpinit(p,xi);

	Lpav = 0;
	for(w = 0; w < nweek; w++){                                // We step through one measurement at a time
		for(p = 0; p < npart; p++){
			mbppart[p]->mbp(w*timestep,(w+1)*timestep,xi);
			mbppart[p]->Lobs(w,mbppart[p]->fev);
		}	
		Lpav += mbpbootstrap(core,ncore,npart,(fediv*(w+1))/nweek,w,backpart,1);
	}
			
	mbpgeneventsample(core,ncore,npart,nweek,backpart,xp);

	// switches over parameter sets
	temp = model.parami; model.parami = model.paramp; model.paramp = temp;
	for(st = 0; st < nsettime; st++){ tmp = model.betai[st]; model.betai[st] = model.betap[st]; model.betap[st] = tmp;}
	
	for(p = 0; p < npart; p++) mbppart[p]->mbpinit(p,xp);
	
	Liav = 0;
	for(w = 0; w < nweek; w++){                                // We step through one measurement at a time
		for(p = 0; p < npart; p++){
			if(core == 0 && p == 0) mbppart[p]->recreate(w*timestep,(w+1)*timestep,xp,xi);
			else mbppart[p]->mbp(w*timestep,(w+1)*timestep,xp);
			mbppart[p]->Lobs(w,mbppart[p]->fev);
		}
		Liav += mbpbootstrap(core,ncore,npart,(fediv*(w+1))/nweek,w,backpart,-1);
	}
	
	//timembp += clock();
	//cout << exp(invT*(Lpav-Liav)) << " " << Lpav << " " << Liav << " al\n";
	
	//Liav = 0; for(w = 0; w < nweek; w++){ mbppart[0]->Lobs(w,xi); Liav += mbppart[0]->Li;} cout << Liav << " Liav\n";
	//Lpav = 0; for(w = 0; w < nweek; w++){ mbppart[0]->Lobs(w,xp); Lpav += mbppart[0]->Li;} cout << Lpav << " Lpav\n";

	return exp(invT*(Lpav-Liav));
}

/// Constructs the event sequence sample by gathering all the pieces from different particles (from the bootstrap function)
static void mbpgeneventsample(short core, short ncore, short npart, short nweek, short *backpart, vector < vector <FEV> > &xp)	
{
	short p, w, co, d, fedivmin, fedivmax, nparttot = npart*ncore;
	int siz;
	
	if(1 == 0){
		for(w = 0; w < nweek; w++){
			cout << w << ": "; for(p = 0; p < nparttot; p++) cout << backpart[w*nparttot+p] << ",";
			cout << " Backpart" << endl;  
		}
	}

	xp.resize(fediv);
	for(w = nweek-1; w >= 0; w--){
		fedivmin = (fediv*w)/nweek;	fedivmax = (fediv*(w+1))/nweek;
		if(w == nweek-1) p = backpart[nparttot*w];  // Picks the first final particle
		else p = backpart[nparttot*w + p];
		
		co = p/npart;
		if(core == 0){  
			if(co == 0){
				for(d = fedivmin; d < fedivmax; d++) xp[d] = mbppart[p%npart]->fev[d];
			}
			else{
				emsg("err");
				
				MPI_Status status;
				MPI_Recv(packbuffer(),MAX_NUMBERS,MPI_DOUBLE,co,0,MPI_COMM_WORLD,&status);
				MPI_Get_count(&status, MPI_DOUBLE, &siz); if(siz >= MAX_NUMBERS) emsg("Buffer not big enough");
		
				packinit();
				unpack(xp,fedivmin,fedivmax);
				if(packsize() != siz) emsg("PMBP: EC10");
			}
		}
		else{
			if(co == core){
				emsg("err");
				packinit();
				pack(mbppart[p%npart]->fev,fedivmin,fedivmax);
				MPI_Send(packbuffer(),packsize(),MPI_DOUBLE,0,0,MPI_COMM_WORLD);	
			}
		}
	}
}

/// This function does the usual bootstrap step of randomly selecting particles based on their observation probability 
static double mbpbootstrap(short core, short ncore, short npart, short fedivmin, short w, short *backpart, short dir)
{
	long p, pp, pmin, pmax, cor, j, jmax, k, kmax, nparttot = npart*ncore, partstep;
	double res; 
	int flag, siz;
	short nreqs = 0, nsendbuf = 0, rec,	npreclist, ncorlist;
	double *buf = packbuffer();
	
	vector <short> corlist;
	vector <short> preclist;
	vector <vector <short> > preclistli;
	
	backpart += w*nparttot;	
	
	if(core == 0){
		double Limax, av, z, Litot[nparttot], sum, sumst[nparttot];
		vector <short> extra;
		
		for(p = 0; p < npart; p++) Litot[p] = mbppart[p]->Li;                    // Gathers together the likelihoods from different cores 
	
		for(cor = 1; cor < ncore; cor++){
			MPI_Recv(Litot+npart*cor,npart,MPI_DOUBLE,cor,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
		}
		
		Limax = -large;	for(p = 0; p < nparttot; p++){ if(Litot[p] > Limax) Limax = Litot[p];}
		
		av = 0;	for(p = 0; p < nparttot; p++) av += exp(Litot[p] - Limax); av /= nparttot;

		res = Limax + log(av);
		
		sum = 0;
		for(p = 0; p < nparttot; p++){
			sum += exp(Litot[p] - Limax);
			sumst[p] = sum;
		}

		for(p = 0; p < nparttot; p++) backpart[p] = -1;
			
		//partstep = nparttot/32; if(partstep < 2) partstep = 2;
		
		for(p = 0; p < nparttot; p++){                                            // Finds "backpart", the particle which is being copied from
			if(p == 0 && dir == -1) pp = 0;  // Used for the reverese transition
			else{
				z = sum*ran();
		
				pp = 0;
				//while(pp < nparttot && z > sumst[pp]) pp += partstep;
				//if(pp > 0) pp -= partstep;
				
				while(pp < nparttot && z > sumst[pp]) pp++;
				if(pp == nparttot) emsg("PMBP: EC1");	
				if(pp > 0){ if(z < sumst[pp-1]) emsg("PMBP: EC1a");}
			}
			
			if(backpart[pp] == -1) backpart[pp] = pp;
			else extra.push_back(pp);
		}
		
		j = 0; jmax = extra.size();
		while(j < jmax){                                                            // Tries to copy to same core (to increase speed)
			pp = extra[j];
			p = pp - pp%npart; pmax = p+npart; while(p < pmax && backpart[p] != -1) p++;

			if(p < pmax){	backpart[p] = pp; jmax--; extra[j] = extra[jmax]; extra.pop_back();}
			else j++;
		}
		
		for(p = 0; p < nparttot; p++){  
			if(backpart[p] == -1){
				pp = extra[extra.size()-1]; extra.pop_back();
				backpart[p] = pp;
			}
		}
		
		//for(p = 0; p < nparttot; p++) cout << exp(Litot[p]) << ","; cout << "Li\n"; 
		//for(p = 0; p < nparttot; p++) cout <<  backpart[p] << ","; cout << "back\n";

		if(extra.size() != 0) emsg("PMBP: EC9");
	}
	else{
		double Li[npart];
		for(p = 0; p < npart; p++) Li[p] = mbppart[p]->Li;
		MPI_Send(Li,npart,MPI_DOUBLE,0,0,MPI_COMM_WORLD);	
		res = 0;
	}
	
	MPI_Bcast(backpart,nparttot,MPI_SHORT,0,MPI_COMM_WORLD);                                  // Broadcasts back part to all cores
	
	pmin = core*npart;                                                                        // Sets up information to be recieved
	npreclist = 0;
	for(p = pmin; p < pmin+npart; p++){
		pp = backpart[p];
		if(p != pp){                   
			cor = pp/npart;
			if(cor != core){
				j = 0; while(j < npreclist && preclist[j] != pp) j++;
				if(j < npreclist) preclistli[j].push_back(p%npart);
				else{
					preclist.push_back(pp); 
					preclistli.push_back(vector <short> ()); 
					preclistli[npreclist].push_back(p%npart);
					
					MPI_Irecv(mbprecibuffer[npreclist],BUFMAX,MPI_DOUBLE,cor,pp,MPI_COMM_WORLD,&mbpreqs[nreqs]); nreqs++;
					npreclist++; if(npreclist == SENDRECMAX) emsg("PMBP: EC11"); 
				}
			}
		}
	}

	for(p = pmin; p < pmin+npart; p++){                                                    // Initiates information to be sent
		if(p == backpart[p]){ 
			corlist.clear(); ncorlist = 0;
			for(pp = 0; pp < nparttot; pp++){
				if(p == backpart[pp]){
					cor = pp/npart;
					if(cor == core){
						if(pp != p) mbppart[pp%npart]->copy(*mbppart[p%npart],fedivmin);
					}
					else{
						j = 0; while(j < ncorlist && corlist[j] != cor) j++;
						if(j == ncorlist){ corlist.push_back(cor); ncorlist++;}
					}
				}
			}
						
			if(ncorlist > 0){
				mbppart[p%npart]->partpack(fedivmin);
				kmax = packsize(); 
				mbpsendbuffer[nsendbuf][0] = kmax; if(kmax >= BUFMAX-1) emsg("PMBP: EC20");
				for(k = 0; k < kmax; k++) mbpsendbuffer[nsendbuf][k+1] = buf[k];
				
				for(j = 0; j < ncorlist; j++){
					MPI_Isend(mbpsendbuffer[nsendbuf],kmax+1,MPI_DOUBLE,corlist[j],p,MPI_COMM_WORLD,&mbpreqs[nreqs]);
					nreqs++;		
				}				
				nsendbuf++; if(nsendbuf == SENDRECMAX) emsg("PMBP: EC21");
			}
		}
	}
	
	if(nreqs > 0){
		if(MPI_Waitall(nreqs,mbpreqs,mbpstats) != MPI_SUCCESS) emsg("PMBP: EC22");
	}
		
	for(rec = 0; rec < npreclist; rec++){                                              	// Unpacks the recieved information
		kmax = mbprecibuffer[rec][0]; for(k = 0; k < kmax; k++) buf[k] = mbprecibuffer[rec][k+1];		
		p = preclistli[rec][0];
		mbppart[p]->partunpack(fedivmin);
		
		for(j = 1; j < preclistli[rec].size(); j++){
			pp = preclistli[rec][j];
			mbppart[pp]->copy(*mbppart[p],fedivmin);
		}
	}
	
	MPI_Barrier(MPI_COMM_WORLD); 
		
	return res;
}

/// Checks that recreate function is working correctly
static void checkrecreate(short nweek, vector < vector <FEV> > &xi)
{
	long p, w, c, cmax, l, i, imax, j, jmax, d;
	double dd;
	vector < vector <FEV> > xp;
		
	MBPPART *p1, *p2;
	
	p1 = mbppart[0]; p2 = mbppart[1];
	p1->mbpinit(0,xi);
	p2->mbpinit(1,xi);
	
	for(w = 0; w < nweek; w++){    
		p1->mbp(w*timestep,(w+1)*timestep,xi);
		xp = p1->fev;
		p2->recreate(w*timestep,(w+1)*timestep,xi,xp);
			
		imax = p1->stati.size();
		for(i = 0; i < imax; i++){
			if(p1->stati[i] !=  p2->stati[i]) emsg("stati");
			if(p1->statp[i] !=  p2->statp[i]) emsg("stati");
		}
		
		cmax = p1->susboth.size();
		for(c = 0; c < cmax; c++){
			dd = p1->susboth[c] - p2->susboth[c]; if(dd*dd > tiny) emsg("susboth");
			dd = p1->susp[c] - p2->susp[c]; if(dd*dd > tiny) emsg("susp");
			
			dd = p1->lami[c] - p2->lami[c]; if(dd*dd > tiny) emsg("lami");
			dd = p1->lamp[c] - p2->lamp[c]; if(dd*dd > tiny) emsg("lamp");
			
			dd = p1->MIi[c] - p2->MIi[c]; if(dd*dd > tiny) emsg("MIi");
			dd = p1->MIp[c] - p2->MIp[c]; if(dd*dd > tiny) emsg("MIp");
		}
		
		for(l = 0; l < p1->poptree.level; l++){
			for(c = 0; c < p1->Rtot[l].size(); c++){
				dd = p1->Rtot[l][c] - p2->Rtot[l][c]; if(dd*dd > tiny) emsg("Rtot");
			}
		}
		
		if(p1->tdnext != p2->tdnext) emsg("tdnext");
		if(p1->tdfnext != p2->tdfnext) emsg("tdfnext");
		
		if(p1->xitdnext != p2->xitdnext) emsg("tdnext");
		if(p1->xitdfnext != p2->xitdfnext) emsg("tdfnext");
	
		if(p1->sett != p2->sett) emsg("sett");
		
		for(d = 0; d < fediv; d++){
			jmax = p1->fev[d].size();
			if(jmax != p2->fev[d].size()){ cout << d << " " << p1->fev[d].size() << " " << p2->fev[d].size() << " d\n"; emsg("fev size");}
			for(j = 0; j < jmax; j++){
				if(p1->fev[d][j].ind != p2->fev[d][j].ind) emsg("fev ind");
				if(p1->fev[d][j].t != p2->fev[d][j].t) emsg("fev t");
				if(p1->fev[d][j].trans != p2->fev[d][j].trans) emsg("fev trans");
			}
		}
		
		for(c = 0; c < p1->N.size(); c++){
			if(p1->N.size() != p2->N.size()) emsg("N");
		}		
	}	
}
