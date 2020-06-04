/// Initialises a particle for a MBP
void MBPPART::mbpinit(long p, vector < vector <FEV> > &xi)
{
	long c, cmax, l, i;
	
	fev.clear(); fev.resize(fediv);
		
	N.resize(comp.size()); for(c = 0; c < comp.size(); c++) N[c] = 0;
 	N[0] = data.popsize;
	
	sett = 0;
		
	MIi.resize(poptree.Cfine); MIp.resize(poptree.Cfine);
	susboth.resize(poptree.Cfine); susp.resize(poptree.Cfine);
	lami.resize(poptree.Cfine); lamp.resize(poptree.Cfine);
	
	l = poptree.level-1;
	for(c = 0; c < poptree.Cfine; c++){
		MIi[c] = 0; MIp[c] = 0;
		susboth[c] = lev[l].node[c].sussum; susp[c] = 0;
	}
	
	stati.resize(data.popsize); statp.resize(data.popsize);
	for(i = 0; i < data.popsize; i++){
		stati[i] = 0;	statp[i] = 0;
	}
	
	tdnext = fediv;
	
	xitdnext = 0; xitdfnext = 0;
	while(xitdnext < fediv && xi[xitdnext].size() == 0) xitdnext++;	
}
	
void MBPPART::mbp(double ti, double tf, vector < vector <FEV> > &xi, long &timegen)
{
	long td, j, c, cmax;
	double t, tpl, betai, phii, betap, phip, sum, dlam, z;
	vector <double> sumst;
	
	NEV n;
	vector <NEV> nev;

	long d; for(d = 0; d < fediv; d++){ for(j = 0; j < xi[d].size(); j++) xi[d][j].done = 0;}
	
	cmax = poptree.Cfine;

	sumst.resize(cmax);
	if(sett == nsettime) emsg("MBP: EC1");

	t = ti; tpl = t;
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

		timegen -= clock();
		sum = 0;		
		betai = model.betai[sett]; betap = model.betap[sett];
		phii = model.parami[model.phiparam]; phip = model.paramp[model.phiparam];	
		for(c = 0; c < cmax; c++){
			lami[c] = betai*MIi[c] + phii;
			lamp[c] = betap*MIp[c] + phip;
				
			dlam = (lamp[c]-lami[c])*susboth[c]; if(dlam < 0) dlam = 0;
			sum += dlam;
			sum += lamp[c]*susp[c];
			sumst[c] = sum;
		}
		//if(sum != 0) emsg("sum not zero");
		timegen += clock();
			
		if(sum < tiny) n.t = tf; else n.t = t - log(ran())/sum;
		n.type = INF_EV;
		nev.push_back(n);
		
		sort(nev.begin(),nev.end(),compNEV);
		
		if(1 == 0){
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
			sett++; if(sett >= nsettime) emsg("MBP: EC2");
			break;
		
		case INF_EV:                 // These are infection events within the system
			z = ran()*sum;
			c = 0; while(c < cmax && z > sumst[c]) c++;
			if(c == cmax) emsg("MBPcode: EC3");
			//if(c > 0) cout <<  sumst[c]-sumst[c-1] << "comp\n";
			
			mbpaddinfc(c,t);	
			break;
			
		case FEV_EV:                 // These correspond to other compartmental transitions (e.g. E->A, E->I etc...)
			//cout << t << " f e\n";
			mbpdofe();
			break;
			
		case XIFEV_EV:
		//cout << t << " xife\n";
			mbpxidofe(xi);
			break;
		
		default: emsg("MBP: EC4"); break;
		}
		
		if(checkon == 1) check(0);
	}while(t < tf);
	
	/*
	for(d = 0; d < fediv; d++){
		if(fev[d].size() != xi[d].size()) emsg("p1");
		for(j = 0; j < fev[d].size(); j++){
			if(fev[d][j].t != xi[d][j].t) emsg("p2");
		}
	}
	*/
}

void MBPPART::mbpdofe()
{
	long i, c, cmax, cc, ccc, j, jmax, h, ii, k, kmax, l, ll, tra;
	double fac, val, num, dd, MInew;
	TRANS tr;

	long **&nMval(poptree.nMval);
	long ***&Mnoderef(poptree.Mnoderef);
	float ***&Mval(poptree.Mval);
	
	i = fev[tdnext][tdfnext].ind; if(fev[tdnext][tdfnext].done != 0){ cout << i << "\n"; emsg("MBP: EC5");}
	fev[tdnext][tdfnext].done = 1;
	c = poptree.ind[i].noderef;

	//cout << i << " " <<  fev[tdnext][tdfnext].t << " FEE\n";
	tra = fev[tdnext][tdfnext].trans;
	tr = trans[tra];
	
	//cout << comp[tr.from].name << "->" << comp[tr.to].name << " frp\n";
	
	N[tr.from]--; if(N[tr.from] < 0) emsg("MBP: EC6"); 
	N[tr.to]++;
	
	fac = poptree.ind[i].inf*(comp[tr.to].infectivity - comp[tr.from].infectivity);

	tdfnext++;
	if(tdfnext == fev[tdnext].size()){
		tdnext++; tdfnext = 0; 
		while(tdnext < fediv && fev[tdnext].size() == 0) tdnext++;
	}
		
	//cout << fev[tdnext][tdfnext].t << " nexttime\n";
	
	if(fac == 0) return;

	for(l = poptree.level-1; l >= 0; l--){
		kmax = nMval[c][l];
		for(k = 0; k < kmax; k++){
			cc = Mnoderef[c][l][k];
			val = fac*Mval[c][l][k];
			
			num = val*sussum[l][cc];
			jmax = lev[l].node[cc].fine.size();
			for(j = 0; j < jmax; j++){
				ccc = lev[l].node[cc].fine[j];
				MInew = MIp[ccc]+val; if(MInew < 0){ if(MInew < -tiny) emsg("MBPcode: EC7"); MInew = 0;}
				MIp[ccc] = MInew;
			}
		}
	}
}

void MBPPART::mbpxidofe(vector < vector <FEV> > &xi)
{
	long i, c, cmax, cc, ccc, j, jmax, h, ii, k, kmax, l, ll, tra;
	double fac, val, num, dd, MInew, t, sus, al;
	TRANS tr;

	long **&nMval(poptree.nMval);
	long ***&Mnoderef(poptree.Mnoderef);
	float ***&Mval(poptree.Mval);
	
	i = xi[xitdnext][xitdfnext].ind; if(xi[xitdnext][xitdfnext].done != 0) emsg("MBP: EC8");
	xi[xitdnext][xitdfnext].done = 1;
	c = poptree.ind[i].noderef;

	tra = xi[xitdnext][xitdfnext].trans;
	tr = trans[tra];
	//cout << comp[tr.from].name << "->" << comp[tr.to].name << " fr\n";
	fac = poptree.ind[i].inf*(comp[tr.to].infectivity - comp[tr.from].infectivity);

	if(fac != 0){
		for(l = poptree.level-1; l >= 0; l--){
			kmax = nMval[c][l];
			for(k = 0; k < kmax; k++){
				cc = Mnoderef[c][l][k];
				val = fac*Mval[c][l][k];
				
				num = val*sussum[l][cc];
				jmax = lev[l].node[cc].fine.size();
				for(j = 0; j < jmax; j++){
					ccc = lev[l].node[cc].fine[j];
					MInew = MIi[ccc]+val; if(MInew < 0){ if(MInew < -tiny) emsg("MBPcode: EC9"); MInew = 0;}
					MIi[ccc] = MInew;
				}
			}
		}
	}
	
	// decides to copy over event or not
	
	t = xi[xitdnext][xitdfnext].t;
	if(tr.type == INFECTION){
		sus = poptree.ind[i].sus;
		//cout << stati[i] << " " << statp[i] << " sta\n";
		stati[i] = 1;
	
		if(statp[i] == 0){
			susboth[c] -= sus;
		
			al = lamp[c]/lami[c]; 
			//if(al != 1) emsg("al not 1");
			if(ran() < al){  // keeps infection event
				statp[i] = 2;
			
				addfev(t,tra,i,1);
			}
			else{
				susp[c] += sus;
			}
		}	
	}
	else{
		if(statp[i] == 2){
			addfev(t,tra,i,0);
			mbpdofe();
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
	//cout << sum << "sum\n";
	
	z = ran()*sum;       
	j = 0; while(j < jmax && z > sumst[j]) j++; if(j == jmax) emsg("MBP: EC11");
	i = poptree.subpop[c][j];

	sus = poptree.ind[i].sus;
	if(stati[i] == 0) susboth[c] -= sus;
	else susp[c] -= sus;

	statp[i] = 1;
	
	simmodel(i,t);
}
