// This header contains functions associated with simulation

// Generates weekly case data (similar to the actual data we currently have from Covid-19)
void simulatedata()
{
	long week, r;
	vector <long> num;
	
	part[0] = new PART;
  npart = 1;
	
	part[0]->partinit(0);
	
	timesim -= clock();
	part[0]->gillespie(0,tmax);
	timesim += clock();
		
	num = part[0]->getnumtrans("I","H",0,tmax);
	cout << "\nTotal number of hospitalised cases:\n";
	for(r = 0; r < nregion; r++) cout <<	"Region " <<  r << ": " << num[r] << "\n";
	cout << "\n";
	
	ofstream regplot("Weekly case data.txt");
	//ofstream regplot("Weekly case data.txt");
	for(week = 0; week < tmax/7; week++){
		regplot << week << " ";
		num = part[0]->getnumtrans("I","H",week*7,week*7+7);
		for(r = 0; r < nregion; r++) regplot << num[r] << " "; regplot << "\n";
	}
}

// Performs the modified Gillespie algorithm between times ti and tf 
void PART::gillespie(double ti, double tf)
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
