//  The data inputs

#include <iostream>
#include <vector>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <math.h> 

using namespace std;

#include "utils.hh"
#include "data.hh"	
#include "consts.hh"
#include "pack.hh"

/// Reads in transition and area data
void DATA::readdata(unsigned int core, unsigned int ncore, unsigned int mod, unsigned int per)
{
	unsigned int t, r, i, c, imax, k, nreg, td, j, jst, jmax, jj, cc, fl, d, dp, a1, a2, a, aa, vi, q, s;
	int dc;
	double v, sum;
	string line, ele, name, regcode, st;
	REGION reg;
	AREA are;
	DEMOCAT dem;
	IND indi;
	vector <unsigned int> count;
	vector <vector <double> > val;
	
	mode = mod;
	period = per;
	fepertime = 100;
	
	settpertime = 7;
	nsettime = settpertime*per;
	settime.resize(nsettime+1);
	for(s = 0; s <= nsettime; s++) settime[s] = double(s*period)/nsettime;
			
	fediv = nsettime*fepertime;
	
	if(core == 0){
		ifstream demoin(democatfile.c_str());                             	// Reads in demographic categories
		if (!demoin.is_open()) {
			emsg("Failed to open file '"+democatfile+"'");
		}
		do{
			getline(demoin,line); line = strip(line);
			if(demoin.eof()) break;
			
			j = 0; jmax = line.length(); while(j < jmax && line.substr(j,1) != ":") j++;
			if(j == jmax) emsg("Problem loading: ",democatfile);
			dem.name = line.substr(0,j);
			j++;
			
			dem.value.clear();
			while(j < jmax){	
				while(j < jmax && line.substr(j,1) == " " ) j++;
				jst = j;
				while(j < jmax && line.substr(j,1) != "," ) j++;
				dem.value.push_back(line.substr(jst,j-jst));
				j++;
			}
			democat.push_back(dem);
		}while(1 == 1);
		ndemocat = democat.size();
		nage = democat[0].value.size();
		
		cout << endl << "Demographic data loaded:" << endl;
		for(d = 0; d < ndemocat; d++){
			cout << democat[d].name << ": "; 
			for(j = 0; j < democat[d].value.size(); j++){
				if(j > 0) cout << ",";
				cout << democat[d].value[j] << " ";
			}
			cout << endl;
		}
		
		count.resize(ndemocat);                                     // Defines all the demographic states
		for(dc = 0; dc < int(ndemocat); dc++) count[dc] = 0;
		
		do{
			democatpos.push_back(count);
			
			dc = ndemocat-1;
			do{
				fl = 0;
				count[dc]++; if(count[dc] == democat[dc].value.size()){ fl = 1; count[dc] = 0; dc--;}
			}while(fl == 1 && dc >= 0);
		}while(dc >= 0);
		ndemocatpos = democatpos.size();
		ndemocatposperage = ndemocatpos/nage;
 
		ifstream regionin(regiondatafile.c_str());      	             // Loads information about data regions
		if(!regionin) emsg("Cannot find the file '"+regiondatafile+"'");
		getline(regionin,line); 
		do{
			getline(regionin,line); line = strip(line);
			if(regionin.eof()) break;
			stringstream ss(line);
			getline(ss,reg.name,'\t');
			getline(ss,reg.code,'\t');
			region.push_back(reg);
		}while(1 == 1);
		nregion = region.size();
		
		cout << endl << "Region data loaded." << endl;
		if(checkon == 1){
			for(r = 0; r < nregion; r++) cout << region[r].code << ", " << region[r].name  << " regionload" << endl;
		}
	
		ifstream areain(areadatafile.c_str());                                        // Loads information about areas
		if(!areain) emsg("Cannot find the file '"+areadatafile+"'");
			
		getline(areain,line);
		do{
			getline(areain,line); line = strip(line);
			if(areain.eof()) break;
			
			stringstream ss(line);	
			
			ss >> are.code >> are.x >> are.y >> regcode >> are.density;
			
			r = 0; while(r < nregion && region[r].code != regcode) r++;
			if(r == nregion) emsg("Region code not recognised: ",regcode);
			are.region = r;
			
			val.resize(democat.size());
			for(d = 0; d < democat.size(); d++){
				jmax = democat[d].value.size();
				val[d].resize(jmax);
				for(j = 0; j < jmax; j++){
					ss >> st;
					val[d][j] = atof(st.c_str());
					if(isnan(val[d][j])) emsg("In file '"+areadatafile+"' "+st+" is not a number");	
				}
			}
			
			are.agepop.resize(nage);
			for(a = 0; a < nage; a++){
				are.agepop[a] = val[0][a];
			}
				
			are.pop.resize(ndemocatpos);
			for(dp = 0; dp < ndemocatpos; dp++){
				for(j = 0; j < democat.size(); j++){
					if(j == 0) v = val[j][democatpos[dp][j]];
					else v *= val[j][democatpos[dp][j]]/100.0;
				}
				are.pop[dp] = (unsigned int)(v+0.5);
			}
		
			area.push_back(are);
		}while(1 == 1);
		narea = area.size();
		
		cout << endl << "Area data loaded." << endl;
		if(checkon == 1){
			for(c = 0; c < narea; c++){
				cout << nregion << " " << area[c].region << "region" << endl;
				cout << area[c].code << " " << region[area[c].region].code << " " <<  area[c].x << " "<<  area[c].y << " " <<  area[c].density << "  ***";
			
				for(j = 0; j < area[c].pop.size(); j++) cout << area[c].pop[j] << ", ";
				cout << endl;	
			}
		}			
		
		ifstream Min(Mdatafile.c_str());                                // Loads information about mixing matrix
		if(!Min) emsg("Cannot find the file '"+Mdatafile+"'");
			
		getline(Min,line);
		
		nM.resize(narea); Mto.resize(narea); Mval.resize(narea);
		do{
			getline(Min,line);
			if(Min.eof()) break;
			
			stringstream ss(line);	
			ss >> a1 >> a2 >> v;
			Mto[a1].push_back(a2); Mval[a1].push_back(v);
			if(a1 != a2){ Mto[a2].push_back(a1); Mval[a2].push_back(v);}
		}while(1 == 1);
		
		for(c = 0; c < narea; c++) nM[c] = Mto[c].size();

		ifstream Nin(Ndatafile.c_str());                                   // Loads information about age mixing
		if(!Nin) emsg("Cannot find the file '"+Ndatafile+"'");
	
		N.resize(nage);
		for(j = 0; j < nage; j++){
			getline(Nin,line);
			stringstream ss(line);	
			for(jj = 0; jj < nage; jj++){
				ss >> v;
				N[j].push_back(v);
			}
		}		
		
		cout << endl << "Age contact structure data loaded:" << endl;
		for(j = 0; j < nage; j++){
			for(jj = 0; jj < nage; jj++) cout << N[j][jj] << ","; cout << " N" << endl;
		}
		
		if(mode != MODE_SIM){        // Loads transition data for inference
			if(transdata.size() == 0) emsg("Transition data must be loaded");
			
			for(td = 0; td < transdata.size(); td++){
				ifstream transin(transdata[td].file.c_str());                            // Loads the transition data
				if(!transin) emsg("Cannot find the file '"+transdata[td].file+"'");
				
				nreg = 0;                                                                // Checks regions names match with area file
				getline(transin,line);
				stringstream ss(line);
				getline(ss,ele,'\t');
				while(!ss.eof()){
					getline(ss,ele,'\t');
					if(ele != region[nreg].name) emsg("Region names in file '"+transdata[td].file+"' do not agree with '"+regiondatafile+"'");
					nreg++;
				}
				
				t = 0;
				transdata[td].num.resize(nregion);
				do{
					getline(transin,line);
					if(transin.eof()) break;
					stringstream ss(line);
					getline(ss,ele,'\t');
					for(r = 0; r < nregion; r++){
						getline(ss,ele,'\t');
						transdata[td].num[r].push_back(atoi(ele.c_str()));
					}
					t++;
				}while(1 == 1);
				if(t != period){
					stringstream ss; ss << "The file 'transdata[td].file' has " << t << "' instead of period='" << period << "' rows";
					emsg(ss.str());
				}
			}
		}
	}

	if(ncore > 1) copydata(core);

	for(c = 0; c < narea; c++){                                              // Adds individuals to the system
		area[c].ind.resize(ndemocatpos);
		for(dp = 0; dp < ndemocatpos; dp++){
			imax = area[c].pop[dp];	
			for(i = 0; i < imax; i++){
				area[c].ind[dp].push_back(ind.size());
				
				indi.area = c;
				indi.dp = dp;
				ind.push_back(indi);
			}
		}
	}
	popsize = ind.size();

	narage = narea*nage;                                              // Generates the mixing matrix between ages/areas
	nardp = narea*ndemocatpos; 
	nsettardp = nsettime*nardp;

	ntimeperiod = 2;
	timeperiod.push_back(period/2); timeperiod.push_back(large);      // lockdown happens half way
	Qnum = 6;
	
	Qcomp.resize(Qnum); Qtimeperiod.resize(Qnum); 
	nQ.resize(Qnum); Qto.resize(Qnum); Qval.resize(Qnum);
	for(q = 0; q < Qnum; q++){
		switch(q){
			case 0: Qcomp[q] = "I"; Qtimeperiod[q] = 0; break;
			case 1: Qcomp[q] = "I"; Qtimeperiod[q] = 1; break;
			case 2: Qcomp[q] = "P"; Qtimeperiod[q] = 0; break;
			case 3: Qcomp[q] = "P"; Qtimeperiod[q] = 1; break;
			case 4: Qcomp[q] = "A"; Qtimeperiod[q] = 0; break;
			case 5: Qcomp[q] = "A"; Qtimeperiod[q] = 1; break;
		}
		
		nQ[q].resize(narage); Qto[q].resize(narage); Qval[q].resize(narage);
		for(c = 0; c < narea; c++){
			for(a = 0; a < nage; a++){
				vi = c*nage + a;
			
				for(j = 0; j < nM[c]; j++){
					cc = Mto[c][j]; v = Mval[c][j]; 
					if(Qtimeperiod[q] == 1) v = sqrt(v);                       // Distribution changed after lockdown
				
					k = Qto[q][vi].size();
					Qto[q][vi].push_back(cc);
					Qval[q][vi].push_back(vector <double> ());
					for(aa = 0; aa < nage; aa++) Qval[q][vi][k].push_back(N[a][aa]*v);
				}
			}
		}
		
		for(vi = 0; vi < narage; vi++) nQ[q][vi] = Qto[q][vi].size();
	                                              
		for(vi = 0; vi < narage; vi++){                                  // Normalisation
			sum = 0; 
			for(j = 0; j < nQ[q][vi]; j++){
				c = Qto[q][vi][j];
				for(a = 0; a < nage; a++) sum += Qval[q][vi][j][a]*area[c].agepop[a];
			}
			for(j = 0; j < nQ[q][vi]; j++){
				for(a = 0; a < nage; a++) Qval[q][vi][j][a] /= sum;
			}
		}
	}
}

/// Copies data from core zero to all the others
void DATA::copydata(unsigned int core)
{
	unsigned int td;
	int si;
	
	if(core == 0){                                  				   // Copies the above information to all the other cores
		packinit();
		pack(ndemocat);
		pack(democat);
		pack(ndemocatpos);
		pack(democatpos);
		pack(nregion);
		pack(region);
		pack(narea);
		pack(area);
		pack(nM);
		pack(Mto);
		pack(Mval);
		pack(N);
		pack(nage);
		pack(ndemocatposperage);
		for(td = 0; td < transdata.size(); td++) pack(transdata[td].num);
		si = packsize();
	}
	
	MPI_Bcast(&si,1,MPI_UNSIGNED,0,MPI_COMM_WORLD);
	MPI_Bcast(packbuffer(),si,MPI_DOUBLE,0,MPI_COMM_WORLD);
		
	if(core != 0){
		packinit();
		unpack(ndemocat);
		unpack(democat);
		unpack(ndemocatpos);
		unpack(democatpos);
		unpack(nregion);
		unpack(region);
		unpack(narea);
		unpack(area);
		unpack(nM);
		unpack(Mto);
		unpack(Mval);
		unpack(N);
		unpack(nage);
		unpack(ndemocatposperage);
		for(td = 0; td < transdata.size(); td++) unpack(transdata[td].num);
		if(si != packsize()) emsg("Data: EC9");
	}
}

/*
/// Calculates density of houses for a given house
void DATA::housedensity()
{
	int L, h, hh, i, j, imin, imax, jmin, jmax, num, c, k;
	double x, y, xx, yy, dd, fac = 0.9999999, val;
	vector< vector <int> > grid;
	
	L = int(scale/(2*rden));

	grid.resize(L*L);
	for(h = 0; h < nhouse; h++){
		c = int(fac*house[h].y*L)*L + int(fac*house[h].x*L); if(c < 0 || c >= L*L) emsg("Data: EC11");
		
		grid[c].push_back(h);
	}
	
	for(h = 0; h < nhouse; h++){
		x = house[h].x; y = house[h].y;
		imin = int(x*L+0.5)-1; imax = imin+1; if(imin < 0) imin = 0; if(imax >= L) imax = L-1;
		jmin = int(y*L+0.5)-1; jmax = jmin+1; if(jmin < 0) jmin = 0; if(jmax >= L) jmax = L-1;
		
		num = 0;
		for(i = imin; i < imax; i++){
			for(j = jmin; j < jmax; j++){
				c = j*L+i;
				for(k = 0; k < grid[c].size(); k++){
					hh = grid[c][k];
					xx = house[hh].x; yy = house[hh].y;
					dd = (xx-x)*(xx-x) + (yy-y)*(yy-y);
					if(dd < rden*rden/(scale*scale)) num++;
				}
			}
		}
		val = num/(2*3.14159*rden*rden); if(val < 0.5) val = 0.5;
		house[h].density = val;
	}
}
*/

string DATA::strip(string line)
{
	int len = line.length();
	if(len > 0){
		if(line.substr(len-1,1) == "\r") line = line.substr(0,len-1);
	}
	return line;
}	

void DATA::sortX(vector <unsigned int> &vec){ sort(vec.begin(),vec.end(),compX);}
void DATA::sortY(vector <unsigned int> &vec){ sort(vec.begin(),vec.end(),compY);}
