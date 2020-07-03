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
	unsigned int t, r, i, c, imax, k, nreg, td, j, jmax, jj, cc, fl, d, dp, a1, a2, a, aa, vi, q, s, row;
	unsigned int namecol, codecol, xcol, ycol, regcol;
	int dc;
	double v;
	string line, ele, name, regcode, st, file;
	REGION reg;
	AREA are;
	DEMOCAT dem;
	IND indi;
	TABLE tab;
	vector <unsigned int> count;
	vector <vector <double> > val;
	vector <double> vec;
	
	mode = mod;
	period = per;
	fepertime = 10;
	
	settpertime = 1;
	nsettime = settpertime*per;
	settime.resize(nsettime+1);
	for(s = 0; s <= nsettime; s++) settime[s] = double(s*period)/nsettime;
			
	fediv = nsettime*fepertime;
	
	if(core == 0){
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
 
		tab = loadtable(datadir+"/"+regiondatafile);
		namecol = findcol(tab,"Name");
		codecol = findcol(tab,"Code");
		
		for(row = 0; row < tab.nrow; row++){
			reg.name = tab.ele[row][namecol];
			reg.code = tab.ele[row][codecol];
			region.push_back(reg);
		}		
		nregion = region.size();
		
		cout << endl << "Region data loaded." << endl;
		if(checkon == 1){
			for(r = 0; r < nregion; r++) cout << region[r].code << ", " << region[r].name  << " regionload" << endl;
		}
	
		file = datadir+"/"+areadatafile;
		tab = loadtable(file);

		codecol = findcol(tab,"area");                                                 // Works out ccolumns for different columns
		xcol = findcol(tab,"x");
		ycol = findcol(tab,"y");
		regcol = findcol(tab,"region");
		
		for(j = 0; j < ncovar; j++) covar[j].col = findcol(tab,covar[j].name);

		for(j = 0; j < ndemocat; j++){
			democat[j].col.resize(democat[j].value.size());
			for(k = 0; k < democat[j].col.size(); k++) democat[j].col[k] = findcol(tab,democat[j].value[k]);  
		}

		for(row = 0; row < tab.nrow; row++){
			are.code = tab.ele[row][codecol];
			are.x = atof(tab.ele[row][xcol].c_str());
			are.y = atof(tab.ele[row][ycol].c_str());
	
			regcode = tab.ele[row][regcol];
			r = 0; while(r < nregion && region[r].code != regcode) r++;
			if(r == nregion) emsg("Region code not recognised: ",regcode);
			are.region = r;
					
			are.covar.resize(ncovar);
			for(j = 0; j < ncovar; j++){
				st = tab.ele[row][covar[j].col];
				v = atof(st.c_str());
				if(std::isnan(are.covar[j])) emsg("In file '"+areadatafile+"' the expression '"+st+"' is not a number");	
			
				if(covar[j].func == "log"){
					if(v <= 0) emsg("Log transformed quantities must be positive.");
					are.covar[j] = log(v);
				}
				else{
					if(covar[j].func == "linear") are.covar[j] = v;
					else emsg("The functional relationship '"+covar[j].func+"' is not recognised.");
				}
			}
			
			val.resize(democat.size());
			for(d = 0; d < democat.size(); d++){
				jmax = democat[d].value.size();
				val[d].resize(jmax);
				for(j = 0; j < jmax; j++){
					st = tab.ele[row][democat[d].col[j]];
					val[d][j] = atof(st.c_str());
					if(std::isnan(val[d][j])) emsg("In file '"+areadatafile+"' the expression '"+st+"' is not a number");	
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
		}
		narea = area.size();
		
		cout << endl << "Area data loaded." << endl;
		if(checkon == 1){
			for(c = 0; c < narea; c++){
				cout << nregion << " " << area[c].region << "region" << endl;
				cout << area[c].code << " " << region[area[c].region].code << " " <<  area[c].x << " " <<  area[c].y << "  ***";
			
				for(j = 0; j < area[c].pop.size(); j++) cout << area[c].pop[j] << ", ";
				cout << endl;	
			}
		}			
		
		if(mode != MODE_SIM){        // Loads transition data for inference
			if(transdata.size() == 0) emsg("Transition data must be loaded");
			
			for(td = 0; td < transdata.size(); td++){
				file = datadir+"/"+transdata[td].file;

				ifstream transin(file.c_str());            // Loads the transition data
				if(!transin) emsg("Cannot open the file '"+file+"'");
				
				nreg = 0;                                                            // Checks regions names match with area file
				getline(transin,line);
				
				if(transdata[td].type == "reg"){
					stringstream ss(line);
					getline(ss,ele,'\t');
					while(!ss.eof()){
						getline(ss,ele,'\t');
						if(ele != region[nreg].name){
							emsg("Region names in file '"+transdata[td].file+"' do not agree with '"+regiondatafile+"'");
						}
						nreg++;
					}
				}
				
				row = 0;
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
					row++;
				}while(1 == 1);
				transdata[td].rows = row;
				
				if(transdata[td].start + row*transdata[td].units > period){
					emsg("The file '"+file+"' has more data than will fit in the time period.");
				}
			}
		}
		
		vec.resize(nage);                                                 // Reads in Q tensors
		for(q = 0; q < Q.size(); q++){
			tab = loadtable(datadir+"/"+Q[q].file);
			
			Q[q].to.resize(narea*nage); Q[q].val.resize(narea*nage);
			for(row = 0; row < tab.nrow; row++){
				c = atoi(tab.ele[row][0].c_str()); cc = atoi(tab.ele[row][1].c_str());
				
				for(a = 0; a < nage; a++){
					for(aa = 0; aa < nage; aa++) vec[aa] = atof(tab.ele[row][2+a*nage+aa].c_str()); 
					
					Q[q].to[c*nage+a].push_back(cc);
					Q[q].val[c*nage+a].push_back(vec);
					
					
					if(c < cc){
						Q[q].to[cc*nage+a].push_back(c);
						Q[q].val[cc*nage+a].push_back(vec);
					}
				}
			}
			
			normaliseQ(q);
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
}

/// Normalises the Q tensors
void DATA::normaliseQ(unsigned int q)
{
	unsigned int c, a, cc, aa, vi, j, jmax;
	double sum, sum2;
	
	for(c = 0; c < narea; c++){                                       // Normalisation
		sum = 0; sum2 = 0;
		for(a = 0; a < nage; a++){
			sum2 += area[c].agepop[a];
			vi = c*nage + a;
		
			jmax = Q[q].to[vi].size();
			for(j = 0; j < jmax; j++){
				cc = Q[q].to[vi][j];
				for(aa = 0; aa < nage; aa++) sum += area[c].agepop[a]*Q[q].val[vi][j][aa]*area[cc].agepop[aa];
			}
		}
		sum /= sum2;

		for(a = 0; a < nage; a++){
			vi = c*nage + a;
			jmax = Q[q].to[vi].size();
			for(j = 0; j < jmax; j++){
				for(aa = 0; aa < nage; aa++) Q[q].val[vi][j][aa] /= sum;
			}
		}
	}
}

/// Adds a time period
void DATA::addtimep(string name, double tend)
{
	TIMEP timep;
	
	timep.name = name;
	timep.tend = tend;
	timeperiod.push_back(timep);
}
	
/// Adds a Q tensor
void DATA::addQtensor(string timep, string comp, string file)
{
	unsigned int tp;
	QTENSOR qten;
	
	tp = 0; while(tp < timeperiod.size() && timeperiod[tp].name != timep) tp++;
	if(tp == timeperiod.size()) emsg("Cannot find '"+timep+"' as a time period");
	
	qten.timep = tp;
	qten.comp = comp;
	qten.file = file;
	Q.push_back(qten);
}
	
/// Loads a table from a file
TABLE DATA::loadtable(string file)
{
	TABLE tab;
	string line, st;
	vector <string> vec;
	
	ifstream in(file.c_str());                             // Loads information about areas
	if(!in) emsg("Cannot open the file '"+file+"'");
		
	getline(in,line);

	stringstream ss(line);
	do{
		getline(ss,st,'\t'); st = strip(st);
		tab.heading.push_back(st);
		if(ss.eof()) break;
	}while(1 == 1);
	tab.ncol = tab.heading.size();
	
	do{
		vec.clear();
		getline(in,line);
		if(in.eof()) break;
				
		stringstream ss(line);
		do{
			getline(ss,st,'\t'); st = strip(st);
			vec.push_back(st);
			if(ss.eof()) break;
		}while(1 == 1);
		if(vec.size() != tab.ncol) emsg("Rows in file '"+file+"' do not all have the same length.");
		
		tab.ele.push_back(vec);
	}while(1 == 1);
	tab.nrow = tab.ele.size();
	
	return tab;
}
		
/// Finds a column in a table
unsigned int DATA::findcol(TABLE &tab, string name)
{
	unsigned int c;
	
	for(c = 0; c < tab.ncol; c++) if(tab.heading[c] == name) break;
	if(c == tab.ncol) emsg("Cannot find the column heading '"+name+"'.");
	return c;
}		
			
/// Copies data from core zero to all the others
void DATA::copydata(unsigned int core)
{
	unsigned int td, q;
	int si;
	
	if(core == 0){                                  				   // Copies the above information to all the other cores
		packinit();
		pack(ndemocatpos);
		pack(democatpos);
		pack(nregion);
		pack(region);
		pack(narea);
		pack(area);
		for(q = 0; q < Q.size(); q++){ pack(Q[q].to); pack(Q[q].val);}
		pack(nage);
		pack(ndemocatposperage);
		for(td = 0; td < transdata.size(); td++){
			pack(transdata[td].num);
			pack(transdata[td].units);
			pack(transdata[td].rows);
		}
		si = packsize();
	}
	
	MPI_Bcast(&si,1,MPI_UNSIGNED,0,MPI_COMM_WORLD);
	MPI_Bcast(packbuffer(),si,MPI_DOUBLE,0,MPI_COMM_WORLD);
		
	if(core != 0){
		packinit();
		unpack(ndemocatpos);
		unpack(democatpos);
		unpack(nregion);
		unpack(region);
		unpack(narea);
		unpack(area);
		for(q = 0; q < Q.size(); q++){ unpack(Q[q].to); unpack(Q[q].val);}
		unpack(nage);
		unpack(ndemocatposperage);
		for(td = 0; td < transdata.size(); td++){
			unpack(transdata[td].num);
			unpack(transdata[td].units);
			unpack(transdata[td].rows);
		}
		if(si != packsize()) emsg("Data: EC9");
	}
}

/// Adds demographic categories
void DATA::adddemocat(string name, vector <string> &st, vector <string> &params)
{
	unsigned int j;
	DEMOCAT dem;
	
	dem.name = name;	
	dem.value = st;
	dem.param = params;
	
	democat.push_back(dem);

	ndemocat = democat.size();
	if(ndemocat == 1) nage = st.size();
}
	
/// Add a covariate for the areas
void DATA::addcovar(string name, string param, string func)
{
	COVAR cov;
	cov.name = name;
	cov.param = param;
	cov.func = func;
	
	covar.push_back(cov);
	ncovar = covar.size();
}
	
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
