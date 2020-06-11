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

/// Reads in transition and house data
void DATA::readdata(int core, int mod, int per)
{
	int t, tt, r, nu, h, hh, si, i, k, nreg, td;
	string line, ele, name;
	HOUSE ho;
	vector <int> regionpop;
	
	mode = mod;
	period = per;
	fediv = 100*per;
	
	if(core == 0){
		int RX, RY;
		
		if(mode == MODE_SIM){      
			if(simtype == "smallsim"){ RX = 2; RY = 2; popsize = 10000; nhouse = 1024;}  
			else{
				if(simtype == "scotsim"){ RX = 4; RY = 4; popsize = 5500000; nhouse = 1500000;}
				else{
					if(simtype == "uksim"){ RX = 10; RY = 10; popsize = 68000000; nhouse = 20000000;}
					else emsg("The property 'simtype' is wrong");
				}
			}
	
			nregion = RX*RY;
				
			for(r = 0; r < nregion; r++){
				stringstream ss; ss << "region" << r%RX << "_" << r/RX; regionname.push_back(ss.str());
			}
	
			for(h = 0; h < nhouse; h++){                            // Randomly distributes houses
				ho.x = ran();
				ho.y = ran();
				ho.region = int(ho.y*RY)*RX + int(ho.x*RX);
				house.push_back(ho);
			}
		}
		
		if(mode != MODE_SIM){            // Loads houses	
			ifstream housein(housefile.c_str());
			if(housein.fail()){
				stringstream ss; ss << "The file '" << housefile << "' does not exist";
				emsg(ss.str());
			}
			
			getline(housein,line); 
			stringstream ssh(line);	
			ssh >> popsize >> nhouse >> nregion;
			for(r = 0; r < nregion; r++){
				getline(housein,name); name = name.substr(0,name.length());
				regionname.push_back(name);
			}

			regionpop.resize(nregion); for(r =0; r < nregion; r++) regionpop[r] = 0;
			for(h = 0; h < nhouse; h++){
				housein >> ho.x >> ho.y >> ho.region;
				house.push_back(ho);
				regionpop[ho.region]++;
			}
			
			for(r = 0; r < nregion; r++) cout << regionname[r] << "  # Houses: " << regionpop[r] << endl;
		}
		
		if(mode != MODE_SIM){        // Loads transition data for inference
			if(transdata.size() == 0) emsg("Transition data must be loaded");
			
			for(td = 0; td < transdata.size(); td++){
				ifstream transin(transdata[td].file.c_str());                            // Loads the transition data
			
				if(transin.fail()){
					stringstream ss; ss << "The file '" << housefile << "' does not exist";
					emsg(ss.str());
				}
			
				nreg = 0;                                                                // Checks regions names match with house file
				getline(transin,line);
				stringstream ss(line);
				getline(ss,ele,'\t');
				while(!ss.eof()){
					getline(ss,ele,'\t');
					if(ele != regionname[nreg]) emsg("Data: EC12");
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
		
		for(h = 0; h < nhouse; h++) house[h].ind.push_back(h);   // Randomly distributes individuals into houses
    
		for(i = nhouse; i < popsize; i++){
			h = int(ran()*nhouse);
			house[h].ind.push_back(i);
		}
	}
	
	if(core == 0){                                  				   // Copies the above information to all the other cores
		packinit();
		for(td = 0; td < transdata.size(); td++) pack(transdata[td].num);
		pack(popsize);	
		pack(nregion);
		pack(regionname);
		pack(nhouse);
		pack(regionpop);
		si = packsize();
	}
	
	MPI_Bcast(&si,1,MPI_INT,0,MPI_COMM_WORLD);
	MPI_Bcast(packbuffer(),si,MPI_DOUBLE,0,MPI_COMM_WORLD);
		
	if(core != 0){
		packinit();
		for(td = 0; td < transdata.size(); td++) unpack(transdata[td].num);
		unpack(popsize);
		unpack(nregion);
		unpack(regionname);
		unpack(nhouse);
		unpack(regionpop);
		if(si != packsize()) emsg("Data: EC9");
	}

	h = 0; 
	do{                                                      // Houses are copied in pieces because they fill up the buffer
		hh = h + 10000; if(hh > nhouse) hh = nhouse;
		if(core == 0){ packinit(); pack(house,h,hh); si = packsize();}
		
		MPI_Bcast(&si,1,MPI_INT,0,MPI_COMM_WORLD);
		MPI_Bcast(packbuffer(),si,MPI_DOUBLE,0,MPI_COMM_WORLD);

		if(core != 0){ packinit(); unpack(house,h,hh);}
		h = hh;
	}while(h < nhouse);
		
	housedensity();
}

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

/// Check that the transition data is correct
void DATA::checktransdata(MODEL &model)
{
	short td, tra;
	string from, to;
	TRANS tr;
	
	for(td = 0; td < transdata.size(); td++){
		from = transdata[td].from; to = transdata[td].to; 
		for(tra = 0; tra < model.trans.size(); tra++){
			tr = model.trans[tra];
			if(model.comp[tr.from].name == from && model.comp[tr.to].name == to) break;
		}
		
		if(tra == model.trans.size()){
			stringstream ss; ss << "Cannot find the transition " << from << "→" << to << ".";
			emsg(ss.str());
		}
	}
}

void DATA::sortX(vector <int> &vec){ sort(vec.begin(),vec.end(),compX);}
void DATA::sortY(vector <int> &vec){ sort(vec.begin(),vec.end(),compY);}
