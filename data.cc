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

/// Reads in case and house data
void DATA::readdata(short core, short siminf)
{
	long week, tt, r, nu, h, hh, si, i, k, nreg;
	string line, ele, name;
	HOUSE ho;
	vector <long> regionpop;
	
	tmax = 17*7;  // The current simulations look at 17 weeks worth of data (similar to the real datasets
	
	if(fediv != tmax*10) emsg("Data: EC0");
	
	if(core == 0){
		long RX, RY;
		
		if(type.substr(0,4) != "real" && siminf == 1){                 // Randomly generates houses and regions
			if(type == "small"){ RX = 2; RY = 1; popsize = 3000; nhouse = 1024;}  
			else{
				if(type == "med"){ RX = 3; RY = 3; popsize = 3*16384; nhouse = 16384;}
				else{
					if(type == "scotsim"){ RX = 4; RY = 4; popsize = 5500000; nhouse = 1500000;}
					else{
						if(type == "uksim"){ RX = 10; RY = 10; popsize = 6800000; nhouse = 20000000;}
						else emsg("Type is wrong");
					}
				}
			}
	
			nregion = RX*RY;
					
			for(r = 0; r < nregion; r++){
				stringstream ss; ss << "region" << r%RX << "_" << r/RX; regionname.push_back(ss.str());
			}
	
			for(h = 0; h < nhouse; h++){                            // Randomly distributes houses
				ho.x = ran();
				ho.y = ran();
				ho.region = short(ho.y*RY)*RX + short(ho.x*RX);
				house.push_back(ho);
			}
		}
		
		if(siminf == 0 || type.substr(0,4) == "real"){            // Loads houses
			stringstream ssw; ssw << "houses_" << type << ".txt";
			ifstream houseout(ssw.str().c_str());
				
			getline(houseout,line); 
			stringstream ssh(line);	
			ssh >> popsize >> nhouse >> nregion;
			for(r = 0; r < nregion; r++){
				getline(houseout,name); name = name.substr(0,name.length());
				regionname.push_back(name);
			}

			regionpop.resize(nregion); for(r =0; r < nregion; r++) regionpop[r] = 0;
			for(h = 0; h < nhouse; h++){
				houseout >> ho.x >> ho.y >> ho.region;
				house.push_back(ho);
				regionpop[ho.region]++;
			}
			
			for(r = 0; r < nregion; r++) cout << regionname[r] << "  # Houses: " << regionpop[r] << endl;
		}
		
		if(siminf == 0){  // Loads cases for inference
			stringstream sst; sst << "cases_" << type << ".txt";    // Loads the cases
			ifstream regplot(sst.str().c_str());
		
			nreg = 0;
			getline(regplot,line);
			stringstream ss(line);
			getline(ss,ele,'\t');
			while(!ss.eof()){
				getline(ss,ele,'\t');
				cout << ele << " " <<  regionname[nreg] << " j\n";
				if(ele != regionname[nreg]) emsg("Data: EC12");
				nreg++;
			}
			
			tmax = 0;
			ncase.resize(nregion);
			do{
				getline(regplot,line);
				if(regplot.eof()) break;
				stringstream ss(line);
				getline(ss,ele,'\t');
				for(r = 0; r < nregion; r++){
					getline(ss,ele,'\t');
					ncase[r].push_back(atoi(ele.c_str()));
				}
				tmax += timestep;
			}while(1 == 1);
		}
		
		for(h = 0; h < nhouse; h++) house[h].ind.push_back(h);   // Randomly distributes individuals into houses
    
		for(i = nhouse; i < popsize; i++){
			h = long(ran()*nhouse);
			house[h].ind.push_back(i);
		}	
	}
	
	// Copies the above information to all the other cores
	
	if(core == 0){
		packinit();
		pack(popsize);
		pack(ncase);
		pack(nregion);
		pack(regionname);
		pack(nhouse);
		pack(regionpop);
		si = packsize();
	}
	
	MPI_Bcast(&si,1,MPI_LONG,0,MPI_COMM_WORLD);
	MPI_Bcast(packbuffer(),si,MPI_DOUBLE,0,MPI_COMM_WORLD);
		
	if(core != 0){
		packinit();
		unpack(popsize);
		unpack(ncase);
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
		
		MPI_Bcast(&si,1,MPI_LONG,0,MPI_COMM_WORLD);
		MPI_Bcast(packbuffer(),si,MPI_DOUBLE,0,MPI_COMM_WORLD);

		if(core != 0){ packinit(); unpack(house,h,hh);}
		h = hh;
	}while(h < nhouse);
		
	if(tmax%timestep != 0) emsg("The time must be a multiple of weeks");
	nweek = tmax/timestep;
	
	housedensity();
}

/// Calculates density of houses for a given house
void DATA::housedensity()
{
	long L, h, hh, i, j, imin, imax, jmin, jmax, num, c, k;
	double x, y, xx, yy, dd, fac = 0.9999999, val;
	vector< vector <long> > grid;
	
	L = short(scale/(2*rden));

	grid.resize(L*L);
	for(h = 0; h < nhouse; h++){
		c = short(fac*house[h].y*L)*L + short(fac*house[h].x*L); if(c < 0 || c >= L*L) emsg("Data: EC11");
		
		grid[c].push_back(h);
	}
	
	for(h = 0; h < nhouse; h++){
		x = house[h].x; y = house[h].y;
		imin = short(x*L+0.5)-1; imax = imin+1; if(imin < 0) imin = 0; if(imax >= L) imax = L-1;
		jmin = short(y*L+0.5)-1; jmax = jmin+1; if(jmin < 0) jmin = 0; if(jmax >= L) jmax = L-1;
		
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

void DATA::sortX(vector <long> &vec){ sort(vec.begin(),vec.end(),compX);}
void DATA::sortY(vector <long> &vec){ sort(vec.begin(),vec.end(),compY);}
