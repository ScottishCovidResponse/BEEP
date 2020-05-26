//  The data inputs

#include <iostream>
#include <vector>
#include <fstream>
#include <sstream>
#include <algorithm>

using namespace std;

#include "utils.hh"
#include "data.hh"	
#include "consts.hh"
	
/// Reads in simulated case data
void DATA::readdata(short core, short siminf)
{
	long week, tt, r, nu, h, RX, RY, si, i, k;
	string line, ele;
	HOUSE ho;
	vector <long> pac;
	
	if(type == "small"){ RX = 2; RY = 2; popsize = 3000; nhouse = 1024;}  // use 1024 areas
	else{
		if(type == "med"){ RX = 3; RY = 3; popsize = 3*16384; nhouse = 16384;} // use 16384
		else{
			if(type == "scotsim"){ RX = 4; RY = 4; popsize = 5500000; nhouse = 1500000;}
			else{
				if(type == "uksim"){ RX = 10; RY = 10; popsize = 6800000; nhouse = 20000000;}
				else emsg("Type is wrong");
			}
		}
	}
		
	tmax = 105;
	
	nregion = RX*RY;
				 
	for(r = 0; r < nregion; r++){
		stringstream ss; ss << "R" << r%RX << "_" << r/RX; regionname.push_back(ss.str());
	}
		
		
	if(siminf == 1){
		
		
		
	}
	else{
		if(core == 0){	
			stringstream sst; sst << "cases_" << type << ".txt";
			ifstream regplot(sst.str().c_str());
		
		/*
			nregion = 0;
			getline(regplot,line);
			stringstream ss(line);
			getline(ss,ele,'\t');
			while(!ss.eof()){
				getline(ss,ele,'\t');
				regionname.push_back(ele); nregion++;
			}
			*/
			
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
		
			pack(pac,ncase);
			si = pac.size();
		}
	
		MPI_Bcast(&si,1,MPI_LONG,0,MPI_COMM_WORLD);
		pac.resize(si);
		MPI_Bcast(&pac[0],si,MPI_LONG,0,MPI_COMM_WORLD);
		if(core != 0){
			k = 0;
			unpack(k,pac,ncase);
		}
	}
	
	
		/*
		stringstream ssw; ssw << "houses_" << type << ".txt";
		ifstream houseout(ssw.str().c_str());
		houseout >> popsize >> nhouse;
		*/
		/*
		for(h = 0; h < nhouse; h++){
			houseout >> ho.x >> ho.y >> ho.region;
			house.push_back(ho);
		}
		*/
	//}
	
	for(h = 0; h < nhouse; h++){                                // Randomly distributes houses
		ho.x = ran();
		ho.y = ran();
		ho.region = short(ho.y*RY)*RX + short(ho.x*RX);
		house.push_back(ho);
	}
	
	if(tmax%timestep != 0) emsg("The time must be a multiple of weeks");
	nweek = tmax/timestep;
}

void DATA::pack(vector <long> &pac, vector< vector <long> > &vec)
{
	long imax, i, jmax, j;
	imax = vec.size(); pac.push_back(imax); 
	for(i = 0; i < imax; i++){
		jmax = vec[i].size(); pac.push_back(jmax); for(j = 0; j < jmax; j++) pac.push_back(vec[i][j]);
	}
}

void DATA::unpack(long &k, vector <long> &pac, vector< vector <long> > &vec)
{
	long imax, i, jmax, j;
	imax = pac[k]; k++; vec.resize(imax);
	for(i = 0; i < imax; i++){
		jmax = pac[k]; k++; vec[i].resize(jmax); for(j = 0; j < jmax; j++){ vec[i][j] = pac[k]; k++;}
	}
}

void DATA::sortX(vector <long> &vec){ sort(vec.begin(),vec.end(),compX);}
void DATA::sortY(vector <long> &vec){ sort(vec.begin(),vec.end(),compY);}
