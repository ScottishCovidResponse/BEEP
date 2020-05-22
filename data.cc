//  The data inputs

#include <vector>
#include <fstream>
#include <sstream>
#include <algorithm>

using namespace std;

#include "utils.hh"
#include "data.hh"	
#include "consts.hh"
	
/// Reads in simulated case data
void DATA::readdata(short siminf)
{
	long week, tt, r, nu, h, RX, RY;
	string line, ele;
	HOUSE ho;
		
	if(siminf == 1){
		tmax = 105;
		
		if(type == "small"){ RX = 1; RY = 1; popsize = 1000; nhouse = 300;}
		else{
			if(type == "scotsim"){ RX = 4; RY = 4; popsize = 5500000; nhouse = 1500000;}
			else{
				if(type == "uksim"){ RX = 10; RY = 10; popsize = 6800000; nhouse = 20000000;}
				else emsg("Type is wrong");
			}
		}
		
		nregion = RX*RY;
		       
		for(r = 0; r < nregion; r++){
			stringstream ss; ss << "R" << r%RX << "_" << r/RX; regionname.push_back(ss.str());
		}
		
		for(h = 0; h < nhouse; h++){                                // Randomly distributes houses
			ho.x = ran();
			ho.y = ran();
			ho.region = short(ho.y*RY)*RX + short(ho.x*RX);
			house.push_back(ho);
		}
	}
	else{
		stringstream sst; sst << "cases_" << type << ".txt";
		ifstream regplot(sst.str().c_str());
	
		nregion = 0;
		getline(regplot,line);
		stringstream ss(line);
		getline(ss,ele,'\t');
		while(!ss.eof()){
			getline(ss,ele,'\t');
			regionname.push_back(ele); nregion++;
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
			tmax += 7;
		}while(1 == 1);
		
		stringstream ssw; ssw << "houses_" << type << ".txt";
		ifstream houseout(ssw.str().c_str());
		houseout >> popsize >> nhouse;
		for(h = 0; h < nhouse; h++){
			houseout >> ho.x >> ho.y >> ho.region;
			house.push_back(ho);
		}
	}
	
	if(tmax%7 != 0) emsg("The time must be a multiple of weeks");
	nweek = tmax/7;
}

void DATA::sortX(vector <long> &vec){ sort(vec.begin(),vec.end(),compX);}
void DATA::sortY(vector <long> &vec){ sort(vec.begin(),vec.end(),compY);}
