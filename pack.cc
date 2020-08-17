// These are function for packing up quantities to be sent by MPI

#include <vector>
#include <iostream>
#include <algorithm>
#include <string>
#include <sstream>  

using namespace std;

#include "consts.hh"
#include "utils.hh"
#include "data.hh"
#include "model.hh"

static unsigned int k;

double buffer[MAX_NUMBERS];

void packinit()
{
	k = 0;
}

int packsize()
{
	return k;
}

double * packbuffer()
{
	return buffer;	
}

void pack(const unsigned int num)
{
	buffer[k] = num; k++; if(k == MAX_NUMBERS) emsgEC("Pack",1);
}

void pack(const unsigned short num)
{
	buffer[k] = num; k++; if(k == MAX_NUMBERS) emsgEC("Pack",2);
}

void pack(const double num)
{
	buffer[k] = num; k++; if(k == MAX_NUMBERS) emsgEC("Pack",3);
}

void pack(const string &vec)
{
	unsigned int jmax, j;
	jmax = vec.length(); buffer[k] = jmax; k++; if(k == MAX_NUMBERS) emsgEC("Pack",4);  
	for(j = 0; j < jmax; j++){ buffer[k] = vec.at(j); k++; if(k == MAX_NUMBERS) emsgEC("Pack",5);}
}

void pack(const vector <unsigned int> &vec)
{
	unsigned int imax, i; 
	imax = vec.size(); buffer[k] = imax; k++; if(k == MAX_NUMBERS) emsgEC("Pack",6); 
	for(i = 0; i < imax; i++){ buffer[k] = vec[i]; k++; if(k == MAX_NUMBERS) emsgEC("Pack",7);}
}

void pack(const vector <unsigned short> &vec)
{
	unsigned short imax, i; 
	imax = vec.size(); buffer[k] = imax; k++; if(k == MAX_NUMBERS) emsgEC("Pack",8); 
	for(i = 0; i < imax; i++){ buffer[k] = vec[i]; k++; if(k == MAX_NUMBERS) emsgEC("Pack",9);}
}

void pack(const vector <int> &vec)
{
	unsigned int imax, i; 
	imax = vec.size(); buffer[k] = imax; k++; if(k == MAX_NUMBERS) emsgEC("Pack",10); 
	for(i = 0; i < imax; i++){ buffer[k] = vec[i]; k++; if(k == MAX_NUMBERS) emsgEC("Pack",11);}
}

void pack(const vector <double> &vec)
{	
	unsigned int imax, i; 
	imax = vec.size(); buffer[k] = imax; k++; if(k == MAX_NUMBERS) emsgEC("Pack",12); 
	for(i = 0; i < imax; i++){ buffer[k] = vec[i]; k++; if(k == MAX_NUMBERS) emsgEC("Pack",13);}
}

void pack(const vector< vector <unsigned int> > &vec)
{
	unsigned int imax, i, jmax, j;
	imax = vec.size(); buffer[k] = imax; k++; if(k == MAX_NUMBERS) emsgEC("Pack",14); 
	for(i = 0; i < imax; i++){
		jmax = vec[i].size(); buffer[k] = jmax; k++; if(k == MAX_NUMBERS) emsgEC("Pack",15);  
		for(j = 0; j < jmax; j++){ buffer[k] = vec[i][j]; k++; if(k == MAX_NUMBERS) emsgEC("Pack",16);}
	}
}

void pack(const vector< vector <double> > &vec)
{
	unsigned int imax, i, jmax, j;
	imax = vec.size(); buffer[k] = imax; k++; if(k == MAX_NUMBERS) emsgEC("Pack",17); 
	for(i = 0; i < imax; i++){
		jmax = vec[i].size(); buffer[k] = jmax; k++; if(k == MAX_NUMBERS) emsgEC("Pack",18);  
		for(j = 0; j < jmax; j++){ buffer[k] = vec[i][j]; k++; if(k == MAX_NUMBERS) emsgEC("Pack",19);}
	}
}

void pack(const vector< vector <float> > &vec)
{
	unsigned int imax, i, jmax, j;
	imax = vec.size(); buffer[k] = imax; k++; if(k == MAX_NUMBERS) emsgEC("Pack",20); 
	for(i = 0; i < imax; i++){
		jmax = vec[i].size(); buffer[k] = jmax; k++; if(k == MAX_NUMBERS) emsgEC("Pack",21);  
		for(j = 0; j < jmax; j++){ buffer[k] = vec[i][j]; k++; if(k == MAX_NUMBERS) emsgEC("Pack",22);}
	}
}

void pack(const vector< vector< vector <double> > > &vec)
{
	unsigned int imax, i, jmax, j, nmax, n;
	imax = vec.size(); buffer[k] = imax; k++; if(k == MAX_NUMBERS) emsgEC("Pack",23); 
	for(i = 0; i < imax; i++){
		jmax = vec[i].size(); buffer[k] = jmax; k++; if(k == MAX_NUMBERS) emsgEC("Pack",24);  
		for(j = 0; j < jmax; j++){ 
			nmax = vec[i][j].size(); buffer[k] = nmax; k++; if(k == MAX_NUMBERS) emsgEC("Pack",25);  
			for(n = 0; n < nmax; n++){ buffer[k] = vec[i][j][n]; k++; if(k == MAX_NUMBERS) emsgEC("Pack",26);}
		}
	}
}

void pack(const vector <string> &vec)
{
	unsigned int imax, i, jmax, j;
	imax = vec.size(); buffer[k] = imax; k++; if(k == MAX_NUMBERS) emsgEC("Pack",27); 
	for(i = 0; i < imax; i++){
		jmax = vec[i].length(); buffer[k] = jmax; k++; if(k == MAX_NUMBERS) emsgEC("Pack",28);  
		for(j = 0; j < jmax; j++){ buffer[k] = vec[i].at(j); k++; if(k == MAX_NUMBERS) emsgEC("Pack",29);}
	}
}

void pack(const vector< vector <FEV> > &vec, unsigned int fedivmin, unsigned int fedivmax)
{
	unsigned int imax, i, jmax, j;
	imax = vec.size(); buffer[k] = imax; k++; if(k == MAX_NUMBERS) emsgEC("Pack",30);  
	for(i = fedivmin; i < fedivmax; i++){
		jmax = vec[i].size(); buffer[k] = jmax; k++; if(k == MAX_NUMBERS) emsgEC("Pack",31); 
		for(j = 0; j < jmax; j++){
			buffer[k] = vec[i][j].trans; k++; if(k == MAX_NUMBERS) emsgEC("Pack",32);  
			buffer[k] = vec[i][j].ind; k++; if(k == MAX_NUMBERS) emsgEC("Pack",33);  
			buffer[k] = vec[i][j].t; k++; if(k == MAX_NUMBERS) emsgEC("Pack",34);  
		}
	}
}

void pack(const vector <AREA> &vec)
{
	unsigned int imax, i;
	imax = vec.size(); buffer[k] = imax; k++; if(k == MAX_NUMBERS) emsgEC("Pack",35);  
	for(i = 0; i < imax; i++){
		pack(vec[i].code);
		pack(vec[i].region);
		pack(vec[i].x);
		pack(vec[i].y);
		pack(vec[i].agepop);
		pack(vec[i].pop);
		pack(vec[i].covar);             
	}
}

void pack(const vector <REGION> &vec)
{
	unsigned int imax, i;
	imax = vec.size(); buffer[k] = imax; k++; if(k == MAX_NUMBERS) emsgEC("Pack",36);  
	for(i = 0; i < imax; i++){
		pack(vec[i].name);
		pack(vec[i].code);
	}
}

void pack(const vector <DEMOCAT> &vec)
{
	unsigned int imax, i;
	imax = vec.size(); buffer[k] = imax; k++; if(k == MAX_NUMBERS) emsgEC("Pack",37);  
	for(i = 0; i < imax; i++){
		pack(vec[i].name);
		pack(vec[i].value);
	}
}

void pack(const vector <vector <EVREF> > &vec)
{
	unsigned int imax, i, jmax, j;
	imax = vec.size(); buffer[k] = imax; k++; if(k == MAX_NUMBERS) emsgEC("Pack",38);  
	for(i = 0; i < imax; i++){
		jmax = vec[i].size(); buffer[k] = jmax; k++; if(k == MAX_NUMBERS) emsgEC("Pack",39); 
		for(j = 0; j < jmax; j++){
			buffer[k] = vec[i][j].ind; k++; if(k == MAX_NUMBERS) emsgEC("Pack",40);  
			buffer[k] = vec[i][j].e; k++; if(k == MAX_NUMBERS) emsgEC("Pack",41);  
		}
	}
}

void pack(const vector < vector <vector <unsigned int> > > &vec)
{
	unsigned int i, imax;
	imax = vec.size(); buffer[k] = imax; k++; if(k == MAX_NUMBERS) emsgEC("Pack",44);  
	for(i = 0; i < imax; i++){
		pack(vec[i]);
	}
}

void unpack(unsigned int &num)
{
	num = buffer[k]; k++;
}

void unpack(unsigned short &num)
{
	num = buffer[k]; k++;
}

void unpack(double &num)
{
	num = buffer[k]; k++;
}

void unpack(string &vec)
{
	unsigned int jmax, j;
	
	jmax = buffer[k]; k++;
	stringstream ss; for(j = 0; j < jmax; j++){ ss << (char) buffer[k]; k++;}
	vec = ss.str();
}

void unpack(vector <unsigned int> &vec)
{
	unsigned int imax, i; 
	imax = buffer[k]; k++; vec.resize(imax);
	for(i = 0; i < imax; i++){ vec[i] = buffer[k]; k++;}
}

void unpack(vector <unsigned short> &vec)
{
	unsigned int imax, i; 
	imax = buffer[k]; k++; vec.resize(imax);
	for(i = 0; i < imax; i++){ vec[i] = buffer[k]; k++;}
}

void unpack(vector <int> &vec)
{
	unsigned int imax, i; 
	imax = buffer[k]; k++; vec.resize(imax);
	for(i = 0; i < imax; i++){ vec[i] = buffer[k]; k++;}
}

void unpack(vector <double> &vec)
{
	unsigned int imax, i;
	imax = buffer[k]; k++; vec.resize(imax); 
	for(i = 0; i < imax; i++){ vec[i] = buffer[k]; k++;}
}

void unpack(vector< vector <unsigned int> > &vec)
{
	unsigned int imax, i, jmax, j;
	imax = buffer[k]; k++; vec.resize(imax);
	for(i = 0; i < imax; i++){
		jmax = buffer[k]; k++; vec[i].resize(jmax);
		for(j = 0; j < jmax; j++){ vec[i][j] = buffer[k]; k++;}
	}
}

void unpack(vector< vector <double> > &vec)
{
	unsigned int imax, i, jmax, j;
	imax = buffer[k]; k++; vec.resize(imax);
	for(i = 0; i < imax; i++){
		jmax = buffer[k]; k++; vec[i].resize(jmax);
		for(j = 0; j < jmax; j++){ vec[i][j] = buffer[k]; k++;}
	}
}

void unpack(vector< vector <float> > &vec)
{
	unsigned int imax, i, jmax, j;
	imax = buffer[k]; k++; vec.resize(imax);
	for(i = 0; i < imax; i++){
		jmax = buffer[k]; k++; vec[i].resize(jmax);
		for(j = 0; j < jmax; j++){ vec[i][j] = buffer[k]; k++;}
	}
}

void unpack(vector< vector< vector <double> > > &vec)
{
	unsigned int imax, i, jmax, j, nmax, n;
	imax = buffer[k]; k++; vec.resize(imax);
	for(i = 0; i < imax; i++){
		jmax = buffer[k]; k++; vec[i].resize(jmax);
		for(j = 0; j < jmax; j++){ 
			nmax = buffer[k]; k++; vec[i][j].resize(nmax);
			for(n = 0; n < nmax; n++){ vec[i][j][n] = buffer[k]; k++;}
		}
	}
}

void unpack(vector <string> &vec)
{
	unsigned int imax, i, jmax, j;
	imax = buffer[k]; k++; vec.resize(imax);
	for(i = 0; i < imax; i++){
		jmax = buffer[k]; k++;
		stringstream ss; for(j = 0; j < jmax; j++){ ss << (char) buffer[k]; k++;}
		vec[i] = ss.str();
	}
}

void unpack(vector< vector <FEV> > &vec, unsigned int fedivmin, unsigned int fedivmax)
{
	unsigned int imax, i, jmax, j;
	imax = buffer[k]; k++; vec.resize(imax);
	for(i = fedivmin; i < fedivmax; i++){
		jmax = buffer[k]; k++; vec[i].resize(jmax);
		for(j = 0; j < jmax; j++){ 
			vec[i][j].trans = buffer[k]; k++;
			vec[i][j].ind = buffer[k]; k++;
			vec[i][j].t = buffer[k]; k++;
		}
	}
}

void unpack(vector <AREA> &vec)
{
	unsigned int imax, i;
	imax = buffer[k]; k++; vec.resize(imax); 
	for(i = 0; i < imax; i++){
		unpack(vec[i].code);
		unpack(vec[i].region);
		unpack(vec[i].x);
		unpack(vec[i].y);
		unpack(vec[i].agepop);
		unpack(vec[i].pop);
		unpack(vec[i].covar);  
	}
}

void unpack(vector <REGION> &vec)
{
	unsigned int imax, i;
	imax = buffer[k]; k++; vec.resize(imax); 
	for(i = 0; i < imax; i++){
		unpack(vec[i].name);
		unpack(vec[i].code);
	}
}

void unpack(vector <DEMOCAT> &vec)
{
	unsigned int imax, i;
	imax = buffer[k]; k++; vec.resize(imax); 
	for(i = 0; i < imax; i++){
		unpack(vec[i].name);
		unpack(vec[i].value);
	}
}

void unpack(vector <vector <EVREF> > &vec)
{
	unsigned int imax, i, jmax, j;
	imax = buffer[k]; k++; vec.resize(imax);
	for(i = 0; i < imax; i++){
		jmax = buffer[k]; k++; vec[i].resize(jmax);
		for(j = 0; j < jmax; j++){ 
			vec[i][j].ind = buffer[k]; k++;
			vec[i][j].e = buffer[k]; k++;
		}
	}
}

void unpack(vector < vector <vector <unsigned int> > > &vec)
{
	unsigned int i, imax;
	imax = buffer[k]; k++; vec.resize(imax);
	for(i = 0; i < imax; i++){
		unpack(vec[i]);
	}
}
