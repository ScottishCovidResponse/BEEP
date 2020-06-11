// These are function for packing up quantities to be sent by MPI

#include <vector>
#include <iostream>
#include <algorithm>
#include <string>
#include <sstream>  

using namespace std;

#include "PART.hh"
#include "consts.hh"
#include "utils.hh"

static int k;

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

void pack(int num)
{
	buffer[k] = num; k++; if(k == MAX_NUMBERS) emsg("Pack: EC1");
}

void pack(double num)
{
	buffer[k] = num; k++; if(k == MAX_NUMBERS) emsg("Pack: EC1");
}

void pack(vector <int> &vec)
{
	int imax, i; 
	imax = vec.size(); buffer[k] = imax; k++; if(k == MAX_NUMBERS) emsg("Pack: EC1"); 
	for(i = 0; i < imax; i++){ buffer[k] = vec[i]; k++; if(k == MAX_NUMBERS) emsg("Pack: EC1");}
}

void pack(vector <double> &vec)
{	
	int imax, i; 
	imax = vec.size(); buffer[k] = imax; k++; if(k == MAX_NUMBERS) emsg("Pack: EC1"); 
	for(i = 0; i < imax; i++){ buffer[k] = vec[i]; k++; if(k == MAX_NUMBERS) emsg("Pack: EC1");}
}

void pack(vector< vector <int> > &vec)
{
	int imax, i, jmax, j;
	imax = vec.size(); buffer[k] = imax; k++; if(k == MAX_NUMBERS) emsg("Pack: EC1"); 
	for(i = 0; i < imax; i++){
		jmax = vec[i].size(); buffer[k] = jmax; k++; if(k == MAX_NUMBERS) emsg("Pack: EC1");  
		for(j = 0; j < jmax; j++){ buffer[k] = vec[i][j]; k++; if(k == MAX_NUMBERS) emsg("Pack: EC1");}
	}
}

void pack(vector< vector <double> > &vec)
{
	int imax, i, jmax, j;
	imax = vec.size(); buffer[k] = imax; k++; if(k == MAX_NUMBERS) emsg("Pack: EC1"); 
	for(i = 0; i < imax; i++){
		jmax = vec[i].size(); buffer[k] = jmax; k++; if(k == MAX_NUMBERS) emsg("Pack: EC1");  
		for(j = 0; j < jmax; j++){ buffer[k] = vec[i][j]; k++; if(k == MAX_NUMBERS) emsg("Pack: EC1");}
	}
}

void pack(vector <string> &vec)
{
	int imax, i, jmax, j;
	imax = vec.size(); buffer[k] = imax; k++; if(k == MAX_NUMBERS) emsg("Pack: EC1"); 
	for(i = 0; i < imax; i++){
		jmax = vec[i].length(); buffer[k] = jmax; k++; if(k == MAX_NUMBERS) emsg("Pack: EC1");  
		for(j = 0; j < jmax; j++){ buffer[k] = vec[i].at(j); k++; if(k == MAX_NUMBERS) emsg("Pack: EC1");}
	}
}

void pack(vector< vector <FEV> > &vec, int fedivmin, int fedivmax)
{
	int imax, i, jmax, j;
	imax = vec.size(); buffer[k] = imax; k++; if(k == MAX_NUMBERS) emsg("Pack: EC1");  
	for(i = fedivmin; i < fedivmax; i++){
		jmax = vec[i].size(); buffer[k] = jmax; k++; if(k == MAX_NUMBERS) emsg("Pack: EC1"); 
		for(j = 0; j < jmax; j++){
			buffer[k] = vec[i][j].trans; k++; if(k == MAX_NUMBERS) emsg("Pack: EC1");  
			buffer[k] = vec[i][j].ind; k++; if(k == MAX_NUMBERS) emsg("Pack: EC1");  
			buffer[k] = vec[i][j].t; k++; if(k == MAX_NUMBERS) emsg("Pack: EC1");  
			buffer[k] = vec[i][j].done; k++; if(k == MAX_NUMBERS) emsg("Pack: EC1");  
		}
	}
}

void pack(vector <HOUSE> &vec, int min, int max)
{
	int imax, i, jmax, j;
	imax = vec.size(); buffer[k] = imax; k++; if(k == MAX_NUMBERS) emsg("Pack: EC1");  
	for(i = min; i < max; i++){
		buffer[k] = vec[i].x; k++; if(k == MAX_NUMBERS) emsg("Pack: EC1");  
		buffer[k] = vec[i].y; k++; if(k == MAX_NUMBERS) emsg("Pack: EC1"); 
		buffer[k] = vec[i].region; k++; if(k == MAX_NUMBERS) emsg("Pack: EC1"); 
		jmax = vec[i].ind.size(); buffer[k] = jmax; k++; if(k == MAX_NUMBERS) emsg("Pack: EC1");  
		for(j = 0; j < jmax; j++){ buffer[k] = vec[i].ind[j]; k++; if(k == MAX_NUMBERS) emsg("Pack: EC1");}  
	}
}

void unpack(int &num)
{
	num = buffer[k]; k++;
}

void unpack(double &num)
{
	num = buffer[k]; k++;
}

void unpack(vector <int> &vec)
{
	int imax, i; 
	imax = buffer[k]; k++; vec.resize(imax);
	for(i = 0; i < imax; i++){ vec[i] = buffer[k]; k++;}
}

void unpack(vector <double> &vec)
{
	int imax, i;
	imax = buffer[k]; k++; vec.resize(imax); 
	for(i = 0; i < imax; i++){ vec[i] = buffer[k]; k++;}
}

void unpack(vector< vector <int> > &vec)
{
	int imax, i, jmax, j;
	imax = buffer[k]; k++; vec.resize(imax);
	for(i = 0; i < imax; i++){
		jmax = buffer[k]; k++; vec[i].resize(jmax);
		for(j = 0; j < jmax; j++){ vec[i][j] = buffer[k]; k++;}
	}
}

void unpack(vector< vector <double> > &vec)
{
	int imax, i, jmax, j;
	imax = buffer[k]; k++; vec.resize(imax);
	for(i = 0; i < imax; i++){
		jmax = buffer[k]; k++; vec[i].resize(jmax);
		for(j = 0; j < jmax; j++){ vec[i][j] = buffer[k]; k++;}
	}
}

void unpack(vector <string> &vec)
{
	int imax, i, jmax, j;
	imax = buffer[k]; k++; vec.resize(imax);
	for(i = 0; i < imax; i++){
		jmax = buffer[k]; k++;
		stringstream ss; for(j = 0; j < jmax; j++){ ss << (char) buffer[k]; k++;}
		vec[i] = ss.str();
	}
}

void unpack(vector< vector <FEV> > &vec, int fedivmin, int fedivmax)
{
	int imax, i, jmax, j;
	imax = buffer[k]; k++; vec.resize(imax);
	for(i = fedivmin; i < fedivmax; i++){
		jmax = buffer[k]; k++; vec[i].resize(jmax);
		for(j = 0; j < jmax; j++){ 
			vec[i][j].trans = buffer[k]; k++;
			vec[i][j].ind = buffer[k]; k++;
			vec[i][j].t = buffer[k]; k++;
			vec[i][j].done = buffer[k]; k++;
		}
	}
}

void unpack(vector <HOUSE> &vec, int min, int max)
{
	int imax, i, jmax, j;
	imax = buffer[k]; k++; vec.resize(imax); 
	for(i = min; i < max; i++){
		vec[i].x = buffer[k]; k++;  
		vec[i].y = buffer[k]; k++; 
		vec[i].region = buffer[k]; k++; 
		jmax = buffer[k]; k++; vec[i].ind.resize(jmax);
		for(j = 0; j < jmax; j++){ vec[i].ind[j] = buffer[k]; k++;}  
	}
}
