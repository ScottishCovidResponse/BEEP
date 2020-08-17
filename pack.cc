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

namespace {
	unsigned int k;

	std::vector<double> buffer;
}

void packinit(size_t size)
{
	k = 0;
	buffer.resize(size);
}

int packsize()
{
	return k;
}

double * packbuffer()
{
	return &buffer[0];
}

void pack(const unsigned int num)
{
	buffer.push_back(num); k++; 
}

void pack(const unsigned short num)
{
	buffer.push_back(num); k++; 
}

void pack(const double num)
{
	buffer.push_back(num); k++; 
}

void pack(const string &vec)
{
	unsigned int jmax, j;
	jmax = vec.length(); buffer.push_back(jmax); k++; 
	for(j = 0; j < jmax; j++){
		buffer.push_back(vec.at(j)); k++;
	}
}

template <class T>
void pack_vector(const vector<T>& vec)
{
	pack(vec.size());
	for (auto& item : vec) {
		pack(item);
	}
}

void pack(const vector <unsigned int> &vec)
{
	unsigned int imax, i;
	imax = vec.size(); buffer.push_back(imax); k++;
	for(i = 0; i < imax; i++) {
		buffer.push_back(vec[i]); k++;
	}
}

void pack(const vector <unsigned short> &vec)
{
	unsigned short imax, i; 
	imax = vec.size(); buffer.push_back(imax); k++;
	for(i = 0; i < imax; i++) {
		buffer.push_back(vec[i]); k++;
	}
}

void pack(const vector <int> &vec)
{
	unsigned int imax, i;
	imax = vec.size(); buffer.push_back(imax); k++;
	for(i = 0; i < imax; i++){
		buffer.push_back(vec[i]); k++;
	}
}

void pack(const vector <double> &vec)
{	
	unsigned int imax, i;
	imax = vec.size(); buffer.push_back(imax); k++;
	for(i = 0; i < imax; i++) {
		buffer.push_back(vec[i]); k++;
	}
}

void pack(const vector< vector <unsigned int> > &vec)
{
	unsigned int imax, i, jmax, j;
	imax = vec.size(); buffer.push_back(imax); k++; 
	for(i = 0; i < imax; i++){
		jmax = vec[i].size(); buffer.push_back(jmax); k++;
		for(j = 0; j < jmax; j++){
			buffer.push_back(vec[i][j]); k++;
		}
	}
}

void pack(const vector< vector <double> > &vec)
{
	unsigned int imax, i, jmax, j;
	imax = vec.size(); buffer.push_back(imax); k++; 
	for(i = 0; i < imax; i++){
		jmax = vec[i].size(); buffer.push_back(jmax); k++;
		for(j = 0; j < jmax; j++) {
			buffer.push_back(vec[i][j]); k++;
		}
	}
}

void pack(const vector< vector <float> > &vec)
{
	unsigned int imax, i, jmax, j;
	imax = vec.size(); buffer.push_back(imax); k++;
	for(i = 0; i < imax; i++){
		jmax = vec[i].size(); buffer.push_back(jmax); k++;
		for(j = 0; j < jmax; j++){
			buffer.push_back(vec[i][j]); k++;
		}
	}
}

void pack(const vector< vector< vector <double> > > &vec)
{
	unsigned int imax, i, jmax, j, nmax, n;
	imax = vec.size(); buffer.push_back(imax); k++;
	for(i = 0; i < imax; i++){
		jmax = vec[i].size(); buffer.push_back(jmax); k++;
		for(j = 0; j < jmax; j++){
			nmax = vec[i][j].size(); buffer.push_back(nmax); k++;
			for(n = 0; n < nmax; n++){
				buffer.push_back(vec[i][j][n]); k++;
			}
		}
	}
}

void pack(const vector <string> &vec)
{
	unsigned int imax, i, jmax, j;
	imax = vec.size(); buffer.push_back(imax); k++;
	for(i = 0; i < imax; i++){
		jmax = vec[i].length(); buffer.push_back(jmax); k++;
		for(j = 0; j < jmax; j++){
			buffer.push_back(vec[i].at(j)); k++;
		}
	}
}

void pack(const vector< vector <FEV> > &vec, unsigned int fedivmin, unsigned int fedivmax)
{
	unsigned int imax, i, jmax, j;
	imax = vec.size(); buffer.push_back(imax); k++;
	for(i = fedivmin; i < fedivmax; i++){
		jmax = vec[i].size(); buffer.push_back(jmax); k++;
		for(j = 0; j < jmax; j++){
			buffer.push_back(vec[i][j].trans); k++;
			buffer.push_back(vec[i][j].ind); k++;
			buffer.push_back(vec[i][j].t); k++;
		}
	}
}

void pack(const vector <AREA> &vec)
{
	unsigned int imax, i;
	imax = vec.size(); buffer.push_back(imax); k++;
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
	imax = vec.size(); buffer.push_back(imax); k++;
	for(i = 0; i < imax; i++){
		pack(vec[i].name);
		pack(vec[i].code);
	}
}

void pack(const vector <DEMOCAT> &vec)
{
	unsigned int imax, i;
	imax = vec.size(); buffer.push_back(imax); k++;
	for(i = 0; i < imax; i++){
		pack(vec[i].name);
		pack(vec[i].value);
	}
}

void pack(const vector <vector <EVREF> > &vec)
{
	unsigned int imax, i, jmax, j;
	imax = vec.size(); buffer.push_back(imax); k++;
	for(i = 0; i < imax; i++){
		jmax = vec[i].size(); buffer.push_back(jmax); k++;
		for(j = 0; j < jmax; j++){
			buffer.push_back(vec[i][j].ind); k++;
			buffer.push_back(vec[i][j].e); k++;
		}
	}
}

void pack(const vector < vector <vector <unsigned int> > > &vec)
{
	unsigned int i, imax;
	imax = vec.size(); buffer.push_back(imax); k++;
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
