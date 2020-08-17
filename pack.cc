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

template <class T>
void pack_item(T t)
{
	buffer.push_back(t); k++;
}

// Provide overloads for cases where default behaviour needs
// to be modified
void pack_item(const string& vec)
{
	unsigned int jmax = vec.length();
	pack_item(jmax);
	for(unsigned int j = 0; j < jmax; j++){
		pack_item(vec.at(j));
	}
}

void pack(const unsigned int num)
{
	pack_item(num); 
}

void pack(const int num)
{
	pack_item(num); 
}

void pack(const unsigned short num)
{
	pack_item(num); 
}

void pack(const double num)
{
	pack_item(num); 
}

void pack(const string &vec)
{
	pack_item(vec);
}

// Use template to share implementation details between cases
template <class T>
void pack_item(const vector<T>& vec)
{
	pack_item(vec.size());
	for (auto& item : vec) {
		pack_item(item);
	}
}

void pack(const vector <unsigned int> &vec)
{
	pack_item(vec);
}

void pack(const vector <unsigned short> &vec)
{
	pack_item(vec);
}

void pack(const vector <int> &vec)
{
	pack_item(vec);
}

void pack(const vector <double> &vec)
{
	pack_item(vec);
}

void pack(const vector< vector <unsigned int> > &vec)
{
	pack_item(vec);
}

void pack(const vector< vector <double> > &vec)
{
	pack_item(vec);
}

void pack(const vector< vector <float> > &vec)
{
	pack_item(vec);
}

void pack(const vector< vector< vector <double> > > &vec)
{
	pack_item(vec);
}

void pack(const vector <string> &vec)
{
	pack_item(vec);
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

void pack_item(const AREA& area)
{
	pack_item(area.code);
	pack_item(area.region);
	pack_item(area.x);
	pack_item(area.y);
	pack_item(area.agepop);
	pack_item(area.pop);
	pack_item(area.covar);
	pack_item(area.ind);
}

void pack(const vector <AREA> &vec)
{
	pack_item(vec);
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
	pack_item(vec);
}

template<class T>
void unpack_item(T &num)
{
	num = buffer[k]; k++;
}
void unpack_item(string &vec)
{
	unsigned int jmax, j;
	
	unpack_item(jmax);
	stringstream ss; for(j = 0; j < jmax; j++){ ss << (char) buffer[k]; k++;}
	vec = ss.str();
}

void unpack(unsigned int &num)
{
	unpack_item(num);
}

void unpack(unsigned short &num)
{
	unpack_item(num);
}

void unpack(double &num)
{
	unpack_item(num);
}

void unpack(string &vec)
{
	unpack_item(vec);
}

template <class T>
void unpack_item(vector<T>& vec)
{
	unsigned int size;
	unpack_item(size);
	vec.resize(size);
	for (auto& item : vec) {
		unpack_item(item);
	}
}

void unpack(vector <unsigned int> &vec)
{
	unpack_item(vec);
}

void unpack(vector <unsigned short> &vec)
{
	unpack_item(vec);
}

void unpack(vector <int> &vec)
{
	unpack_item(vec);
}

void unpack(vector <double> &vec)
{
	unpack_item(vec);
}

void unpack(vector< vector <unsigned int> > &vec)
{
	unpack_item(vec);
}

void unpack(vector< vector <double> > &vec)
{
	unpack_item(vec);
}

void unpack(vector< vector <float> > &vec)
{
	unpack_item(vec);
}

void unpack(vector< vector< vector <double> > > &vec)
{
	unpack_item(vec);
}

void unpack(vector <string> &vec)
{
	unpack_item(vec);
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
	unpack_item(imax); vec.resize(imax); 
	for(i = 0; i < imax; i++){
		unpack(vec[i].code);
		unpack(vec[i].region);
		unpack(vec[i].x);
		unpack(vec[i].y);
		unpack(vec[i].agepop);
		unpack(vec[i].pop);
		unpack(vec[i].covar);  
		unpack(vec[i].ind);  
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
	unpack_item(vec);
}
