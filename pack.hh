#pragma once

using namespace std;

#include "PART.hh"

void packinit();
int packsize();
double * packbuffer();

void pack(short num);
void pack(long num);
void pack(vector <long> &vec);
void pack(vector <double> &vec);
void pack(vector< vector <long> > &vec);
void pack(vector< vector <double> > &vec);
void pack(vector <string> &vec);
void pack(vector< vector <FEV> > &vec, short fedivmin, short fedivmax);
void pack(vector <HOUSE> &vec, long min, long max);

void unpack(short &num);
void unpack(long &num);
void unpack(vector <long> &vec);
void unpack(vector <double> &vec);
void unpack(vector< vector <long> > &vec);
void unpack(vector< vector <double> > &vec);
void unpack(vector <string> &vec);
void unpack(vector< vector <FEV> > &vec, short fedivmin, short fedivmax);
void unpack(vector <HOUSE> &vec, long min, long max);
	
	