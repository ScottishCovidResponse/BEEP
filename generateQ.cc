// This code is used to generate the Q tensor from the census flow data and age mixing matrices.
// This Q tensor is then used as an input into the full code.

// Compile: g++ generateQ.cc -O3
// Run: ./a.out "Name of data directory" 

// E.g. ./a.out Data_ScotlandMSOA   
// E.g. ./a.out Data_example

#include <iostream>
#include <sstream>
#include <math.h>
#include <map>
#include <algorithm>
#include <vector>
#include <fstream>
#include <cmath>

using namespace std;

#include "generateQ.hh"
#include "utils.hh"

#ifdef USE_DATA_PIPELINE
#include "datapipeline.hh"
#include "array.hh"
#endif

unsigned int nage;                         // The number of age categories used 
//const short normon = 1;                    // Determines if matrix normalised
//const short symetric = 1;                  // Set to 1 if Q matrix symetric in area

struct MATRIX {                            // Loads a matrix
	unsigned int N;													 // The size of the matrix
	vector <vector <double> > ele;           // The elements of the matrix
};

struct SPARSEMATRIX {                      // Loads a matrix
	unsigned int N;													 // The size of the matrix
	vector <unsigned int> i;
	vector <unsigned int> j;
	vector <double> val;                     // The elements of the matrix
};

TABLE loadtable(string file, string head);
TABLE loadarray(string file, string dir);
unsigned int findcol(TABLE &tab, string name);
vector <AREA> loadarea(TABLE tab);
MATRIX matfromtable(TABLE tab, unsigned int N);
SPARSEMATRIX loadsparse(string file, string dir, unsigned int N);
SPARSEMATRIX identity(unsigned int N);
void plotmat(MATRIX mat, string title);
void generateQten(SPARSEMATRIX &M, MATRIX &N, string name, GENQ &genQ, vector <AREA> &area);
TABLE loadarrayfromdatapipeline(string file);

string strip(string line);

DataPipeline *datapipeline;             // DataPipeline object

void generateQ(unsigned int nage, string datadir, GENQ &genQ, vector <AREA> &area,
							 DataPipeline *dp)
{
	TABLE tab;
	MATRIX N_all, N_home, N_other, N_school, N_work;
	SPARSEMATRIX M, I;

	datapipeline = dp;

	cout << "Generating Q tensors." << endl;
	
	tab = loadarray(genQ.Nall, datadir);        // Loads age stratified mixing matrices for different activities
	N_all = matfromtable(tab,nage);
	
	tab = loadarray(genQ.Nhome, datadir);
	N_home = matfromtable(tab,nage);
	
	tab = loadarray(genQ.Nother, datadir);
	N_other = matfromtable(tab,nage);
	
	tab = loadarray(genQ.Nschool, datadir);
	N_school = matfromtable(tab,nage);
	
	tab = loadarray(genQ.Nwork, datadir);
	N_work = matfromtable(tab,nage);

	//tab = loadtable(datadir+"/"+genQ.areadata,"head");           // Loads information about the areas
	//area = loadarea(tab);
	
	M = loadsparse(genQ.M,datadir,area.size());
	
	I = identity(area.size());                                 // Generates the identity matrix
	
	//plotmat(N_all,"'All' matrix");                             // Outputs matrices to the terminal
	//plotmat(N_home,"'Home' matrix");
		
	generateQten(M,N_all,genQ.flowall,genQ,area);      
	
	generateQten(I,N_home,genQ.localhome,genQ,area);      
	//generateQten(M,N_all,genQ.localhome,genQ,area);      

	cout << endl;
}

/// Outputs a matrix
void plotmat(MATRIX mat, string title)
{
	unsigned int i, j;
	
	cout << title << ":" << endl;
	for(j = 0; j < mat.N; j++){
		for(i = 0; i < mat.N; i++) cout << mat.ele[j][i] << "\t";
		cout << endl;
	}
	cout << endl;
}

/// Loads age stratified population data for areas
vector <AREA> loadarea(TABLE tab)
{
	unsigned int c, a;
	vector <int> agecol;
	AREA are;
	vector <AREA> area;
	
	switch(nage){
	case 1:
		agecol.push_back(findcol(tab,"age0-14"));
		break;
		
	case 3:		
		agecol.push_back(findcol(tab,"age0-19"));
		agecol.push_back(findcol(tab,"age20-64"));
		agecol.push_back(findcol(tab,"age65+"));
		break;
		
	case 4:	
		agecol.push_back(findcol(tab,"age0-14"));
		agecol.push_back(findcol(tab,"age15-44"));
		agecol.push_back(findcol(tab,"age45-64"));
		agecol.push_back(findcol(tab,"age65+"));
		break;
	
	default: emsg("nage not recognised"); break;	
	}
	
	for(c = 0; c < tab.nrow; c++){
		are.agepop.resize(nage);
		for(a = 0; a < nage; a++) are.agepop[a] = atoi(tab.ele[c][agecol[a]].c_str());
		area.push_back(are);
	}
	
	return area;
}

void generateQten(SPARSEMATRIX &M, MATRIX &N, string name, GENQ &genQ, vector <AREA> &area)
{
	unsigned int nage = N.N, narea = M.N, k, c, cc, a, aa, vi, j, jmax, q;
	double v, sum, sum2;
	vector <float> vec;
	vector <vector <unsigned short> > to;      // Stores the mixing matrix between areas and ages at different times
	vector <vector< vector <float> > > val;
	
	vec.resize(nage);
	to.resize(nage*narea);
	val.resize(nage*narea);

	long num = 0;
	for(k = 0; k < M.val.size(); k++){
		c = M.i[k]; cc = M.j[k]; v = M.val[k];
		for(a = 0; a < nage; a++){
			for(aa = 0; aa < nage; aa++){
				vec[aa] = v*N.ele[a][aa];
				if(std::isnan(vec[aa])) emsg("Raw value in '"+name+"' not a number");
			}

			to[c*nage+a].push_back(cc); num++;
			val[c*nage+a].push_back(vec);
					
			if(c < cc){
				to[cc*nage+a].push_back(c);
				val[cc*nage+a].push_back(vec);
			}
		}
	}

	for(c = 0; c < narea; c++){                       // Normalises the tensor
		sum = 0; sum2 = 0;
		for(a = 0; a < nage; a++){
			sum2 += area[c].agepop[a];
			vi = c*nage + a;
		
			jmax = to[vi].size();
			for(j = 0; j < jmax; j++){
				cc = to[vi][j];
				for(aa = 0; aa < nage; aa++) sum += area[c].agepop[a]*val[vi][j][aa]*area[cc].agepop[aa];
			}
		}
		sum /= sum2;
		
		for(a = 0; a < nage; a++){
			vi = c*nage + a;
			jmax = to[vi].size();
			for(j = 0; j < jmax; j++){
				for(aa = 0; aa < nage; aa++){
					val[vi][j][aa] /= sum;
					if(std::isnan(val[vi][j][aa])) emsg("Value in '"+name+"' not a number");
				}
			}
		}
	}
	
	q = genQ.Qten.size();
	genQ.Qten.push_back(SPARSETENSOR ());

	genQ.Qten[q].name = name;
	
	genQ.Qten[q].tof = new unsigned short*[to.size()];
	genQ.Qten[q].ntof = new unsigned short[to.size()];
	for(vi = 0; vi < to.size(); vi++){
		genQ.Qten[q].ntof[vi] = to[vi].size();
		genQ.Qten[q].tof[vi] = new unsigned short[to[vi].size()];
		for(j = 0; j < to[vi].size(); j++){
			genQ.Qten[q].tof[vi][j] = to[vi][j];
		}
	}		

	genQ.Qten[q].valf = new float**[val.size()];
	for(vi = 0; vi < val.size(); vi++){
		genQ.Qten[q].valf[vi] = new float*[val[vi].size()];
		for(j = 0; j < val[vi].size(); j++){
			genQ.Qten[q].valf[vi][j] = new float[val[vi][j].size()];
			for(a = 0; a < nage; a++) genQ.Qten[q].valf[vi][j][a] = val[vi][j][a];
		}
	}
}

/// Generates the identity matrix
SPARSEMATRIX identity(unsigned int N)
{
	unsigned int k;
	SPARSEMATRIX M;
	
	M.N = N;
	M.i.resize(N); M.j.resize(N); M.val.resize(N);
	for(k = 0; k < N; k++){
		M.i[k] = k; M.j[k] = k; M.val[k] = 1; 
	}
	
	return M;
}

/// Creates a matrix from a table
MATRIX matfromtable(TABLE tab, unsigned int N)
{
	unsigned int i, j, ii, jj;
	double sum, val;
	MATRIX mat;
	vector <unsigned int> div;
	
	switch(N){
	case 1:
		div.push_back(0); div.push_back(16);
		break;
		
	case 3:
	  div.push_back(0); div.push_back(4); div.push_back(13); div.push_back(16);
		break;
		
	case 4:
	  div.push_back(0); div.push_back(3); div.push_back(9); div.push_back(13); div.push_back(16);
		break;
		
	default: emsg("nage not recognised"); break;
	}

	if(tab.nrow != 16 || tab.ncol != 16) emsg("Table size not right");
	
	mat.N = N;
	mat.ele.resize(N);
	for(j = 0; j < N; j++){
		mat.ele[j].resize(N);
		for(i = 0; i < N; i++){
			sum = 0;
			for(jj = div[j]; jj < div[j+1]; jj++){
				for(ii = div[i]; ii < div[i+1]; ii++){
					val = atof(tab.ele[jj][ii].c_str());
					if(std::isnan(val)) emsg("Not a number!");
					sum += val;
				}
			}
			mat.ele[j][i] = sum/((div[j+1]-div[j])*(div[i+1]-div[i]));
		}		
	}
	
	return mat;
}

/// Uses 'from' and 'to' columns to generate a sparse matrix
SPARSEMATRIX loadsparsefromdatapipeline(string file, unsigned int N)
{
	SPARSEMATRIX mat;
#ifdef USE_DATA_PIPELINE

	unsigned int a1, a2;
	double v;
	string line;
	
	mat.N = N;	

	Table  dptable = datapipeline->read_table(file, "default");

	int nrows = dptable.get_column_size();
	
	vector<long> area1 = dptable.get_column<long>("area1");
	vector<long> area2 = dptable.get_column<long>("area2");
	vector<double> contact = dptable.get_column<double>("contact");

	for (int i = 0; i < nrows; i++) {
		a1 = area1.at(i);
		a2 = area2.at(i);
		v = contact.at(i);
		mat.i.push_back(a1);
		mat.j.push_back(a2);
		if(std::isnan(v)) emsg("Value in file '"+file+"' is not a number");
		mat.val.push_back(v);
	}

	// Array<double> dparray = datapipeline->read_array(file,"default");
	// vector<int>   dims = dparray.size();

	// for (int i = 0; i < dims.at(0); i++) {
	// 	a1 = dparray(i,0);
	// 	a2 = dparray(i,1);
	// 	v = dparray(i,2);
	// 	mat.i.push_back(a1);
	// 	mat.j.push_back(a2);
	// 	if(std::isnan(v)) emsg("Value in file '"+file+"' is not a number");
	// 	mat.val.push_back(v);

	// }

	cout << "Loaded sparse matrix " << file << " from data pipeline" << endl;

#else
	N = N;
	emsg("loadsparsefromdatapipeline for "+file+" cannot be called as data pipeline is not compiled in");
#endif

	return mat;
}


/// Uses 'from' and 'to' columns to generate a sparse matrix
SPARSEMATRIX loadsparsefromfile(string file, unsigned int N)
{
	unsigned int a1, a2;
	double v;
	string line;
	SPARSEMATRIX mat;
	
	ifstream in(file.c_str());                             // Loads information about areas
	if(!in) emsg("Cannot open the file '"+file+"'");
	
	mat.N = N;	
	getline(in,line);
	do{
		getline(in,line);
		if(in.eof()) break;
				
		stringstream ss(line);
		ss >> a1 >> a2 >> v;
		mat.i.push_back(a1);
		mat.j.push_back(a2);
		if(std::isnan(v)) emsg("Value in file '"+file+"' is not a number");
		mat.val.push_back(v);
	}while(1 == 1);

	cout << "Loaded sparse matrix " << file << " from file" << endl;

	return mat;
}

/// Uses 'from' and 'to' columns to generate a sparse matrix
SPARSEMATRIX loadsparse(string file, string dir, unsigned int N)
{
	if (stringhasending(file, ".txt")) {
		return loadsparsefromfile(dir+"/"+file, N);
	} else {
		return loadsparsefromdatapipeline(file, N);
	}
}

/// Loads a table from the data pipeline
TABLE loadarrayfromdatapipeline(string file)
{
	TABLE tab;

#ifdef USE_DATA_PIPELINE

	Array<double> dparray = datapipeline->read_array(file,"default");
	vector<int>   dims = dparray.size();


	tab.file = file;
//	tab.heading = ????;
	tab.ncol = dims.at(1);
	
	for (size_t i = 0; i < dims.at(0); i++) {
		vector<string> vec;

		for (size_t j = 0; j < dims.at(1); j++) {
			vec.push_back(to_string(dparray(i,j)));
		}
		tab.ele.push_back(vec);
	}

	tab.nrow = tab.ele.size();

	cout << "Loaded array " << file << " from data pipeline" << endl;
#else
	emsg("loadarrayfromdatapipeline for "+file+" cannot be called as data pipeline is not compiled in");
#endif

	return tab;
}


/// Loads a table from a file
TABLE loadarray(string file, string dir)
{
	if (stringhasending(file, ".txt")) {
		return loadtable(dir+"/"+file, "nohead");
	} else {
		return loadarrayfromdatapipeline(file);
	}
}



/// Loads a table from a file
TABLE loadtable(string file, string head)
{
	TABLE tab;
	string line, st;
	vector <string> vec;
	
	ifstream in(file.c_str());                             // Loads information about areas
	if(!in) emsg("Cannot open the file '"+file+"'");
		
	if(head == "head"){	
		getline(in,line);

		stringstream ss(line);
		do{
			getline(ss,st,'\t'); st = strip(st);
			tab.heading.push_back(st);
			if(ss.eof()) break;
		}while(1 == 1);
		tab.ncol = tab.heading.size();
	}
	
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
		
		if(tab.ele.size() == 0 && head == "nohead") tab.ncol = vec.size();
		if(vec.size() != tab.ncol) emsg("Rows in file '"+file+"' do not all have the same length.");
		
		tab.ele.push_back(vec);
	}while(1 == 1);
	tab.nrow = tab.ele.size();
	
	cout << "Loaded array " << file << " from file " << file << endl;

	return tab;
}

/// Finds a column in a table
unsigned int findcol(TABLE &tab, string name)
{
	unsigned int c;
	
	for(c = 0; c < tab.ncol; c++) if(tab.heading[c] == name) break;
	if(c == tab.ncol) emsg("Cannot find the column heading '"+name+"'.");
	return c;
}		

/// Strips off '\\r' character if necessary
string strip(string line)
{
	unsigned int len = line.length();
	if(len > 0){
		if(line.substr(len-1,1) == "\r") line = line.substr(0,len-1);
	}
	return line;
}	
