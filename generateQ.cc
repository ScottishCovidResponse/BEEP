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

using namespace std;

const short nage = 3;                      // The number of age categories used 
const short normon = 0;                    // Determines if matrix normalised
const short symetric = 1;                  // Set to 1 if Q matrix symetric in area
	
struct AREA {                              // Stores information about areas
	vector <unsigned int> agepop;            // The populations in different age groups
};

struct TABLE {                             // Loads a table
	unsigned int ncol;                       // The number of columns
	unsigned int nrow;                       // The number if rows
	vector <string> heading;                 // The headings for the columns
	vector <vector <string> > ele;        	 // The elements of the table
};

struct MATRIX {                            // Loads a matrix
	unsigned int N;													 // The size of the matrix
	vector <vector <double> > ele;           // The elements of the matrix
};

struct TENSOR {                            // Loads a tensor
	vector <vector <vector <vector <double> > > > ele;  // The elements of the tensor
};

TABLE loadtable(string file, string head);
unsigned int findcol(TABLE &tab, string name);
vector <AREA> loadarea(TABLE tab);
MATRIX matfromtable(TABLE tab, unsigned int N);
MATRIX matfromtable_sparse(TABLE tab, unsigned int N);
TENSOR generateQ(MATRIX M, MATRIX N);
MATRIX identity(unsigned int N);
void normalise(TENSOR &Q, vector <AREA> &area);
void outputQ(string file, TENSOR Q);
void plotmat(MATRIX mat, string title);

string strip(string line);
void emsg(string msg);

int main(int argc, char** argv)
{
	string datadir;
	vector <AREA> area;
	TABLE tab;
	MATRIX N_all, N_home, N_other, N_school, N_work, M, I;
	TENSOR Q;
	
	if(argc != 2) emsg("Incorrect number of arguments");
	
	datadir = argv[1];                                         // The data directory
	
	tab = loadtable(datadir+"/Ndata_all.txt","nohead");        // Loads age stratified mixing matrices for different activities
	N_all = matfromtable(tab,nage);
	
	tab = loadtable(datadir+"/Ndata_home.txt","nohead");
	N_home = matfromtable(tab,nage);
	
	tab = loadtable(datadir+"/Ndata_other.txt","nohead");
	N_other = matfromtable(tab,nage);
	
	tab = loadtable(datadir+"/Ndata_school.txt","nohead");
	N_school = matfromtable(tab,nage);
	
	tab = loadtable(datadir+"/Ndata_work.txt","nohead");
	N_work = matfromtable(tab,nage);

	tab = loadtable(datadir+"/areadata.txt","head");           // Loads information about the areas
	area = loadarea(tab);
	
	tab = loadtable(datadir+"/Mdata.txt","head");              // Loads the census flow data
	M = matfromtable_sparse(tab,area.size());

	I = identity(area.size());                                 // Generates the identity matrix
	
	plotmat(N_all,"'All' matrix");                             // Outputs matrices to the terminal
	plotmat(N_home,"'Home' matrix");
		
	Q = generateQ(M,N_all);                                    // Outputs a Q matrix representative of 'normal' mixing
	if(normon == 1) normalise(Q,area);
	outputQ(datadir+"/Q_flow_all_data.txt",Q);

	Q = generateQ(I,N_home);                                   // Outputs a Q matrix representative of someone at home 
	if(normon == 1) normalise(Q,area);
	outputQ(datadir+"/Q_local_home_data.txt",Q);
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

/// Outputs the non-zero elements of Q
void outputQ(string file, TENSOR Q)
{
	unsigned int c, cc, ccmin, a, aa, narea;
	
	ofstream op(file.c_str());
	if(!op) emsg("Cannot write file");
	
	cout << "Outputing file '" << file << "'." << endl;
	
	op << "area A\tarea B";
	for(a = 0; a < nage; a++){
		for(aa = 0; aa < nage; aa++){
			op << "\tAG" << a << "/" << "AG" << aa; 
		}
	}
	op << endl;

	narea = Q.ele.size();
	for(c = 0; c < narea; c++){
		if(symetric == 1) ccmin = c; else ccmin = 0;
		for(cc = ccmin; cc < narea; cc++){
			for(a = 0; a < nage; a++){
				for(aa = 0; aa < nage; aa++) if(Q.ele[c][a][cc][aa] != 0) break;
				if(aa != nage) break;
			}
			
			if(a != nage){
				op << c << "\t" << cc;
				for(a = 0; a < nage; a++){
					for(aa = 0; aa < nage; aa++){
						op << "\t" << Q.ele[c][a][cc][aa];
					}
				}
				op << endl;
			}
		}
	}
}

/// Loads age stratified population data for areas
vector <AREA> loadarea(TABLE tab)
{
	unsigned int c, a;
	vector <int> agecol;
	AREA are;
	vector <AREA> area;
		
	agecol.push_back(findcol(tab,"0-18"));
	agecol.push_back(findcol(tab,"19-64"));
	agecol.push_back(findcol(tab,"65+"));

	for(c = 0; c < tab.nrow; c++){
		are.agepop.resize(nage);
		for(a = 0; a < nage; a++) are.agepop[a] = atoi(tab.ele[c][agecol[a]].c_str());
		area.push_back(are);
	}
	
	return area;
}

/// Generates a Q tensor from a movement matrix and an age contact matrix
TENSOR generateQ(MATRIX M, MATRIX N)
{
	unsigned int c, a, cc, aa;
	TENSOR Q;
	
	Q.ele.resize(M.N);
	for(c = 0; c < M.N; c++){
		Q.ele[c].resize(N.N);
		for(a = 0; a < N.N; a++){
			Q.ele[c][a].resize(M.N);
			for(cc = 0; cc < M.N; cc++){
				Q.ele[c][a][cc].resize(N.N);
				for(aa = 0; aa < N.N; aa++){
					Q.ele[c][a][cc][aa] = M.ele[c][cc]*N.ele[a][aa];
				}
			}
		}
	}	
	
	return Q;
}

/// Generates the identity matrix
MATRIX identity(unsigned int N)
{
	unsigned int i, j;
	MATRIX M;
	
	M.N = N;
	M.ele.resize(N);
	for(j = 0; j < N; j++){
		M.ele[j].resize(N);
		for(i = 0; i < N; i++){
			if(i == j) M.ele[j][i] = 1; else M.ele[j][i] = 0;
		}
	}	
	
	return M;
}
	
/// Normalises the matrix by area
void normalise(TENSOR &Q, vector <AREA> &area)
{
	unsigned int c, a, cc, aa, narea;
	double sum, sum2;
	
	narea = area.size();
	
	for(c = 0; c < narea; c++){ 
		sum = 0; sum2 = 0;
		for(a = 0; a < nage; a++){
			sum2 += area[c].agepop[a];
			for(cc = 0; cc < narea; cc++){  
				for(aa = 0; aa < nage; aa++){ 
					sum += area[c].agepop[a]*Q.ele[c][a][cc][aa]*area[cc].agepop[aa];
				}
			}
		}
		sum /= sum2;

		for(a = 0; a < nage; a++){
			for(cc = 0; cc < narea; cc++){  
				for(aa = 0; aa < nage; aa++){ 
					Q.ele[c][a][cc][aa] /= sum;
				}
			}
		}
	}
}

/// Creates a matrix from a table
MATRIX matfromtable(TABLE tab, unsigned int N)
{
	unsigned int i, j, ii, jj;
	double sum, val;
	MATRIX mat;
	vector <unsigned int> div;
	
	if(N == 3){
	  div.push_back(0); div.push_back(4); div.push_back(13); div.push_back(16);
	}
	else emsg("N not recognised");

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
					if(isnan(val)) emsg("Not a number!");
					sum += val;
				}
			}
			mat.ele[j][i] = sum/((div[j+1]-div[j])*(div[i+1]-div[i]));
		}		
	}
	
	return mat;
}

/// Uses 'from' and 'to' columns to generate a matrix
MATRIX matfromtable_sparse(TABLE tab, unsigned int N)
{
	unsigned int i, j, r, a1, a2;
	double v;
	MATRIX mat;
	
	mat.N = N;
	mat.ele.resize(N);
	for(j = 0; j < N; j++){
		mat.ele[j].resize(N);
		for(i = 0; i < N; i++){
			mat.ele[j][i] = 0;
		}
	}
	
	for(r = 0; r < tab.nrow; r++){
		a1 = atoi(tab.ele[r][0].c_str()); a2 = atoi(tab.ele[r][1].c_str()); v = atof(tab.ele[r][2].c_str());
		mat.ele[a1][a2] = v; mat.ele[a2][a1] = v;
	}
	
	return mat;
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

/// Strips off '\r' character if necessary
string strip(string line)
{
	unsigned int len = line.length();
	if(len > 0){
		if(line.substr(len-1,1) == "\r") line = line.substr(0,len-1);
	}
	return line;
}	

/// Displays an error message
void emsg(string msg)
{
	cout << msg << endl;
	exit (EXIT_FAILURE);
}
