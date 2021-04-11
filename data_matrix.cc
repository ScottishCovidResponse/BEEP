#include <iostream>
#include <vector>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <math.h> 
#include <iomanip>

using namespace std;

#include "data.hh"	

/// Initialises all the quantities in genQ
void Data::generate_matrices()
{
	for(auto i = 0u; i < genQ.N_name.size(); i++){
		Matrix N;
		if(genQ.N_name[i] == "UNIT MATRIX"){
			N.N = 1; N.ele.resize(1); N.ele[0].resize(1); N.ele[0][0] = 1;
		}
		else{
			Table tab = load_table(genQ.N_name[i],data_directory,false);      
		
			Matrix N = matfromtable(tab,nage);
		}
		agematrix_normalise(N);
		genQ.N.push_back(N);
	}

	if(genQ.M_name == ""){
		if(area.size() != 1) emsg("'geo_mixing_matrix' needs to be specified to specify geograpical mixing.");
		SparseMatrix M;
		M.N = 1;
		M.diag.resize(1);
		M.diag[0] = 1;
		M.to.resize(1);
		M.val.resize(1);
	
		genQ.M = M;
	}
	else{
		genQ.M = loadsparse(data_directory+"/"+genQ.M_name,area.size(),genQ);
	}
	geo_normalise(genQ.M);
}


/// Normalises the geographical matrix
void Data::geo_normalise(SparseMatrix &mat)
{
	if(mat.N != area.size()) emsgEC("data_matrix",10);

	vector <unsigned int> pop(area.size());
	for(auto c = 0u; c < area.size(); c++){    
		auto sum = 0u; for(auto dp = 0u; dp < area[c].pop.size(); dp++) sum += area[c].pop[dp];
		pop[c] = sum;
	}
	
	genQ.onlydiag.resize(area.size());         
	for(auto c = 0u; c < area.size(); c++) genQ.onlydiag[c] = 1.0/pop[c];

	for(auto c = 0u; c < area.size(); c++){
		if(mat.diag[c] == 0) emsg("No interaction");
		
		auto sum = mat.diag[c]*pop[c];
		
		auto jmax = mat.to[c].size();
		for(auto j = 0u; j < jmax; j++) sum += mat.val[c][j]*pop[mat.to[c][j]];

		mat.diag[c] /= sum;
		for(auto j = 0u; j < jmax; j++) mat.val[c][j] /= sum;
	}
	
	genQ.factor.resize(area.size());
	for(auto c = 0u; c < area.size(); c++){
		genQ.factor[c] = genQ.onlydiag[c]/genQ.M.diag[c];
	}
	
	if(false){
		for(auto c = 0u; c < area.size(); c++){
			cout << pop[c] << " " << area[c].code << " " << c << " " <<	genQ.onlydiag[c] << " " << genQ.M.diag[c] << " " << 	genQ.M.diag[c] / genQ.onlydiag[c] <<  " area\n";
		}
	}
	
	if(false){
		for(auto c = 0u; c < area.size(); c++){
			cout << c << " " << mat.diag[c] << " : ";
			auto jmax = mat.to[c].size();
			for(auto j = 0u; j < jmax; j++) cout << mat.to[c][j] << "," << mat.val[c][j] << "  ";
			cout << "value\n";		
		}
	}
}
	

/// Normalises an age mixing matrix
void Data::agematrix_normalise(Matrix &mat)
{
	auto nage = mat.N;
			
	vector <unsigned int> pop(nage);
	for(auto a = 0u; a < nage; a++) pop[a] = 0;
	
	for(auto c = 0u; c < area.size(); c++){   
		for(auto dp = 0u; dp < area[c].pop.size(); dp++){
			pop[democatpos[dp][0]] += area[c].pop[dp];
		}
	}

	auto sum = 0.0, sum2 = 0.0;
	for(auto a = 0u; a < nage; a++){
		for(auto aa = 0u; aa < nage; aa++){
			sum += pop[aa]*mat.ele[aa][a];
			sum2 += pop[aa];
		}
	}
	sum /= sum2;
	
	for(auto a = 0u; a < nage; a++){         
		for(auto aa = 0u; aa < nage; aa++) mat.ele[aa][a] /= sum;
	}
	
	if(false){
		for(auto a = 0u; a < nage; a++){    
			for(auto aa = 0u; aa < nage; aa++) cout <<  mat.ele[a][aa] << ", "; cout << "m\n";
		}
	}
}


/// Creates a matrix from a table
Matrix Data::matfromtable(const Table& tab, unsigned int N)
{		
	vector <unsigned int> pop(N);
	for(auto a = 0u; a < N; a++) pop[a] = 0;
	
	auto poptot = 0.0;
	for(auto c = 0u; c < area.size(); c++){   
		for(auto dp = 0u; dp < area[c].pop.size(); dp++){
			pop[democatpos[dp][0]] += area[c].pop[dp];
			poptot += area[c].pop[dp];
		}
	}

	Matrix mat;
	mat.N = N;
	
	if(tab.nrow != N || tab.ncol != N) emsg("For the file '"+tab.file+"' the table size is not right");

	mat.ele.resize(N);
	for(auto j = 0u; j < N; j++){
		mat.ele[j].resize(N);
		for(auto i = 0u; i < N; i++){
			mat.ele[j][i] = atof(tab.ele[j][i].c_str());
		}
	}
				 
	for(auto j = 0u; j < N; j++){
		for(auto i = 0u; i < N; i++){
			mat.ele[j][i] /= double(pop[j])/poptot;
		}
	}
	
	if(false){       // Equal mixing of age groups
		for(auto j = 0u; j < N; j++){
			for(auto i = 0u; i < N; i++){
				mat.ele[j][i] = 1.0;
			}
		}
	}
		
	auto mat_st = mat;  // Makes the matrix symetric
	for(auto j = 0u; j < N; j++){
		for(auto i = 0u; i < N; i++){
			mat.ele[j][i] = 0.5*(mat_st.ele[j][i]+mat_st.ele[i][j]);
		}
	}
	
	if(false){
		for(auto j = 0u; j < N; j++){
			for(auto i = 0u; i < N; i++){
				cout << int(mat.ele[j][i]*10) << ", ";
			}
			cout << "\n";
		}
		emsg("Age Matrix");
	}
	
	return mat;
}


/// Uses 'from' and 'to' columns to generate a sparse matrix
SparseMatrix Data::loadsparse(const string& file, unsigned int N, GenerateQ &genQ)
{
	if(N != narea) emsg("N area problem");
	
	ifstream in(file.c_str());                             // Loads information about areas
	if(!in) emsg("Cannot open the file '"+file+"'");
	
	vector <unsigned int> mati;
	vector <unsigned int> matj;
	vector <double> val;  
	
	string line;
	getline(in,line);
	do{
		getline(in,line);
		if(in.eof()) break;
				
		stringstream ss(line);
		unsigned int a1, a2;
		double v;
		ss >> a1 >> a2 >> v;
	
		a1 = genQ.area_filter_ref[a1];
		a2 = genQ.area_filter_ref[a2];
	
		if(a1 != UNSET && a2 != UNSET){
			if(a1 > a2) emsg("Matrix is not triangular");
			mati.push_back(a1);
			matj.push_back(a2);
			if(std::isnan(v)) emsg("The value '"+to_string(v)+"' in file '"+file+"' is not a number");
			val.push_back(v);
		}
	}while(true);

	SparseMatrix mat;
	mat.N = N;	
	mat.diag.resize(N); for(auto c = 0u; c < N; c++) mat.diag[c] = 0; 
	mat.to.resize(N);
	mat.val.resize(N);
	
	for(auto k = 0u; k < val.size(); k++){
		if(mati[k] == matj[k]) mat.diag[mati[k]] = val[k];
	}
	
	for(auto k = 0u; k < val.size(); k++){
		if(mati[k] != matj[k]){
			mat.to[mati[k]].push_back(matj[k]);
			mat.val[mati[k]].push_back(val[k]);
			
			mat.to[matj[k]].push_back(mati[k]);
			mat.val[matj[k]].push_back(val[k]);
		}
	}
	
	cout << "Loaded sparse matrix " << file << "." << endl;

	return mat;
}


/// Outputs a matrix
void Data::plotmat(const Matrix& mat, const string& title)
{
	cout << title << ":" << endl;
	for(auto j = 0u; j < mat.N; j++){
		for(auto i = 0u; i < mat.N; i++) cout << mat.ele[j][i] << "\t";
		cout << endl;
	}
	cout << endl;
}
