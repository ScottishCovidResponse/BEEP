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
			N = age_mixing_matrix(tab);
		}
		agematrix_normalise(N);
		genQ.N.push_back(N);
	}

	if(genQ.M_name == ""){
		if(area.size() != 1){
			emsg("A value for 'geo_mixing_matrix' must be set (to specify the geograpical mixing of individuals).");
		}

		SparseMatrix M;
		M.N = 1;
		M.diag.resize(1);
		M.diag[0] = 1;
		M.to.resize(1);
		M.val.resize(1);
	
		genQ.M = M;
	}
	else{
		Table tab = load_table(genQ.M_name,data_directory,false);   
		genQ.M = load_geo_mixing_matrix(tab);
	}

	geo_normalise(genQ.M);

	genQ.treenode = area_split(genQ.M);
}


/// Normalises the geographical matrix (this ensures that individuals in each area 
void Data::geo_normalise(SparseMatrix &mat)
{
	if(mat.N != area.size()) emsgEC("data_matrix",1);

	vector <unsigned int> pop(area.size());
	for(auto c = 0u; c < area.size(); c++) pop[c] = area[c].total_pop;
	
	genQ.onlydiag.resize(area.size());         
	for(auto c = 0u; c < area.size(); c++) genQ.onlydiag[c] = 1.0/pop[c];

	for(auto c = 0u; c < area.size(); c++){
		if(mat.diag[c] == 0) emsg("The matrix specified in 'geo_mixing_matrix' cannot have zero on the diagonal");
		
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
			cout << pop[c] << " " << area[c].code << " " << c << " " <<	genQ.onlydiag[c] << " " << genQ.M.diag[c] << " " << 	genQ.M.diag[c] / genQ.onlydiag[c] <<  " area" << endl;
		}
	}
	
	if(false){
		for(auto c = 0u; c < area.size(); c++){
			cout << c << " " << mat.diag[c] << " : ";
			auto jmax = mat.to[c].size();
			for(auto j = 0u; j < jmax; j++) cout << mat.to[c][j] << "," << mat.val[c][j] << "  ";
			cout << "value" << endl;		
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
		for(auto co = 0u; co < area[c].pop_init.size(); co++){
			for(auto dp = 0u; dp < ndemocatpos; dp++) pop[democatpos[dp][0]] += area[c].pop_init[co][dp];
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
	
	mat.norm_factor = 1.0/sum;
	
	for(auto a = 0u; a < nage; a++){         
		for(auto aa = 0u; aa < nage; aa++) mat.ele[aa][a] /= sum;
	}
	
	if(false){
		for(auto a = 0u; a < nage; a++){    
			for(auto aa = 0u; aa < nage; aa++) cout <<  mat.ele[a][aa] << ", "; cout << "m" << endl;
		}
	}
}


/// Creates a matrix from a table
Matrix Data::age_mixing_matrix(const Table &tab) const
{	
	auto N = nage;
	if(tab.nrow != N || tab.ncol != N) emsg("In the file '"+tab.file+"' the table size is not right (it should be a "+to_string(N)+"x"+to_string(N)+" matrix because there are "+to_string(N)+" age groups.");

	vector <unsigned int> pop(N);
	for(auto a = 0u; a < N; a++) pop[a] = 0;
	
	auto poptot = 0.0;
	for(auto c = 0u; c < area.size(); c++){
		for(auto co = 0u; co < area[c].pop_init.size(); co++){
			for(auto dp = 0u; dp < ndemocatpos; dp++){
				auto po = area[c].pop_init[co][dp];
				pop[democatpos[dp][0]] += po;
				poptot += po;
			}
		}
	}

	Matrix mat;
	mat.N = N;
	
	mat.ele.resize(N);
	for(auto j = 0u; j < N; j++){
		mat.ele[j].resize(N);
		for(auto i = 0u; i < N; i++){
			mat.ele[j][i] = get_double(tab.ele[j][i],"In the file '"+tab.file+"'");
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
	
	if(true){
		auto mat_st = mat;  // Makes the matrix symetric
		for(auto j = 0u; j < N; j++){
			for(auto i = 0u; i < N; i++){
				mat.ele[j][i] = 0.5*(mat_st.ele[j][i]+mat_st.ele[i][j]);
			}
		}
	}
	
	if(false){
		for(auto j = 0u; j < N; j++){
			for(auto i = 0u; i < N; i++){
				cout << int(mat.ele[j][i]*10) << ", ";
			}
			cout  << endl;
		}
		emsg("Done");
	}
	
	return mat;
}


/// Loads a square matrix to give the geograpical mixing
SparseMatrix Data::load_geo_mixing_matrix(const Table &tab) const
{
	auto N = narea;
	if(tab.nrow != N || tab.ncol != N) emsg("In the file '"+tab.file+"' the table size is not right (it should be a "+to_string(N)+"x"+to_string(N)+" matrix because there are "+to_string(N)+" areas.");

	SparseMatrix M;
	M.N = N;
	M.to.resize(N);
	M.val.resize(N);
	for(auto j = 0u; j < N; j++){
		auto v = get_double(tab.ele[j][j],"In the file '"+tab.file+"'");		
		M.diag.push_back(v);
		for(auto i = 0u; i < N; i++){
			auto val = get_double(tab.ele[j][i],"In the file '"+tab.file+"'");
			if(i != j && val != 0){
				M.to[j].push_back(i);
				M.val[j].push_back(val);
			}
		}
	}
	
	return M;
}
	
	
/// Uses 'from' and 'to' columns to generate a sparse matrix (not currently implemented)
SparseMatrix Data::load_geo_mixing_matrix_sparse() const
{
	auto N = narea;
	
	SparseMatrix M;
	M.N = N;

	string file = data_directory+"/"+genQ.M_name;
	ifstream in(file.c_str());   
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
	
		if(a1 != UNSET && a2 != UNSET){
			if(a1 > a2) emsg("In the file '"+file+"' the matrix is not triangular");
			mati.push_back(a1);
			matj.push_back(a2);
			if(std::isnan(v)) emsg("The value '"+to_string(v)+"' in file '"+file+"' is not a number");
			val.push_back(v);
		}
	}while(true);

	for(auto k = 0u; k < val.size(); k++){
		if(mati[k] == matj[k]) M.diag[mati[k]] = val[k];
	}
	
	for(auto k = 0u; k < val.size(); k++){
		if(mati[k] != matj[k]){
			M.to[mati[k]].push_back(matj[k]);
			M.val[mati[k]].push_back(val[k]);
			
			M.to[matj[k]].push_back(mati[k]);
			M.val[matj[k]].push_back(val[k]);
		}
	}
	
	cout << "Loaded sparse matrix " << file << "." << endl;

	return M;
}


/// Outputs a matrix
void Data::plotmat(const Matrix &mat, const string &title)
{
	cout << title << ":" << endl;
	for(auto j = 0u; j < mat.N; j++){
		for(auto i = 0u; i < mat.N; i++) cout << mat.ele[j][i] << "\t";
		cout << endl;
	}
	cout << endl;
}


/// Splits area into small pieces (this is used for FixedTree MBPs)
vector <TreeNode> Data::area_split(SparseMatrix M) const
{
	if(M.N != narea) emsgEC("Data",2);
	
	vector <TreeNode> treenode;   
		
	vector <unsigned int> area_pop(narea);
	for(auto c = 0u; c < narea; c++) area_pop[c] = area[c].total_pop;
	
	vector< vector <double> > G;
	G.resize(narea);
	for(auto j = 0u; j < narea; j++){
		G[j].resize(narea);
		for(auto i = 0u; i < narea; i++) G[j][i] = 0;
			
		G[j][j] = M.diag[j]*area_pop[j];
		for(auto k = 0u; k < M.to[j].size(); k++){
			auto i = M.to[j][k]; if(i == j) emsgEC("Data",3);
		
			G[j][i] = M.val[j][k]*sqrt(area_pop[j]*area_pop[i]);
		}
	}	

	if(false){                                                        // Used for testing if orders correctly
		vector <double> posx(narea), posy(narea);
		for(auto j = 0u; j < narea; j++){ posx[j] = j%8; posy[j] = j/8;}
		
		for(auto j = 0u; j < narea; j++){
			for(auto i = 0u; i < narea; i++){
				G[j][i] = exp(-sqrt((posx[i]-posx[j])*(posx[i]-posx[j]) + (posy[i]-posy[j])*(posy[i]-posy[j])));
			}		
		}
	}
	
	if(false){
		for(auto j = 0u; j < narea; j++){
			for(auto i = 0u; i < narea; i++) cout << int(100*G[j][i])<< " ";
			cout << "G" << endl;
		}
		emsg("Done");
	}
	
	TreeNode no;
	for(auto i = 0u; i < narea; i++) no.arearef.push_back(i);      
	treenode.push_back(no);
	
	for(auto n = 0u; n <  treenode.size(); n++){
		if(treenode[n].arearef.size() > 1){
			TreeNode child1, child2;
			split_in_two(treenode[n].arearef,child1.arearef,child2.arearef,G);
			treenode[n].child.push_back(treenode.size());
			treenode.push_back(child1);
						
			treenode[n].child.push_back(treenode.size());
			treenode.push_back(child2);
		}
	}
	
	if(false){
		for(auto n = 0u; n < treenode.size(); n++){
			cout << "Node: " << n << endl;
			cout << "Area: "; for(auto i : treenode[n].arearef) cout << area[i].code << " "; cout << endl;
			cout << "Child: "; for(auto i : treenode[n].child) cout << i << " "; cout << endl;
			cout << endl;
		}
		emsg("Done");
	}
	
	return treenode;
}


/// Splits a vector of areas into two groups based on geographical matrix G
void Data::split_in_two(const vector <unsigned int> &arearef, vector <unsigned int> &ch1,  vector <unsigned int> &ch2, const vector <vector <double> > &G) const
{
	auto si = arearef.size();
	vector <bool> gr(si);                                                 // Starts with a random allocation into groups
	for(auto i = 0u; i < si; i++){ if(i < si/2) gr[i] = true; else gr[i] = false;}
	
	auto fit = get_split_fit(arearef,gr,G);

	auto fail = 0u;
	do{			
		unsigned int i, j;                                                  // Randomly swaps two allocations
		do{
			i = (unsigned int)(ran()*si); j = (unsigned int)(ran()*si); 
		}while(i == j || gr[i] == gr[j]);

		auto fit_new = fit;
		fit_new += get_split_fit_dif(arearef,gr,G,i);	
		fit_new += get_split_fit_dif(arearef,gr,G,j);
		
		if(fit_new > fit){ fail = 0; fit = fit_new;}
		else{	fail++; bool temp = gr[i]; gr[i] = gr[j]; gr[j] = temp;}
	}while(fail < 50);
	
	auto d = fit - get_split_fit(arearef,gr,G);
	if(d*d > TINY) emsg("prob");
		
	for(auto i = 0u; i < si; i++){
		if(gr[i] == true) ch1.push_back(arearef[i]);  
		else ch2.push_back(arearef[i]);  
	}
}

			
/// Returns a value which is higher if the groups are better divided 
double Data::get_split_fit(const vector <unsigned int> &area_list, const vector <bool> &gr, const vector <vector <double> > &G) const
{
	auto N = area_list.size();
	auto sum = 0.0;
	for(auto i = 0u; i < N; i++){
		for(auto j = 0u; j < N; j++){
			auto val = G[area_list[j]][area_list[i]];
			if(val != 0 && i != j){
				if(gr[i] == gr[j]) sum += val; else sum -= val;
			}
		}
	}
	return sum;
}

/// Calculates the change by swithing one gr 
double Data::get_split_fit_dif(const vector <unsigned int> &area_list, vector <bool> &gr, const vector <vector <double> > &G, unsigned int isel) const
{
	auto N = area_list.size();
	auto sum = 0.0;
	for(auto i = 0u; i < N; i++){
		if(i != isel){
			auto val = G[area_list[isel]][area_list[i]] + G[area_list[i]][area_list[isel]];
			if(val != 0){
				if(gr[isel] == gr[i]) sum -= 2*val; else sum += 2*val; 
			}
		}
	}
	if(gr[isel] == false) gr[isel] = true; else gr[isel] = false;
	return sum;
}
