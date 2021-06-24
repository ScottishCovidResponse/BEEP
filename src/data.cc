//  The data inputs

#include <iostream>
#include <vector>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <math.h> 
#include <iomanip>

using namespace std;

#include "data.hh"	
#include "inputs.hh"
#include "mpi.hh"

#ifdef USE_Data_PIPELINE
#include "datapipeline.hh"
#include "table.hh"
#endif

/// Initialises data
Data::Data(Inputs &inputs, const Details &details, Mpi &mpi, DataPipeline *dp) : datapipeline(dp), data_directory(inputs.find_string("datadir","UNSET")), details(details), inputs(inputs)
{
	// The data directory
	if(data_directory == "UNSET") emsgroot("'data_directory' must be set.");

	threshold_str = inputs.find_string("threshold_str","*");         // The threshold string (if specified)
	
	nodata_str = inputs.find_string("nodata_str",".");               // The string which represents no data (if specified) 
	
	democat = inputs.find_demographic_category(strain);              // Gets information about demographic categories
	ndemocat = democat.size();	
	nstrain = strain.size();
	nage = democat[0].value.size();
		
	democat_change = inputs.find_democat_change();
	
	calc_democatpos();                                               // Calcualtes to demographic possibilities

	covar = inputs.find_covariates();                                // Finds information on covariates
	ncovar = covar.size();
	
	init_pop = inputs.find_string("init_pop","");                    // Finds file giving initial population
	
	level_effect = inputs.find_level_effect();                       // Finds information about level effects
	
	inputs.find_genQ(genQ,details,democat[0].value);                 // Finds information about geographical and age mixing

	read_data_files(inputs,mpi);                                     // Reads the data files
	
	area_effect = inputs.find_area_effect(area);                     // Finds information about area effect
}


/// Based on the different demographic categories, this calculates all the possible combinations
void Data::calc_democatpos()
{
	vector <unsigned int> count(ndemocat);                           // Defines all the demographic states
	for(auto &co : count) co = 0;
	
	unsigned int dc;
	do{
		democatpos.push_back(count);
		
		auto fl=0u;
		dc = 0;
		do{
			fl = 0;
			count[dc]++; if(count[dc] == democat[dc].value.size()){ fl = 1; count[dc] = 0; dc++;}
		}while(fl == 1 && dc < ndemocat);
	}while(dc < ndemocat);
	ndemocatpos = democatpos.size();

	ndemocatpos_per_age = ndemocatpos/nage;
	ndemocatpos_per_strain = ndemocatpos/nstrain;
	
	democatpos_name.resize(ndemocatpos);
	for(auto dp = 0u; dp < ndemocatpos; dp++){
		bool fl = false;
		stringstream ss;
		for(auto i = 0u; i < ndemocat; i++){
			if(democat[i].value.size() > 1){
				if(fl == true) ss << ",";
				fl = true;
				ss << democat[i].value[democatpos[dp][i]];
			}
		}
		democatpos_name[dp] = ss.str();
	}
	
	if(false){
		for(auto dp = 0u; dp < ndemocatpos; dp++){
			cout << dp << ": ";
			for(auto i = 0u; i < ndemocat; i++){
				cout << democat[i].value[democatpos[dp][i]] << ", ";
			}
			cout << endl;
		}
		emsg("Done");
	}
}


/// Reads in transition and area data
void Data::read_data_files(Inputs &inputs, Mpi &mpi)
{
	datatable = inputs.find_datatable(details);
	if(datatable.size() == 0) emsgroot("'data_tables' must be set.");	

	if(mpi.core == 0){
		check_datatable();

		string file = inputs.find_string("areas","UNSET");
		if(file == "UNSET") emsgroot("'areas' must be set");
			
		Table tab = load_table(file);

		read_initial_population(tab,file,inputs);
	
		read_covars();
		
		read_level_effect();
		
		load_datatable(tab);	
	
		load_democat_change(tab);	
	
		load_modification(inputs,tab);
		
		generate_matrices();
	}
	
	if(mpi.ncore > 1) mpi.copy_data(narea,area,nobs,obs,genQ,modification,covar,level_effect,democat_change);

	if(details.siminf == INFERENCE) set_datatable_weights();

	agedist.resize(nage); for(auto &aged : agedist) aged = 0;
	democatpos_dist.resize(ndemocatpos_per_strain); for(auto &dpd : democatpos_dist) dpd = 0;
	democat_dist.resize(ndemocat-1); 
	for(auto d = 0u; d < ndemocat-1; d++){
		democat_dist[d].resize(democat[d].value.size()); for(auto &dvd : democat_dist[d]) dvd = 0;
	}
		
	popsize = 0;
	for(auto c = 0u; c < narea; c++){ 
		for(auto co = 0u; co < area[c].pop_init.size(); co++){
			for(auto dp = 0u; dp < ndemocatpos; dp++){
				auto pop = area[c].pop_init[co][dp];

				democatpos_dist[dp%ndemocatpos_per_strain] += pop;
		
				agedist[democatpos[dp][0]] += pop;
			
				for(auto d = 0u; d < ndemocat-1; d++){
					democat_dist[d][democatpos[dp][d]] += pop; 
				}
			
				popsize += pop;
			}
		}
	}
	
	if(false){
		for(auto &aged : agedist) cout << aged << endl; cout << "agedist" << endl; emsg("Done");
	}
	
	for(auto &aged : agedist) aged /= popsize;
	for(auto &dpd : democatpos_dist) dpd /= popsize;
	for(auto d = 0u; d < ndemocat-1; d++){
		for(auto &dvd : democat_dist[d]) dvd /= popsize;
	}
	
	if(false){
		for(auto d = 0u; d < democat_dist.size(); d++){
			cout << " Demo :" << democat[d].name << endl;
			for(auto i = 0u; i < democat_dist[d].size(); i++) cout << d << " " << i << " " << democat_dist[d][i] << "," << endl;
		}
		
		for(auto a = 0u; a < agedist.size(); a++) cout << agedist[a] << ","; cout << "h" << endl;
		emsg("Done");	
	}
	
	narage = narea*nage;
	
	raw();                                                          // Does raw data analysis
}


/// Checks for any errors in the definition of the datatables
void Data::check_datatable()
{
	if(details.mode == SIM){                                        // If simulating checks file names are different
		for(auto dt = 0u; dt < datatable.size(); dt++){
			for(auto dt2 = 0u; dt2 < dt; dt2++){
				if(datatable[dt].file != "" && datatable[dt].file == datatable[dt2].file){
					stringstream ss; ss << "The file name '" << datatable[dt].file << "' is used more than once.";
					emsg(ss.str());
				}				
			}
		}
	}		
}		


/// If a heading doesn't exist in the table then this attempts to create one based on number ranges in head
void Data::check_or_create_column(Table &tab, string head, unsigned int d) const 
{
	auto col = 0u; while(col < tab.ncol && tab.heading[col] != head) col++;
	if(col < tab.ncol) return;
		
	vector <unsigned int> cols;
	
	auto spl = split(head,'*');
	if(spl.size() > 1){
		if(spl.size() > 2){
			emsg("The expression '"+head+"' from '"+democat[d].name+ "' can contain no more than one wildcard character '*'");
		}
		
		for(auto col = 0u; col < tab.ncol; col++){
			auto st = tab.heading[col];
			if(st.length() >= spl[0].length()+spl[1].length()){
				if(st.substr(0,spl[0].length()) == spl[0] && st.substr(st.length()-spl[1].length(),spl[1].length()) == spl[1]){
					cols.push_back(col);
				}				
			}
		}
	}
	else{
		auto i = find_char(head,"0123456789");
	
		if(i != UNSET){
			auto root = head.substr(0,i);
			auto rest = head.substr(i,head.length()-i);
			
			unsigned int from = UNSET, to = LARGE;
			if(rest.substr(rest.length()-1,1) == "+"){
				from = get_int_error(rest.substr(0,rest.length()-1));
			}
			else{
				auto spl = split(rest,'-');
				if(spl.size() == 2){
					from = get_int_error(spl[0]);
					to = get_int_error(spl[1]);
				}
			}
			
			if(from != UNSET && to != UNSET){
				for(auto col = 0u; col < tab.ncol; col++){
					auto st = tab.heading[col];
					if(st.length() > root.length()){
						if(st.substr(0,root.length()) == root){
							auto num = get_int_error(st.substr(root.length(),st.length()-root.length()));
							if(num != UNSET && num >= from && num <= to) cols.push_back(col);
						}
					}
				}
			}
		}
	}
	
	if(cols.size() == 0){
		if(d > 0) return;
		emsgroot("In file '"+tab.file+"' the column '"+head+"' as specified in '"+democat[d].name+"' does not exist");
	}
	
	cout << "  Generating '"+head+"' from columns: ";
	for(auto c = 0u; c < cols.size(); c++){ if(c != 0) cout << ", "; cout << "'" << tab.heading[cols[c]] << "'";}
	cout << endl;
	
	table_create_column(head,cols,tab);
}


/// Reads in the initial population
void Data::read_initial_population(Table &tab, string file, Inputs &inputs)
{
	auto comps = inputs.find_compartments();
	auto co_sus = inputs.find_susceptible_compartment();
	
	auto codecol = find_column(tab,"area");          
	for(auto row = 0u; row < tab.nrow; row++){                  // Sets up areas with zero population
		Area are;
		are.code = tab.ele[row][codecol];
		are.pop_init.resize(comps.size());
		for(auto co = 0u; co < comps.size(); co++){
			are.pop_init[co].resize(ndemocatpos);
			for(auto dp = 0u; dp < ndemocatpos; dp++) are.pop_init[co][dp] = 0;
		}
		area.push_back(are);			
	}
	narea = area.size();
	
	if(init_pop != "") read_init_pop_file(comps);              // Uses an 'init_pop' file to load initial population
	else read_initial_population_areas(co_sus,tab,file);       // Uses columns in 'areas' to load initial population
	
	for(auto &are : area){                                     // Sets the area population total
		auto total = 0.0; for(const auto &vec : are.pop_init){ for(auto val : vec) total += val;}
		are.total_pop = total;
	}
	
	if(false){
		for(const auto &are : area){
			cout << are.code << " " << are.total_pop << ": " << endl;
			for(auto co = 0u; co < comps.size(); co++){
				cout << " " << comps[co] << ": ";
				for(const auto &pop : are.pop_init[co]) cout << pop << ", ";
				cout << endl;	
			}
			cout << endl;	
		}
		emsg("Done");
	}	
}


/// Reads in information about the initial population from the 'init_pop' file
void Data::read_init_pop_file(const vector <string> comps)
{		
	Table tab = load_table(init_pop);
	
	auto area_col = find_column(tab,"area");
	auto pop_col = find_column(tab,"population");
	auto compartment_col = find_column(tab,"compartment");
	
	vector <unsigned int> demo_col(ndemocat);
	for(auto d = 0u; d < ndemocat; d++){
		if(democat[d].value.size() > 1) demo_col[d] = find_column(tab,democat[d].name);
		else demo_col[d] = UNSET;
	}
	
	for(auto row = 0u; row < tab.nrow; row++){
		auto areacode = tab.ele[row][area_col];
		auto c = 0u; while(c < narea && area[c].code != areacode) c++;
		if(c == narea) emsg("In 'init_pop' file '"+init_pop+"' the area '"+areacode+"' is not recognised");
		
		vector <unsigned int> index(ndemocat);
		for(auto d = 0u; d < ndemocat; d++){
			if(democat[d].value.size() == 1) index[d] = 0;
			else{
				auto demo = tab.ele[row][demo_col[d]];	
				auto imax = democat[d].value.size(); 
				auto i = 0u; while(i < imax && demo != democat[d].value[i]) i++;
				if(i == imax){
					emsg("In 'init_pop' file '"+init_pop+"' the value '"+demo+"' is not recognised as a demographic category");
				}
				index[d] = i;
			}
		}
		
		auto pop = get_double_positive(tab.ele[row][pop_col],"In file '"+init_pop+"'");
		
		auto comp = tab.ele[row][compartment_col];
		auto co = find_in(comps,comp);
		if(co == UNSET) emsg("In 'init_pop' file '"+init_pop+"' the compartment '"+comp+"' is not recognised.");
			
		unsigned int dp;
		for(dp = 0u; dp < ndemocatpos; dp++){
			auto d = 0u; while(d < ndemocat && index[d] == democatpos[dp][d]) d++; 
			if(d == ndemocat){ area[c].pop_init[co][dp] = pop; break;}
		}
		if(dp == ndemocatpos) emsgEC("Data",23);
	}
}


/// Reads in information about the initial population from the 'areas' file
void Data::read_initial_population_areas(const unsigned int co_sus, Table &tab, string file)
{
	for(auto d = 0u; d < democat.size()-1; d++){                    // Creates new columns in the table if needed
		for(auto val : democat[d].value){
			check_or_create_column(tab,val,d);
		}
	}

	vector <unsigned int> age_col(nage); 
	for(auto a = 0u; a < nage; a++) age_col[a] = find_column(tab,democat[0].value[a]);  
	
	vector <double> total_pop(narea);                               // Calculates the populations for different ages
	vector < vector <double> > pop_age; pop_age.resize(narea);
	for(auto c = 0u; c < narea; c++){
		pop_age[c].resize(nage);
		auto sum = 0.0;
		for(auto a = 0u; a < nage; a++){
			pop_age[c][a] = get_double_positive(tab.ele[c][age_col[a]],"In file '"+file+"'");
			sum += pop_age[c][a];
		}
		total_pop[c] = sum;
	}

	vector < vector < vector <double> > > pop_demo_frac; pop_demo_frac.resize(ndemocat-1);
	for(auto d = 1u; d < ndemocat-1; d++) pop_demo_frac[d] = get_demo_frac(d,tab,total_pop);
	
	for(auto c = 0u; c < narea; c++){
		for(auto dp = 0u; dp < ndemocatpos_per_strain; dp++){
			auto v = 0.0;
			for(auto j = 0u; j < democat.size()-1; j++){
				if(j == 0) v = pop_age[c][democatpos[dp][j]];
				else v *= pop_demo_frac[j][c][democatpos[dp][j]];
			}
			area[c].pop_init[co_sus][dp] = (unsigned int)(v+0.5);
		}
	}
}


/// Reads files which gives information about covaraites
void Data::read_covars()
{
	for(auto &cov : covar){
		cov.value.resize(narea);
		for(auto c = 0u; c < narea; c++) cov.value[c].resize(details.period);
		
		auto tab = load_table(cov.file);
		switch(cov.type){
			case AREA_COVAR:
				{
					auto col = find_column(tab,cov.name);
					
					if(tab.nrow != narea) emsgEC("Data",14);
					for(auto c = 0u; c < narea; c++){
						auto v = get_double(tab.ele[c][col],"In file '"+cov.file+"'");
			
						for(auto t = 0u; t < details.period; t++) cov.value[c][t] = v;
					}
				}
				break;
				
			case TV_COVAR:
				{
					auto date_col = find_column(tab,"Date");
					auto col = find_column(tab,cov.name);
					
					for(auto r = 0u; r < tab.nrow-1; r++){
						int ti = details.gettime(tab.ele[r][date_col],"For 'tv_covar' in file '"+tab.file+"'") - details.start;
						int tf = details.gettime(tab.ele[r+1][date_col],"For 'tv_covar' in file '"+tab.file+"'") - details.start;
				
						if(r == 0 && ti > 0){
							emsg("For 'tv_covar' in file '"+tab.file+"' the dates must start at or before the analysis period.");
						}
						
						if(r == tab.nrow-2 && tf < (int)details.period){
							emsg("For 'tv_covar' in file '"+tab.file+"' the dates must end at or after the analysis period.");
						}
						
						if(tf < ti) emsg("For 'tv_covar' in file '"+tab.file+"' the ordering of times is not right");
						if(ti == tf) emsg("For 'tv_covar' in file '"+tab.file+"' two rows have equal time");
						
						for(auto t = ti; t <= tf; t++){
							if(t >= 0 && t < (int)details.period){
								auto v = get_double(tab.ele[r][col],"In file '"+cov.file+"'");
								for(auto c = 0u; c < narea; c++) cov.value[c][t] = v;
							}
						}						
					}
				}
				break;
				
			case AREA_TV_COVAR:
				break;
		}
		
		if(cov.func == LOG_TRANS){
			for(auto c = 0u; c < narea; c++){
				for(auto t = 0u; t < details.period; t++){ 
					auto v = cov.value[c][t];
					if(v <= 0) emsg("In file '"+cov.file+"' the values must be positive for log transformation.");
					cov.value[c][t] = log(v);
				}
			}
		}
		
		auto av = 0.0;                                                         // Shifts covariates so average is zero
		for(auto c = 0u; c < narea; c++){
			for(auto t = 0u; t < details.period; t++) av += cov.value[c][t];
		}
		av /= narea*details.period;
		
		for(auto c = 0u; c < narea; c++){
			for(auto t = 0u; t < details.period; t++) cov.value[c][t] -= av;
		}
	}	
}


/// Reads in a file giving level effects
void Data::read_level_effect()
{
	if(level_effect.on == true){
		auto tab = load_table(level_effect.file);
		auto date_col = find_column(tab,"Date");
		vector <unsigned int> area_col(narea);
		for(auto c = 0u; c < narea; c++) area_col[c] = find_column(tab,area[c].code);
			
		level_effect.param_map.resize(details.period);
		for(auto t = 0u; t < details.period; t++) level_effect.param_map[t].resize(narea);
		
		vector <unsigned int> param_area(narea);
		
		auto imax = level_effect.ps.size();
		auto t = 0u;
		for(auto r = 0u; r <= tab.nrow; r++){
			unsigned int t_new;
			if(r == tab.nrow) t_new = details.period;
			else{
				auto ti = details.gettime(tab.ele[r][date_col],"For 'level_effect' in file '"+tab.file+"'");
				t_new = (ti - details.start);
			}
			
			if(t_new < t || t_new > details.period){
				emsg("In the file '"+level_effect.file+"' the dates are out of the analysis range");
			}				
			
			while(t < t_new){
				level_effect.param_map[t] = param_area;
				t++;
			}
			if(r == tab.nrow) break;
			
			for(auto c = 0u; c < narea; c++){
				auto val = tab.ele[r][area_col[c]];
				auto i = 0u; while(i < imax && val != level_effect.ps[i].name) i++;
				if(i == imax) emsg("In 'level_effect' the value '"+val+"' in file '"+level_effect.file+"' is not contained in 'param'");
				
				param_area[c] = i;
			}
		}
		
		level_effect.frac.resize(imax); for(auto i = 0u; i < imax; i++) level_effect.frac[i] = 0;
		auto N = 0u;	
		for(auto t = 0u; t < details.period; t++){
			for(auto c = 0u; c < narea; c++){
				level_effect.frac[level_effect.param_map[t][c]]++;
				N++;
			}
		}
		
		for(auto i = 0u; i < imax; i++) level_effect.frac[i] /= N;
		
		if(false){
			for(auto t = 0u; t < details.period; t++){
				cout << t << " ";
				for(auto c = 0u; c < narea; c++) cout << level_effect.param_map[t][c] << " "; cout << endl;
			}
			
			for(auto i = 0u; i < imax; i++) cout << level_effect.frac[i] << " Parameter Fraction" << endl;
			emsg("Done");	
		}
	}
}


/// Chops the root directory from a file, it is present
void Data::chop_dir(string &file, const string dir) const
{
	if(file.size() > dir.size()){
		if(file.substr(0,dir.size()) == dir) file = file.substr(dir.size()+1,file.size()-(dir.size()+1));
	}
}


/// Gets a column from a table
vector <string> Data::get_table_column_str(const unsigned int col, string file, const string dir) const
{
	vector <string> result;
	chop_dir(file,dir);
	auto tab = load_table(file,dir,true,true);
	for(auto r = 0u; r < tab.nrow; r++) result.push_back(tab.ele[r][col]);
	
	return result;
}


/// Gets a column from a table
vector <string> Data::get_table_column(const string col_name, string file, const string dir) const
{
	vector <string> result;
	chop_dir(file,dir);
	auto tab = load_table(file,dir,true,true);
	auto col = find_column(tab,col_name);
	for(auto r = 0u; r < tab.nrow; r++) result.push_back(tab.ele[r][col]);
	
	return result;
}


/// Loads a column of numbers from a table
vector <double> Data::get_table_column(const unsigned int col, string file, const string dir) const
{
	vector <double> result;	
	chop_dir(file,dir);
	auto tab = load_table(file,dir,true,true);
	if(col >= tab.ncol) emsg("The file '"+file+"' does not have column "+to_string(col));
	for(auto row = 0u; row < tab.nrow; row++) result.push_back(get_double(tab.ele[row][col],"In file '"+tab.file+"'"));
	return result;
}


/// Gets an array of data from a file and encodes using JSON (used to encode maps)
string Data::get_array_JSON(const string file, const string dir) const
{
	stringstream ss; 
	auto tab = load_table(file,dir,true,true); 
	ss << "[";
	for(auto row = 0u; row < tab.nrow; row++){
		if(row > 0) ss << ",";
		ss << "[";
		for(auto col = 1u; col < tab.ncol; col++){
			if(col > 1) ss << ",";
			ss << tab.ele[row][col];
		}		
		ss << "]";
	}
	ss << "]";
	return ss.str();
}


/// Gets a table from a file and encodes using JSON (used to encode maps)
string Data::get_table_JSON(const string file, const string dir) const
{
	stringstream ss; 
	auto tab = load_table(file,dir,true,true); 
	ss << "{ \"heading\":[";
	for(auto col = 0u; col < tab.ncol; col++){
		if(col > 0) ss << ",";
		ss << "\"" << tab.heading[col] << "\"";
	}
	ss << "]";
	
	
	ss << ",\"ele\":[";
	for(auto row = 0u; row < tab.nrow; row++){
		if(row > 0) ss << ",";
		ss << "[";
		for(auto col = 0u; col < tab.ncol; col++){
			if(col > 0) ss << ",";
			ss << "\"" << replace(tab.ele[row][col],"|",",") << "\"";
		}		
		ss << "]";
	}
	ss << "]}";
	
	return ss.str();
}


/// Gets the fraction of individuals in different demographic values of democat d 
vector < vector <double> > Data::get_demo_frac(const unsigned int d, const Table &tab, const vector <double> &pop) const
{
	vector < vector <double> > frac; frac.resize(tab.nrow);

	auto nval = democat[d].value.size();
	
	vector <unsigned int> v_st, unset_st, col_st;
	for(auto i = 0u; i < democat[d].value.size(); i++){
	auto c = find_column_noerror(tab,democat[d].value[i]);
		if(c != UNSET){ v_st.push_back(i); col_st.push_back(c);}
		else unset_st.push_back(i);
	}
	
	if(v_st.size() < nval-1){
		emsg("File '"+tab.file+"' must have columns for "+to_string(nval-1)+" of the "+to_string(nval)+" categories in '"+democat[d].name+"'");
	}
	
	vector <double> row_sum(tab.nrow);
	vector < vector <double> > row_val; row_val.resize(tab.nrow);
	for(auto row = 0u; row < tab.nrow; row++){
		auto sum = 0.0; 
		row_val[row].resize(col_st.size());
		for(auto i = 0u; i < col_st.size(); i++){
			row_val[row][i] = get_double_positive(tab.ele[row][col_st[i]],"In file '"+tab.file+"'");
			if(row_val[row][i] < 0){
				emsg("In file '"+tab.file+"' cannot have a negative value in column '"+tab.heading[col_st[i]]);
			}
			sum += row_val[row][i];
		}
		row_sum[row] = sum;
	}
	
	if(v_st.size() == nval){                                          // We have columns for each value in category
		auto error_fraction = 0.0, error_percent = 0.0, error_pop = 0.0;
		for(auto row = 0u; row < tab.nrow; row++){
			auto sum = row_sum[row]; 
			
			frac[row].resize(nval); for(auto i = 0u; i < nval; i++) frac[row][v_st[i]] = row_val[row][i]/sum;
				
			error_fraction += (1.0-sum)*(1.0-sum);
			error_percent += (100.0-sum)*(100.0-sum);
			error_pop += (pop[row]-sum)*(pop[row]-sum);
		}
		
		NumberType nt;
		if(error_fraction <= error_percent && error_fraction <= error_pop) nt = FRACTION;
		else{
			if(error_percent <= error_fraction && error_percent <= error_pop) nt = PERCENT;
			else{
				if(error_pop <= error_fraction && error_pop <= error_percent) nt = POPULATION;
				else emsgEC("Data",1);
			}
		}
				
		for(auto row = 0u; row < tab.nrow; row++){
			double target;
			switch(nt){
				case FRACTION: target = 1.0; break;
				case PERCENT: target = 100.0; break;
				case POPULATION: target = pop[row]; break;
			}
			
			if(row_sum[row] < 0.99*target || row_sum[row] > 1.01*target){
				stringstream ss; ss << "On row " << row+1 << " of file '"+tab.file+"' the columns ";
				for(auto i = 0u; i < nval; i++){
					if(i != 0){
						if(i < nval-1) ss << ", ";
						else ss << " and ";
					}
					ss << "'" << tab.heading[col_st[i]] << "'";
				}
				ss << " add up to " << row_sum[row] << " but should add up to " << target;
				warning(ss.str());
			}
		}
	}		
	else{                                                               // One column is missing (this must be imputted)
		auto max = 0.0;
		for(auto row = 0u; row < tab.nrow; row++){
			if(row_sum[row] > max) max = row_sum[row];
		}
		
		NumberType nt;
		if(max > 100) nt = POPULATION;
		else{
			if(max > 1){ 
				nt = PERCENT; 
				cout << "  The values are assumed to be percentages (add a column for '" << democat[d].value[unset_st[0]] << "' to disambiguate)." << endl;
			}
			else{
				nt = FRACTION; 
				cout << "  The values are assumed to be fractions (add a column for '" << democat[d].value[unset_st[0]] << "' to disambiguate)." << endl;
			}
		}
	
		for(auto row = 0u; row < tab.nrow; row++){
			double target;
			switch(nt){
				case FRACTION: target = 1.0; break;
				case PERCENT: target = 100.0; break;
				case POPULATION: target = pop[row]; break;
			}
			
			if(row_sum[row] > target){
				stringstream ss; ss << "On row " << row+1 << " of file '"+tab.file+"' ";
				if(nval-1 > 1){
					ss << "the columns ";
					for(auto i = 0u; i < nval-1; i++){
						if(i != 0){
							if(i < nval-2) ss << ", ";
							else ss << " and ";
						}
						ss << "'" << tab.heading[col_st[i]] << "'";
					}
					ss << " add up to " << row_sum[row];
				}
				else{
					ss << "the column '" << tab.heading[col_st[0]] << "' is " << row_sum[row];
				}
				ss << ", which is more than the total population size " << target;
				emsgroot(ss.str()); 
			}
			
			frac[row].resize(nval);
			for(auto i = 0u; i < nval-1; i++) frac[row][v_st[i]] = row_val[row][i]/target;
			frac[row][unset_st[0]] = 1-row_sum[row]/target;
		}
	}
	
	return frac;
}


/// Loads files giving demographic changes
void Data::load_democat_change(const Table &tabarea)
{
	for(auto &dcc : democat_change){   
		dcc.area = create_area_sel(tabarea,dcc.geo_filt);
		
		auto tab = load_table(dcc.file);
		
		auto d = 0u; while(d < ndemocat && democat[d].name != dcc.name) d++;
		if(d == ndemocat) emsg("In 'democat_change' the name '"+dcc.name+"' is not in 'democats'");
		dcc.d = d;
		
		auto date_col = find_column(tab,details.time_format_str);
		
		vector <int> times(tab.nrow);
		for(auto row = 0u; row < tab.nrow; row++){
			times[row] = details.gettime(tab.ele[row][date_col],"In file '"+tab.file+"'") + dcc.shift - details.start;
		}
		
		auto nval = democat[d].value.size();
	
		auto total_pop = 0.0; for(auto c : dcc.area) total_pop += area[c].total_pop;
		vector <double> pop_vec(tab.nrow); for(auto row = 0u; row < tab.nrow; row++) pop_vec[row] = total_pop;
		
		auto frac = get_demo_frac(d,tab,pop_vec);
		
		if(false){
			for(auto row = 0u; row < tab.nrow; row++){
				cout << row << " " << times[row] << " "; 
				for(auto i = 0u; i < nval; i++) cout << frac[row][i] << " "; 
				cout << endl;
			}
			emsg("Done");
		}
		
		if(times[0] > 0){
			cout << "  Demographic proportions assumed the same before '" << tab.ele[0][date_col] << "'" << endl;
		}
		if(tab.nrow < 2) emsg("File '"+dcc.file+"' must contain at least two rows of data.");
		
		if(times[tab.nrow-1] < (int)details.period-1){
			cout << "  Demographic proportions assumed the same after '" << tab.ele[tab.nrow-1][date_col] << "'" << endl;
		}
		
		dcc.frac.resize(details.ndivision);
		auto row = 0u;
		for(auto sett = 0u; sett < details.ndivision; sett++){
			int t = sett/details.division_per_time;
			while(row < tab.nrow-2 && times[row+1] <= t) row++;
		
			auto fr = 0.0;
			if(t < times[row]) fr = 0;
			else{	
				if(t > times[row+1]) fr = 1;
				else fr = double(t-times[row])/(times[row+1]-times[row]);
			}
			dcc.frac[sett].resize(nval);
			for(auto i = 0u; i < nval; i++) dcc.frac[sett][i] = frac[row][i]*(1-fr) + frac[row+1][i]*fr;
		}
		
		if(false){
			for(auto sett = 0u; sett < details.ndivision; sett++){
				cout << sett << " "; for(auto i = 0u; i < nval; i++) cout << dcc.frac[sett][i] << " "; cout << endl;
			}
			emsg("Done");
		}
		
		auto dp_sel = create_dp_sel(dcc.democats_filt);
		
		for(auto v = 0u; v < nval; v++){   // Groups dp into groups which vary in the selected demographic classification
			for(auto dp : dp_sel){
				if(democatpos[dp][d] == v && dp < ndemocatpos_per_strain){
					auto g = 0u; 
					while(g < dcc.dp_group.size()){
						auto dd = 0u; 
						while(dd < ndemocat && (dd == d || democatpos[dp][dd] == democatpos[dcc.dp_group[g][0]][dd])) dd++;
						if(dd == ndemocat) break;
						g++;
					}
				
					if(g < dcc.dp_group.size()) dcc.dp_group[g].push_back(dp);
					else{ 
						vector <unsigned int> vec; vec.push_back(dp);
						dcc.dp_group.push_back(vec);
					}
				}
			}
		}
		
		for(auto g = 0u; g < dcc.dp_group.size(); g++){
			if(dcc.dp_group[g].size() != nval) emsgEC("Data",2);
		}
		
		if(false){
			for(auto g = 0u; g < dcc.dp_group.size(); g++){
				cout << " Group " << g << endl;
				for(auto dp : dcc.dp_group[g]){
					for(auto dd = 0u; dd < ndemocat; dd++) cout << democat[dd].value[democatpos[dp][dd]] << ",";
					cout << "   ";
				}
				cout << endl;
			}
			emsg("Done");
		}
	}
}

	
/// Loads datatables
void Data::load_datatable(const Table &tabarea)
{
	bool sim;
	if(details.mode == SIM) sim = true; else sim = false;     

	for(auto i = 0u; i < datatable.size(); i++){ 
		switch(datatable[i].type){
			case POP: case POPFRAC: case TRANS: load_timeseries_datatable(tabarea,i,sim); break;	
			case MARGINAL: load_marginal_datatable(tabarea,i,sim); break;	
		}
	}
	
	nobs = obs.size();
	if(false){ cout << nobs <<" Number of observations" << endl; emsg("Done");}
}


/// Loads a table giving time series data
void Data::load_timeseries_datatable(const Table &tabarea, const unsigned int i, const bool sim)
{
	auto &dt = datatable[i];
		
	Observation ob;
	ob.factor_spline = UNSET;
	ob.w = UNSET;
	ob.obsmodel = dt.obsmodel; ob.shape = dt.shape; ob.invT = dt.invT;
	ob.datatable = i;
	
	auto datafilter = get_datafilter(tabarea,dt.geo_dep,dt.democats_dep,dt.geo_filt,dt.democats_filt,dt.type,dt.observation);

	vector <int> obs_times;
	
	bool table_loaded = false;
	Table tab;
	if(dt.optype == DATA){
		if(sim == false){
			table_loaded = true; 
			tab = load_table(dt.file);
			table_select_dates(dt.start,dt.end,dt.timestep,tab,"trans",dt.shift,obs_times);
		}
		else{
			for(auto t = dt.start; t <= dt.end; t += dt.timestep) obs_times.push_back(t+dt.shift);
		}
	}

	auto flag = false;
	for(const auto &df : datafilter){
		auto col=0u; if(table_loaded == true) col = find_column(tab,df.colname);
		
		ob.dp_sel = df.dp_sel;
		
		flag = true;
		ob.area = df.area;
		ob.graph = graph.size();
		
		string fulldesc;
		
		Graph gr;
		if(dt.file != ""){ 
			gr.tab = "Data Tables"; gr.tab2 = dt.file; gr.tab3 = df.colname;
			
			
			if(details.siminf == SIMULATE){
				fulldesc = "Data table "+dt.file+": This shows the temporal variation corresponding to the observations "+dt.observation+". ";
				fulldesc += " The red line shows this variation for the underlying state.";
				fulldesc += " The black line shows data generated from this state (for example if data is measured weekly, this stepped curve will have a period of 7 days). This data is saved in the file "+dt.file+" in the simulated data directory.";
			}
			else{
				fulldesc = "Data table "+dt.file+": This gives the temporal variation corresponding to the observations "+dt.observation+". ";
				fulldesc += " The black line shows the raw data (from the input file "+dt.file+").";
				fulldesc += " The red line (with dashed 95% credible intervals) shows the posterior distribution for the inferred underlying system state. Under a suitable model, and with accurate inferece, it would be expected that the black and red curves exhibit the same temporal behaviour (barring variation coming from a weak observation model).";	
				fulldesc += " Large deviations between the two lines indicate that either the model is not able to properly account for the actual observed system dynamics, or that inference has not converged on the true posterior distribution.";
			}
		}
		else{ 
			gr.tab = "State Outputs"; gr.tab2 = dt.plot_name; gr.tab3 = df.colname;
			fulldesc = "State Output: This shows the plot *"+dt.plot_name+"*, as specifed in the input TOML file.";
		}
		gr.fulldesc = fulldesc;
		
		gr.type = GRAPH_TIMESERIES;
		gr.file = df.file;
		gr.name = df.name;
		gr.desc = df.desc;
		gr.colname = df.colname;
		gr.factor_spline = UNSET;
		gr.datatable = i;
		gr.dp_sel = ob.dp_sel;
		gr.area = ob.area;
		
		ob.factor = dt.factor;
		if(dt.type == POPFRAC){
			auto total_pop = 0.0;
			if(dt.type == POPFRAC){
				for(auto c :	gr.area){
					for(auto dp : gr.dp_sel){
						for(const auto &vec : area[c].pop_init) total_pop += vec[dp];
					}
				}
			}
			ob.factor *= 1.0/total_pop;	
		}
		
		gr.factor = ob.factor;
			
		switch(dt.type){
			case POP: case POPFRAC:		
				for(auto row = 0u; row < obs_times.size(); row++){
					int ti = obs_times[row]*details.division_per_time;
					if(ti < 0 || ti > int(details.ndivision-1)){
						emsg("In 'data_tables' the file '"+dt.file+"' provides information which is out of the specified time range");
					}
											
					GraphPoint gp;
					gp.xi = ti/details.division_per_time;
					gp.xf = (ti+1)/details.division_per_time;
					gp.obs = obs.size();
					gr.point.push_back(gp);
					
					if(table_loaded == true) ob.value = get_data(tab.ele[row][col],"In file '"+df.file+"'",threshold_str,nodata_str); 
					else ob.value = UNSET;
					ob.sett_i = ti;
					ob.sett_f = ti+1;
				
					obs.push_back(ob);
				}
				break;
			
			case TRANS:
				{
					auto val = 0u;
					auto ti = obs_times[0]*details.division_per_time;
			
					for(auto row = 0u; row < obs_times.size(); row++){
						if(table_loaded == true){
							auto v = get_data(tab.ele[row][col],"In file '"+df.file+"'",threshold_str,nodata_str); 
							if(details.trans_combine != UNSET){
								if(v == UNSET){
									emsg("In file '"+df.file+"' a value for 'trans_combine' cannot be set in conjunction with unknown data");
								}
								if(v == THRESH){
									emsg("In file '"+df.file+"' a value for 'trans_combine' cannot be set in conjunction with thresholded data");
								}
							}
							val += v;
						}
						else val = UNSET;
										
						GraphPoint gp;
						unsigned int tf = (obs_times[row]+dt.timestep)*details.division_per_time;

						if(details.obs_section == true){
							for(auto t = ti+1; t < tf; t++){
								if(details.sec_define[t] == true) emsg("Observations cross particle filtering time points.");
							}
						}
						
						if(sim == true || row == obs_times.size()-1 || details.trans_combine == UNSET  
     						 || val >= details.trans_combine || (details.obs_section == true && details.sec_define[tf] == true)){
							gp.xi = ti/details.division_per_time;
							gp.xf = tf/details.division_per_time;
						
							gp.obs = obs.size();
							gr.point.push_back(gp);
						
							ob.sett_i = ti;
							ob.sett_f = tf;
							ob.value = val;
						
							obs.push_back(ob);
							val = 0;
							ti = tf;
						}
					}
				}
				break;
			
			case MARGINAL: emsgEC("Data",2); break;
		}
		
		dt.graph_ref.push_back(graph.size());
		graph.push_back(gr);
	} 
	if(flag == false) emsg("In 'data_tables' the data in '"+dt.file+"' was not used");
}


/// Loads a table giving marginal data
void Data::load_marginal_datatable(const Table &tabarea, const unsigned int i, const bool sim)
{
	auto &dt = datatable[i];
	
	if(dt.geo_dep == "" && dt.democats_dep == ""){
		emsg("For marginal distributions in 'data_tables' either 'democats_dep' or 'geo_dep' must be set");
	}

	if(dt.geo_dep != "" && dt.democats_dep != ""){
		emsg("For marginal distributions in 'data_tables', 'democats_dep' and 'geo_dep' cannot both be set");
	}
		
	string dep;
	if(dt.democats_dep != "") dep = dt.democats_dep; else dep = dt.geo_dep;
	
	Observation ob;
	ob.factor_spline = UNSET;
	ob.datatable = i;
	ob.w = UNSET;
	ob.obsmodel = dt.obsmodel; ob.shape = dt.shape; ob.invT = dt.invT; 
	ob.graph = graph.size();
	ob.factor = dt.factor;
	
	auto datafilter = get_datafilter(tabarea,dt.geo_dep,dt.democats_dep,dt.geo_filt,dt.democats_filt,dt.type,dt.observation);

	if(datafilter.size() == 0){ 
		emsg("In 'data_tables' the marginal distribution for '"+dt.observation+"' does not contain any valid entries");
	}
	
	Table tab; if(sim == false) tab = load_table(dt.file);
	
	auto df = datafilter[0];
	Graph gr;
	gr.type = GRAPH_MARGINAL;
	gr.file = df.file;
	gr.name = df.name;
	gr.desc = df.desc;
	
	gr.colname = "";
	gr.factor_spline = UNSET;
	gr.datatable = i;
	
	auto row = 0u;
	for(const auto &df : datafilter){
		if(sim == false){
			for(row = 0; row < tab.nrow; row++) if(tab.ele[row][0] == df.colname) break;
			if(row == tab.nrow) emsg("In 'data_tables' cannot find the value '"+df.colname+"' in the first column of '"+dt.file+"'");
		}
		
		dt.demolist.push_back(df.colname);
	
		ob.dp_sel = df.dp_sel;
		ob.area = df.area;
		
		GraphPoint gp;
		gp.xi = row;
		gp.xf = row;
		gp.obs = obs.size();
		gr.point.push_back(gp);
	
		int ti = dt.start*details.division_per_time;
		int tf = dt.end*details.division_per_time;
	
		if(ti < 0 || tf > int(details.ndivision)) emsg("In 'data_tables' the file '"+dt.file+"' contains information which is outside the inference time range");
		
		ob.sett_i = ti;
		ob.sett_f = tf;
		
		if(sim == true) ob.value = UNSET;
		else ob.value = get_data(tab.ele[row][1],"In file '"+df.file+"'",threshold_str,nodata_str); 
					
		obs.push_back(ob);
	
		if(sim == true) row++;
	}

	dt.graph_ref.push_back(graph.size());
	graph.push_back(gr);
}


/// Based on expressions for geo and democats in the data file this generates information for all the data columns
vector <DataFilter> Data::get_datafilter(const Table &tabarea, const string geo_dep, const string democats_dep, const string geo_filt, const string democats_filt, const DataType type, const string observation) const
{ 
	string file, name, desc;
	
	switch(type){
		case POP: desc = "Population"; break;
		case POPFRAC: desc = "Population fraction"; break;
		case TRANS: desc = "Transitions";  break;
		case MARGINAL: desc = "Marginal"; break;
	}
	
	file = desc+"_"+replace(observation,"->","-");
	name += observation;
	desc += " in "+observation;
	
	if(geo_filt != ""){
		file += "-"+ replace(geo_filt,":","=");
		name += " "+geo_filt;
		desc += " "+geo_filt;
	}
	
	if(democats_filt != ""){
		file += "-"+replace(democats_filt,":","=");
		name += " "+democats_filt;
		desc += " "+democats_filt;
	}

	vector <DataFilter> datafilter;
	
	if(democats_dep != ""){
		auto dc = 0u; while(dc < democat.size() && democat[dc].name != democats_dep) dc++;
		if(dc == democat.size()) emsg("In the expression 'democats_dep="+democats_dep+"' as a name in 'democats'");
		
		if(type != MARGINAL){
			file += "-"+democats_dep;
			name += " "+democats_dep;
			desc += " "+democats_dep;
			file += "="; name += ":"; desc += ":";
		}
	
		for(const auto &val : democat[dc].value){
			auto valst = val; if(type == MARGINAL) valst="";
			DataFilter df; 
			df.file = file+valst+".csv";
			df.colname = val;
			df.name = name+valst;
			df.desc = desc+valst;
			df.area = create_area_sel(tabarea,geo_filt);
			auto state = democat[dc].name+":"+val; if(democats_filt != "") state += ","+democats_filt;
			df.dp_sel = create_dp_sel(state);
			
			if(df.area.size() > 0 && df.dp_sel.size() > 0) datafilter.push_back(df);	
		}		
	}
	else{
		if(geo_dep != ""){
			auto geomap = create_geomap(tabarea,geo_dep);
			
			file += "-"+geo_dep;
			name += " "+geo_dep;
			desc += " "+geo_dep;
			if(type != MARGINAL){ file += "="; name += ":"; desc += ":";}
		
			for(const auto &gm : geomap){
				auto valst = gm.region; if(type == MARGINAL) valst="";
				DataFilter df; 
				df.file = file+valst+".csv";
				df.colname = gm.region;
				df.name = name+valst;
				df.desc = desc+valst;
				auto state = geo_dep+":"+gm.region; if(geo_filt != "") state += ","+geo_filt;
				df.area = create_area_sel(tabarea,state);
				df.dp_sel = create_dp_sel(democats_filt);
				if(df.area.size() > 0 && df.dp_sel.size() > 0) datafilter.push_back(df);		
			}		
		}
		else{
			DataFilter df; 
			df.file = file+".csv";
			df.colname = "Data";
			df.name = name;
			df.desc = desc;
			df.area = create_area_sel(tabarea,geo_filt);
			df.dp_sel = create_dp_sel(democats_filt);
			if(df.area.size() > 0 && df.dp_sel.size() > 0) datafilter.push_back(df);	
		}
	}
	
	if(datafilter.size() == 0){                                        // Checks that the filter generates valid columns
		stringstream ss;
		ss << "In 'data_tables' the filter ";
		if(democats_dep != "") ss << "'democats_dep=\"" << democats_dep << "\"' ";
		if(geo_dep != "") ss << "'geo_dep=\"" << geo_dep << "\"' ";
		if(democats_filt != "") ss << "'democats_filt=\"" << democats_filt << "\"' ";
		if(geo_filt != "") ss << "'geo_filt=\"" << geo_filt << "\"' ";
		ss << "does not generate any valid states";
		emsg(ss.str());
	}
	
	if(false){
		cout << endl;
		cout << "geo_dep: " << geo_dep << "   democats_dep: " << democats_dep <<  "      geo_filt: " <<  geo_filt << "   democats: " << democats_filt << endl;
		for(auto  df : datafilter){
			cout << df.name << "  Datafilter: ";
			cout << "Areas: "; for(auto c : df.area) cout << area[c].code << ","; cout << "  ";
			cout << "DP: "; 
			for(auto dp : df.dp_sel){ 
				for(auto i = 0u; i < ndemocat; i++) cout << democat[i].value[democatpos[dp][i]] << ","; cout << "  ";
			}
			cout << endl;
		}
		emsg("Done");
	}
	
	return datafilter;
}


/// Sets the weights used for all the observations
void Data::set_datatable_weights()
{
	vector <unsigned int> num_obs(datatable.size());
	vector <double> value_max(datatable.size());

	auto ngraph = 0u;
	for(auto &ob : obs){ if(ob.graph > ngraph) ngraph = ob.graph;}
	ngraph++;

	vector <double> graph_value_max(ngraph);
	for(auto gr = 0u; gr < ngraph; gr++){
		graph_value_max[gr] = 0;
	}
		
	for(auto dt = 0u; dt < datatable.size(); dt++){		
		num_obs[dt] = 0;
		value_max[dt] = 0;
	}
	
	for(auto &ob : obs){		
		auto dt = ob.datatable;
		auto gr = ob.graph;
		num_obs[dt]++;
		auto val = ob.value;
		if(val != UNKNOWN && val != THRESH){
			if(val > value_max[dt]) value_max[dt] = val;
			if(val > graph_value_max[gr]) graph_value_max[gr] = val;
		}
	}
	
	for(auto &ob : obs){		
		auto dt = ob.datatable;
		auto gr = ob.graph;
		if(gr >= ngraph) emsgEC("data",3);
		
		ob.w = datatable[dt].weight/num_obs[dt];
	}
}


/// Creates all the demographic categories consistent with a string 
vector <unsigned int> Data::create_dp_sel(const string dp_str) const
{
	vector <unsigned int> dp_sel;

	if(dp_str == ""){			
		for(auto dp = 0u; dp < ndemocatpos; dp++) dp_sel.push_back(dp);
	}
	else{
		vector < vector <unsigned int> > flag;
		flag.resize(ndemocat);
		for(auto d = 0u; d < ndemocat; d++){
			flag[d].resize(democat[d].value.size());
			for(auto i = 0u; i < democat[d].value.size(); i++) flag[d][i] = 0;
		}
			
		auto st = split(dp_str,',');
		
		for(auto val : st){
			auto valsp = split(val,':');
			if(valsp.size() != 2) emsg("The expression '"+dp_str+"' must specify the demographic category and its value");
			
			auto d = 0u; while(d < ndemocat && democat[d].name != valsp[0]) d++;
			if(d == ndemocat) emsg("Cannot find the value '"+valsp[0]+"' in 'ages' or 'democats'");
			
			auto val_list = split(valsp[1],'|');
			
			for(auto i = 0u; i < val_list.size(); i++){                 // Checks for repeated values
				for(auto j = i+1; j < val_list.size(); j++){
					if(val_list[i] == val_list[j]) emsg("In 'data_tables' the value '"+val_list[i]+"' must not be repeated.");
				}
			}
			
			for(auto i = 0u; i < val_list.size(); i++){ 
				auto pos = 0u; while(pos < democat[d].value.size() && democat[d].value[pos] != val_list[i]) pos++;
				if(pos == democat[d].value.size()){
					emsg("In 'data_tables' cannot find the value '"+val_list[i]+"' in 'name=\""+valsp[0]+"\"'");
				}
			
				flag[d][pos]++;
			}

			for(auto dd = 0u; dd < ndemocat; dd++){
				if(dd != d){ for(auto i = 0u; i < democat[dd].value.size(); i++) flag[dd][i]++;}
			}
		}
		
		for(auto dp = 0u; dp < ndemocatpos; dp++){
			auto d = 0u; while(d < ndemocat && flag[d][democatpos[dp][d]] == st.size()) d++;
			if(d == ndemocat) dp_sel.push_back(dp);
		}
		
		if(false){
			for(auto d = 0u; d < ndemocat; d++){
				cout << d << ": ";
				for(auto i = 0u; i < democat[d].value.size(); i++) cout << flag[d][i] <<","; cout << "flag" << endl;
			}
		}
	}

	return dp_sel;
}
 
 
/// Creates all the areas consistent with a string 
vector <unsigned int> Data::create_area_sel(const Table &tabarea, const string str) const
{
	vector <unsigned int> area_sel; 
	
	if(str == ""){
		for(auto c = 0u; c < narea; c++) area_sel.push_back(c);
	}
	else{
		vector <unsigned int> flag(narea);
		for(auto c = 0u; c < narea; c++) flag[c] = 0;
		
		auto spl = split(str,',');
		for(const auto &sp : spl){
			auto geosp = split(sp,':');
			if(geosp.size() != 2) emsg("There was a problem with the expression '"+str+"'");
		
			auto geomap = create_geomap(tabarea,geosp[0]);
			auto reg_list = split(geosp[1],'|');
			
			for(auto i = 0u; i < reg_list.size(); i++){                  // Checks for repeated values
				for(auto j = i+1; j < reg_list.size(); j++){
					if(reg_list[i] == reg_list[j]) emsg("In 'data_tables' the value '"+reg_list[i]+"' must not be repeated.");
				}
			}
			
			for(auto reg_str : reg_list){
				bool fl = false;
				for(auto gm : geomap){
					if(gm.region == reg_str){
						fl = true;
						for(auto c : gm.area) flag[c]++;
					}
				}
				if(fl == false)  emsg("In 'data_tables' the value '"+reg_str+"' is not found in '"+geosp[0]+"'."); 
			} 
		}
		
		for(auto c = 0u; c < narea; c++){
			if(flag[c] == spl.size()) area_sel.push_back(c);
		}
	}
	
	return area_sel;
}
	

/// Based on a column in the area file a geogaphical mapping is created 
vector <GeographicMap> Data::create_geomap(const Table &tab, const string geo) const 
{
	vector <GeographicMap> geomap;
	
	if(geo == "all"){
		GeographicMap geoadd;
		geoadd.region = "all";
		for(auto row = 0u; row < tab.nrow; row++) geoadd.area.push_back(row);
		geomap.push_back(geoadd);	
	}
	else{
		auto col = find_column(tab,geo);
		for(auto row = 0u; row < tab.nrow; row++){
			auto st = tab.ele[row][col];
			if(st != "NA"){
				auto g = 0u; while(g < geomap.size() && geomap[g].region != st) g++;
				if(g < geomap.size()) geomap[g].area.push_back(row);
				else{
					GeographicMap geoadd;
					geoadd.region = st;
					geoadd.area.push_back(row);
					geomap.push_back(geoadd);
				}
			}
		}
	}
	
	return geomap;
}

		
/// Loads a table from the data pipeline
Table Data::load_table_from_datapipeline(const string file) const
{
	Table tab;
#ifdef USE_Data_PIPELINE
	Table dptable = datapipeline->read_table(file,"default");

	tab.file = filebasename(file);
	tab.heading = dptable.get_column_names();
	tab.ncol = tab.heading.size();
	
	vector< vector <string> > cols;
	
	// Load each column as a string into cols
	for (size_t j = 0; j < tab.ncol; j++) {
		cols.push_back(dptable.get_column_as_string(tab.heading.at(j)));
	}

	for (size_t i = 0; i < dptable.get_column_size(); i++) {
		vector<string> vec;

		for (size_t j = 0; j < tab.ncol; j++) {
			vec.push_back(cols.at(j).at(i));
		}
		tab.ele.push_back(vec);
	}

	tab.nrow = tab.ele.size();

	cout << "Loaded table '" << file << "' from data pipeline" << endl;
#else
	emsg("load_tablefromdatapipeline for '"+file+"' cannot be called as data pipeline is not compiled in");
#endif

	return tab;
}


/// Loads a table from a file
Table Data::load_table_from_file(const string file, const string dir, const bool heading, const bool supop, const char sep) const
{
	ifstream in;

	string used_file;
	if(dir != ""){
    used_file = dir+"/"+file;
		in.open(used_file.c_str());
		if(!in) emsg("Cannot open the file '"+dir+"/"+file+"'.");
	}
	else{
    used_file = details.output_directory+"/Simulated_data/"+file;
		in.open(used_file.c_str());
		if(!in){
      used_file = data_directory+"/"+file;
			in.open(used_file.c_str());
			if(!in) emsg("Cannot open the file '"+used_file+"'");
		}
	}
	
	if(supop == false) cout << "Loaded table '" << file << "'." << endl;

	Table tab;
	tab.file = file;
	
	string line;
	if(heading == true){
		do{
			getline(in,line);
		}while(line.substr(0,1) == "#");
		
		stringstream ss(line);
		do{
			string st;
			getline(ss,st,sep); strip(st);
			tab.heading.push_back(st);
			if(ss.eof()) break;
		}while(true);
		tab.ncol = tab.heading.size();
	}
	else tab.ncol = 0;
		
	do{
		vector <string> vec;
		getline(in,line);
		if(in.eof()) break;
				
		stringstream ss(line);
		do{
			string st;
			getline(ss,st,sep); strip(st);
			vec.push_back(st);
			if(ss.eof()) break;
		}while(true);
		if(tab.ncol == 0) tab.ncol = vec.size();
		else{
			if(vec.size() != tab.ncol) emsg("Rows in the file '"+file+"' do not all share the same number of columns.");
		}
		
		tab.ele.push_back(vec);
	}while(true);
	tab.nrow = tab.ele.size();
	
	return tab;
}


/// Loads a table from a file (if dir is specified then this directory is used)
Table Data::load_table(const string file, const string dir, const bool heading, const bool supop) const
{
	if(stringhasending(file,".txt")){ return load_table_from_file(file,dir,heading,supop,'\t');} 
	else {
		if(stringhasending(file,".csv")){ return load_table_from_file(file,dir,heading,supop,',');} 
		else return load_table_from_datapipeline(file);
	}
}


/// Creates a new column by adding together existing columns		
void Data::table_create_column(const string head, const vector <unsigned int> &cols, Table &tab) const
{
	tab.heading.push_back(head);
	for(auto &trow : tab.ele){
		auto sum = 0u; for(auto i = 0u; i < cols.size(); i++) sum += get_int(trow[cols[i]].c_str(),"In file '"+tab.file+"'");
		stringstream ss; ss << sum;
		trow.push_back(ss.str());
	}
	tab.ncol++;
}


/// Selects dates as specified in the TOML file
void Data::table_select_dates(int &start, int &end, unsigned int &timestep, Table &tab, const string type, const unsigned int shift, vector <int> &times) const 
{
	auto c = find_column(tab,details.time_format_str);
	
	times.resize(tab.nrow);
	for(auto row = 0u; row < tab.nrow; row++){
		times[row] = details.gettime(tab.ele[row][c],"In file '"+tab.file+"'") + shift - details.start;
	}
	
	int row_start = 0u; 
	while(row_start < int(tab.nrow) && (times[row_start] < 0 || (start != UNSET && times[row_start] < start))) row_start++;
	if(row_start == (int)tab.nrow) emsgroot("In file '"+tab.file+"' there is no data in the specified time period");
	if(start != UNSET && start != times[row_start]){
		emsg("In data_tables' the value for 'start' does not agree with file '"+tab.file+"'");
	}
	start = times[row_start];
	
	int row_end = tab.nrow-1; 
	while(row_end >= 0 && (times[row_end] >= (int)details.period || (end != UNSET && times[row_end] >= int(end)))) row_end--;
	
	if(row_end == -1) emsgroot("In file '"+tab.file+"' there is no data in the specified time period");
	if(end != UNSET && end != times[row_end]){
		emsg("In data_tables' the value for 'end' does not agree with file '"+tab.file+"'");
	}
	end = times[row_end];
	
	if(type != "trans") timestep = 1;
	else{
		if(row_start == row_end){
			if(timestep == UNSET) emsg("In file '"+tab.file+"' with only one row of data a value for 'timestep' must be set.");
		}
		else{
			auto tstep = UNSET;
			for(auto row = row_start; row < row_end; row++){
				if(times[row+1] < times[row]) emsg("In file '"+tab.file+"' the measurements must be time ordered.");
				if(times[row+1] == times[row]) emsg("In file '"+tab.file+"' each measurement must be at a different time.");
				
				unsigned int ts = times[row+1]-times[row]; 
				if(tstep == UNSET) tstep = ts;
				else{ if(ts != tstep) emsg("In file '"+tab.file+"' measurements should have an equal timestep.");}
			}
			
			if(timestep == UNSET) timestep = tstep;
			else{
				if(timestep != tstep){
					emsg("In 'data_tables' the value for 'timestep' does not agree with the timestep in file '"+tab.file+"'.");
				}
			}
			
			if(times[row_end]+timestep > details.period){ row_end--; end = times[row_end];}
		}
	}

	if(row_start != 0){
		if(row_start == 1) cout << "  The first row lies outside the analysis period and is not used.";
		else cout << "  The first " << row_start << " rows lie outside the analysis period and are not used.";
		cout << endl;
	}

	if(row_end != (int)tab.nrow-1){
		auto num = (int)tab.nrow-1 - row_end;
		if(num == 1) cout << "  The last row lies outside the analysis period and is not used.";
		else cout << "  The last " << num << " rows lie outside the analysis period and are not used.";
		cout << endl;
	}

	for(int row = tab.nrow-1; row > row_end; row--){
		tab.ele.erase(tab.ele.begin()+row);
		times.erase(times.begin()+row);
		tab.nrow--;
	}

	for(int row = row_start-1; row >= 0; row--){
		tab.ele.erase(tab.ele.begin()+row);
		times.erase(times.begin()+row);
		tab.nrow--;
	}

	if(false){
		for(auto row = 0u; row < tab.ele.size(); row++){
			cout << times[row] << ": ";
			for(const auto &ele : tab.ele[row]) cout << ele << " ";
			cout << "TABLE" << endl;
		}
		emsg("Done");
	}
}				


/// Finds a column in a table
unsigned int Data::find_column(const Table &tab, const string name) const
{
	auto c = find_column_noerror(tab,name); 
	if(c == UNSET) emsg("Cannot find the column heading '"+name+"' in file '"+tab.file+"'.");
	return c;
}		


// Finds a column in a table (but doesn't return an error if can't find)
unsigned int Data::find_column_noerror(const Table &tab, const string name) const
{
	for(auto c = 0u; c < tab.ncol; c++){
		if(toLower(tab.heading[c]) == toLower(name)) return c;
	}
	return UNSET;	
}
			
			
/// Determines if a vector contains a given element
bool Data::vector_contains(const vector <unsigned int> &vec, const unsigned int num) const
{
	if(find(vec.begin(),vec.end(),num) != vec.end()) return true;
	return false;
}	


/// Determines if a string vector contains a given element
bool Data::vector_contains(const vector <string> &vec, const string num) const
{
	if(find(vec.begin(),vec.end(),num) != vec.end()) return true;
	return false;
}	


/// Removes an element from a vector
void Data::vector_remove(vector <unsigned int> &vec, const unsigned int num) const
{
	auto i = find(vec.begin(),vec.end(),num);
	if(i == vec.end()) emsgEC("Data",4);
	vec.erase(i);
}	


/// Loads up model modification information  
void Data::load_modification(Inputs &inputs, const Table &tabarea)  
{
	if(details.mode == PREDICTION){
		inputs.find_modification(details,modification);
		
		for(auto &cf : modification){
			cf.dp_sel = create_dp_sel(cf.democats_filt);
			cf.area = create_area_sel(tabarea,cf.geo_filt);  
		}
	}
}


/// Outputs properties of data to the terminal
string Data::print() const
{
	stringstream ss;
	ss << endl;
	
	ss << "Age categories: " << endl << "  ";
	ss << democat[0].name << ": ";
	for(auto j = 0u; j < democat[0].value.size(); j++){
		if(j != 0) ss << ", ";
		ss << democat[0].value[j];
		if(democat[0].sus_vari == true){
			ss << " sus='" << democat[0].ps[j].name << "'";
			ss << " value='" << democat[0].ps[j].value << "'";
			ss << " prior='" << democat[0].ps[j].prior << "'";
		}
	}
	ss << endl << endl;

	if(covar.size() > 0){
		ss << "Area covariates: " << endl;
		ss << "  ";
		for(const auto &cov : covar){
			ss << cov.name;
			ss << "   param='" << cov.ps.name << "'";
			if(cov.ps.value != "") ss << "   value='" << cov.ps.value << "'"; 
			if(cov.ps.prior != "") ss << "   prior='" << cov.ps.prior << "'";
			ss << endl;
		}
		ss << endl;
	}	
	
	if(democat.size() > 1){
		ss << "Demographic categories: " << endl;
		for(auto k = 1u; k < democat.size(); k++){
			ss << "  ";
			ss << democat[k].name << ": ";
			for(auto j = 0u; j < democat[k].value.size(); j++){
				if(j != 0) ss << ", ";
				ss << democat[k].value[j];
				if(democat[k].sus_vari == true){
					ss << " sus='" << democat[k].ps[j].name << "'";
					if(democat[k].ps[j].value != "") ss << " value='" << democat[k].ps[j].value << "'";
					if(democat[k].ps[j].prior != "") ss << " prior='" << democat[k].ps[j].prior << "'";
				}
			}	
			ss << endl;
		}
		ss << endl;
	}
	
	return ss.str();
}


/// Prints all the observations made
void Data::print_obs() const
{
	for(auto &ob : obs){
		cout << "Datatable: " << ob.datatable << "  ";
		cout << "Graph: " << ob.graph << "  ";
		cout << "Value: " << ob.value << "  ";
		cout << "ti: " << ob.sett_i << "  ";
		cout << "tf: " << ob.sett_f << "  ";
		cout << "w: " << ob.w << "  ";
		cout << "factor: " << ob.factor << "  ";
		cout << "area: "; for(auto c : ob.area) cout <<c << ",  ";
		cout << "dp_sel: "; for(auto dp : ob.dp_sel) cout << dp << ",  ";
		cout << endl << endl;
	}
}
