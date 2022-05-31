//  The data inputs

#include <iostream>
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
		emsgroot("Done");
	}
}


/// Reads in transition and area data
void Data::read_data_files(Inputs &inputs, Mpi &mpi)
{
	datatable = inputs.find_datatable(details);
	if(datatable.size() == 0) emsgroot("'data_tables' must be set.");	

	if(mpi.core == 0){
		areas_file = inputs.find_string("areas","UNSET");
		if(areas_file == "UNSET") emsgroot("'areas' must be set");
			
		Table tab = load_table(areas_file);

		read_initial_population(tab,inputs);
	
		check_datatable();

		read_covars();
		
		read_level_effect();
		
		load_datatable(tab);	

		load_democat_change(tab);	
	
		load_modification(inputs,tab);

		generate_matrices();
	}
	
	if(mpi.ncore > 1) mpi.copy_data(narea,area,nobs,obs,genQ,modification,covar,level_effect,democat_change);

	if(details.siminf == INFERENCE){
		set_obs_sd();
		normal_approximation();
	}
	
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
		for(auto &aged : agedist) cout << aged << endl; 
		cout << "agedist" << endl; 
		emsgroot("Done");
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
		
		for(auto a = 0u; a < agedist.size(); a++) cout << agedist[a] << ","; 
		cout << "h" << endl;
		emsgroot("Done");	
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
					emsgroot(ss.str());
				}				
			}
		}
	}		
	
	for(const auto &dt : datatable){
		if(narea == 1 && dt.geo_dep != ""){
			emsgroot("The data table with file '"+dt.file+"' cannot have 'geo_dep' set because there is only one area.");
		}
	}
}		


/// If a heading doesn't exist in the table then this attempts to create one based on number ranges in head
void Data::check_or_create_column(Table &tab, string head, unsigned int d) const 
{
	auto col = find_column_noerror(tab,head);  
	if(col != UNSET) return;
		
	vector <unsigned int> cols;
	
	auto spl = split(head,'*');
	if(spl.size() > 1){
		if(spl.size() > 2){
			emsgroot("The expression '"+head+"' from '"+democat[d].name+ "' can contain no more than one wildcard character '*'");
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
	for(auto c = 0u; c < cols.size(); c++){ if(c != 0) cout << ","; cout << "'" << tab.heading[cols[c]] << "'";}
	cout << endl;
	
	table_create_column(head,cols,tab);
}


/// Reads in the initial population
void Data::read_initial_population(Table &tab, Inputs &inputs)
{
	unsigned int co_sus;
	auto comps = inputs.find_compartment_names(co_sus,democat,democatpos);

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
	
	if(init_pop != "") read_init_pop_file(comps);         // Uses an 'init_pop' file to load initial population
	else read_initial_population_areas(co_sus,tab);       // Uses columns in 'areas' to load initial population
	
	for(auto &are : area){                                // Sets the area population total
		auto total = 0.0; for(const auto &vec : are.pop_init){ for(auto val : vec) total += val;}
		are.total_pop = total;
		
		are.pop_dp.resize(ndemocatpos);
		for(auto dp = 0u; dp < ndemocatpos; dp++){
			auto sum = 0.0; for(auto co = 0u; co < comps.size(); co++) sum += are.pop_init[co][dp];
			are.pop_dp[dp] = sum;
		}
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
		emsgroot("Done");
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
	
	for(auto c = 0u; c < narea; c++){
		for(auto co = 0u; co < comps.size(); co++){
			for(auto dp = 0u; dp < ndemocatpos; dp++){
				area[c].pop_init[co][dp] = UNSET;
			}
		}
	}
	
	for(auto row = 0u; row < tab.nrow; row++){
		auto areacode = tab.ele[row][area_col];
		auto c = 0u; while(c < narea && area[c].code != areacode) c++;
		if(c == narea) emsgroot("In 'init_pop' file '"+init_pop+"' the area '"+areacode+"' is not recognised");
		
		vector <unsigned int> index(ndemocat);
		for(auto d = 0u; d < ndemocat; d++){
			if(democat[d].value.size() == 1) index[d] = 0;
			else{
				auto demo = tab.ele[row][demo_col[d]];	
				auto imax = democat[d].value.size(); 
				auto i = 0u; while(i < imax && demo != democat[d].value[i]) i++;
				if(i == imax){
					emsgroot("In 'init_pop' file '"+init_pop+"' the value '"+demo+"' is not recognised as a demographic category");
				}
				index[d] = i;
			}
		}
		
		auto pop = get_double_positive(tab.ele[row][pop_col],"In file '"+init_pop+"'");
		
		auto comp = tab.ele[row][compartment_col];
		
		vector <unsigned int> co_list;
		for(auto co = 0u; co < comps.size(); co++){
			if(comps[co] == comp) co_list.push_back(co);
		}		
		if(co_list.size() == 0) emsgroot("In 'init_pop' file '"+init_pop+"' the compartment '"+comp+"' is not recognised.");

		unsigned int dp;
		for(dp = 0u; dp < ndemocatpos; dp++){
			auto d = 0u; while(d < ndemocat && index[d] == democatpos[dp][d]) d++; 
			if(d == ndemocat){ 
				for(auto i = 0u; i < co_list.size(); i++){   // This places all the population in the first compartment with the specified name
					auto co = co_list[i];
					
					if(i == 0) area[c].pop_init[co][dp] = pop;
					else area[c].pop_init[co][dp] = 0;
				}
				break;
			}
		}
		if(dp == ndemocatpos) emsgEC("Data",23);
	}
	
	for(auto co = 0u; co < comps.size(); co++){
		auto fl = false, fl2 = false;
		for(auto c = 0u; c < narea; c++){
			for(auto dp = 0u; dp < ndemocatpos; dp++){
				if(area[c].pop_init[co][dp] == UNSET) fl = true;
				else fl2 = true;
					
			}
		}
		if(fl == true){
			if(fl2 == false) emsgroot("The initial populations for '"+comps[co]+"' are not specified");
			else emsgroot("The initial populations for '"+comps[co]+"' are not all specified");
		}
	}
}


/// Reads in information about the initial population from the 'areas' file
void Data::read_initial_population_areas(const unsigned int co_sus, Table &tab)
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
			pop_age[c][a] = get_double_positive(tab.ele[c][age_col[a]],"In file '"+areas_file+"'");
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
					auto date_col = find_column(tab,details.time_format_str);
					auto col = find_column(tab,cov.name);
					
					for(auto r = 0u; r < tab.nrow-1; r++){
						int ti = details.gettime(tab.ele[r][date_col],"For 'tv_covar' in file '"+tab.file+"'") - details.start;
						int tf = details.gettime(tab.ele[r+1][date_col],"For 'tv_covar' in file '"+tab.file+"'") - details.start;
				
						if(r == 0 && ti > 0){
							emsgroot("For 'tv_covar' in file '"+tab.file+"' the dates must start at or before the analysis period.");
						}
						
						if(r == tab.nrow-2 && tf < (int)details.period-1){
							emsgroot("For 'tv_covar' in file '"+tab.file+"' the dates must end at or after the analysis period.");
						}
						
						if(tf < ti) emsgroot("For 'tv_covar' in file '"+tab.file+"' the ordering of times is not right");
						if(ti == tf) emsgroot("For 'tv_covar' in file '"+tab.file+"' two rows have equal time");
						
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
				{
					auto date_col = find_column(tab,details.time_format_str);
					vector <unsigned int> cols(narea);
					for(auto c = 0u; c < narea; c++) cols[c] = find_column(tab,area[c].code);
			
					for(auto r = 0u; r < tab.nrow-1; r++){
						int ti = details.gettime(tab.ele[r][date_col],"For 'area_tv_covar' in file '"+tab.file+"'");
						int tf = details.gettime(tab.ele[r+1][date_col],"For 'area_tv_covar' in file '"+tab.file+"'");
						ti -= details.start; tf -= details.start;
				
						if(r == 0 && ti > 0){
							emsgroot("For 'area_tv_covar' in file '"+tab.file+"' the dates must start at or before the analysis period.");
						}
						
						if(r == tab.nrow-2 && tf < (int)details.period-1){
							emsgroot("For 'area_tv_covar' in file '"+tab.file+"' the dates must end at or after the analysis period.");
						}
						
						if(tf < ti) emsgroot("For 'area_tv_covar' in file '"+tab.file+"' the ordering of times is not right");
						if(ti == tf) emsgroot("For 'area_tv_covar' in file '"+tab.file+"' two rows have equal time");
						
						for(auto t = ti; t <= tf; t++){
							if(t >= 0 && t < (int)details.period){
								for(auto c = 0u; c < narea; c++){
									auto v = get_double(tab.ele[r][cols[c]],"In file '"+cov.file+"'");
									cov.value[c][t] = v;
								}
							}
						}			
					}
				}
				break;
		}
		
		if(cov.func == LOG_TRANS){
			for(auto c = 0u; c < narea; c++){
				for(auto t = 0u; t < details.period; t++){ 
					auto v = cov.value[c][t];
					if(v <= 0) emsgroot("In file '"+cov.file+"' the values must be positive for log transformation.");
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
		auto date_col = find_column(tab,details.time_format_str);
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
				emsgroot("In the file '"+level_effect.file+"' the dates are out of the analysis range");
			}				
			
			while(t < t_new){
				level_effect.param_map[t] = param_area;
				t++;
			}
			if(r == tab.nrow) break;
			
			for(auto c = 0u; c < narea; c++){
				auto val = tab.ele[r][area_col[c]];
				auto i = 0u; while(i < imax && val != level_effect.ps[i].name) i++;
				if(i == imax) emsgroot("In 'level_effect' the value '"+val+"' in file '"+level_effect.file+"' is not contained in 'param'");
				
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
				for(auto c = 0u; c < narea; c++) cout << level_effect.param_map[t][c] << " ";
				cout << endl;
			}
			
			for(auto i = 0u; i < imax; i++) cout << level_effect.frac[i] << " Parameter Fraction" << endl;
			emsgroot("Done");	
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


/// Gets a table (used in output)
Table Data::get_table(const string file, const string dir) const 
{
	string file_copy = file;
	chop_dir(file_copy,dir);
	return load_table(file_copy,dir,true,true);
}


/// Gets a column from a table
vector <string> Data::get_table_column_str(const unsigned int col, const Table &tab) const
{
	if(tab.file == "") emsgEC("Data",31);
	
	vector <string> result;
	for(auto r = 0u; r < tab.nrow; r++) result.push_back(tab.ele[r][col]);
	
	return result;
}


/// Gets a column from a table
vector <string> Data::get_table_column(const string col_name, const Table &tab) const
{
	if(tab.file == "") emsgEC("Data",32);
	
	vector <string> result;
	auto col = find_column(tab,col_name);
	for(auto r = 0u; r < tab.nrow; r++) result.push_back(tab.ele[r][col]);
	
	return result;
}


/// Loads a column of numbers from a table
vector <double> Data::get_table_column(const unsigned int col, const Table &tab) const
{
	vector <double> result;	
	
	if(tab.file == "") emsgEC("Data",33);
	if(col >= tab.ncol) emsgroot("The file '"+tab.file+"' does not have column "+to_string(col));
	for(auto row = 0u; row < tab.nrow; row++){
		result.push_back(get_double(tab.ele[row][col],"In file '"+tab.file+"'"));
	}
	
	return result;
}


// Gets a matrix from a file
Matrix Data::get_matrix(const string file, const string dir) const
{
	Matrix mat;
	
	auto tab = load_table(file,dir,true,true);

	mat.N = tab.nrow;
	mat.ele.resize(tab.nrow);
	for(auto j = 0u; j < tab.nrow; j++){
		mat.ele[j].resize(tab.nrow);
		for(auto i = 0u; i < tab.nrow; i++){
			mat.ele[j][i] = get_double(tab.ele[j][i+1],"Error loading file'"+file+"'");
		}		
	}
	
	return mat;
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
			//ss << "\"" << replace(tab.ele[row][col],"|",",") << "\"";
			ss << "\"" << tab.ele[row][col] << "\"";
		}		
		ss << "]";
	}
	ss << "]}";
	
	return ss.str();
}


/// Gets table colums from a file and encodes using JSON (used to encode maps)
string Data::get_table_cols_JSON(const string file, const string dir, const vector <unsigned int> cols) const
{
	stringstream ss; 
	auto tab = load_table(file,dir,true,true); 
	ss << "{ \"heading\":[";
	bool fl = false;
	for(auto col : cols){
		if(col >= tab.ncol) emsgEC("Data",34);
		
		if(fl == true) ss << ","; 
		fl = true;
		ss << "\"" << tab.heading[col] << "\"";
	}
	ss << "]";
	
	
	ss << ",\"ele\":[";
	for(auto row = 0u; row < tab.nrow; row++){
		if(row > 0) ss << ",";
		ss << "[";
		bool fl = false;
		for(auto col : cols){
			if(fl == true) ss << ","; 
			fl = true;
			ss << "\"" << replace(tab.ele[row][col],"|",",") << "\"";
		}		
		ss << "]";
	}
	ss << "]}";
	
	return ss.str();
}


/// Converts data file to have a column for time instead of data (to allow it to be plotted)
void Data::make_table_with_time(const string file_in, const string file_out, const string col) const
{
	auto filefull = details.output_directory+"/"+file_out;
	Table tab = load_table(file_in,data_directory,true,true);
	
	auto date_col = find_column(tab,details.time_format_str);
	auto c = find_column(tab,col);	

	ofstream fout(filefull.c_str());
	if(!fout) emsgroot("Cannot open the file '"+filefull+"'");
	fout << "Time," << tab.heading[c] << endl;
	for(auto r = 0u; r < tab.nrow; r++){
		fout << details.gettime(tab.ele[r][date_col],"In file '"+tab.file+"'")-details.start;
		fout << "," << tab.ele[r][c] << endl;
	}
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
		emsgroot("File '"+tab.file+"' must have columns for "+to_string(nval-1)+" of the "+to_string(nval)+" categories in '"+democat[d].name+"'");
	}
	
	vector <double> row_sum(tab.nrow);
	vector < vector <double> > row_val; row_val.resize(tab.nrow);
	for(auto row = 0u; row < tab.nrow; row++){
		auto sum = 0.0; 
		row_val[row].resize(col_st.size());
		for(auto i = 0u; i < col_st.size(); i++){
			row_val[row][i] = get_double_positive(tab.ele[row][col_st[i]],"In file '"+tab.file+"'");
			if(row_val[row][i] < 0){
				emsgroot("In file '"+tab.file+"' cannot have a negative value in column '"+tab.heading[col_st[i]]);
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
		if(max > 100){
			nt = POPULATION;
			cout << "  The values are assumed to be population sizes (add a column for '" << democat[d].value[unset_st[0]] << "' to disambiguate)." << endl;
		}
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
				ss << ", which is more than the expected value " << target;
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
		dcc.area = create_area_sel(tabarea,dcc.geo_filt,"In 'geo_filt' for 'democat_change'");
		
		auto tab = load_table(dcc.file);
		
		auto d = 0u; while(d < ndemocat && democat[d].name != dcc.name) d++;
		if(d == ndemocat) emsgroot("In 'democat_change' the name '"+dcc.name+"' is not in 'democats'");
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
			emsgroot("Done");
		}
		
		if(times[0] > 0){
			cout << "  Demographic proportions assumed the same before '" << tab.ele[0][date_col] << "'" << endl;
		}
		if(tab.nrow < 2) emsgroot("File '"+dcc.file+"' must contain at least two rows of data.");
		
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
			emsgroot("Done");
		}
		
		auto dp_sel = create_dp_sel(dcc.democats_filt,"In 'democats_filt' for 'democat_change'");
		
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
			emsgroot("Done");
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

	if(false){ cout << nobs <<" Number of observations" << endl; emsgroot("Done");}
}


/// Generates a normal approximation to the observation model (used in approximate methods)
void Data::normal_approximation()
{
	for(auto &ob : obs){
		auto value = ob.value;
		
		switch(ob.obsmodel){
		case NORMAL_OBSMODEL: case NORMAL_PERCENT_OBSMODEL:
			{
				ob.var_approx = ob.sd*ob.sd;
			}
			break;
		
		case POISSON_OBSMODEL:
			{
				auto lam = value; if(lam < 1) lam = 1;
				ob.var_approx = lam;
				break;
			}
			
		case NEGBINO_OBSMODEL:
			{
				auto m = value; if(m < 1) m = 1;
				ob.var_approx = m + m*m/ob.shape;
			}
			break;
		}
		
		if(std::isinf(ob.var_approx)) emsgroot("Cannot have infinite variance for observation model");
		if(ob.var_approx == 0) emsgroot("Cannot have zero variance for observation model");
	}		
}


/// Loads a table giving time series data
void Data::load_timeseries_datatable(const Table &tabarea, const unsigned int i, const bool sim)
{
	auto &dt = datatable[i];
		
	Observation ob;
	ob.factor_spline = UNSET;
	ob.obsmodel = dt.obsmodel; ob.shape = dt.shape;
	ob.datatable = i;
	
	auto datafilter = get_datafilter(tabarea,dt.geo_dep,dt.democats_dep,dt.geo_filt,dt.democats_filt,dt.type,dt.observation,dt.file);

	vector <int> obs_times;
	
	bool table_loaded = false;
	Table tab;
	if(dt.optype == DATA){
		if(sim == false){
			table_loaded = true; 
			tab = load_table(dt.file);
			table_select_dates(dt.start,dt.end,dt.timestep,tab,dt.type,dt.shift,obs_times);
		}
		else{
			for(auto t = dt.start; t <= dt.end; t += dt.timestep) obs_times.push_back(t+dt.shift);
		}
	}

	vector <string> no_column;
	
	auto flag = false;
	for(const auto &df : datafilter){
		auto col=0u; if(table_loaded == true) col = find_column_noerror(tab,df.colname);
							
		if(col == UNSET){
			no_column.push_back(df.colname);
		}
		else{
			auto colSD = UNSET; 
			if(table_loaded == true && dt.load_sd == true) colSD = find_column(tab,df.colname+" SD");
	
			ob.dp_sel = df.dp_sel;
			
			flag = true;
			ob.area = df.area;
			ob.graph = graph.size();
			
			string fulldesc, fulldesc_data;
			
			Graph gr;
			if(dt.file != ""){ 
				gr.tab = details.analysis_type; 
				gr.tab2 = "Data Table";
				gr.tab3 = dt.file; 
				if(df.colname != "Data") gr.tab4 = df.colname; else gr.tab4 = "";
				
				if(details.siminf == SIMULATE){
					fulldesc = "Data table "+dt.file+": This shows the temporal variation in "+observation_description(dt.type,dt.observation);
					if(datafilter.size() > 1) fulldesc += " for "+df.colname;
					fulldesc += ".";
					fulldesc += "  The red line shows the underlying state";
					if(details.mode != SIM) fulldesc += " (with the dashed lines giving 95% of the stochastic variation across simulations)";
					fulldesc += ".";
					
					if(details.mode == SIM)	fulldesc += "  The black line shows data generated from this state (for example if data is measured weekly, this stepped curve would have a periodicity of 7 days). This data is saved in the file "+dt.file+" in the simulated data directory.";
					else fulldesc += "  The black line shows the raw data (from the input file "+dt.file+").";
				}
				else{
					fulldesc_data = "Data table "+dt.file+": This visualises the file "+dt.file+" which gives the temporal variation in "+observation_description(dt.type,dt.observation);
					if(datafilter.size() > 1) fulldesc_data += " for "+df.colname;
					fulldesc_data += ".";
					
					fulldesc = "Data table "+dt.file+": This shows the temporal variation in "+observation_description(dt.type,dt.observation);
					if(datafilter.size() > 1) fulldesc += " for "+df.colname;
					fulldesc += ".";
					fulldesc += "  The black line shows the raw data (from the input file "+dt.file+").";
					fulldesc += "  The red line (with dashed 95% credible intervals) shows the posterior distribution for the inferred underlying system state.";
					fulldesc += "  Under a suitable model, and with accurate inferece, it would be expected that the black and red curves exhibit the same temporal behaviour (barring variation coming from the observation model).";	
					fulldesc += " Large deviations between the two lines indicate that either the model is not able to properly account for the actual observed system dynamics, or that inference has not converged on the true posterior distribution.";
				}
			}
			else{ 
				gr.tab = details.analysis_type;
				gr.tab2 = "State Output";
				gr.tab3 = dt.plot_name; 
				if(df.colname != "Data") gr.tab4 = df.colname; else gr.tab4 = "";
				
				fulldesc = "State Output: This shows the plot *"+dt.plot_name+"*, as specifed in the input TOML file. ";
				if(details.siminf == SIMULATE){
					if(details.mode != SIM) fulldesc += "The dashed lines giving 95% of the stochastic variation across simulations.";
				}
				else{
					fulldesc += "The dashed lines give 95% credible intervals.";
				}
			}
			gr.fulldesc = fulldesc;
			gr.fulldesc_data = fulldesc_data;
			
			gr.type = GRAPH_TIMESERIES;
			gr.file = df.file;
			gr.name = df.name; if(dt.label != "" && datafilter.size() == 1) gr.name = dt.label;
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
							emsgroot("In 'data_tables' the file '"+dt.file+"' provides information which is out of the specified time range");
						}
												
						GraphPoint gp;
						gp.xi = ti/details.division_per_time;
						gp.xf = (ti+1)/details.division_per_time;
						gp.obs = obs.size();
						gr.point.push_back(gp);
						
						if(table_loaded == true) ob.value = get_data(tab.ele[row][col],"In file '"+df.file+"'",threshold_str,nodata_str); 
						else ob.value = UNSET;
						
						ob.sd = UNSET; 
						if(dt.load_sd == true && colSD != UNSET){							
							if(table_loaded == true) ob.sd = get_double(tab.ele[row][colSD],"In file '"+df.file+"'");
						}
						
						ob.sett_i = ti;
						ob.sett_f = ti+1;
					
						obs.push_back(ob);
					}
					break;
				
				case TRANS:
					{
						auto val = 0.0;
						auto ti = obs_times[0]*details.division_per_time;
				
						for(auto row = 0u; row < obs_times.size(); row++){
							if(table_loaded == true){
								auto v = get_data(tab.ele[row][col],"In file '"+df.file+"'",threshold_str,nodata_str); 
								if(details.trans_combine != UNSET){
									if(v == UNSET){
										emsgroot("In file '"+df.file+"' a value for 'trans_combine' cannot be set in conjunction with unknown data");
									}
									if(v == THRESH){
										emsgroot("In file '"+df.file+"' a value for 'trans_combine' cannot be set in conjunction with thresholded data");
									}
								}
								val += v;
							}
							else val = UNSET;
											
							ob.sd = UNSET; 
							if(dt.load_sd == true && colSD != UNSET){							
								if(details.trans_combine != UNSET) emsg("'loadSD' cannot be used in conjuction with 'trans_combine'"); 
								if(table_loaded == true) ob.sd = get_double(tab.ele[row][colSD],"In file '"+df.file+"'");
							}

							GraphPoint gp;
							unsigned int tf = (obs_times[row]+dt.timestep)*details.division_per_time;

							if(details.obs_section == true){
								for(auto t = ti+1; t < tf; t++){
									if(details.sec_define[t] == true) emsgroot("Observations cross particle filtering time points.");
								}
							}
							
							//load_datatable zz
							// TEMP turn off combine
							if(1 == 1 || sim == true || row == obs_times.size()-1 || details.trans_combine == UNSET  
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
	} 
	
	if(flag == false){
		string em = "In 'data_tables' the data in '"+dt.file+"' was not used. ";
		em += "Expected columns such as '" +no_column[0]+"' etc...";
		emsgroot(em);
	}
	
	if(no_column.size() > 0){
		string em = "WARNING! In the file '"+dt.file+"' there ";
		if(no_column.size() == 1) em += "was no column '"+no_column[0]+"'.";
		else{
			em += "were no columns :";
			for(auto i = 0u; i < no_column.size(); i++){
				if(i > 0) em += ", ";
				em += "'"+no_column[i]+"'";
			}
			em += ".";
		}
		cout << em << endl;
	}
}


/// Loads a table giving marginal data
void Data::load_marginal_datatable(const Table &tabarea, const unsigned int i, const bool sim)
{
	auto &dt = datatable[i];
	
	if(dt.geo_dep == "" && dt.democats_dep == ""){
		emsgroot("For marginal distributions in 'data_tables' either 'democats_dep' or 'geo_dep' must be set");
	}

	if(dt.geo_dep != "" && dt.democats_dep != ""){
		emsgroot("For marginal distributions in 'data_tables', 'democats_dep' and 'geo_dep' cannot both be set");
	}
		
	string dep;
	if(dt.democats_dep != "") dep = dt.democats_dep; else dep = dt.geo_dep;
	
	Observation ob;
	ob.factor_spline = UNSET;
	ob.datatable = i;
	ob.obsmodel = dt.obsmodel; ob.shape = dt.shape;
	ob.graph = graph.size();
	ob.factor = dt.factor;
	
	auto datafilter = get_datafilter(tabarea,dt.geo_dep,dt.democats_dep,dt.geo_filt,dt.democats_filt,dt.type,dt.observation,dt.file);

	if(datafilter.size() == 0){ 
		emsgroot("In 'data_tables' the marginal distribution for '"+dt.observation+"' does not contain any valid entries");
	}
	
	Table tab; if(sim == false) tab = load_table(dt.file);
	
	auto df = datafilter[0];
	Graph gr;
	gr.tab = details.analysis_type;
	gr.tab2 = "Data Table";
	gr.tab3 = dt.file; 
	gr.tab4 = "";
			
	string fulldesc = "Data table "+dt.file+": ";
	string fulldesc_data = "Data table "+dt.file+": ";
	if(details.siminf == SIMULATE){
		fulldesc += "This histogram shows the marginal distribution for "+observation_description(dt.type,dt.observation)+" between "+details.getdate(dt.start)+" and "+ details.getdate(dt.end)+". This is output in the data file *"+dt.file+"*.";
	}
	else{
		fulldesc_data += "This visualises the file "+dt.file+" which gives the marginal distribution for "+observation_description(dt.type,dt.observation)+" between "+details.getdate(dt.start)+" and "+ details.getdate(dt.end)+".";
		
		fulldesc += "This histogram shows the marginal distribution for "+observation_description(dt.type,dt.observation)+" between "+details.getdate(dt.start)+" and "+ details.getdate(dt.end)+".";

		fulldesc += "  The red bars represent the posterior mean with the error bars giving 95% credible intervals.";
		fulldesc += "  The horizontal black lines show the data from *"+dt.file+"*.";
	}
		
	gr.fulldesc = fulldesc;
	gr.fulldesc_data = fulldesc_data;
		
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
			if(row == tab.nrow) emsgroot("In 'data_tables' cannot find the value '"+df.colname+"' in the first column of '"+dt.file+"'");
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
	
		if(ti < 0 || tf > int(details.ndivision)) emsgroot("In 'data_tables' the file '"+dt.file+"' contains information which is outside the inference time range");
		
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
vector <DataFilter> Data::get_datafilter(const Table &tabarea, const string geo_dep, const string democats_dep, const string geo_filt, const string democats_filt, const DataType type, const string observation, const string datafile) const
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
		auto filt = democats_filt;
		if(democats_filt.substr(0,7) == "age:age") filt = filt.substr(7); // TO DO
		file += "-"+replace(filt,":","=");
		name += " "+filt;
		desc += " "+filt;
	}

	vector <DataFilter> datafilter;
	
	if(democats_dep != ""){
		auto dc = 0u; while(dc < democat.size() && toLower(democat[dc].name) != toLower(democats_dep)) dc++;
		if(dc == democat.size()){
			emsgroot("For file '"+datafile+"', the expression 'democats_dep="+democats_dep+"' is not a value in 'democats'");
		}
		
		if(democat[dc].value.size() == 1){
			emsgroot("For file '"+datafile+"', the expression 'democats_dep=\""+democats_dep+"\"' can only be used if '"+democats_dep+"' has more than one value.");
		}
		
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
			df.area = create_area_sel(tabarea,geo_filt,"In 'geo_filt' for '"+datafile+"'");
			auto state = democat[dc].name+":"+val; if(democats_filt != "") state += ","+democats_filt;
			df.dp_sel = create_dp_sel(state,"In 'democats_filt' for '"+datafile+"'");
			
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
				df.area = create_area_sel(tabarea,state,"In 'geo_filt' for '"+file+"'");
				df.dp_sel = create_dp_sel(democats_filt,"In 'democats_filt' for '"+datafile+"'");
				if(df.area.size() > 0 && df.dp_sel.size() > 0) datafilter.push_back(df);		
			}		
		}
		else{
			DataFilter df; 
			df.file = file+".csv";
			df.colname = "Data";
			df.name = name;
			df.desc = desc;
			df.area = create_area_sel(tabarea,geo_filt,"In 'geo_filt' for '"+datafile+"'");
			df.dp_sel = create_dp_sel(democats_filt,"In 'democats_filt' for '"+datafile+"'");
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
		emsgroot(ss.str());
	}
	
	if(false){
		cout << endl;
		cout << "geo_dep: " << geo_dep << "   democats_dep: " << democats_dep <<  "      geo_filt: " <<  geo_filt << "   democats: " << democats_filt << endl;
		for(auto  df : datafilter){
			cout << df.name << "  Datafilter: ";
			cout << "Areas: "; for(auto c : df.area) cout << area[c].code << ","; 
			cout << "  ";
			cout << "DP: "; 
			for(auto dp : df.dp_sel){ 
				for(auto i = 0u; i < ndemocat; i++) cout << democat[i].value[democatpos[dp][i]] << ",";
				cout << "  ";
			}
			cout << endl;
		}
		emsgroot("Done");
	}
	
	return datafilter;
}


/// Sets the standard deviations used for the observations
void Data::set_obs_sd()
{
	auto ngraph = 0u;
	for(auto &ob : obs){ if(ob.graph > ngraph) ngraph = ob.graph;}
	ngraph++;
	
	vector <double> graph_value_max(ngraph);
	for(auto gr = 0u; gr < ngraph; gr++){
		graph_value_max[gr] = 0;
	}
	
	for(auto &ob : obs){		
		auto val = ob.value;
		auto gr = ob.graph;
		if(val != UNKNOWN && val != THRESH){
			if(val > graph_value_max[gr]) graph_value_max[gr] = val;
		}
	}
	
	for(auto &ob : obs){		
		if(ob.sd == UNSET){
			auto dt = ob.datatable;
			auto gr = ob.graph;
		
			if(ob.obsmodel == NORMAL_OBSMODEL){
				ob.sd = datatable[dt].sd;
			}
			
			if(ob.obsmodel == NORMAL_PERCENT_OBSMODEL){
				auto ef = datatable[dt].epsilon_factor;
				//cout << ef << "ef\n";
				ob.sd = (ob.value*(1-ef) + graph_value_max[gr]*ef)*0.01*datatable[dt].percent;
			}
		}			
	}
	
	/*
	vector <unsigned int> num_obs(datatable.size());
	vector <unsigned int> num_obs_cut(datatable.size());
	vector <double> value_max(datatable.size());

	auto ngraph = 0u;
	for(auto &ob : obs){ if(ob.graph > ngraph) ngraph = ob.graph;}
	ngraph++;

	vector <double> graph_value_max(ngraph);
	vector <double> sd_max(ngraph);
	vector <unsigned int> num_obs_graph_cut(ngraph);
	for(auto gr = 0u; gr < ngraph; gr++){
		graph_value_max[gr] = 0;
		sd_max[gr] = 0;
		num_obs_graph_cut[gr] = 0;
	}
		
	for(auto dt = 0u; dt < datatable.size(); dt++){		
		num_obs[dt] = 0;
		value_max[dt] = 0;
		num_obs_cut[dt] = 0;
	}
	
	for(auto &ob : obs){		
		auto dt = ob.datatable;
		auto gr = ob.graph;
		num_obs[dt]++;
		auto val = ob.value;
		if(val != UNKNOWN && val != THRESH){
			if(val > value_max[dt]) value_max[dt] = val;
			if(val > graph_value_max[gr]) graph_value_max[gr] = val;
			
			if(ob.sd > sd_max[gr]) sd_max[gr] = ob.sd;
		}
	}
	
	for(auto &ob : obs){	
		auto dt = ob.datatable;	
		auto gr = ob.graph;
		auto val = ob.value;
		if(val != UNKNOWN && val != THRESH){
			if(val > 0.1*graph_value_max[gr]){
				num_obs_cut[dt]++;
				num_obs_graph_cut[gr]++;
			}
		}
	}
	
	
	
	if(false){  // Provides an reassignment of weights based on matching different data sources
		// Accounts for many observations on some graphs compared to others
		vector <double> graph_weightsum(ngraph);		
		for(auto gr = 0u; gr < ngraph; gr++) graph_weightsum[gr] = 0;
		for(auto &ob : obs) graph_weightsum[ob.graph] += ob.w;	
		for(auto &ob : obs) ob.w /= graph_weightsum[ob.graph];
	
		// Accounts for the fact that some graphs within a datatable have fewer point (combining transition events)
		vector <double> ngraph_point(ngraph);
		for(auto gr = 0u; gr < ngraph; gr++) ngraph_point[gr] = 0;
		for(auto &ob : obs) ngraph_point[ob.graph]++;
		
		
		auto ndatatable = datatable.size();
		vector <double> ndatatable_point(ndatatable);
		for(auto dt = 0u; dt < ndatatable; dt++) ndatatable_point[dt] = 0;
		
		for(auto &ob : obs){
			auto dt = ob.datatable;	
			auto gr = ob.graph;
			if(ngraph_point[gr] > ndatatable_point[dt]) ndatatable_point[dt] = ngraph_point[gr];
		}
		
		for(auto &ob : obs){
			auto dt = ob.datatable;	
			auto gr = ob.graph;
			auto fac = sqrt(ngraph_point[gr]/ndatatable_point[dt]);
			ob.w *= fac;
		}
	}
	
	if(false){
		for(auto &ob : obs){	
			auto dt = ob.datatable;
			cout << ob.value << " " << ob.sd << " " << datatable[dt].file << " weight" << endl;
		}
	}
	*/
}


/// Creates all the demographic categories consistent with a string 
vector <unsigned int> Data::create_dp_sel(const string dp_str, const string em) const
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
			if(valsp.size() != 2) emsgroot(em+" the expression '"+dp_str+"' must specify the demographic category and its value");
			
			auto d = 0u; while(d < ndemocat && democat[d].name != valsp[0]) d++;
			if(d == ndemocat) emsgroot(em+" cannot find the value '"+valsp[0]+"' in 'ages' or 'democats'");
			
			auto val_list = split(valsp[1],'|');
			
			for(auto i = 0u; i < val_list.size(); i++){                 // Checks for repeated values
				for(auto j = i+1; j < val_list.size(); j++){
					if(val_list[i] == val_list[j]) emsgroot(em+" the value '"+val_list[i]+"' must not be repeated.");
				}
			}
			
			for(auto i = 0u; i < val_list.size(); i++){ 
				auto pos = 0u; while(pos < democat[d].value.size() && democat[d].value[pos] != val_list[i]) pos++;
				if(pos == democat[d].value.size()){
					emsgroot(em+" cannot find the value '"+val_list[i]+"' in 'name=\""+valsp[0]+"\"'");
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
				for(auto i = 0u; i < democat[d].value.size(); i++) cout << flag[d][i] << ","; 
				cout << "flag" << endl;
			}
		}
	}

	return dp_sel;
}
 
 
/// Creates all the areas consistent with a string 
vector <unsigned int> Data::create_area_sel(const Table &tabarea, const string str, const string em) const
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
			if(geosp.size() != 2) emsgroot(em+" there was a problem with the expression '"+str+"'");
		
			auto geomap = create_geomap(tabarea,geosp[0]);
			auto reg_list = split(geosp[1],'|');
			
			for(auto i = 0u; i < reg_list.size(); i++){                  // Checks for repeated values
				for(auto j = i+1; j < reg_list.size(); j++){
					if(reg_list[i] == reg_list[j]) emsgroot(em+" the value '"+reg_list[i]+"' must not be repeated.");
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
				if(fl == false)  emsgroot(em+" the value '"+reg_str+"' is not found in '"+geosp[0]+"'."); 
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
	emsgroot("load_tablefromdatapipeline for '"+file+"' cannot be called as data pipeline is not compiled in");
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
		if(!in) emsgroot("Cannot open the file '"+dir+"/"+file+"'.");
	}
	else{
    used_file = details.output_directory+"/Simulated_data/"+file;
		in.open(used_file.c_str());
		if(!in){
      used_file = data_directory+"/"+file;
			in.open(used_file.c_str());
			if(!in) emsgroot("Cannot open the file '"+used_file+"'");
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
		
		tab.heading = split(line,sep);
	
		tab.ncol = tab.heading.size();
	}
	else tab.ncol = 0;
		
	do{
		getline(in,line);
		if(in.eof()) break;

		auto vec = split(line,sep);
	
		if(tab.ncol == 0) tab.ncol = vec.size();
		else{
			if(vec.size() != tab.ncol) emsgroot("Rows in the file '"+file+"' do not all share the same number of columns.");
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
void Data::table_select_dates(int &start, int &end, unsigned int &timestep, Table &tab, const DataType type, const unsigned int shift, vector <int> &times) const 
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
		emsgroot("In data_tables' the value for 'start' does not agree with file '"+tab.file+"'");
	}
	start = times[row_start];
	
	int row_end = tab.nrow-1; 
	while(row_end >= 0 && (times[row_end] >= (int)details.period || (end != UNSET && times[row_end] >= int(end)))) row_end--;
	
	if(row_end == -1) emsgroot("In file '"+tab.file+"' there is no data in the specified time period");
	if(end != UNSET && end != times[row_end]){
		emsgroot("In data_tables' the value for 'end' does not agree with file '"+tab.file+"'");
	}
	end = times[row_end];
	
	if(type != TRANS) timestep = 1;
	else{
		if(row_start == row_end){
			if(timestep == UNSET) emsgroot("In file '"+tab.file+"' with only one row of data a value for 'timestep' must be set.");
		}
		else{
			auto tstep = UNSET;
			for(auto row = row_start; row < row_end; row++){
				if(times[row+1] < times[row]) emsgroot("In file '"+tab.file+"' the measurements must be time ordered.");
				if(times[row+1] == times[row]) emsgroot("In file '"+tab.file+"' each measurement must be at a different time.");
				
				unsigned int ts = times[row+1]-times[row]; 
				if(tstep == UNSET) tstep = ts;
				else{ if(ts != tstep) emsgroot("In file '"+tab.file+"' measurements should have an equal timestep.");}
			}
			
			if(timestep == UNSET) timestep = tstep;
			else{
				if(timestep != tstep){
					emsgroot("In 'data_tables' the value for 'timestep' does not agree with the timestep in file '"+tab.file+"'.");
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
		emsgroot("Done");
	}
}				


/// Finds a column in a table
unsigned int Data::find_column(const Table &tab, const string name) const
{
	auto c = find_column_noerror(tab,name); 
	if(c == UNSET) emsgroot("Cannot find the column heading '"+name+"' in file '"+tab.file+"'.");
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
			
			
/// Generates a file by selected a certain set of colums from a file
void Data::generate_file_from_table(const string file_data, const string file, const vector <string> &cols) const
{
	Table tab = load_table(file,"",true,true); 
	
	auto ncol = cols.size();
	
	string sep;
	if(stringhasending(file,".txt")) sep = "\t"; else sep = ",";
		
	vector <unsigned int> co;
	for(auto c = 0u; c < ncol; c++) co.push_back(find_column(tab,cols[c])); 
	
	auto filefull = details.output_directory+"/"+file_data;
	ofstream fout(filefull);
	if(!fout) emsgroot("Cannot open the file '"+file+"'");
	
	for(auto c = 0u; c < ncol; c++){ fout << tab.heading[co[c]]; if(c < ncol-1) fout << sep;}
	fout << endl;
	
	for(auto r = 0u; r < tab.nrow; r++){
		for(auto c = 0u; c < ncol; c++){ fout << tab.ele[r][co[c]]; if(c < ncol-1) fout << sep;}
		fout << endl;
	}
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
			cf.dp_sel = create_dp_sel(cf.democats_filt,"In 'democats_filt' for 'modification'");
			cf.area = create_area_sel(tabarea,cf.geo_filt,"In 'geo_filt' for 'modification'");  
		}
	}
}


/// Converst from an obervation to a description
string Data::observation_description(const DataType type, const string obs, const unsigned int timestep) const
{
	stringstream ss;
	auto spl = split(obs,',');
	
	switch(type){
		case TRANS:
			ss << "the number of individual transitions ";
			if(timestep == 1) ss << "per unit time";
			else ss << "per " << timestep << " time units";
			ss <<	" moving from the ";
			for(auto i = 0u; i < spl.size(); i++){
				auto st = spl[i];
				if(i > 0) ss << " plus from the ";
				
				auto j = 0u; 
				bool fl = false;
				if(st.length() < 2) fl = true;
				else{				
					while(j < st.length()-2 && st.substr(j,2) != "->") j++;
					if(j == st.length()-2) fl = true;
				}
				
				if(fl == true){
					emsgroot("In 'data_tables' the transition observation '"+st+"' does not contain '->'");
				}
				
				ss << st.substr(0,j) << " to the " <<  st.substr(j+2,st.length()-(j+2)) << " compartment";
			}			
			break;
		
		case POP:	case POPFRAC:
			if(spl.size() == 1){
				ss << "the population";
				if(type == POPFRAC) ss << " fraction";
				ss <<" in the " << spl[0] << " compartment";
			}
			else{
				ss << "the sum of the population";
				if(type == POPFRAC) ss << " fraction";
				ss << "s in the ";
				for(auto i = 0u; i < spl.size(); i++){
					if(i > 0){
						if(i < spl.size()-1) ss << ", "; else ss << " and ";
					}
					ss << spl[i];
				}
				ss << " compartments";
			}
			break;
			
		case MARGINAL:
			ss << "the number of individual transitions from the ";
			for(auto i = 0u; i < spl.size(); i++){
				auto st = spl[i];
				if(i > 0) ss << " plus from the ";
				
				bool fl = false;
				if(st.length() < 2) fl = true;
				else{
					auto j = 0u; while(j < st.length()-2 && st.substr(j,2) != "->") j++;
					if(j < st.length()-2){
						ss << st.substr(0,j) << " to the " <<  st.substr(j+2,st.length()-(j+2)) << " compartment";
					}
					else fl = true;
				}
				
				if(fl == true){
					emsgroot("In 'data_tables' the marginal observation '"+st+"' does not contain '->'");
				}
			}			
			break;
	}
	
	return ss.str();
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
		cout << "factor: " << ob.factor << "  ";
		cout << "area: "; for(auto c : ob.area) cout <<c << ",  ";
		cout << "dp_sel: "; for(auto dp : ob.dp_sel) cout << dp << ",  ";
		cout << endl << endl;
	}
}
