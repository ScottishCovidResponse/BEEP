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

const string filter = "";                                       // Used to filter certain regions
//const string filter = "northing50";                                       // Used to filter certain regions
//const string filter = "scotland";  
//const string filter = "notscotland";
//const string filter = "england";  

#ifdef USE_Data_PIPELINE
#include "datapipeline.hh"
#include "table.hh"
#endif

/// Initialises data
Data::Data(const Inputs &inputs, const Details &details, Mpi &mpi, DataPipeline *dp) :
	datapipeline(dp), data_directory(inputs.find_string("datadir","UNSET")), details(details)
{
	// The data directory
	if(data_directory == "UNSET") emsgroot("The 'data_directory' must be set.");

	threshold = inputs.find_integer("threshold",UNSET);                                       // The threshold (if specified)
	
	democat = inputs.find_demographic_category();
	ndemocat = democat.size();	
	nage = democat[0].value.size();
	calc_democatpos();

	covar = inputs.find_covariate();
	ncovar = covar.size();
	
	inputs.find_genQ(genQ,details,democat[0].value);
	
	read_data_files(inputs,mpi);                                  // Reads the data files
}


/// Based on the different demographic categories, this calculates all the possible combinations
void Data::calc_democatpos()
{
	vector <unsigned int> count(ndemocat);                        // Defines all the demographic states
	for(auto& co : count) co = 0;
	
	int dc;
	do{
		democatpos.push_back(count);
		
		auto fl=0u;
		dc = ndemocat-1;
		do{
			fl = 0;
			count[dc]++; if(count[dc] == democat[dc].value.size()){ fl = 1; count[dc] = 0; dc--;}
		}while(fl == 1 && dc >= 0);
	}while(dc >= 0);
	ndemocatpos = democatpos.size();

	ndemocatposperage = ndemocatpos/nage;
}


/// Reads in transition and area data
void Data::read_data_files(const Inputs &inputs, Mpi &mpi)
{
	datatable = inputs.find_datatable(details);
	if(datatable.size() == 0) emsgroot("'datatables' must be set.");	

	if(mpi.core == 0){
		check_datatable();

		string file = inputs.find_string("areas","UNSET");
		if(file == "UNSET") emsgroot("A 'areas' file must be specified");
			
		Table tab = load_table(file);
	
		read_areas(tab,file);
		
		load_region_effect(inputs,tab);

		load_datatable(tab);	
	
		load_counterfactual(inputs,tab);
		
		generate_matrices();
		
		coarse_grain(tab,inputs.find_string("geography","UNSET"));
	}
	if(mpi.ncore > 1) mpi.copy_data(narea,area,region_effect,nobs,obs,genQ,counterfact);

	if(details.siminf == INFERENCE) set_datatable_weights();

	agedist.resize(nage); for(auto& aged : agedist) aged = 0;
	democatpos_dist.resize(ndemocatpos); for(auto& dpd : democatpos_dist) dpd = 0;
	democat_dist.resize(ndemocat); 
	for(auto d = 0u; d < ndemocat; d++){
		democat_dist[d].resize(democat[d].value.size()); for(auto& dvd : democat_dist[d]) dvd = 0;
	}
		
	popsize = 0;
	for(auto c = 0u; c < narea; c++){                                    // Adds individuals to the system
		area[c].ind.resize(ndemocatpos);
		for(auto dp = 0u; dp < ndemocatpos; dp++){
			auto imax = area[c].pop[dp];
			if(modeltype == IND_MODEL){
				for(auto i = 0u; i < imax; i++){
					area[c].ind[dp].push_back(ind.size());
				
					Individual indi;
					indi.area = c;
					indi.dp = dp;
					ind.push_back(indi);
				}
			}
			democatpos_dist[dp] += imax;
			
			auto a = democatpos[dp][0];
			agedist[a] += imax;
			
			for(auto d = 0u; d < ndemocat; d++){
				democat_dist[d][democatpos[dp][d]] += imax; 
			}
			
			popsize += imax;
		}
	}
	
	if(false){
		for(auto& aged : agedist) cout << aged << "\n"; cout << "agedist\n"; emsg("O");
	}
	
	for(auto& aged : agedist) aged /= popsize;
	for(auto& dpd : democatpos_dist) dpd /= popsize;
	for(auto d = 0u; d < ndemocat; d++){
		for(auto& dvd : democat_dist[d]) dvd /= popsize;
	}
	
	if(false){
		for(auto d = 0u; d < ndemocat; d++){
			for(auto i = 0u; i < democat_dist[d].size(); i++) cout << d << " " << i << " " << democat_dist[d][i] << ",\n";
		}
		
		for(auto a = 0u; a < agedist.size(); a++) cout << agedist[a] << ","; cout << "h\n";
		emsg("P");	
	}
	
	narage = narea*nage;
	nardp = narea*ndemocatpos; 
	nsettardp = details.ndivision*nardp;
	
	raw();  // Does raw data analysis
}


/// Checks for any errors in the definition of the datatables
void Data::check_datatable()
{
	if(details.mode == SIM){ // If siulting checks that all the file names are different
		for(auto dt = 0u; dt < datatable.size(); dt++){
			for(auto dt2 = 0u; dt2 < dt; dt2++){
				if(datatable[dt].file == datatable[dt2].file){
					stringstream ss; ss << "The file name '" << datatable[dt].file << "' is used more than once.";
					emsg(ss.str());
				}				
			}
		}
	}		
}		


/// Reads the file which gives information about areas
void Data::read_areas(Table &tab, string file)
{
	// If only one age group then combines all columns with "age" to generate an "all" column
	
	age_from.resize(nage); age_to.resize(nage);
	for(auto a = 0u; a < nage; a++){
		string name = democat[0].value[a];
		
		auto col = 0u; while(col < tab.ncol && tab.heading[col] != name) col++;

		if(col == tab.ncol){		
			if(name == "all"){ age_from[a] = 0; age_to[a] = 100;} 
			else{
				auto len = name.length();
				if(name.substr(len-1,1)== "+"){ age_from[a] = atoi(name.substr(0,len-1).c_str()); age_to[a] = 100;} 
				else{
					auto st = split(name,'-');
					if(st.size() != 2) emsg("Age problem");
					age_from[a] = atoi(st[0].c_str());
					age_to[a] = atoi(st[1].c_str());
				}
			}
			table_add_age(name,age_from[a],age_to[a],tab);
		}
	}
	
	if(false){ // Multiplies all the populations by a factor
		auto fac = 13u;
		for(auto a = 0u; a < nage; a++){
			auto col = find_column(tab,democat[0].value[a]);
			for(auto row = 0u; row < tab.nrow; row++){
				stringstream ss; ss << atoi(tab.ele[row][col].c_str())*fac;
				tab.ele[row][col] = ss.str();
			}
		}
	}
	
	auto codecol = find_column(tab,"area");            // Works out columns for different quantities
	//auto xcol = find_column(tab,"easting");
	//auto ycol = find_column(tab,"northing");
	
	vector <unsigned int> cov_col(ncovar); 
	for(auto cov = 0u; cov < ncovar; cov++) cov_col[cov] = find_column(tab,covar[cov].name);

	vector < vector <unsigned int> > democat_col(ndemocat); 
	for(auto democ = 0u; democ < ndemocat; democ++){
		democat_col[democ].resize(democat[democ].value.size());
		for(auto k = 0u; k < democat[democ].value.size(); k++){
			democat_col[democ][k] = find_column(tab,democat[democ].value[k]);  
		}
	}

	filter_areas(tab);
	
	for(const auto& trow : tab.ele){
		Area are;
		are.code = trow[codecol];
		
		are.x = ran();
		//if(trow[xcol] == "NA") emsg("Easting");
		//else are.x = atof(trow[xcol].c_str());
		//if(std::isnan(are.x)) emsg("Problem with eastings");
		
		are.y = ran();
		//if(trow[ycol] == "NA") emsg("Northing");
		//else are.y = atof(trow[ycol].c_str());
		//if(std::isnan(are.y)) emsg("Problem with northings");
		
		are.region_effect = UNSET;
			
		are.covar.resize(ncovar);
		for(auto j = 0u; j < ncovar; j++){
			auto st = trow[cov_col[j]];
			auto v = atof(st.c_str());
			if(std::isnan(are.covar[j])) emsg("In file '"+file+"' the expression '"+st+"' is not a number");	
		
			if(covar[j].func == "log"){
				if(v == 0) v = 0.01;
				if(v <= 0) emsg("In file '"+file+"' the log transformed quantities from 'covar' must be positive.");
				are.covar[j] = log(v);
			}
			else{
				if(covar[j].func == "linear") are.covar[j] = v;
				else emsg("n file '"+file+"' the functional relationship '"+covar[j].func+"' is not recognised.");
			}
		}
		
		vector <vector <double> > val;
		val.resize(democat.size());
		for(auto d = 0u; d < democat.size(); d++){
			auto jmax = democat[d].value.size();
			val[d].resize(jmax);
			for(auto j = 0u; j < jmax; j++){
				auto st = trow[democat_col[d][j]];
				val[d][j] = atof(st.c_str());
				if(std::isnan(val[d][j])) emsg("In file '"+file+"' the expression '"+st+"' is not a number");	
			}
		}
		
		for(auto d = 1u; d < democat.size(); d++){  // Normalises non-age categories 
			auto sum = 0.0;
			for(auto j = 0u; j < democat[d].value.size(); j++) sum += val[d][j];
			for(auto j = 0u; j < democat[d].value.size(); j++) val[d][j] /= sum;
		}

		are.pop.resize(ndemocatpos);
		for(auto dp = 0u; dp < ndemocatpos; dp++){
			auto v = 0.0;
			for(auto j = 0u; j < democat.size(); j++){
				if(j == 0) v = val[j][democatpos[dp][j]];
				else v *= val[j][democatpos[dp][j]];
			}
			are.pop[dp] = (unsigned int)(v+0.5);
		}

		area.push_back(are);			
	}
	
	narea = area.size();

	for(auto j = 0u; j < ncovar; j++){            // Shifts covariates so average is zero
		auto sum = 0.0; for(auto c = 0u; c < narea; c++) sum += area[c].covar[j];
		sum /= narea;
	
		for(auto c = 0u; c < narea; c++) area[c].covar[j] -= sum;
	}		
	
	if(checkon == true){
		for(const auto& are : area){
			cout << are.code << are.x << " " <<  are.y << "  ***";
			for(const auto& pop : are.pop) cout << pop << ", ";
			cout << endl;	
		}
	}	
}


/// Potentially reduces the number of areas via a filter
void Data::filter_areas(Table &tab)
{
	auto codecol = find_column(tab,"area"); 
		
	genQ.area_filter_ref.resize(tab.nrow);
	if(filter == ""){
		for(auto c = 0u; c < tab.nrow; c++) genQ.area_filter_ref[c] = c;
	}
	else{             // Generates a mapping between raw and filters so that M.txt can be loaded
		auto narea_old = tab.nrow; 
		vector <string> code_store(narea_old);
		for(auto c = 0u; c < narea_old; c++){ 
			code_store[c] = tab.ele[c][codecol]; 
			genQ.area_filter_ref[c] = UNSET;
		}
	
		filter_table(filter,tab); 
		if(tab.nrow == 0) emsg("All areas have been filtered out");
		
		for(auto c = 0u; c < tab.nrow; c++){ 
			auto cc = 0u; while(cc < narea_old && code_store[cc] != tab.ele[c][codecol]) cc++;
			if(cc == narea_old) emsgEC("Data",65);
			genQ.area_filter_ref[cc] = c;
		}
	}
}


/// Deletes rows from a table based on a filter
void Data::filter_table(const string st, Table &tab) const
{		
	if(st.substr(0,8) == "northing" || st == "scotland" || st == "notscotland" || st == "england"){
		if(st == "scotland" || st == "notscotland" || st == "england"){
			if(st == "england"){
				for(auto row = 0u; row < tab.nrow; row++){
					if(tab.ele[row][0].substr(0,1) != "E") tab.ele[row][0] = "";
				}
			}
			
			if(st == "scotland"){
				for(auto row = 0u; row < tab.nrow; row++){
					if(tab.ele[row][0].substr(0,1) != "S") tab.ele[row][0] = "";
				}
			}
			
			if(st == "notscotland"){
				for(auto row = 0u; row < tab.nrow; row++){
					if(tab.ele[row][0].substr(0,1) == "S") tab.ele[row][0] = "";
				}
			}
		}
		else{
			unsigned int num = atoi(st.substr(8,st.length()-8).c_str());
		
			Table tab2 = load_table("areadata.txt");
			auto colno = find_column(tab2,"northing");
		
			vector <double> northing;
			for(auto row = 0u; row < tab2.nrow; row++){
				if(tab2.ele[row][0].substr(0,1) == "S") tab2.ele[row][colno] = to_string(-LARGE);
				
				northing.push_back(atof(tab2.ele[row][colno].c_str()));
			}
			
			sort(northing.begin(),northing.end());
			auto cut = (northing[northing.size()-num] + northing[northing.size()-num-1])/2;
			
			for(auto row = 0u; row < tab2.nrow; row++){
				if(atof(tab2.ele[row][colno].c_str()) < cut) tab.ele[row][0] = "";
			}
		}
		
		remove_empty_rows(tab);
		
		if(false){ for(auto row = 0u; row < tab.nrow; row++) cout << tab.ele[row][0] << " keep" << endl;}
	}
	else emsg("filter prob");
	
	if(false){
		for(auto r = 0u; r < tab.nrow; r++){
			for(auto c = 0u; c < tab.ncol; c++){
				cout << tab.ele[r][c] << " ";
			}
			cout << endl;
		}
	}
}


/// Removes any empty rows
void Data::remove_empty_rows(Table& tab) const
{
	auto row = 0u;
	while(row < tab.nrow){
		if(tab.ele[row][0] != ""){
			row++;
		}
		else{
			tab.ele.erase(tab.ele.begin()+row);
			tab.nrow--;
		}
	}
}


/// Loads datatables
void Data::load_datatable(const Table &tabarea)
{
	bool sim;
	if(details.mode == SIM) sim = true; else sim = false;     

	for(auto i = 0u; i < datatable.size(); i++){        // Loads data tables for inference
		switch(datatable[i].type){
			case POP: case TRANS: load_timeseries_datatable(tabarea,i,sim); break;	
			case MARGINAL: load_marginal_datatable(tabarea,i,sim); break;	
		}
	}
	
	nobs = obs.size();
	if(false){ cout << nobs <<" Number of observations\n"; emsg("Number of observations");}
}


/// Loads a table giving time series data
void Data::load_timeseries_datatable(const Table &tabarea, unsigned int i, bool sim)
{
	auto& dt = datatable[i];
		
	Observation ob;
	ob.fraction_spline_ref = UNSET;
	ob.w = UNSET;
	ob.logdif_offset = UNSET;
	ob.obsmodel = dt.obsmodel; ob.shape = dt.shape; ob.invT = dt.invT;
	
	ob.datatable = i;

	auto file = dt.file; 
	Table tab;
	if(sim == false && dt.optype == DATA){
		tab = load_table(file);
		
		table_select_dates(dt.start,dt.units,tab,"trans",dt.shift);
	}
	
	auto datafilter = get_datafilter(tabarea,dt.geo,dt.democats);
	
	string filename, name, desc;
	switch(dt.type){
		case POP: desc = "Population"; break;
		case TRANS: desc += "Transitions";  break;
		case MARGINAL: emsgEC("Data",66); break;
	}
	filename += dt.observation;
	name += dt.observation;
	desc += " in "+dt.observation;
	
	if(dt.geo != "all"){
		filename += "-"+dt.geo;
		name += " "+dt.geo;
		desc += " in "+dt.geo+": ";
	}
	if(dt.democats != "all"){
		desc += " in "+dt.democats+": ";
	}
	
	auto nrow = 0u;
	if(dt.optype == DATA){
		if(sim == false) nrow = tab.nrow;
		else nrow = (unsigned int)((details.period - dt.start)/dt.units);
	}
	
	auto flag = false;
	for(const auto& df : datafilter){
		auto col=0u;
		if(sim == false && dt.optype == DATA){
			for(col = 0; col < tab.ncol; col++) if(tab.heading[col] == df.name) break;
			if(col == tab.ncol) col = UNSET;
		}
		
		if(col != UNSET){
			ob.dp_sel = df.dp_sel;
			
			flag = true;
			
			if(dt.optype == DATA && dt.start + (nrow-1)*dt.units > details.period){
				emsg("The file '"+file+"' has more data than will fit in the defined 'start' and 'end' time period.");
			}
			
			ob.area = df.area;
			ob.graph = graph.size();
			
			Graph gr;
			gr.type = GRAPH_TIMESERIES;
			
			gr.filename = filename;
			gr.name = name;
			gr.desc = desc;
			if(df.name != "all"){
				gr.filename += "-"+df.name;
				gr.name += " "+df.name;
				gr.desc += " "+df.name;
			}
			gr.colname = df.name;
			gr.fraction_spline_ref = UNSET;
			gr.datatable = i;
			gr.dp_sel = ob.dp_sel;
			gr.area = ob.area;
			
			switch(dt.type){
				case POP: 				
					for(auto row = 0u; row < nrow; row++){
						int ti = (dt.start+row*dt.units)*details.division_per_time;
						if(ti < 0 || ti > int(details.ndivision-1)) emsg("Data table '"+dt.file+"' is out of the inference time range");
												
						GraphPoint gp;
						gp.xi = ti/details.division_per_time;
						gp.xf = (ti+1)/details.division_per_time;
						gp.obs.push_back(obs.size());
						gr.point.push_back(gp);
						
						if(sim == true) ob.value = UNSET;
						else ob.value = get_integer(tab.ele[row][col],file);
						ob.sett_i = ti;
						ob.sett_f = ti+1;
					
						obs.push_back(ob);
					}
				break;
				
				case TRANS:
					{
						auto val = 0u;
						auto ti = dt.start;
				
						for(auto row = 0u; row < nrow; row++){
							if(sim == true) val = UNSET;
							else val += get_integer(tab.ele[row][col],file);
											
							GraphPoint gp;
							unsigned int tf = (dt.start+(row+1)*dt.units)*details.division_per_time;

							if(details.obs_combine_type == COMBINE_SECTION){
								for(auto t = ti+1; t < tf; t++){
									if(details.sec_define[t] == true) emsg("Observations cross particle filtering time points.");
								}
							}
							
							if(sim == true || row == nrow-1 || (details.obs_combine_type == COMBINE_VALUE && val >= VALUE_LIMIT) 
					    //   || (details.obs_combine_type == COMBINE_SECTION && details.sec_define[tf] == true)){
					      || (details.obs_combine_type == COMBINE_SECTION && (val >= VALUE_LIMIT || details.sec_define[tf] == true))){
								gp.xi = ti/details.division_per_time;
								gp.xf = tf/details.division_per_time;
							
								gp.obs.push_back(obs.size());
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
				
				case MARGINAL:
					emsgEC("Data",57);
					break;
			}
			
			dt.graph_ref.push_back(graph.size());
			graph.push_back(gr);
		}
	} 
	if(flag == false) emsg("The data in '"+file+"' was not used");
}


/// Based on expressions for geo and democats in the data file this generates information for all the data columns
vector <DataFilter> Data::get_datafilter(const Table &tabarea, string geo, string democats)
{ 
	vector < vector<unsigned int> > dp_sel_list;    // Finds possible demographic possibilities                
	vector <string> colname;

	if(democats == "all"){
		colname.push_back("all");
		dp_sel_list.push_back(create_dp_sel("all"));
	}
	else{
		auto dem = split(democats,',');
		for(auto i = 0u; i < dem.size(); i++){
			auto dsp = split(dem[i],':');
			if(dsp.size() == 1){
				auto dc = 0u; while(dc < democat.size() && democat[dc].name != dsp[0]) dc++;
				if(dc == democat.size()) emsg("Cannot find '"+dsp[0]+"'");
				
				colname.clear();
				for(const auto &val : democat[dc].value){
					string st = "";
					for(auto ii = 0u; ii < i; ii++){ if(st != "") st += ","; st += dem[i];}
					if(st != "") st += ","; st += democat[dc].name+":"+val;
					for(auto ii = i+1; ii < dem.size(); ii++){ if(st != "") st += ","; st += dem[i];}
					
					colname.push_back(val);
					dp_sel_list.push_back(create_dp_sel(st));
				}
			}
			else{
				emsg("TO DO");
			}
		}
		
		if(dp_sel_list.size() == 0){
			dp_sel_list.push_back(create_dp_sel(democats));
		}
	}
		
	vector < vector<unsigned int> > area_list; 
	auto geosp = split(geo,':');
	
	auto geomap = create_geomap(tabarea,geosp[0]);
		
	switch(geosp.size()){
		case 1:
			if(geo != "all") colname.clear();
			for(const auto &geom : geomap){
				if(geo != "all") colname.push_back(geom.region);
				area_list.push_back(geom.area);
			}
			break;

		case 2:
			{
				auto areas = split(geosp[1],'|');
				auto used = 0u;
			
				vector <unsigned int> area;
				for(const auto &gm : geomap){
					if(vector_contains(areas,gm.region) == true){
						for(auto c : gm.area) area.push_back(c);
						used++;
					}
				}
				
				if(used != areas.size()) emsg("'"+geosp[1]+"' are not all used");
				area_list.push_back(area);
			}
			break;
			
		default: emsg("'"+geo+"' cannot be interpreted."); break;
	}
	
	if(dp_sel_list.size() > 1 && area_list.size() > 1){
		emsg("With '"+geo+"' and '"+democats+"' cannot have variation in both geograph and demographics");
	}
	
	vector <DataFilter> datafilter;
	
	if(area_list.size() > 1){
		for(auto i = 0u; i < area_list.size(); i++){
			DataFilter df; 
			df.name = colname[i];
			df.area = area_list[i];
			df.dp_sel = dp_sel_list[0];
			datafilter.push_back(df);
		}		
	}
	else{
		for(auto i = 0u; i < dp_sel_list.size(); i++){
			DataFilter df; 
			df.name = colname[i];
			df.area = area_list[0];
			df.dp_sel = dp_sel_list[i];
			datafilter.push_back(df);
		}		
	}
	
	if(false){
		cout << "\n";
		cout << "geo: " <<  geo << "   democats: " << democats << "\n";
		for(auto  df : datafilter){
			cout << df.name << "  Datafilter\n";
			for(auto c : df.area) cout << area[c].code << ","; cout << " areas\n";
			for(auto dp : df.dp_sel) cout << dp << ","; cout << " dp_sel\n";
		}
		emsg("P");
	}

	return datafilter;
}


// Loads up any regional effects
void Data::load_region_effect(const Inputs &inputs, const Table& tab) 
{
	string geography;

	inputs.regional_effects(geography,sigma);	
	if(geography != ""){	
		region_effect = create_geomap(tab,geography);
		for(auto r = 0u; r < region_effect.size(); r++){
			for(auto c : region_effect[r].area) area[c].region_effect = r;
		}
	}
}
	

/// Loads a table giving marginal data
void Data::load_marginal_datatable(const Table &tabarea, unsigned int i, bool sim)
{
	auto& dt = datatable[i];
		
	Observation ob;
	ob.fraction_spline_ref = UNSET;
	ob.datatable = i;
	ob.w = UNSET;
	ob.logdif_offset = UNSET;
	ob.obsmodel = dt.obsmodel; ob.shape = dt.shape; ob.invT = dt.invT; 
	ob.graph = graph.size();
	
	auto datafilter = get_datafilter(tabarea,dt.geo,dt.democats);

	auto file = dt.file; 
	Table tab;
	if(sim == false){
		tab = load_table(file);
	}
	
	string filename = dt.observation, name = dt.observation, desc = dt.observation;

	if(dt.geo != "all"){
		filename += "-"+dt.geo;
		name += " "+dt.geo;
		desc += " in "+dt.geo;
	}
	if(dt.democats != "all"){
		filename += "-"+dt.democats;
		name += "-"+dt.democats;
		desc += " in "+dt.democats;
	}
	
	Graph gr;
	gr.type = GRAPH_MARGINAL;
	
	gr.filename = filename;
	gr.name = name;
	gr.desc = desc;
	
	gr.colname = "";
	gr.fraction_spline_ref = UNSET;
	
	auto flag = false;
	auto row=0u;
	for(const auto& df : datafilter){
		if(sim == false){
			for(row = 0; row < tab.nrow; row++) if(tab.ele[row][0] == df.name) break;
			if(row == tab.nrow) row = UNSET;
			else dt.demolist.push_back(df.name);
		}
		else{
			dt.demolist.push_back(df.name);
		}
	
		if(row != UNSET){
			flag = true;
			
			ob.dp_sel = df.dp_sel;
			ob.area = df.area;
			
			GraphPoint gp;
			gp.xi = row;
			gp.xf = row;
			gp.obs.push_back(obs.size());
			gr.point.push_back(gp);
		
			int ti = dt.start*details.division_per_time;
			int tf = dt.end*details.division_per_time;
		
			if(ti < 0 || tf > int(details.ndivision)) emsg("Data table '"+dt.file+"' is out of the inference time range");
			
			ob.sett_i = ti;
			ob.sett_f = tf;
			
			if(sim == true) ob.value = UNSET;
			else ob.value = get_integer(tab.ele[row][1],file);
						
			obs.push_back(ob);
		}
		if(sim == true) row++;
	}
	if(flag == false) emsg("The data in '"+file+"' was not used");

	dt.graph_ref.push_back(graph.size());
	graph.push_back(gr);
}


/// Sets the offset and weight used for all the observations
void Data::set_datatable_weights()
{
	vector <unsigned int> num_obs(datatable.size());
	vector <double> value_max(datatable.size());

	auto ngraph = 0u;
	for(auto& ob : obs){ if(ob.graph > ngraph) ngraph = ob.graph;}
	ngraph++;

	vector <double> graph_value_max(ngraph);
	for(auto gr = 0u; gr < ngraph; gr++){
		graph_value_max[gr] = 0;
	}
		
	for(auto dt = 0u; dt < datatable.size(); dt++){		
		num_obs[dt] = 0;
		value_max[dt] = 0;
	}
	
	for(auto& ob : obs){		
		auto dt = ob.datatable;
		auto gr = ob.graph;
		num_obs[dt]++;
		auto val = ob.value;
		if(val != UNKNOWN && val != THRESH){
			if(val > value_max[dt]) value_max[dt] = val;
			if(val > graph_value_max[gr]) graph_value_max[gr] = val;
		}
	}
	
	vector <double> graph_offset_store(ngraph);
	for(auto gr = 0u; gr < ngraph; gr++) graph_offset_store[gr] = LARGE;
	
	auto frac = 0.8;
	for(auto& ob : obs){		
		auto dt = ob.datatable;
		auto gr = ob.graph;
		if(gr >= ngraph) emsgEC("data",56);
		
		ob.w = datatable[dt].weight/num_obs[dt];
		
		ob.logdif_offset = 0.05*(frac*graph_value_max[gr] + (1-frac)*value_max[dt]);
		if(ob.logdif_offset < 3) ob.logdif_offset = 3;
		
		if(graph_offset_store[gr] == LARGE) graph_offset_store[gr] = ob.logdif_offset;
		else{
			if(graph_offset_store[gr] != ob.logdif_offset) emsg("offset");
		}
	}
	
	if(false){
		for(auto gr = 0u; gr < graph.size(); gr++){
			cout << graph[gr].name << " offset: " << graph_offset_store[gr] << "   graphmax:" << graph_value_max[gr] << "  tablemax:" << value_max[graph[gr].datatable] << endl;
		}
		emsg("done");
	}
}


/// Creates all the demographic categories consistent with a string 
vector <unsigned int> Data::create_dp_sel(string dp_str) const
{
	vector <unsigned int> dp_sel;

	if(dp_str == "all"){			
		for(auto dp = 0u; dp < ndemocatpos; dp++) dp_sel.push_back(dp);
	}
	else{
		auto st = split(dp_str,',');
		
		vector <unsigned int> demo(ndemocat);
		for(auto d = 0u; d < ndemocat; d++) demo[d] = UNSET;
			
		for(auto val : st){
			auto valsp = split(val,':');
			
			if(valsp.size() != 2) emsg("'"+dp_str+"' must specify the demographic category and its value");
			
			auto d = 0u; while(d < ndemocat && democat[d].name != valsp[0]) d++;
			if(d == ndemocat) emsg("Cannot find '"+valsp[0]+"' for the demographic selection");
			
			auto pos = 0u; while(pos < democat[d].value.size() && democat[d].value[pos] != valsp[1]) pos++;
			if(pos == democat[d].value.size()) emsg("Cannot find '"+valsp[1]+"'");
			
			demo[d] = pos; 
		}
		
		for(auto dp = 0u; dp < ndemocatpos; dp++){
			auto d = 0u; while(d < ndemocat && (demo[d] == UNSET || demo[d] == democatpos[dp][d])) d++;
			
			if(d == ndemocat) dp_sel.push_back(dp);
		}
	}
	
	return dp_sel;
}
 
 
/// Based on a column in the area file a geogaphical mapping is created 
vector <GeographicMap> Data::create_geomap(const Table &tab, string geo) const 
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
 
 
/// Gets a positive integer from a string
unsigned int Data::get_integer(const string& st, const string& file) const
{
	try {
		return ::get_integer(st,threshold);
	} catch (const std::runtime_error& e) {
		emsg("In file '"+file+"', "+e.what());
	}
}


/// Adds age columns together to get a new combined age column
void Data::table_add_age(string name, unsigned int ti, unsigned int tf, Table &tab)
{
	vector <unsigned int> agecols;
	
	if(name == "all"){
		for(auto col = 0u; col < tab.heading.size(); col++){
			if(tab.heading[col] == "population"){
				agecols.push_back(col); break;
			}
		}
	}
	
	if(agecols.size() == 0){
		for(auto col = 0u; col < tab.heading.size(); col++){
			auto head = tab.heading[col];
			if(head.length() > 3){
				if(head.substr(0,3) == "age"){
					unsigned int num = atoi(head.substr(3,head.length()-3).c_str());
					if(num >= ti && num <= tf){
						agecols.push_back(col);
					}
				}
			}
		}
	}
	
	if(false){
		cout << "Add " << agecols.size() << " age columns: ";
		for(auto col : agecols) cout << tab.heading[col] << " ";
		cout << "to generate " << name << endl;
	}
	
	if(agecols.size() == 0) emsgroot("No age columns");
	table_create_column(name,agecols,tab);
}
		
		
/// Loads a table from the data pipeline
Table Data::load_table_from_datapipeline(string file) const
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
Table Data::load_table_from_file(string file, string dir, bool heading, char sep) const
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
	
	cout << "Loaded table " << file << "." << endl;

	Table tab;
	tab.file = file;
	
	string line;
	if(heading == true){
		getline(in,line);

		stringstream ss(line);
		do{
			string st;
			getline(ss,st,sep); st = strip(st);
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
			getline(ss,st,sep); st = strip(st);
			vec.push_back(st);
			if(ss.eof()) break;
		}while(true);
		if(tab.ncol == 0) tab.ncol = vec.size();
		else{
			if(vec.size() != tab.ncol) emsg("Rows in file '"+file+"' do not all share the same number of columns.");
		}
		
		tab.ele.push_back(vec);
	}while(true);
	tab.nrow = tab.ele.size();
	
	return tab;
}


/// Loads a table from a file (if dir is specified then this directory is used
Table Data::load_table(string file, string dir, bool heading) const
{
	if(stringhasending(file, ".txt")){ return load_table_from_file(file,dir,heading,'\t');} 
	else {
		if(stringhasending(file, ".csv")){ return load_table_from_file(file,dir,heading,',');} 
		else return load_table_from_datapipeline(file);
	}
}


/// Creates a new column by adding together existing columns		
void Data::table_create_column(string head, const vector <unsigned int> &cols, Table &tab) const
{
	tab.heading.push_back(head);
	for(auto& trow : tab.ele){
		auto sum = 0u; for(auto i = 0u; i < cols.size(); i++) sum += atoi(trow[cols[i]].c_str());
		stringstream ss; ss << sum;
		trow.push_back(ss.str());
	}
	tab.ncol++;
}


/// Selects dates as specified in the TOML file
void Data::table_select_dates(unsigned int t, unsigned int units, Table &tab, string type, unsigned int shift) const 
{
	auto row = 0u;
	while(row < tab.nrow){
		auto tt = details.gettime(tab.ele[row][0]) + shift - details.start;
		tab.ele[row][0] = details.getdate(tt);
		
		if(tt < t){
			if(type == "trans" && row > 0){ // In the case of transitions adds up the contributions from other days to make a week 
				for(auto i = 1u; i < tab.ncol; i++){
					auto num1 = get_integer(tab.ele[row-1][i],tab.file);
					auto num2 = get_integer(tab.ele[row][i],tab.file);
					if(num1 == THRESH || num2 == THRESH || num1 == UNKNOWN || num2 == UNKNOWN){
						emsg("In the file '"+tab.file+"', amalgamation of data with unknown values is not possible.");
					}
					tab.ele[row-1][i] = to_string(num1+num2);
				}				
			}
			tab.ele.erase(tab.ele.begin()+row);
			tab.nrow--;
		}
		else{
			if(tt > t) emsg("In file '"+tab.file+"' there is no observed data at time '"+details.getdate(t)+"'."); 
			t += units;
			row++;
		}
	}
	if(tab.nrow == 0) emsg("The file '"+tab.file+"' does not contain any information.");
	
	if(type=="trans"){                         // Removes the last line if incomplete
		while(details.gettime(tab.ele[tab.nrow-1][0]) + shift > details.end-units){
			tab.ele.erase(tab.ele.begin()+tab.nrow-1);
			tab.nrow--;
		}
	}	
	
	if(checkon == true){
		for(const auto& trow : tab.ele){
			for(const auto& ele : trow) cout << ele << " ";
			cout << "TABLE" << endl;
		}
	}
}				


/// Finds a column in a table
unsigned int Data::find_column(const Table &tab, string name) const
{
	unsigned int c;
	for(c = 0; c < tab.ncol; c++) if(tab.heading[c] == name) break;
	if(c == tab.ncol) emsg("Cannot find the column heading '"+name+"' in file '"+tab.file+"'.");
	return c;
}		
			

/// Coarse grains the areas to a rougher geography
void Data::coarse_grain(const Table &tab, string coarse)
{		
	if(coarse == "UNSET") return;
	
	vector <double> pop_c(narea);
	auto pop_total = 0.0;
	for(auto c = 0u; c < narea; c++){
		pop_c[c] = 0; for(auto p : area[c].pop) pop_c[c] += p;
		pop_total += pop_c[c];
	}
	
	auto geomap = create_geomap(tab,coarse);
	vector <Area> areanew;
	vector <unsigned int> areaold_ref(narea);
	
	for(auto cnew = 0u; cnew < geomap.size(); cnew++){
		Area are;
		are.code = geomap[cnew].region; 
		are.region_effect = UNSET;
		
		are.x = 0; are.y = 0;
		auto sum = 0.0;
		for(auto c : geomap[cnew].area){
			areaold_ref[c] = cnew;
			 
			are.x += area[c].x*pop_c[c];
			are.y += area[c].y*pop_c[c];
			sum += pop_c[c];
		}
		are.x /= sum; are.y /= sum;
		
		are.pop.resize(ndemocatpos);
		for(auto dp = 0u; dp < ndemocatpos; dp++){
			are.pop[dp] = 0; for(auto c : geomap[cnew].area) are.pop[dp] += area[c].pop[dp]; 
		}
		
		are.covar.resize(ncovar);
		for(auto co = 0u; co < ncovar; co++){
			are.covar[co] = 0;
			auto sum = 0.0;
			for(auto c : geomap[cnew].area){
				are.covar[co] += area[c].covar[co]*pop_c[c];
				sum += pop_c[c];
			}
			are.covar[co] /= sum;
		}
		
		areanew.push_back(are);
	}
	auto areaold = area;
	
	if(false){
		for(auto v: area) cout <<v.code << " codes before coarse graining\n";
	}
	
	area = areanew;
	narea = area.size();
	
	if(false){
		for(auto v: area) cout << v.code << " codes after coarse graining\n";
	}
	
	vector <double> pop_cnew(narea);
	auto pop_totalnew = 0.0;
	for(auto c = 0u; c < narea; c++){
		pop_cnew[c] = 0; for(auto p : area[c].pop) pop_cnew[c] += p;
		pop_totalnew += pop_cnew[c];
	}
	
	SparseMatrix M;                      // Recalibrates matrix of interactions
	M.N = narea;
	M.diag.resize(narea);
	M.to.resize(narea);
	M.val.resize(narea);
	
	for(auto c = 0u; c < narea; c++) M.diag[c] = 0;
	
	for(auto cold = 0u; cold < areaold.size(); cold++){
		auto c = areaold_ref[cold];
		M.diag[c] += genQ.M.diag[c]*pop_c[c]*pop_c[c];
		
		for(auto i = 0u; i < genQ.M.to[cold].size(); i++){
			auto coldto = genQ.M.to[cold][i];
			auto val = genQ.M.val[cold][i];
			
			auto cto = areaold_ref[coldto];
			if(cto == c){
				M.diag[c] += val*pop_c[c]*pop_c[c]; 
			}
			else{
				auto j = 0u; while(j < M.to[c].size() && M.to[c][j] != cto) j++;
				if(j == M.to[c].size()){
					M.to[c].push_back(cto);
					M.val[c].push_back(val*pop_c[c]*pop_c[cto]);
				}
				else{
					M.val[c][j] += val*pop_c[c]*pop_c[cto];
				}
			}
		}
	}
	
	for(auto c = 0u; c < narea; c++){
		M.diag[c] /= (pop_cnew[c]*pop_cnew[c]);
		for(auto j = 0u; j < M.to[c].size(); j++) M.val[c][j] /= (pop_cnew[c]*pop_cnew[M.to[c][j]]);
	}

	genQ.M = M;
	
	geo_normalise(genQ.M);
 
 	if(false){
		//for(auto c = 0u; c < areaold_ref.size(); c++) cout << c << " " << areaold_ref[c] << " arearef\n";
		
		for(auto c = 0u; c < narea; c++){
			cout << c << " diag " << genQ.M.diag[c] << ": ";
			for(auto j = 0u; j < genQ.M.to[c].size(); j++){
				cout << genQ.M.to[c][j] << " " << genQ.M.val[c][j] << ", ";
			}
			cout << "M\n";
		}
		emsg("do");
	}
	
	
	for(auto& ob : obs)	ob.area = convert_areas("observation",ob.area,geomap);                 // Recalibrates observations
	 cout << "beg\n";
	for(auto& gr : graph){
		//for(auto c: gr.area) cout << areaold[c].code << ", "; cout << "code\n";
		gr.area = convert_areas("graph",gr.area,geomap);                     // Recalibrates graphs
	}
	 cout << "end\n";
	for(auto& re : region_effect) re.area = convert_areas("region effect",re.area,geomap);     // Recalibrates regional effects

	for(auto r = 0u; r < region_effect.size(); r++){
		for(auto c : region_effect[r].area) area[c].region_effect = r;
	}
	
	if(false) print_obs();
	
	if(checkon == true){
		for(const auto& are : area){
			cout << are.code << are.x << " " <<  are.y << "  ***";
			for(const auto& pop : are.pop) cout << pop << ", ";
			cout << endl;	
		}
	}	
}


/// Determines if a vector contains a given element
bool Data::vector_contains(const vector <unsigned int> &vec, unsigned int num) const
{
	if(find(vec.begin(),vec.end(),num) != vec.end()) return true;
	return false;
}	


/// Determines if a string vector contains a given element
bool Data::vector_contains(const vector <string> &vec, string num) const
{
	if(find(vec.begin(),vec.end(),num) != vec.end()) return true;
	return false;
}	


/// Removes an element from a vector
void Data::vector_remove(vector <unsigned int> &vec, unsigned int num) const
{
	auto i = find(vec.begin(),vec.end(),num);
	if(i == vec.end()) emsgEC("Data",101);
	vec.erase(i);
}	


/// Converts a set of areas to a new set provides by a geographical map
vector <unsigned int> Data::convert_areas(string type, vector <unsigned int> &vec, const vector <GeographicMap> &map) const
{
	vector <unsigned int> vec_new;
	for(auto cnew = 0u; cnew < map.size(); cnew++){
		auto i = 0u; while(i < map[cnew].area.size() && vector_contains(vec,map[cnew].area[i]) == true) i++;
		if(i == map[cnew].area.size()){
			vec_new.push_back(cnew);
			for(auto c: map[cnew].area) vector_remove(vec,c); 
		}
	}	

	if(vec.size() != 0) emsg("Cannot course grain "+type);
	
	return vec_new;
}


/// Loads up counterfactual information  
void Data::load_counterfactual(const Inputs &inputs, const Table &tabarea)  
{
	if(details.mode == COUNTER){
		inputs.find_counterfact(details,counterfact);
		for(auto &cf : counterfact){
			cf.dp_sel = create_dp_sel(cf.democats);
			
			auto geosp = split(cf.geo,':');
			auto geomap = create_geomap(tabarea,geosp[0]);
			if(geosp.size() == 2) emsg("This is not implemented yet");
			if(geomap.size() != 1) emsg("This is not implemented yet");
			cf.area = geomap[0].area;
		}
	}
	
	if(details.mode == PPC){
		inputs.find_ppc_start(details,counterfact);
	}
}


/// Outputs properties of data to the terminal
string Data::print() const
{
	stringstream ss;
	ss << endl;
	
	ss << "Age categories: " << endl << "  ";
	for(auto j = 0u; j < democat[0].value.size(); j++){
		if(j != 0) ss << " | ";
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
		for(const auto& cov : covar){
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
	for(auto& ob : obs){
		cout << "Datatable: " << ob.datatable << "  ";
		cout << "Graph: " << ob.graph << "  ";
		cout << "Value: " << ob.value << "  ";
		cout << "ti: " << ob.sett_i << "  ";
		cout << "tf: " << ob.sett_f << "  ";
		cout << "offset: " << ob.logdif_offset << "  ";
		cout << "w: " << ob.w << "  ";
		cout << "area: "; for(auto c : ob.area) cout <<c << ",  ";
		cout << "dp_sel: "; for(auto dp : ob.dp_sel) cout << dp << ",  ";
		cout << "\n\n";
	}
}
