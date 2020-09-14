//  The data inputs

#include <iostream>
#include <vector>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <math.h> 

using namespace std;

#include "utils.hh"
#include "data.hh"	
#include "consts.hh"
#include "pack.hh"
#include "generateQ.hh"
#include "output.hh"
#include "details.hh"
#include "inputs.hh"

#ifdef USE_Data_PIPELINE
#include "datapipeline.hh"
#include "table.hh"
#endif

/// Initialises data
Data::Data(const Inputs &inputs, const Details &details, const Mpi &mpi, DataPipeline *dp) :
	datapipeline(dp), data_directory(inputs.find_string("datadir","UNSET")),
	compX(area), compY(area), details(details)
{
	// The data directory
	if(data_directory == "UNSET") emsgroot("The 'data_directory' must be set.");

	threshold = inputs.find_integer("threshold",UNSET);                                       // The threshold (if specified)
	if(threshold != UNSET) thres_h = log(1.0/(threshold + 0.5*sqrt(2*M_PI*MINIMUM_VARIANCE)));
	else thres_h = UNSET;
	
	democat = inputs.find_demographic_category(details);
	ndemocat = democat.size();	
	nage = democat[0].value.size();
	calc_democatpos();

	covar = inputs.find_covariate(details);
	ncovar = covar.size();
	
	time_period = inputs.find_time_period(details);
	
	inputs.find_genQ(genQ,details);
	inputs.find_Q(Q,time_period,details);
	
	read_data_files(inputs,mpi);                                  // Reads the data files
}

/// Outputs properties of data to the terminal
void Data::print_to_terminal() const
{
	if(covar.size() > 0){
		cout << "Area covariates: " << endl;
		cout << "  ";
		for(const auto& cov : covar) cout << cov.name << "   param='" << cov.param << "'" << endl; 
		cout << endl;
	}	
	
	cout << "Age categories: " << endl << "  ";
	for(auto j = 0u; j < democat[0].value.size(); j++){
		if(j != 0) cout << ", ";
		cout << democat[0].value[j] << " sus='" << democat[0].param[j] << "'";
	}
	cout << endl << endl;

	if(democat.size() > 1){
		cout << "Demographic categories: " << endl;
		for(auto k = 1u; k < democat.size(); k++){
			cout << "  ";
			for(auto j = 0u; j < democat[k].value.size(); j++){
				if(j != 0) cout << ", ";
				cout << democat[k].value[j] << " sus='" <<  democat[k].param[j] << "'";
			}	
			cout << endl;
		}
		cout << endl;
	}
	
	cout << "Time periods defined:" << endl;
	for(auto j = 0u; j < time_period.size(); j++){
		cout << "  ";
		cout << time_period[j].name << ": ";
		if(j == 0) cout << "0"; else cout << time_period[j-1].tend;
		cout << " - " <<  time_period[j].tend << endl;
	}
	cout << endl;
	
	cout << "Q tensors loaded:" << endl;
	for(const auto& QQ : Q){
		cout << "    ";
		cout << "timep: " << time_period[QQ.timep].name << "  ";
		cout << "compartment: " << QQ.comp << "  ";
		cout << "name: " << QQ.name << "  ";
		cout << endl;
	}
	cout << endl;
}

/// Based to the different demographic categories, this calculates all the possible combinations
void Data::calc_democatpos()
{
	vector <unsigned int> count(ndemocat);       // Defines all the demographic states
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

/// Loads the region data file
void Data::load_region_file(const Inputs &inputs)
{
	string file = inputs.find_string("regions","UNSET");
	if(file == "UNSET") emsgroot("A 'regions' file must be specified");
	
	Table tab = load_table(file);
	
	auto namecol = find_column(tab,"name");
	auto codecol = find_column(tab,"code");
	
	for(const auto& trow : tab.ele){
		DataRegion reg;
		reg.name = trow[namecol];
		reg.code = trow[codecol];
		region.push_back(reg);
	}		
	nregion = region.size();
}

/// Reads in transition and area data
void Data::read_data_files(const Inputs &inputs, const Mpi &mpi)
{
	transdata = inputs.find_transition_data(details);	                                          // Loads data
	
	popdata = inputs.find_population_data(details);	                                       
	
	margdata = inputs.find_marginal_data(details,democat);	 

	if(transdata.size() == 0 && popdata.size() == 0)  emsgroot("'transdata' and/or 'popdata' must be set.");	

	if(mpi.core == 0){
		load_region_file(inputs);
	
		string file = inputs.find_string("areas","UNSET");
		if(file == "UNSET") emsgroot("A 'areas' file must be specified");
		Table tab = load_table(file);
		
		// If onle one age group then combines all columns with "age" to generate an "all" column
		if(nage == 1){
			vector <unsigned int> agecols;
			for(auto col = 0u; col < tab.heading.size(); col++){
				if(tab.heading[col].length() > 3){
					if(tab.heading[col].substr(0,3) == "age") agecols.push_back(col);
				}
			}
			table_create_column("all",agecols,tab);
		}
		
		auto codecol = find_column(tab,"area");                                                 // Works out ccolumns for different columns
		auto xcol = find_column(tab,"easting");
		auto ycol = find_column(tab,"northing");
		auto regcol = find_column(tab,"region");
		
		for(auto& cov : covar) cov.col = find_column(tab,cov.name);

		for(auto& democ : democat){
			democ.col.resize(democ.value.size());
			for(auto k = 0u; k < democ.col.size(); k++) democ.col[k] = find_column(tab,democ.value[k]);  
		}

		for(const auto& trow : tab.ele){
			Area are;
			are.code = trow[codecol];
			are.x = atof(trow[xcol].c_str());
			are.y = atof(trow[ycol].c_str());
	
			auto regcode = trow[regcol];
			auto r = 0u; while(r < nregion && region[r].code != regcode) r++;
			if(r == nregion) emsg("In file '"+file+"' the region code '"+regcode+"' is not recognised.");
			are.region = r;
					
			are.covar.resize(ncovar);
			for(auto j = 0u; j < ncovar; j++){
				auto st = trow[covar[j].col];
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
					auto st = trow[democat[d].col[j]];
					val[d][j] = atof(st.c_str());
					if(std::isnan(val[d][j])) emsg("In file '"+file+"' the expression '"+st+"' is not a number");	
				}
			}
			
			are.agepop.resize(nage);
			for(auto a = 0u; a < nage; a++){
				are.agepop[a] = val[0][a];
			}
				
			are.pop.resize(ndemocatpos);
			for(auto dp = 0u; dp < ndemocatpos; dp++){
				auto v = 0.0;
				for(auto j = 0u; j < democat.size(); j++){
					if(j == 0) v = val[j][democatpos[dp][j]];
					else v *= val[j][democatpos[dp][j]]/100.0;
				}
				are.pop[dp] = (unsigned int)(v+0.5);
			}

			area.push_back(are);			
		}
		narea = area.size();
		
		if(false){  // Averages covariates across regions
			vector <double> av, nav;
			av.resize(region.size()); nav.resize(region.size());
			for(auto j = 0u; j < ncovar; j++){  
				for(auto r = 0u; r < region.size(); r++){ av[r] = 0; nav[r] = 0;}
				
				for(auto c = 0u; c < narea; c++){
					auto r = area[c].region;
					if(covar[j].func == "log") av[r] += exp(area[c].covar[j]);
					else av[r] += area[c].covar[j];
					nav[r]++;
				}
				
				for(auto c = 0u; c < narea; c++){
					auto r = area[c].region;
					if(covar[j].func == "log") area[c].covar[j] = log(av[r]/nav[r]);
					else area[c].covar[j] = av[r]/nav[r];
				}
				
				for(auto r = 0u; r < region.size(); r++) cout << region[r].name << " " << av[r]/nav[r] << " average density" << endl;
			}
		}
		
		for(auto j = 0u; j < ncovar; j++){            // Shifts covariates so average is zero
			auto sum = 0.0; for(auto c = 0u; c < narea; c++) sum += area[c].covar[j];
			sum /= narea;
			
			for(auto c = 0u; c < narea; c++) area[c].covar[j] -= sum;
		}		
		
		if(checkon == true){
			for(const auto& are : area){
				cout << nregion << " " << are.region << "region" << endl;
				cout << are.code << " " << region[are.region].code << " " << are.x << " " <<  are.y << "  ***";
			
				for(const auto& pop : are.pop) cout << pop << ", ";
				cout << endl;	
			}
		}

		if(details.mode != SIM && details.mode != MULTISIM){         
			for(auto& tdata : transdata){                                // Loads transition data for inference
				file = tdata.file; 
				Table tab = load_table(file);
				table_select_dates(tdata.start,tdata.units,tab,"trans");
				
				vector <unsigned int> rcol;
				if(tdata.type == "reg"){	for(const auto& reg : region) rcol.push_back(find_column(tab,reg.code));}
				else{ rcol.push_back(find_column(tab,"all"));}
				
				tdata.num.resize(rcol.size());
				for(auto r = 0u; r < rcol.size(); r++){
					
					tdata.num[r].resize(tab.nrow);
					for(auto row = 0u; row < tab.nrow; row++){	
						tdata.num[r][row] = get_integer(tab.ele[row][rcol[r]],file);
					}
				}
				tdata.rows = tab.nrow;
				if(tdata.start + (tab.nrow-1)*tdata.units > details.period){
					emsg("The file '"+file+"' has more data than will fit in the defined 'start' and 'end' time period.");
				}
			}
	                                     
			for(auto& pdata : popdata){                            // Loads population data for inference
				file = pdata.file;
				Table tab = load_table(file);
				table_select_dates(pdata.start,pdata.units,tab,"pop");
			
				vector <unsigned int> rcol;
				if(pdata.type == "reg"){	for(const auto& reg : region) rcol.push_back(find_column(tab,reg.code));}
				else{ rcol.push_back(find_column(tab,"all"));}
				
				pdata.num.resize(rcol.size());
				for(auto r = 0u; r < rcol.size(); r++){
					pdata.num[r].resize(tab.nrow);
					for(auto row = 0u; row < tab.nrow; row++) pdata.num[r][row] = get_integer(tab.ele[row][rcol[r]],file);
				}
				pdata.rows = tab.nrow;
				
				if(pdata.start + (tab.nrow-1)*pdata.units > details.period){
					emsg("The file '"+file+"' has more data than will fit in the defined 'start' and 'end' time period.");
				}
			}
		                  
			for(auto& mdata : margdata){                                        // Loads marginal data for inference
				file = mdata.file;
				Table tab = load_table(file);
	
				vector <unsigned int> rcol;
				if(mdata.type == "reg"){	for(const auto& reg : region) rcol.push_back(find_column(tab,reg.code));}
				else{ rcol.push_back(find_column(tab,"all"));}
				
				mdata.percent.resize(rcol.size());
				for(auto r = 0u; r < rcol.size(); r++){
					mdata.percent[r].resize(tab.nrow);
					for(auto row = 0u; row < tab.nrow; row++){	
						mdata.percent[r][row] = atof(tab.ele[row][rcol[r]].c_str());
					}
				}
			}
		}
		
		generateQ(nage,data_directory,genQ,area,datapipeline);
	}

	if(mpi.ncore > 1) copy_data(mpi.core);
		
	vector <double> vec(nage);                                                           // Reads in Q tensors
	for(auto& QQ : Q){
		auto j = 0u; while(j < genQ.Qten.size() && genQ.Qten[j].name != QQ.name) j++;
		if(j == genQ.Qten.size()) emsg("Cannot find the reference to '"+QQ.name+"' in the input TOML file.");
		QQ.Qtenref = j;
	}
		
	agedist.resize(nage); for(auto& aged : agedist) aged = 0;
	
	for(auto c = 0u; c < narea; c++){                                              // Adds individuals to the system
		area[c].ind.resize(ndemocatpos);
		for(auto dp = 0u; dp < ndemocatpos; dp++){
			auto imax = area[c].pop[dp];
			for(auto i = 0u; i < imax; i++){
				area[c].ind[dp].push_back(ind.size());
				
				Individual indi;
				indi.area = c;
				indi.dp = dp;
				ind.push_back(indi);
			}
			auto a = democatpos[dp][0];
			agedist[a] += imax;
		}
	}
	popsize = ind.size();
	for(auto& aged : agedist) aged /= popsize;
	
	narage = narea*nage;                                              // Generates the mixing matrix between ages/areas
	nardp = narea*ndemocatpos; 
	nsettardp = details.ndivision*nardp;
	
	//plotrawdata(); emsg("done");
	//generatedeathdata(); emsg("done");
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

/// Loads a table from the data pipeline
Table Data::load_table_from_datapipeline(string file) const
{
	Table tab;
#ifdef USE_Data_PIPELINE
	Table dptable = datapipeline->read_table(file,"default");

	tab.file = filebasename(file);
	tab.heading = dptable.get_column_names();
	tab.ncol = tab.heading.size();
	
	vector<vector<string>> cols;
	
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
Table Data::load_table_from_file(string file, string dir) const
{
	ifstream in;

	string used_file;
	if(dir != ""){
    used_file = dir+"/"+file;
		in.open(used_file.c_str());
		if(!in) emsg("Cannot open the file '"+dir+"/"+file+"'.");
	}
	else{
    used_file = details.output_directory+"/"+file;
		in.open(used_file.c_str());
		if(!in){
      used_file = data_directory+"/"+file;
			in.open(used_file.c_str());
			if(!in) emsg("Cannot open the file '"+used_file+"'");
		}
	}
	
	cout << "Loaded table " << file << " from file " << used_file << endl;

	Table tab;
	tab.file = file;
	
	string line;
	getline(in,line);

	stringstream ss(line);
	do{
		string st;
		getline(ss,st,'\t'); st = strip(st);
		tab.heading.push_back(st);
		if(ss.eof()) break;
	}while(true);
	tab.ncol = tab.heading.size();
	
	do{
		vector <string> vec;
		getline(in,line);
		if(in.eof()) break;
				
		stringstream ss(line);
		do{
			string st;
			getline(ss,st,'\t'); st = strip(st);
			vec.push_back(st);
			if(ss.eof()) break;
		}while(true);
		if(vec.size() != tab.ncol) emsg("Rows in file '"+file+"' do not all share the same number of columns.");
		
		tab.ele.push_back(vec);
	}while(true);
	tab.nrow = tab.ele.size();
	
	return tab;
}

/// Loads a table from a file (if dir is specified then this directory is used
Table Data::load_table(string file, string dir) const
{
	if (stringhasending(file, ".txt")) {
		return load_table_from_file(file, dir);
	} else {
		return load_table_from_datapipeline(file);
	}
}

/// Creates a new column by adding together existing columns		
void Data::table_create_column(string head,vector <unsigned int> cols, Table &tab) const
{
	tab.heading.push_back(head);
	for(auto& trow : tab.ele){
		auto sum = 0u; for(auto i = 0u; i < cols.size(); i++) sum += atoi(trow[cols[i]].c_str());
		stringstream ss; ss << sum;
		trow.push_back(ss.str());
	}
	tab.ncol++;
}

/// Selects dates as specified in the TOLM file
void Data::table_select_dates(unsigned int t, unsigned int units, Table &tab, string type) const 
{
	auto row = 0u;
	while(row < tab.nrow){
		auto tt = details.gettime(tab.ele[row][0]) - details.start;
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
		if(details.gettime(tab.ele[tab.nrow-1][0]) > details.end-units){
			tab.ele.erase(tab.ele.begin()+tab.nrow-1);
			tab.nrow--;
		}
	}	
	
	if(checkon == true){
		for(const auto& trow : tab.ele){
			for(const auto& ele : trow) cout << ele << " ";
			cout << endl;
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
			
/// Copies data from core zero to all the others
void Data::copy_data(unsigned int core)
{
	unsigned int si;

	if(core == 0){                                  				   // Copies the above information to all the other cores
		pack_initialise(0);
		pack(nregion);
		pack(region);
		pack(narea);
		pack(area);
		for(const auto& tdata : transdata){
			pack(tdata.num);
			pack(tdata.rows);
		}
		for(const auto& pdata : popdata){
			pack(pdata.num);
			pack(pdata.rows);
		}
		for(const auto& mdata : margdata){
			pack(mdata.percent);
		}
		unsigned int kmax = genQ.Qten.size();
		pack(kmax);
		for(const auto& Qten : genQ.Qten){
			pack(Qten.name);
		}
		si = packsize();
	}

	pack_mpi_bcast();

	if(core != 0){
		unpack(nregion);
		unpack(region);
		unpack(narea);
		unpack(area);
		for(auto& tdata : transdata){
			unpack(tdata.num);
			unpack(tdata.rows);
		}
		for(auto& pdata : popdata){
			unpack(pdata.num);
			unpack(pdata.rows);
		}
		for(auto& mdata : margdata){
			unpack(mdata.percent);
		}
		unsigned int kmax;
		unpack(kmax);
		genQ.Qten.resize(kmax);
		for(auto& Qten : genQ.Qten){
			unpack(Qten.name);
		}
		if(si != packsize()) emsgEC("Data",1);
	}

	for(auto& Qten : genQ.Qten){                                                   // Copies the Q matrices
		auto num = narea*nage;
		mpi_bcast(num);
		
		if(core != 0){
			Qten.ntof.resize(num);
			Qten.tof.resize(num);
			Qten.valf.resize(num);
		}
		
		auto vmin = 0u;
		do{
			auto vmax = vmin+1000; if(vmax > num) vmax = num;
			
			if(core == 0){
				pack_initialise(0);
				for(auto v = vmin; v < vmax; v++){
					pack(Qten.ntof[v]);
					pack(Qten.tof[v]);
					pack(Qten.valf[v]);
				}
				si = packsize();
			}

			pack_mpi_bcast();

			if(core != 0){
				for(auto v = vmin; v < vmax; v++){
					unpack(Qten.ntof[v]);
					unpack(Qten.tof[v]);
					unpack(Qten.valf[v]);
				}
				if(si != packsize()) emsgEC("Data",2);
			}
	
			vmin = vmax;
		}while(vmin < num);
	}	
}

/// Removes 'r' and quotations from a string
string Data::strip(string line) const 
{
	unsigned int len;
	
	len = line.length();
	if(len > 0 && line.substr(len-1,1) == "\r") line = line.substr(0,len-1);
	len = line.length();
	if(len > 0 && line.substr(0,1) == "\"") line = line.substr(1,len-1);
	len = line.length();
	if(len > 0 && line.substr(len-1,1) == "\"") line = line.substr(0,len-1);
	
	return line;
}	

void Data::sortX(vector <unsigned int> &vec){ sort(vec.begin(),vec.end(),compX);}
void Data::sortY(vector <unsigned int> &vec){ sort(vec.begin(),vec.end(),compY);}

