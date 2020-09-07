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

#ifdef USE_DATA_PIPELINE
#include "datapipeline.hh"
#include "table.hh"
#endif

/// Initialises data
DATA::DATA(const Inputs &inputs, const Details &details, const Mpi &mpi, DataPipeline *dp) :
	datapipeline(dp), datadir(inputs.find_string("datadir","UNSET")),
	compX(area), compY(area), details(details)
{
	// The data directory
	if(datadir == "UNSET") emsgroot("The 'datadir' must be set.");

	threshold = inputs.find_int("threshold",UNSET);                                       // The threshold (if specified)
	if(threshold != UNSET) thres_h = log(1.0/(threshold + 0.5*sqrt(2*M_PI*minvar)));
	else thres_h = UNSET;
	
	democat = inputs.find_democat(details);
	ndemocat = democat.size();	
	nage = democat[0].value.size();
	calc_democatpos();

	covar = inputs.find_covar(details);
	ncovar = covar.size();
	
	timeperiod = inputs.find_timeperiod(details);
	
	inputs.find_genQ(genQ,details);
	inputs.find_Q(Q,timeperiod,details);
	
	read_data_files(inputs,mpi);                                  // Reads the data files
}

/// Outputs properties of data to the terminal
void DATA::print_to_terminal() const
{
	if(covar.size() > 0){
		cout << "Area covariates: " << endl;
		cout << "  ";
		for(auto j = 0u; j < covar.size(); j++) cout << covar[j].name << "   param='" << covar[j].param << "'" << endl; 
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
	for(auto j = 0u; j < timeperiod.size(); j++){
		cout << "  ";
		cout << timeperiod[j].name << ": ";
		if(j == 0) cout << "0"; else cout << timeperiod[j-1].tend;
		cout << " - " <<  timeperiod[j].tend << endl;
	}
	cout << endl;
	
	cout << "Q tensors loaded:" << endl;
	for(auto j = 0u; j < Q.size(); j++){
		cout << "    ";
		cout << "timep: " << timeperiod[Q[j].timep].name << "  ";
		cout << "compartment: " << Q[j].comp << "  ";
		cout << "name: " << Q[j].name << "  ";
		cout << endl;
	}
	cout << endl;
}

/// Based to the different demographic categories, this calculates all the possible combinations
void DATA::calc_democatpos()
{
	vector <unsigned int> count;
	count.resize(ndemocat);                                     // Defines all the demographic states
	
	for(auto dc = 0u; dc < ndemocat; dc++) count[dc] = 0;
	
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
void DATA::load_region_file(const Inputs &inputs)
{
	string file = inputs.find_string("regions","UNSET");
	if(file == "UNSET") emsgroot("A 'regions' file must be specified");
	
	TABLE tab = loadtable(file);
	
	auto namecol = findcol(tab,"name");
	auto codecol = findcol(tab,"code");
	
	for(auto row = 0u; row < tab.nrow; row++){
		REGION reg;
		reg.name = tab.ele[row][namecol];
		reg.code = tab.ele[row][codecol];
		region.push_back(reg);
	}		
	nregion = region.size();
}

/// Reads in transition and area data
void DATA::read_data_files(const Inputs &inputs, const Mpi &mpi)
{
	transdata = inputs.find_transdata(details);	                                          // Loads data
	
	popdata = inputs.find_popdata(details);	                                       
	
	margdata = inputs.find_margdata(details,democat);	 

	if(transdata.size() == 0 && popdata.size() == 0)  emsgroot("'transdata' and/or 'popdata' must be set.");	

	if(mpi.core == 0){
		load_region_file(inputs);
	
		string file = inputs.find_string("areas","UNSET");
		if(file == "UNSET") emsgroot("A 'areas' file must be specified");
		TABLE tab = loadtable(file);
		
		// If onle one age group then combines all columns with "age" to generate an "all" column
		if(nage == 1){
			vector <unsigned int> agecols;
			for(auto col = 0u; col < tab.heading.size(); col++){
				if(tab.heading[col].length() > 3){
					if(tab.heading[col].substr(0,3) == "age") agecols.push_back(col);
				}
			}
			table_createcol("all",agecols,tab);
		}
		
		auto codecol = findcol(tab,"area");                                                 // Works out ccolumns for different columns
		auto xcol = findcol(tab,"easting");
		auto ycol = findcol(tab,"northing");
		auto regcol = findcol(tab,"region");
		
		for(auto j = 0u; j < ncovar; j++) covar[j].col = findcol(tab,covar[j].name);

		for(auto j = 0u; j < ndemocat; j++){
			democat[j].col.resize(democat[j].value.size());
			for(auto k = 0u; k < democat[j].col.size(); k++) democat[j].col[k] = findcol(tab,democat[j].value[k]);  
		}

		for(auto row = 0u; row < tab.nrow; row++){
			AREA are;
			are.code = tab.ele[row][codecol];
			are.x = atof(tab.ele[row][xcol].c_str());
			are.y = atof(tab.ele[row][ycol].c_str());
	
			auto regcode = tab.ele[row][regcol];
			auto r = 0u; while(r < nregion && region[r].code != regcode) r++;
			if(r == nregion) emsg("In file '"+file+"' the region code '"+regcode+"' is not recognised.");
			are.region = r;
					
			are.covar.resize(ncovar);
			for(auto j = 0u; j < ncovar; j++){
				auto st = tab.ele[row][covar[j].col];
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
					auto st = tab.ele[row][democat[d].col[j]];
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
		
		if(checkon == 1){
			for(auto c = 0u; c < narea; c++){
				cout << nregion << " " << area[c].region << "region" << endl;
				cout << area[c].code << " " << region[area[c].region].code << " " <<  area[c].x << " " <<  area[c].y << "  ***";
			
				for(auto j = 0u; j < area[c].pop.size(); j++) cout << area[c].pop[j] << ", ";
				cout << endl;	
			}
		}
		
		//convertOAtoM(); emsg("done");
		//convertRegion_M(); emsg("done");

		if(details.mode != sim && details.mode != multisim){                                                    // Loads transition data for inference
			for(auto td = 0u; td < transdata.size(); td++){
				file = transdata[td].file; 
				TABLE tab = loadtable(file);
				table_selectdates(transdata[td].start,transdata[td].units,tab,"trans");
				
				vector <unsigned int> rcol;
				if(transdata[td].type == "reg"){	for(auto k = 0u; k < region.size(); k++) rcol.push_back(findcol(tab,region[k].code));}
				else{ rcol.push_back(findcol(tab,"all"));}
				
				transdata[td].num.resize(rcol.size());
				for(auto r = 0u; r < rcol.size(); r++){
					
					transdata[td].num[r].resize(tab.nrow);
					for(auto row = 0u; row < tab.nrow; row++){	
						transdata[td].num[r][row] = getint(tab.ele[row][rcol[r]],file);
					}
				}
				transdata[td].rows = tab.nrow;
				if(transdata[td].start + (tab.nrow-1)*transdata[td].units > details.period){
					emsg("The file '"+file+"' has more data than will fit in the defined 'start' and 'end' time period.");
				}
			}
		}
		
		if(details.mode != sim && details.mode != multisim){                                                    // Loads population data for inference
			for(auto pd = 0u; pd < popdata.size(); pd++){
				file = popdata[pd].file;
				TABLE tab = loadtable(file);
				table_selectdates(popdata[pd].start,popdata[pd].units,tab,"pop");
			
				vector <unsigned int> rcol;
				if(popdata[pd].type == "reg"){	for(auto k = 0u; k < region.size(); k++) rcol.push_back(findcol(tab,region[k].code));}
				else{ rcol.push_back(findcol(tab,"all"));}
				
				popdata[pd].num.resize(rcol.size());
				for(auto r = 0u; r < rcol.size(); r++){
					popdata[pd].num[r].resize(tab.nrow);
					for(auto row = 0u; row < tab.nrow; row++) popdata[pd].num[r][row] = getint(tab.ele[row][rcol[r]],file);
				}
				popdata[pd].rows = tab.nrow;
				
				if(popdata[pd].start + (tab.nrow-1)*popdata[pd].units > details.period){
					emsg("The file '"+file+"' has more data than will fit in the defined 'start' and 'end' time period.");
				}
			}
		}
		
		if(details.mode != sim && details.mode != multisim){                                                    // Loads marginal data for inference
			for(auto md = 0u; md < margdata.size(); md++){
				file = margdata[md].file;
				TABLE tab = loadtable(file);
	
				vector <unsigned int> rcol;
				if(margdata[md].type == "reg"){	for(auto k = 0u; k < region.size(); k++) rcol.push_back(findcol(tab,region[k].code));}
				else{ rcol.push_back(findcol(tab,"all"));}
				
				margdata[md].percent.resize(rcol.size());
				for(auto r = 0u; r < rcol.size(); r++){
					margdata[md].percent[r].resize(tab.nrow);
					for(auto row = 0u; row < tab.nrow; row++){	
						margdata[md].percent[r][row] = atof(tab.ele[row][rcol[r]].c_str());
					}
				}
			}
		}
		
		generateQ(nage,datadir,genQ,area,datapipeline);
	}

	if(mpi.ncore > 1) copydata(mpi.core);
	
	vector <double> vec(nage);                                                           // Reads in Q tensors
	for(auto q = 0u; q < Q.size(); q++){
		auto j = 0u; while(j < genQ.Qten.size() && genQ.Qten[j].name != Q[q].name) j++;
		if(j == genQ.Qten.size()) emsg("Cannot find the reference to '"+Q[q].name+"' in the input TOML file.");
		Q[q].Qtenref = j;
	}
		
	agedist.resize(nage); for(auto a = 0u; a < nage; a++) agedist[a] = 0;
	
	for(auto c = 0u; c < narea; c++){                                              // Adds individuals to the system
		area[c].ind.resize(ndemocatpos);
		for(auto dp = 0u; dp < ndemocatpos; dp++){
			auto imax = area[c].pop[dp];
			for(auto i = 0u; i < imax; i++){
				area[c].ind[dp].push_back(ind.size());
				
				IND indi;
				indi.area = c;
				indi.dp = dp;
				ind.push_back(indi);
			}
			auto a = democatpos[dp][0];
			agedist[a] += imax;
		}
	}
	popsize = ind.size();
	for(auto a = 0u; a < nage; a++) agedist[a] /= popsize;
	
	narage = narea*nage;                                              // Generates the mixing matrix between ages/areas
	nardp = narea*ndemocatpos; 
	nsettardp = details.nsettime*nardp;
	
	//plotrawdata(); emsg("done");
	//generatedeathdata(); emsg("done");
}

/// Gets a positive integer from a string
unsigned int DATA::getint(const string& st, const string& file) const
{
	try {
		return ::getint(st,threshold);
	} catch (const std::runtime_error& e) {
		emsg("In file '"+file+"', "+e.what());
	}
}

/// Loads a table from the data pipeline
TABLE DATA::loadtablefromdatapipeline(string file) const
{
	TABLE tab;
#ifdef USE_DATA_PIPELINE
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
	emsg("loadtablefromdatapipeline for '"+file+"' cannot be called as data pipeline is not compiled in");
#endif

	return tab;
}

/// Loads a table from a file
TABLE DATA::loadtablefromfile(string file, string dir) const
{
	ifstream in;

	string used_file;
	if(dir != ""){
    used_file = dir+"/"+file;
		in.open(used_file.c_str());
		if(!in) emsg("Cannot open the file '"+dir+"/"+file+"'.");
	}
	else{
    used_file = details.outputdir+"/"+file;
		in.open(used_file.c_str());
		if(!in){
      used_file = datadir+"/"+file;
			in.open(used_file.c_str());
			if(!in) emsg("Cannot open the file '"+used_file+"'");
		}
	}
	
	cout << "Loaded table " << file << " from file " << used_file << endl;

	TABLE tab;
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
TABLE DATA::loadtable(string file, string dir) const
{
	if (stringhasending(file, ".txt")) {
		return loadtablefromfile(file, dir);
	} else {
		return loadtablefromdatapipeline(file);
	}
}

/// Creates a new column by adding together existing columns		
void DATA::table_createcol(string head,vector <unsigned int> cols, TABLE &tab) const
{
	tab.heading.push_back(head);
	for(auto row = 0u; row < tab.nrow; row++){
		auto sum = 0u; for(auto i = 0u; i < cols.size(); i++) sum += atoi(tab.ele[row][cols[i]].c_str());
		stringstream ss; ss << sum;
		tab.ele[row].push_back(ss.str());
	}
	tab.ncol++;
}

/// Selects dates as specified in the TOLM file
void DATA::table_selectdates(unsigned int t, unsigned int units, TABLE &tab, string type) const 
{
	auto row = 0u;
	while(row < tab.nrow){
		auto tt = details.gettime(tab.ele[row][0]) - details.start;
		if(tt < t){
			if(type == "trans" && row > 0){ // In the case of transitions adds up the contributions from other days to make a week 
				for(auto i = 1u; i < tab.ncol; i++){
					auto num1 = getint(tab.ele[row-1][i],tab.file);
					auto num2 = getint(tab.ele[row][i],tab.file);
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
	
	if(checkon == 1){
		for(auto row = 0u; row < tab.nrow; row++){
			for(auto i = 0u; i < tab.ncol; i++) cout << tab.ele[row][i] << " ";
			cout << endl;
		}
	}
}				

/// Finds a column in a table
unsigned int DATA::findcol(const TABLE &tab, string name) const
{
	unsigned int c;
	
	for(c = 0; c < tab.ncol; c++) if(tab.heading[c] == name) break;
	if(c == tab.ncol) emsg("Cannot find the column heading '"+name+"' in file '"+tab.file+"'.");
	return c;
}		
			
/// Copies data from core zero to all the others
void DATA::copydata(unsigned int core)
{
	size_t si;
	
	if(core == 0){                                  				   // Copies the above information to all the other cores
		packinit(0);
		pack(nregion);
		pack(region);
		pack(narea);
		pack(area);
		for(auto td = 0u; td < transdata.size(); td++){
			pack(transdata[td].num);
			pack(transdata[td].rows);
		}
		for(auto pd = 0u; pd < popdata.size(); pd++){
			pack(popdata[pd].num);
			pack(popdata[pd].rows);
		}
		for(auto md = 0u; md < margdata.size(); md++){
			pack(margdata[md].percent);
		}
		unsigned int kmax = genQ.Qten.size();
		pack(kmax);
		for(auto k = 0u; k < kmax; k++){
			pack(genQ.Qten[k].name);
		}
		si = packsize();
	}

	MPI_Bcast(&si,1,MPI_UNSIGNED,0,MPI_COMM_WORLD);
	if(core != 0){
		packinit(si);
	}
	MPI_Bcast(packbuffer(),si,MPI_DOUBLE,0,MPI_COMM_WORLD);

	if(core != 0){
		unpack(nregion);
		unpack(region);
		unpack(narea);
		unpack(area);
		for(auto td = 0u; td < transdata.size(); td++){
			unpack(transdata[td].num);
			unpack(transdata[td].rows);
		}
		for(auto pd = 0u; pd < popdata.size(); pd++){
			unpack(popdata[pd].num);
			unpack(popdata[pd].rows);
		}
		for(auto md = 0u; md < margdata.size(); md++){
			unpack(margdata[md].percent);
		}
		unsigned int kmax;
		unpack(kmax);
		genQ.Qten.resize(kmax);
		for(auto k = 0u; k < kmax; k++){
			unpack(genQ.Qten[k].name);
		}
		if(si != packsize()) emsgEC("Data",1);
	}

	for(auto k = 0u; k < genQ.Qten.size(); k++){                                                   // Copies the Q matrices
		auto num = narea*nage;
		MPI_Bcast(&num,1,MPI_UNSIGNED,0,MPI_COMM_WORLD);
		if(core != 0){
			genQ.Qten[k].ntof.resize(num);
			genQ.Qten[k].tof.resize(num);
			genQ.Qten[k].valf.resize(num);
		}
		
		auto vmin = 0u;
		do{
			auto vmax = vmin+1000; if(vmax > num) vmax = num;
			
			if(core == 0){
				packinit(0);
				for(auto v = vmin; v < vmax; v++){
					pack(genQ.Qten[k].ntof[v]);
					pack(genQ.Qten[k].tof[v]);
					pack(genQ.Qten[k].valf[v]);
				}
				si = packsize();
			}

			MPI_Bcast(&si,1,MPI_UNSIGNED,0,MPI_COMM_WORLD);
			if(core != 0){
				packinit(si);
			}
			MPI_Bcast(packbuffer(),si,MPI_DOUBLE,0,MPI_COMM_WORLD);

			if(core != 0){
				for(auto v = vmin; v < vmax; v++){
					unpack(genQ.Qten[k].ntof[v]);
					unpack(genQ.Qten[k].tof[v]);
					unpack(genQ.Qten[k].valf[v]);
				}
				if(si != packsize()) emsgEC("Data",2);
			}
	
			vmin = vmax;
		}while(vmin < num);
	}	
}

/// Removes 'r' and quotations from a string
string DATA::strip(string line) const 
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

void DATA::sortX(vector <unsigned int> &vec){ sort(vec.begin(),vec.end(),compX);}
void DATA::sortY(vector <unsigned int> &vec){ sort(vec.begin(),vec.end(),compY);}

