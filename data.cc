//  The data inputs

#include <iostream>
#include <vector>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <math.h> 
#include <string.h>

using namespace std;

#include "utils.hh"
#include "data.hh"	
#include "consts.hh"
#include "pack.hh"

/// Reads in transition and area data
void DATA::readdata(unsigned int core, unsigned int ncore, unsigned int mod)
{
	unsigned int r, i, c, imax, k, td, pd, md, j, jmax, cc, fl, d, dp, a, aa, q, s, row;
	unsigned int namecol, codecol, xcol, ycol, regcol;
	int dc;
	double v=0, sum;
	string line, ele, name, regcode, st, file;
	REGION reg;
	AREA are;
	DEMOCAT dem;
	IND indi;
	TABLE tab;
	vector <unsigned int> count;
	vector <vector <double> > val;
	vector <double> vec;
	vector <unsigned int> rcol;
	
	mode = mod;
	fepertime = 10;
	
	settpertime = 1;
	nsettime = settpertime*period;
	settime.resize(nsettime+1);
	for(s = 0; s <= nsettime; s++) settime[s] = double(s*period)/nsettime;
			
	fediv = nsettime*fepertime;
	
	if(core == 0){
		count.resize(ndemocat);                                     // Defines all the demographic states
		for(dc = 0; dc < int(ndemocat); dc++) count[dc] = 0;
		
		do{
			democatpos.push_back(count);
			
			dc = ndemocat-1;
			do{
				fl = 0;
				count[dc]++; if(count[dc] == democat[dc].value.size()){ fl = 1; count[dc] = 0; dc--;}
			}while(fl == 1 && dc >= 0);
		}while(dc >= 0);
		ndemocatpos = democatpos.size();
		ndemocatposperage = ndemocatpos/nage;
 
		tab = loadtable(regiondatafile);
		namecol = findcol(tab,"name");
		codecol = findcol(tab,"code");
		
		for(row = 0; row < tab.nrow; row++){
			reg.name = tab.ele[row][namecol];
			reg.code = tab.ele[row][codecol];
			region.push_back(reg);
		}		
		nregion = region.size();
		
		cout << "Region data loaded." << endl;
		if(checkon == 1){
			for(r = 0; r < nregion; r++) cout << region[r].code << ", " << region[r].name  << " regionload" << endl;
		}
	
		file = areadatafile;
		tab = loadtable(file);
		
		for(c = 0; c < tab.ncol; c++) if(tab.heading[c] == "age0-14") break;
		if(c < tab.ncol){
			vector <unsigned int> agecols;
			agecols.push_back(findcol(tab,"age0-14"));
			agecols.push_back(findcol(tab,"age15-44"));
			agecols.push_back(findcol(tab,"age45-64"));
			agecols.push_back(findcol(tab,"age65+"));	
			table_createcol("all",agecols,tab);
		}
		
		codecol = findcol(tab,"area");                                                 // Works out ccolumns for different columns
		xcol = findcol(tab,"easting");
		ycol = findcol(tab,"northing");
		regcol = findcol(tab,"region");
		
		for(j = 0; j < ncovar; j++) covar[j].col = findcol(tab,covar[j].name);

		for(j = 0; j < ndemocat; j++){
			democat[j].col.resize(democat[j].value.size());
			for(k = 0; k < democat[j].col.size(); k++) democat[j].col[k] = findcol(tab,democat[j].value[k]);  
		}

		for(row = 0; row < tab.nrow; row++){
			are.code = tab.ele[row][codecol];
			are.x = atof(tab.ele[row][xcol].c_str());
			are.y = atof(tab.ele[row][ycol].c_str());
	
			regcode = tab.ele[row][regcol];
			r = 0; while(r < nregion && region[r].code != regcode) r++;
			if(r == nregion) emsg("Region code not recognised: "+regcode);
			are.region = r;
					
			are.covar.resize(ncovar);
			for(j = 0; j < ncovar; j++){
				st = tab.ele[row][covar[j].col];
				v = atof(st.c_str());
				if(std::isnan(are.covar[j])) emsg("In file '"+areadatafile+"' the expression '"+st+"' is not a number");	
			
				if(covar[j].func == "log"){
					if(v <= 0) emsg("Log transformed quantities must be positive.");
					are.covar[j] = log(v);
				}
				else{
					if(covar[j].func == "linear") are.covar[j] = v;
					else emsg("The functional relationship '"+covar[j].func+"' is not recognised.");
				}
			}
			
			val.resize(democat.size());
			for(d = 0; d < democat.size(); d++){
				jmax = democat[d].value.size();
				val[d].resize(jmax);
				for(j = 0; j < jmax; j++){
					st = tab.ele[row][democat[d].col[j]];
					val[d][j] = atof(st.c_str());
					if(std::isnan(val[d][j])) emsg("In file '"+areadatafile+"' the expression '"+st+"' is not a number");	
				}
			}
			
			are.agepop.resize(nage);
			for(a = 0; a < nage; a++){
				are.agepop[a] = val[0][a];
			}
				
			are.pop.resize(ndemocatpos);
			for(dp = 0; dp < ndemocatpos; dp++){
				for(j = 0; j < democat.size(); j++){
					if(j == 0) v = val[j][democatpos[dp][j]];
					else v *= val[j][democatpos[dp][j]]/100.0;
				}
				are.pop[dp] = (unsigned int)(v+0.5);
			}

			area.push_back(are);			
		}
		narea = area.size();
		
		for(j = 0; j < ncovar; j++){            // Shifts covariates so average is zero
			sum = 0; for(c = 0; c < narea; c++) sum += area[c].covar[j];
			sum /= narea;
			
			for(c = 0; c < narea; c++) area[c].covar[j] -= sum;
		}		
		
		cout << endl << "Area data loaded." << endl;
		if(checkon == 1){
			for(c = 0; c < narea; c++){
				cout << nregion << " " << area[c].region << "region" << endl;
				cout << area[c].code << " " << region[area[c].region].code << " " <<  area[c].x << " " <<  area[c].y << "  ***";
			
				for(j = 0; j < area[c].pop.size(); j++) cout << area[c].pop[j] << ", ";
				cout << endl;	
			}
		}
		
		if(mode != MODE_SIM){                                                    // Loads transition data for inference
			for(td = 0; td < transdata.size(); td++){
				file = transdata[td].file;
				tab = loadtable(file);
				table_selectdates(transdata[td].start,transdata[td].units,tab,"trans");
				
				rcol.clear();
				if(transdata[td].type == "reg"){	for(k = 0; k < region.size(); k++) rcol.push_back(findcol(tab,region[k].code));}
				else{ rcol.push_back(findcol(tab,"all"));}
				
				transdata[td].num.resize(rcol.size());
				for(r = 0; r < rcol.size(); r++){
					transdata[td].num[r].resize(tab.nrow);
					for(row = 0; row < tab.nrow; row++){	
						transdata[td].num[r][row] = getint(tab.ele[row][rcol[r]],file);
					}
				}
				transdata[td].rows = tab.nrow;
				
				if(transdata[td].start + (tab.nrow-1)*transdata[td].units > period){
					emsg("The file '"+file+"' has more data than will fit in the time period.");
				}
			}
		}
		
		if(mode != MODE_SIM){                                                    // Loads population data for inference
			for(pd = 0; pd < popdata.size(); pd++){
				file = popdata[pd].file;
				tab = loadtable(file);
				table_selectdates(popdata[pd].start,popdata[pd].units,tab,"pop");
			
				rcol.clear();
				if(popdata[pd].type == "reg"){	for(k = 0; k < region.size(); k++) rcol.push_back(findcol(tab,region[k].code));}
				else{ rcol.push_back(findcol(tab,"all"));}
				
				popdata[pd].num.resize(rcol.size());
				for(r = 0; r < rcol.size(); r++){
					popdata[pd].num[r].resize(tab.nrow);
					for(row = 0; row < tab.nrow; row++) popdata[pd].num[r][row] = getint(tab.ele[row][rcol[r]],file);
				}
				popdata[pd].rows = tab.nrow;
				
				if(popdata[pd].start + (tab.nrow-1)*popdata[pd].units > period){
					emsg("The file '"+file+"' has more data than will fit in the time period.");
				}
			}
		}
		
		if(mode != MODE_SIM){                                                    // Loads marginal data for inference
			for(md = 0; md < margdata.size(); md++){
				file = margdata[md].file;
				tab = loadtable(file);
	
				rcol.clear();
				if(margdata[md].type == "reg"){	for(k = 0; k < region.size(); k++) rcol.push_back(findcol(tab,region[k].code));}
				else{ rcol.push_back(findcol(tab,"all"));}
				
				margdata[md].percent.resize(rcol.size());
				for(r = 0; r < rcol.size(); r++){
					margdata[md].percent[r].resize(tab.nrow);
					for(row = 0; row < tab.nrow; row++){	
						margdata[md].percent[r][row] = atof(tab.ele[row][rcol[r]].c_str());
					}
				}
			}
		}

		vec.resize(nage);                                                 // Reads in Q tensors
		for(q = 0; q < Q.size(); q++){
			tab = loadtable(Q[q].file);
			
			Q[q].to.resize(narea*nage); Q[q].val.resize(narea*nage);
			for(row = 0; row < tab.nrow; row++){
				c = atoi(tab.ele[row][0].c_str()); cc = atoi(tab.ele[row][1].c_str());
				
				for(a = 0; a < nage; a++){
					for(aa = 0; aa < nage; aa++) vec[aa] = atof(tab.ele[row][2+a*nage+aa].c_str()); 
					
					Q[q].to[c*nage+a].push_back(cc);
					Q[q].val[c*nage+a].push_back(vec);
					
					if(c < cc){
						Q[q].to[cc*nage+a].push_back(c);
						Q[q].val[cc*nage+a].push_back(vec);
					}
				}
			}
			
			normaliseQ(q);
		}
	}

	if(ncore > 1) copydata(core);

	agedist.resize(nage); for(a = 0; a < nage; a++) agedist[a] = 0;
	
	for(c = 0; c < narea; c++){                                              // Adds individuals to the system
		area[c].ind.resize(ndemocatpos);
		for(dp = 0; dp < ndemocatpos; dp++){
			imax = area[c].pop[dp];
			for(i = 0; i < imax; i++){
				area[c].ind[dp].push_back(ind.size());
				
				indi.area = c;
				indi.dp = dp;
				ind.push_back(indi);
			}
			a = democatpos[dp][0];
			agedist[a] += imax;
		}
	}
	popsize = ind.size();
	for(a = 0; a < nage; a++) agedist[a] /= popsize;
	
	narage = narea*nage;                                              // Generates the mixing matrix between ages/areas
	nardp = narea*ndemocatpos; 
	nsettardp = nsettime*nardp;
	
	//plotrawdata(); emsg("done");
	//generatedeathdata(); emsg("done");
}


/// Gets a positive integer from a string
unsigned int DATA::getint(string st, string file)
{
	unsigned int i, j;
	int n = st.length();  
  char str[n + 1]; 
	strcpy(str, st.c_str()); 
		
	for(j = 0; j < st.size(); j++) if(!isdigit(str[j])) break;
	if(j < st.size()){
		if(st == "NA") i = UNKNOWN;
		else{
			if(st == "*"){
				i = THRESH;
				if(threshold == UNSET) emsg("If 'NA' is used, there must be a threshold set with the 'threshold' command.");
			}
			else emsg("In file '"+file+"' the quantity '"+st+"' is not a number");
		}
	}
  else i = atoi(str);
	return i;
}

/// Normalises the Q tensors
void DATA::normaliseQ(unsigned int q)
{
	unsigned int c, a, cc, aa, vi, j, jmax;
	double sum, sum2;
	
	for(c = 0; c < narea; c++){
		sum = 0; sum2 = 0;
		for(a = 0; a < nage; a++){
			sum2 += area[c].agepop[a];
			vi = c*nage + a;
		
			jmax = Q[q].to[vi].size();
			for(j = 0; j < jmax; j++){
				cc = Q[q].to[vi][j];
				for(aa = 0; aa < nage; aa++) sum += area[c].agepop[a]*Q[q].val[vi][j][aa]*area[cc].agepop[aa];
			}
		}
		sum /= sum2;

		for(a = 0; a < nage; a++){
			vi = c*nage + a;
			jmax = Q[q].to[vi].size();
			for(j = 0; j < jmax; j++){
				for(aa = 0; aa < nage; aa++) Q[q].val[vi][j][aa] /= sum;
			}
		}
	}
}

/// Adds a time period
void DATA::addtimep(string name, double tend)
{
	TIMEP timep;
	
	timep.name = name;
	timep.tend = tend;
	timeperiod.push_back(timep);
}
	
/// Adds a Q tensor
void DATA::addQtensor(string timep, string comp, string file)
{
	unsigned int tp;
	QTENSOR qten;
	
	tp = 0; while(tp < timeperiod.size() && timeperiod[tp].name != timep) tp++;
	if(tp == timeperiod.size()) emsg("Cannot find '"+timep+"' as a time period");
	
	qten.timep = tp;
	qten.comp = comp;
	qten.file = file;
	Q.push_back(qten);
}

/// Loads a table from a file
TABLE DATA::loadtable(string file)
{
	TABLE tab;
	string line, st;
	vector <string> vec;
	ifstream in;
	
	in.open((datadir+"/"+file).c_str());
	if(!in){
		in.open((outputdir+"/"+file).c_str());
		if(!in) emsg("Cannot open the file '"+file+"'");
	}
	
	tab.file = file;
	
	getline(in,line);

	stringstream ss(line);
	do{
		getline(ss,st,'\t'); st = strip(st);
		tab.heading.push_back(st);
		if(ss.eof()) break;
	}while(1 == 1);
	tab.ncol = tab.heading.size();
	
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
		if(vec.size() != tab.ncol) emsg("Rows in file '"+file+"' do not all have the same length.");
		
		tab.ele.push_back(vec);
	}while(1 == 1);
	tab.nrow = tab.ele.size();
	
	return tab;
}

/// Creates a new column by adding together existing columns		
void DATA::table_createcol(string head,vector <unsigned int> cols, TABLE &tab)
{
	unsigned int row, i, sum;
	
	tab.heading.push_back(head);
	for(row = 0; row < tab.nrow; row++){
		sum = 0; for(i = 0; i < cols.size(); i++) sum += atoi(tab.ele[row][cols[i]].c_str());
		stringstream ss; ss << sum;
		tab.ele[row].push_back(ss.str());
	}
	tab.ncol++;
}

/// Selects dates as specified in the TOLM file
void DATA::table_selectdates(unsigned int t, unsigned int units, TABLE &tab, string type)
{
	unsigned int row, tt, num1, num2, i;
	
	row = 0;
	while(row < tab.nrow){
		tt = gettime(tab.ele[row][0]) - start;
		if(tt < t){
			if(type == "trans" && row > 0){ // In the case of transitions adds up the contributions from other days to make a week 
				for(i = 1; i < tab.ncol; i++){
					num1 = getint(tab.ele[row-1][i],tab.file);
					num2 = getint(tab.ele[row][i],tab.file);
					if(num1 == THRESH || num2 == THRESH || num1 == UNKNOWN || num2 == UNKNOWN) emsg("not done yet");
					tab.ele[row-1][i] = to_string(num1+num2);
				}				
			}
			tab.ele.erase(tab.ele.begin()+row);
			tab.nrow--;
		}
		else{
			if(tt > t) emsg("In file '"+tab.file+"' there is no observed data at time '"+getdate(t)+"'."); 
			t += units;
			row++;
		}
	}
	if(tab.nrow == 0) emsg("The file '"+tab.file+"' does not contain any informations.");
	
	if(type=="trans"){                         // Removes the last line if incomplete
		if(gettime(tab.ele[tab.nrow-1][0]) > end-units){
			tab.ele.erase(tab.ele.begin()+tab.nrow-1);
			tab.nrow--;
		}
	}	
	
	if(checkon == 1){
		for(row = 0; row < tab.nrow; row++){
			for(i = 0; i < tab.ncol; i++) cout << tab.ele[row][i] << " ";
			cout << endl;
		}
	}
}				

/// Finds a column in a table
unsigned int DATA::findcol(TABLE &tab, string name)
{
	unsigned int c;
	
	for(c = 0; c < tab.ncol; c++) if(tab.heading[c] == name) break;
	if(c == tab.ncol) emsg("Cannot find the column heading '"+name+"' in file '"+tab.file+"'.");
	return c;
}		
			
/// Copies data from core zero to all the others
void DATA::copydata(unsigned int core)
{
	unsigned int td, pd, md, q, v, vmin, vmax, num;
	int si;

	if(core == 0){                                  				   // Copies the above information to all the other cores
		packinit();
		pack(ndemocatpos);
		pack(democatpos);
		pack(nregion);
		pack(region);
		pack(narea);
		pack(area);
		pack(nage);
		pack(ndemocatposperage);
		for(td = 0; td < transdata.size(); td++){
			pack(transdata[td].num);
			pack(transdata[td].rows);
		}
		for(pd = 0; pd < popdata.size(); pd++){
			pack(popdata[pd].num);
			pack(popdata[pd].rows);
		}
		for(md = 0; md < margdata.size(); md++){
			pack(margdata[md].percent);
		}
		si = packsize();
	}

	MPI_Bcast(&si,1,MPI_UNSIGNED,0,MPI_COMM_WORLD);
	MPI_Bcast(packbuffer(),si,MPI_DOUBLE,0,MPI_COMM_WORLD);
		
	if(core != 0){
		packinit();
		unpack(ndemocatpos);
		unpack(democatpos);
		unpack(nregion);
		unpack(region);
		unpack(narea);
		unpack(area);
		unpack(nage);
		unpack(ndemocatposperage);
		for(td = 0; td < transdata.size(); td++){
			unpack(transdata[td].num);
			unpack(transdata[td].rows);
		}
		for(pd = 0; pd < popdata.size(); pd++){
			unpack(popdata[pd].num);
			unpack(popdata[pd].rows);
		}
		for(md = 0; md < margdata.size(); md++){
			unpack(margdata[md].percent);
		}
		if(si != packsize()) emsg("Data: EC9");
	}
	
	for(q = 0; q < Q.size(); q++){                                                   // Copies the Q matrices
		num = Q[q].to.size();
		MPI_Bcast(&num,1,MPI_UNSIGNED,0,MPI_COMM_WORLD);
		if(core != 0){ Q[q].to.resize(num); Q[q].val.resize(num);}
		 
		vmin = 0;
		do{
			vmax = vmin+1000; if(vmax > num) vmax = num;
			
			if(core == 0){
				packinit();
				for(v = vmin; v < vmax; v++){
					pack(Q[q].to[v]);
					pack(Q[q].val[v]);
				}
				si = packsize();
			}
				
			MPI_Bcast(&si,1,MPI_UNSIGNED,0,MPI_COMM_WORLD);
			MPI_Bcast(packbuffer(),si,MPI_DOUBLE,0,MPI_COMM_WORLD);
	
			if(core != 0){
				packinit();	
				for(v = vmin; v < vmax; v++){
					unpack(Q[q].to[v]);
					unpack(Q[q].val[v]);
				}
				if(si != packsize()) emsg("Data: EC9");
			}
		
			vmin = vmax;
		}while(vmin < num);
	}		
}

/// Adds demographic categories
void DATA::adddemocat(string name, vector <string> &st, vector <string> &params)
{
	DEMOCAT dem;
	
	dem.name = name;	
	dem.value = st;
	dem.param = params;
	
	democat.push_back(dem);

	ndemocat = democat.size();
	if(ndemocat == 1) nage = st.size();
}
	
/// Add a covariate for the areas
void DATA::addcovar(string name, string param, string func)
{
	COVAR cov;
	cov.name = name;
	cov.param = param;
	cov.func = func;
	
	covar.push_back(cov);
	ncovar = covar.size();
}
	
string DATA::strip(string line)
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

/// Gets the time from a string
unsigned int DATA::gettime(string st)
{
	unsigned int t;
	const char *buf = st.c_str();
	struct tm result;
			
	if(st == "start") return start;
	if(st == "end") return end;
			
	switch(tform){
	case TFORM_NUM:
		t = atoi(buf);
		if(std::isnan(t)) emsg("Time '"+st+"' is not a number");
		break;

	case TFORM_YMD:
		memset(&result, 0, sizeof(result));
		if(strptime(buf,"%Y-%m-%d",&result) != NULL){
			time_t tt = mktime(&result);
			t = tt/(60*60*24);
		}
		else{ 
			emsg("'"+st+"' not regonised as Year-Month-Day format.");
			t = 0;
		}
		break;
		
	default:
		emsg("Do not recognise time format");
		t = 0;
		break;
	}

	return t;	
}

/// Returns a date from a time
string DATA::getdate(unsigned int t)
{
	time_t tt;
	string st;
	stringstream ss; 
	char buffer[80];
  struct tm *timeinfo;
	
	t += start;

	switch(tform){
	case TFORM_NUM:
		ss << t;
		break;
		
	case TFORM_YMD:
		tt = (t + 0.5)*(60*60*24);	
		timeinfo = localtime(&tt);
		strftime(buffer,80,"%Y-%m-%d",timeinfo);
		ss << buffer;
		break;
	}
	
	return ss.str();
	
}

void DATA::sortX(vector <unsigned int> &vec){ sort(vec.begin(),vec.end(),compX);}
void DATA::sortY(vector <unsigned int> &vec){ sort(vec.begin(),vec.end(),compY);}


/************ Code below this is used for diagnostic purposes and will not be in the final version ***********/

/// This is used to plot raw data files (this is used for diagnoistic purposed, and not in the analysis)
void DATA::plotrawdata()
{
	unsigned int row, tt, num, sum;
	TABLE tab;
	
	
	tab = loadtable("DailyDeathsConfirmedCovid.txt");
	ofstream dout(outputdir+"/deathraw.txt");
	sum = 0;
	for(row = 0; row < tab.nrow; row++){
		tt = gettime(tab.ele[row][0]) - start;
		num = getint(tab.ele[row][1],"");
		if(num != UNKNOWN) dout << tt << " " << num << endl;
		if(num != UNKNOWN) sum += num;
	}
	cout << sum << " Deaths" << endl;
	
	tab = loadtable("hospital_admissions_number_per_day.txt");
	ofstream dout2(outputdir+"/hospadminraw.txt");
	sum = 0;
	for(row = 0; row < tab.nrow; row++){
		tt = gettime(tab.ele[row][0]) - start;
		num = getint(tab.ele[row][1],"");
		if(num != UNKNOWN) dout2 << tt << " " << num << endl;
		if(num != UNKNOWN) sum += num;
	}
	cout << sum << " Admissions" << endl;
		
	tab = loadtable("NHS_and_UKG_national_daily_confirmed_cases.txt");
	ofstream dout3(outputdir+"/nhsothercasesraw.txt");
	sum = 0;
	for(row = 0; row < tab.nrow; row++){
		tt = gettime(tab.ele[row][0]) - start;
		num = getint(tab.ele[row][1],"");
		if(num != UNKNOWN) dout3 << tt << " " << num << endl;
		if(num != UNKNOWN) sum += num;
	}
	cout << sum << " NHS+other cases" << endl;
	
	tab = loadtable("NHS_only_national_daily_confirmed_cases.txt");
	ofstream dout4(outputdir+"/nhsonlycasesraw.txt");
	sum = 0;
	for(row = 0; row < tab.nrow; row++){
		tt = gettime(tab.ele[row][0]) - start;
		num = getint(tab.ele[row][1],"");
		if(num != UNKNOWN) dout4 << tt << " " << num << endl;
		if(num != UNKNOWN) sum += num;
	}
	cout << sum << " NHS cases" << endl;
	
	tab = loadtable("HD_NRS.txt");
	ofstream dout5(outputdir+"/HDraw.txt");
	sum = 0;
	for(row = 0; row < tab.nrow; row++){
		tt = gettime(tab.ele[row][0]) - start;
		num = getint(tab.ele[row][1],"");
		if(num != UNKNOWN) dout5 << tt << " " <<  double(num)/7 << endl;
		if(num != UNKNOWN) sum += num;
	}
	cout << sum << "HD Deaths" << endl;
	
	tab = loadtable("ID_NRS.txt");
	ofstream dout6(outputdir+"/IDraw.txt");
	sum = 0;
	for(row = 0; row < tab.nrow; row++){
		tt = gettime(tab.ele[row][0]) - start;
		num = getint(tab.ele[row][1],"");
		if(num != UNKNOWN) dout6 << tt << " " << double(num)/7 << endl;
		if(num != UNKNOWN) sum += num;
	}
	cout << sum << "HI Deaths" << endl;
	
	
	// estimate distribution
	
	vector <unsigned int> numst, numst2;
	int t, tmax = 140, j;
	double dt, mean_ns, sd_ns, sd, mean;
	numst.resize(tmax); numst2.resize(tmax);
	for(t = 0; t < tmax; t++){ numst[t] = 0, numst[t] = 0;}
	
	tab = loadtable("hospital_admissions_number_per_day.txt");
	sum = 0;
	for(row = 0; row < tab.nrow; row++){
		t = gettime(tab.ele[row][0]) - start;
		num = getint(tab.ele[row][1],"");
		if(num != UNKNOWN){
			for(j = 0; j < num; j++){
				mean_ns = 20; sd_ns = 20;
				sd = sqrt(log((1+sd_ns*sd_ns/(mean_ns*mean_ns)))); mean = log(mean_ns) - sd*sd/2;
				dt = exp(normal(mean,sd));
				tt = int(t + dt+ran());
				if(tt < tmax) numst[tt]++;
				tt = int(dt+ran());
				if(tt < tmax) numst2[tt]++;
			}
		}
	}
	
	ofstream deathpred(outputdir+"/deathpred.txt");
	for(t = 0; t < tmax; t++){
		deathpred << t << " " << numst[t] <<  " " << numst2[t] << endl;
	}
}
	
/// This is used to plot raw data files (this is used for diagnoistic purposed, and not in the analysis)
void DATA::generatedeathdata()
{
	unsigned int row, t, nweek = 40, w, ww, r, n, sum, sumH, sumI;
	vector< vector <unsigned int> > num;
	vector <unsigned int> HDnum;
	vector <unsigned int> IDnum;
	vector <string> linedate;
	TABLE tab;
	string reg, date, file, sex, age, cause, loc;
	
	tab = loadtable("deathdata.txt");
	
	num.resize(nweek); linedate.resize(nweek); HDnum.resize(nweek); IDnum.resize(nweek);
	for(w = 0; w < nweek; w++){
		HDnum[w] = 0; IDnum[w] = 0;
				 
		num[w].resize(region.size());
		for(r = 0; r < region.size(); r++) num[w][r] = 0;
	}

	sum = 0; sumI = 0; sumH = 0;
	for(row = 0; row < tab.nrow; row++){
		reg = tab.ele[row][0];
		date = tab.ele[row][1];
		n = getint(tab.ele[row][4],"");
		sex = tab.ele[row][5];
		age = tab.ele[row][6];
		cause = tab.ele[row][7];
		loc = tab.ele[row][8];
		
		if(date.substr(0,4) == "w/c "){
			date = date.substr(4,date.length()-4);
			t = gettime(date);
			if(t >=start){
				w = (t-start)/7; if(w > nweek) emsg("out");
						
				if(sex == "All" && age == "All" && cause == "COVID-19 related"){
					r = 0; while(r < region.size() && region[r].code != reg) r++;
					if(r < region.size()){
						if(loc == "All"){				
							linedate[w] = date;
							num[w][r] += n;
							sum += n;
						}
					}
					else{
						if(reg == "S92000003" && loc != "All"){
							if(loc == "Hospital"){ HDnum[w] += n; sumH += n;}
							else{ IDnum[w] += n; sumI += n;}
						}
					}
				}
			}
		}
	}
	cout << "Total in regions: " << sum << endl;
	cout << "Total from hospitals: " << sumH << endl;
	cout << "Total from community: " << sumI << endl;

	file = outputdir+"/D_NRS_reg.txt";
	ofstream deathout(file.c_str());
	deathout << "date"; for(r = 0; r < region.size(); r++) deathout << "\t" << region[r].code;
	deathout << endl;
	for(w = 0; w < nweek; w++){
		deathout << linedate[w];
		for(r = 0; r < region.size(); r++){
			sum = 0; for(ww = 0; ww < w; ww++) sum += num[ww][r];
			deathout  << "\t" << sum;
		}
		deathout << endl;
	}
	
	file = outputdir+"/HD_NRS.txt";
	ofstream HDout(file.c_str());
	HDout << "date\tall" << endl;
	for(w = 0; w < nweek; w++){
		HDout << linedate[w] << "\t" << HDnum[w] << endl;
	}

	file = outputdir+"/ID_NRS.txt";
	ofstream IDout(file.c_str());
	IDout << "date\tall" << endl;
	for(w = 0; w < nweek; w++){
		IDout << linedate[w] << "\t" << IDnum[w] << endl;
	}
}
