///////////////////////////RAW DATA ANALYSIS //////////////////////////

#include <iostream>
#include <vector>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <math.h> 
#include <iomanip>
#include "json.hpp"

using namespace std;
using json = nlohmann::json;

#include "data.hh"	

/// Selects which raw data analysis to perform
void Data::raw()
{
	//fractional_change("IR_England_age");
	//fractional_change("ca_England_age");
	//fractional_change("IH_England_age");
	//fractional_change("D_England_wales_age"); emsg("done");
	//fractional_change2("HD_England_age_data"); emsg("done");
	
	//calculate_IR();  emsg("done");
	//cases_age(); emsg("done");
	//split_deaths_age_data(); emsg("done");
	//deaths_age(); emsg("done");
	//generate_admission_age(); emsg("done");
	//generatehospdata(); emsg("done");
	//generatedeathdata_england_wales(); emsg("generate death data");
	//generatedeathdata_scotland(); emsg("done");
	//death_age_distribution(); emsg("done");
	//case_age_distribution(); emsg("done");
	
	//generatedeathdata_england_wales(); emsg("done");
	//plotrawdata(); emsg("done");
	//IFR(); emsg("done");
	//deaths_hospital_england(); emsg("done");
	//generate_initpop(); emsg("done");
	//generate_tvcovar(); emsg("done");
	//generate_level_effect(); emsg("done");
	//reducesize_geojson("areas2.geojson"); emsg("done");
	
	//generate_age_mixing_perturb(); emsg("done");
	//generate_geo_mixing_matrix(); emsg("done");
}

void Data::deaths_hospital_england()
{
	Table tab = load_table("deaths_hospital_England_age.txt");

	auto files = details.output_directory+"/HD_England_age_data.txt";
	ofstream IGageout(files.c_str());
	
	
	IGageout << "date\tday\t0-19\t20-39\t40-59\t60-79\t80+\ttotal\n";
	
	for(auto col = 1u; col < tab.ncol; col++){
		auto t = details.gettime(tab.heading[col],"Error") - details.start;
		
		IGageout << tab.heading[col] << "\t" << t;
		auto sum = 0u;
		for(auto row = 0u; row < tab.nrow; row++){
			IGageout << "\t" << tab.ele[row][col];
			sum += get_int( tab.ele[row][col],"erro");
		}
		IGageout << "\t" << sum << "\n";
	}
}

void Data::case_age_distribution()
{
	//Table tab = load_table("ca_England_age.txt");
	Table tab = load_table("IH_England_age.txt");
		
	const unsigned int nage = 4;
	vector <unsigned int> num(nage);
	for(auto a = 0u; a < nage; a++) num[a] = 0;
		
	for(auto row = 0u; row < tab.nrow; row++){
		if(get_int(tab.ele[row][1],"error") < 217){
			cout << tab.ele[row][0] << "\n";
			for(auto a = 0u; a < nage; a++){
				num[a] += get_int(tab.ele[row][2+a],"error");
			}
		}
	}
	for(auto a = 0u; a < nage; a++) cout << a << " " << num[a] << " age dist\n";
}

void Data::death_age_distribution()
{	
	Table tabHD = load_table("HD_England_wales_age.txt");
	
	Table tabCD = load_table("CD_England_wales_age.txt");

	const unsigned int nage = 4;
	vector <unsigned int> num(nage);
	for(auto a = 0u; a < nage; a++) num[a] = 0;
		
	for(auto row = 0u; row < tabHD.nrow; row++){
		if(get_int(tabHD.ele[row][1],"error") < 217){
			cout << tabHD.ele[row][0] << "\n";
			for(auto a = 0u; a < nage; a++){
				num[a] += get_int(tabHD.ele[row][2+a],"error");
				num[a] += get_int(tabCD.ele[row][2+a],"error");
			}
		}
	}
	for(auto a = 0u; a < nage; a++) cout << a << " " << num[a] << " age dist\n";
}

/// Finds the infection fatality rate for a given set of ages
void Data::IFR()	
{
	const int nage = 19;
	
	// From https://www.imperial.ac.uk/media/imperial-college/medicine/mrc-gida/2020-10-29-COVID19-Report-34.pdf
	double IFR[nage] = {0.0,0.01,0.01,0.02,0.02,0.04,0.06,0.09,0.15,0.23,0.36,0.57,0.89,1.39,2.17,3.39,5.3,9.28,16.19};
	
	// From https:www.ons.gov.uk/peoplepopulationandcommunity/populationandmigration/populationestimates/datasets/populationestimatesforukenglandandwalesscotlandandnorthernireland
	
	double pop[nage] = {3465179,3721990,3535065,3262613,3690265,4009669,4000908,3918562,3583853,3914884,4127941,3888131,3304688,2978882,2958612,2067471,1529682,933656,547789};
	
	double sum, sum2;
	
	sum = 0.0; sum2 = 0.0; for(auto i = 0u; i < 4; i++){ sum += IFR[i]*pop[i]; sum2 += pop[i];}	
	cout << sum/sum2 << " ifr\n";
	
	sum = 0.0; sum2 = 0.0; for(auto i = 4u; i < 13; i++){ sum += IFR[i]*pop[i]; sum2 += pop[i];}	
	cout << sum/sum2 << " ifr\n";
	
	sum = 0.0; sum2 = 0.0; for(auto i = 13u; i < 17; i++){ sum += IFR[i]*pop[i]; sum2 += pop[i];}	
	cout << sum/sum2 << " ifr\n";
	
	sum = 0.0; sum2 = 0.0; for(auto i = 17u; i < 19; i++){ sum += IFR[i]*pop[i]; sum2 += pop[i];}	
	cout << sum/sum2 << " ifr\n";
	
	double inf[4];
	inf[0] = 100.0*15/0.00985515;
	inf[1] = 100.0*5555/0.258672;
	inf[2] = 100.0*24541/2.69301;
	inf[3] = 100.0*21984/11.8351;
	
	double suminf = inf[0]+inf[1]+inf[2]+inf[3];
	double sumdeath = 15+ 5555+24541+ 21984;
	//double fac = sumdeath/(0.0115*suminf);
	
	cout << "\n" << sumdeath/suminf  << " IFR\n";
	
	for(auto a = 0u; a < 4; a++){
		cout << inf[a] << " inf\n";
	}
	
	//0.00985515 ifr
//0.258672 ifr
//2.69301 ifr
//11.8351 ifr

//0 15 age dist
//1 5555 age dist
//2 24541 age dist
//3 21984 age dist

}

// Estimates the weekly number of infected recovery (i.e. not counting those hospitalised
void Data::calculate_IR()
{
	Table cases = load_table("ca_England_age.txt");
	
	Table admi = load_table("IH_England_age.txt");
	
	auto files = details.output_directory+"/IR_England_age.txt";
	ofstream IGageout(files.c_str());
	
	auto rowst = 0u; while(rowst < cases.nrow && cases.ele[rowst][0] != admi.ele[0][0]) rowst++;
	cout << rowst << " " << cases.nrow << "rowst\n";
	
	IGageout << "date\tday\t0-19\t20-64\t65-84\t85+\ttotal\n";
	
	for(auto row = 0u; row < admi.nrow-6; row += 7){
		IGageout << admi.ele[row][0] << "\t" << admi.ele[row][1];
		for(auto col = 2u; col < admi.ncol; col++){
			int sum = 0;
			for(auto day = 0u; day < 7; day++){
				sum += get_int(cases.ele[rowst+row+day][col],"error") - get_int(admi.ele[row+day][col],"error");
			}
			if(sum < 0) sum = 0;
			IGageout << "\t" << sum;
		}
		IGageout << "\n";
	}
}


/// Calculates fractional change as a function of time
void Data::fractional_change(string file)
{
	Table tab = load_table(file+".txt");
	
	auto files = details.output_directory+"/"+file+"_frac.txt";
	ofstream IGageout(files.c_str());
	
	IGageout << "date\tday\t0\t20\t65\t85\ttotal\n";
	
	for(auto row = 0u; row < tab.nrow; row++){
		double valst[4];
		double sum = 0;
		for(auto col = 0; col < 4; col++){
			valst[col] = sum;
			sum += get_int(tab.ele[row][col+2],"error");
		}
		
		if(sum > 10){
			IGageout << tab.ele[row][0] << "\t" << tab.ele[row][1];
			for(auto col = 0; col < 4; col++) IGageout << "\t" << (valst[col]/sum);
			IGageout << "\t" << tab.ele[row][6];
			IGageout << "\n";
		}
	}
}

/// Calculates fractional change as a function of time
void Data::fractional_change2(string file)
{
	Table tab = load_table(file+".txt");
	
	auto files = details.output_directory+"/"+file+"_frac.txt";
	ofstream IGageout(files.c_str());
	
	IGageout << "date\tday\t0\t20\t40\t60\t80\ttotal\n";
	
	for(auto row = 0u; row < tab.nrow; row++){
		double valst[5];
		double sum = 0;
		for(auto col = 0; col < 5; col++){
			valst[col] = sum;
			sum += get_int(tab.ele[row][col+2],"error");
		}
		
		if(sum > 10){
			IGageout << tab.ele[row][0] << "\t" << tab.ele[row][1];
			for(auto col = 0; col < 5; col++) IGageout << "\t" << (valst[col]/sum);
			IGageout << "\t" << tab.ele[row][7];
			IGageout << "\n";
		}
	}
}

void Data::cases_age()
{
	Table tab = load_table("cases_England_age.txt");
	
	auto files = details.output_directory+"/ca_England_age.txt";
	ofstream IGageout(files.c_str());
	
	IGageout << "date\tday\t0-19\t20-64\t65-84\t85+\ttotal\n";
	
	auto typecol = find_column(tab,"areaType");
	auto datecol = find_column(tab,"date");
	
	double sum;
	double sum0to19 = 0u, sum20to64 = 0u, sum65to84 = 0u, sum85 = 0u;
	for(auto row = 0u; row < tab.nrow; row++){
		if(tab.ele[row][typecol] == "overview"){
			auto t = details.gettime(tab.ele[row][datecol],"Error") - details.start;
			
			IGageout << tab.ele[row][datecol] << "\t" << t;
			auto sumtot = 0.0;	
			sum = 0; 
			for(auto col = 0u; col < tab.ncol; col++){
				if(tab.heading[col] == "newCasesBySpecimenDate-0_4") sum += get_int(tab.ele[row][col],"error");
				if(tab.heading[col] == "newCasesBySpecimenDate-5_9") sum += get_int(tab.ele[row][col],"error");
				if(tab.heading[col] == "newCasesBySpecimenDate-10_14") sum += get_int(tab.ele[row][col],"error");
				if(tab.heading[col] == "newCasesBySpecimenDate-15_19") sum += get_int(tab.ele[row][col],"error");
			}
			IGageout << "\t" << sum;
			sum0to19 += sum;
			sumtot += sum;
			
			sum = 0; 
			for(auto col = 0u; col < tab.ncol; col++){
				if(tab.heading[col] == "newCasesBySpecimenDate-20_24") sum += get_int(tab.ele[row][col],"error");
				if(tab.heading[col] == "newCasesBySpecimenDate-25_29") sum += get_int(tab.ele[row][col],"error");
				if(tab.heading[col] == "newCasesBySpecimenDate-30_34") sum += get_int(tab.ele[row][col],"error");
				if(tab.heading[col] == "newCasesBySpecimenDate-35_39") sum += get_int(tab.ele[row][col],"error");
				if(tab.heading[col] == "newCasesBySpecimenDate-40_44") sum += get_int(tab.ele[row][col],"error");
				if(tab.heading[col] == "newCasesBySpecimenDate-45_49") sum += get_int(tab.ele[row][col],"error");
				if(tab.heading[col] == "newCasesBySpecimenDate-50_54") sum += get_int(tab.ele[row][col],"error");
				if(tab.heading[col] == "newCasesBySpecimenDate-55_59") sum += get_int(tab.ele[row][col],"error");
				if(tab.heading[col] == "newCasesBySpecimenDate-60_64") sum += get_int(tab.ele[row][col],"error");
			}
			IGageout << "\t" << sum;
			sum20to64 += sum;
			sumtot += sum;
			
			sum = 0; 
			for(auto col = 0u; col < tab.ncol; col++){
				if(tab.heading[col] == "newCasesBySpecimenDate-65_69") sum += get_int(tab.ele[row][col],"error");
				if(tab.heading[col] == "newCasesBySpecimenDate-70_74") sum += get_int(tab.ele[row][col],"error");
				if(tab.heading[col] == "newCasesBySpecimenDate-75_79") sum += get_int(tab.ele[row][col],"error");
				if(tab.heading[col] == "newCasesBySpecimenDate-80_84") sum += get_int(tab.ele[row][col],"error");
			}
			IGageout << "\t" << sum;
			sum65to84 += sum;
			sumtot += sum;
			
			sum = 0; 
			for(auto col = 0u; col < tab.ncol; col++){
				if(tab.heading[col] == "newCasesBySpecimenDate-85_89") sum += get_int(tab.ele[row][col],"error");
				if(tab.heading[col] == "newCasesBySpecimenDate-90+") sum += get_int(tab.ele[row][col],"error");
			}
			IGageout << "\t" << sum;
			sum85 += sum;
			sumtot += sum;
			
			IGageout << "\t" << sumtot;
			
			IGageout << "\n";
		}
	}
	
	double total = sum0to19 + sum20to64 + sum65to84 + sum85;
  
	auto fac= 1.0;
	cout << fac*sum0to19 / total << " " << fac*sum20to64 / total <<  " " << fac*sum65to84 / total << " " << fac*sum85 / total << "dist\n";

	
	fac= 3000000;
	cout << fac*sum0to19 / total << " " << fac*sum20to64 / total <<  " " << fac*sum65to84 / total << " " << fac*sum85 / total << "dist\n";
	
}


string Data::remove_comma(string st)
{
	auto i = 0u;
	while(i < st.length()){
		if(st.substr(i,1) ==",") st = st.substr(0,i)+st.substr(i+1,st.length()-(i+1));
		else i++;
	}
	return st;
}

void Data::deaths_age()
{
	Table tab = load_table("England_wales_deaths_age_occurance.txt");
	//Table tab = load_table("England_wales_deaths_age.txt");
	
	//auto files = details.output_directory+"/D_England_wales_age.txt";
	auto files = details.output_directory+"/D_England_wales_age_occurance.txt";
	ofstream IGageout(files.c_str());
	
	IGageout << "date\tday\t0-19\t20-64\t65-84\t85+\ttotal\tcum total\n";
	
	unsigned int sum, cum = 0;
	for(auto col = 1u; col < tab.ncol; col++){	
		IGageout << details.getdate((col-1)*7) << "\t" << (col-1)*7;
		sum = 0; for(auto row = 0u; row < 5; row++) sum += get_int(remove_comma(tab.ele[row][col]),"error");
		IGageout << "\t" << sum;
		
		sum = 0; for(auto row = 5u; row < 14; row++) sum += get_int(remove_comma(tab.ele[row][col]),"error");
		IGageout << "\t" << sum;
		
		sum = 0; for(auto row = 14u; row < 18; row++) sum += get_int(remove_comma(tab.ele[row][col]),"error");
		IGageout << "\t" << sum;
		
		sum = 0; for(auto row = 18u; row < 20; row++) sum += get_int(remove_comma(tab.ele[row][col]),"error");
		IGageout << "\t" << sum;
		
		auto sumtot = 0u; for(auto row = 0u; row < 20; row++) sumtot += get_int(remove_comma(tab.ele[row][col]),"error");
		IGageout << "\t" << sumtot;
		
		cum += sumtot;
		IGageout << "\t" << cum;
		
		IGageout << "\n";
	}	
	cout << tab.nrow << "nrow\n";
}


/// Splits up the deaths age data into community 
void Data::split_deaths_age_data()
{
	Table CDtab = load_table("CD_region_england_wales.txt");
	auto CDcol = find_column(CDtab,"all");
	
	Table HDtab = load_table("HD_region_england_wales.txt");
	auto HDcol = find_column(HDtab,"all");
	
	
	//Table tab = load_table("D_England_wales_age.txt");
	Table tab = load_table("D_England_wales_age_occurance.txt");

	auto files = details.output_directory+"/agree.txt";
	cout << tab.nrow << " " << CDtab.nrow << " " << HDtab.nrow << "y\n";
	
	auto nrow = CDtab.nrow;
	
	/*
	ofstream IGageout(files.c_str());
	for(auto row = 0u; row < nrow; row++){
		IGageout << row << " " << get_int(tab.ele[row][6],"error") << " " <<  get_int(tab2.ele[row][6],"error") << " " << get_int(CDtab.ele[row][CDcol],"error") + get_int(HDtab.ele[row][HDcol],"error") << " " <<  get_int(CDtab.ele[row][CDcol],"error") << " " << get_int(HDtab.ele[row][HDcol],"error") << "\n";
	}
	return;
	*/
	
	for(auto loop = 0u; loop < 2; loop++){
		string st;
		if(loop == 0) st = "CD"; else st = "HD";
		
		auto files = details.output_directory+"/"+st+"_England_wales_age.txt";
		ofstream IGageout(files.c_str());
		
		IGageout << "date\tday\t0-19\t20-64\t65-84\t85+\ttotal\tcum total\n";
	
		auto cum = 0.0;
		for(auto row = 0u; row < nrow; row++){
			auto numCD = get_int(CDtab.ele[row][CDcol],"error");
			auto numHD = get_int(HDtab.ele[row][HDcol],"error");
			
			double fac;
			if(loop == 0) fac = double(numCD)/(numCD+numHD+TINY);
			else fac = double(numHD)/(numCD+numHD+TINY);
			
			IGageout << tab.ele[row][0] << "\t" <<  tab.ele[row][1];
			auto sum = 0.0;
			for(auto a = 0; a < 4; a++){
				auto val = (unsigned int)(0.5+fac*get_int(tab.ele[row][2+a],"error"));
				//auto val = (double)(0.5+fac*get_int(tab.ele[row][2+a],"error"));
				IGageout << "\t" << val;
				sum += val;
			}
			cum += sum;
			
			IGageout << "\t" << sum << "\t" << cum << "\n";
		}
	}
}


void Data::generate_admission_age()
{
	const unsigned int nage = 5;
	const string agegroup[5] = {"0_to_5","6_to_17","18_to_64","65_to_84","85+"};
	
	const unsigned int tmax = 350;
	vector <string> date(tmax);
	
	vector <vector <unsigned int> > val;
	val.resize(tmax);
	for(auto t = 0u; t < tmax; t++){
		val[t].resize(nage); for(auto a = 0u; a < nage; a++) val[t][a] = UNSET; 
	}
	
	Table tab = load_table("hospital admissions_england_age2.txt");

	auto datecol = find_column(tab,"date");
	auto agecol = find_column(tab,"age");
	auto valuecol = find_column(tab,"value");
		
	//for(auto row = 0u; row < tab.nrow; row++){
		//auto value = tab.ele[row][valuecol];
		//auto t = details.gettime(tab.ele[row][datecol],"Error") - details.start;
		//cout << value << " " << t << " h\n";
	//}
		
	for(auto row = 0u; row < tab.nrow; row++){
		auto value = tab.ele[row][valuecol];
		auto age = tab.ele[row][agecol];
		
		if(value != ""){
			auto t = details.gettime(tab.ele[row][datecol],"Error") - details.start;
			date[t] = tab.ele[row][datecol];
			cout << t << " " << age << "age\n";
			auto a = 0u; while(a < nage && agegroup[a] != age) a++;
			if(a == nage) emsg("Prob");
			val[t][a] = get_int(value,"error");
		}
	}
	
	auto files = details.output_directory+"/IH_England_age.txt";
	ofstream IGageout(files.c_str());
	
	IGageout << "date\tday\t0-19\t20-64\t65-84\t85+\ttotal\n";
	
	for(auto t = 0u; t < tmax-1; t++){
		if(val[t][0] != UNSET && val[t+1][0] != UNSET){
			IGageout << date[t] << "\t" << t; //for(auto a = 0u; a < nage; a++) IGageout << "\t" << val[t+1][a]-val[t][a];
			auto sum = 0u; for(auto a = 0u; a < 2; a++) sum +=  val[t+1][a]-val[t][a];
			IGageout << "\t" << sum;
			for(auto a = 2u; a < nage; a++) IGageout << "\t" << val[t+1][a]-val[t][a];
			
			auto sumtot=0u; for(auto a = 0u; a < nage; a++) sumtot += val[t+1][a]-val[t][a];
			IGageout << "\t" <<sumtot;
		
			IGageout << "\n";
		}
	}
	
	/*
	IGageout << "date"; for(auto a = 0u; a < nage; a++) IGageout << "\t" << agegroup[a];
	IGageout << "\n";
	
	for(auto t = 0u; t < tmax-1; t++){
		if(val[t][0] != UNSET && val[t+1][0] != UNSET){
			IGageout << date[t]; for(auto a = 0u; a < nage; a++) IGageout << "\t" << val[t+1][a]-val[t][a];
			IGageout << "\n";
		}
	}
	*/
}

void Data::generatehospdata()
{
	Table dea = load_table("deaths.txt");
	
	auto files = details.output_directory+"/scottish_deaths.txt";
	ofstream scot(files.c_str());
	cout <<  dea.ncol << " " <<  dea.nrow << "nro\n";
	scot << "date";
	for(auto col = 1u; col < dea.ncol; col++){
		string st = dea.heading[col];
		if(st.substr(0,1) == "S"){
			scot << "\t" << st;
		}
	
	}
	scot << "\n";
		
	for(auto row = 0u; row < dea.nrow; row++){
		scot << dea.ele[row][0];
		for(auto col = 1u; col < dea.ncol; col++){
			string st = dea.heading[col];
			if(st.substr(0,1) == "S"){
				scot << "\t" << dea.ele[row][col];
			}
		}
		scot << "\n";
	}
	
	return;
	
	
	
	//Table tab = load_table("Enland_region_admissions.txt","",false);
	Table tab = load_table("Enland_region_discharge.txt","",false);
	cout << tab.nrow <<" nrow\n";
	
	auto t = details.gettime("2020-03-19","Error") - details.start;
	
	auto file = details.output_directory+"/HR_region.txt";
	//auto file = details.output_directory+"/IH_region.txt";
	ofstream data(file.c_str());
	
	//auto filetot = details.output_directory+"/IH_tot.txt";
	auto filetot = details.output_directory+"/HR_tot.txt";
	ofstream datatot(filetot.c_str());
	
	
	unsigned int step = 7;
	
	data << "date";
	for(auto row = 0u; row < tab.nrow; row++){
		data << "\t" << tab.ele[row][0];
	}	
	data << "\n";
	
	unsigned int nweek = (tab.ncol-1)/step; 
	for(auto w = 0u; w < nweek; w++){
		data << details.getdate(t+w*step);
		
		auto sumtot = 0.0;
		for(auto row = 0u; row < tab.nrow; row++){
			unsigned int sum = 0; for(auto col = 1+w*step; col < 1+w*step+step; col++) sum += get_int(tab.ele[row][col],"error");
			data << "\t" << sum;
			sumtot += sum;
		}	
		data << "\n";
		
		datatot << t+w*step << " " << sumtot << "\n";
	}
	
	
	
	//load_table_from_file(string file,"", bool heading) 
}

/// This is used to plot raw data files (this is used for diagnoistic purposed, and not in the analysis)
void Data::generatedeathdata_england_wales()
{
	const unsigned int nloc = 6;
	const string loc[nloc] = {"Hospital","Care home","Hospice","Home","Other communal establishment", "Elsewhere"};
	
	Table area = load_table("areadata_Scotland_National.txt");
	//Table area = load_table("areadata.txt");
	
	//string file = "deaths england and wales_new.txt";
	string file = "deaths england and wales_March2021.txt";
	//Table tab = load_table(file.c_str(),details.output_directory);
	Table tab = load_table(file.c_str());
	
	vector <unsigned int> areaused(area.nrow);
	
	vector < vector <unsigned int> > HDnum;
	vector < vector <unsigned int> > CDnum;
		
	const int nweek = 56;
	HDnum.resize(area.nrow); CDnum.resize(area.nrow);
	for(auto c = 0u; c < area.nrow; c++){
		areaused[c] = 0;
		HDnum[c].resize(nweek); CDnum[c].resize(nweek);
		for(auto w = 0u; w < nweek; w++){
			HDnum[c][w] = 0; CDnum[c][w] = 0;
		}
	}
	
	vector < vector <unsigned int> > locnum;
	locnum.resize(nloc);
	for(auto l = 0u; l < nloc; l++){
		locnum[l].resize(nweek);
		for(auto w = 0u; w < nweek; w++){
			locnum[l][w] = 0;
		}
	}
	
	for(auto row = 0u; row < tab.nrow; row++){
		cout << row << " " <<  tab.nrow << " row\n";
		if(tab.ele[row][7] == "covid-19" && tab.ele[row][12] == "Occurrences"){	
			auto st = tab.ele[row][5];
			auto w = get_int(st.substr(5,st.length()-5),"error")-1;
			
			auto code = tab.ele[row][3];
			
			auto num = get_int(tab.ele[row][0],"error");
			
			auto c = 0u; while(c < area.nrow && code != area.ele[c][0]) c++;
			if(c ==  area.nrow){
				cout << "cannot find\n";
			}
			else{
				if(tab.ele[row][10] == "Hospital"){
					if(HDnum[c][w] != 0) cout << num << "pr\n";
					HDnum[c][w] += num;
				}
				else CDnum[c][w] += num;
			}
			
			auto lo = 0u; while(lo < nloc && loc[lo] != tab.ele[row][10]) lo++;
			if(lo == nloc) emsg("error");
			
			locnum[lo][w] += num;
			
			 areaused[c] = 1;
			//cout << code  << " " << w << " "<< num << " " << c<< "he \n";
			//emsg("P");
		}
	}
	
	
	auto fileb = details.output_directory+"/HD_region_england_wales.txt";
	ofstream data(fileb.c_str());
	data << "date\tday";
	for(auto c = 0u; c < area.nrow; c++){
		if(areaused[c] == 1) data << "\t" << area.ele[c][0];
	}
	data << "\tall\tcum\n";
	
	auto cum = 0u;
	for(auto w = 0u; w < nweek; w++){
		data << details.getdate(w*7) << "\t" << w*7;
		auto sum = 0u;
		for(auto c = 0u; c < area.nrow; c++){
			if(areaused[c] == 1){
				data << "\t" << HDnum[c][w];
				sum += HDnum[c][w];
			}
		}
		cum += sum;
		data << "\t" << sum << "\t" << cum;
		data << "\n";
	}

	auto fileb2 = details.output_directory+"/CD_region_england_wales.txt";
	ofstream data2(fileb2.c_str());
	data2 << "date\tday";
	for(auto c = 0u; c < area.nrow; c++){
		if(areaused[c] == 1) data2 << "\t" << area.ele[c][0];
	}
	data2 << "\tall\tcum";
	data2 << "\n";
	
	cum = 0;
	for(auto w = 0u; w < nweek; w++){
		data2 << details.getdate(w*7) << "\t" << w*7;
		auto sum = 0u;
		for(auto c = 0u; c < area.nrow; c++){
			if(areaused[c] == 1){
				data2 << "\t" << CDnum[c][w];
				sum += CDnum[c][w];
			}
		}
		cum += sum;
		data2 << "\t" << sum << "\t" << cum;
		data2 << "\n";
	}
	
	
	auto filec = details.output_directory+"/D_locationofdeath_england_wales.txt";
	ofstream data3(filec.c_str());
	data3 << "date\tday";
	for(auto l = 0u; l < nloc; l++){
		data3 << "\t" << loc[l];
	}
	data3 << "\tall\tcum\n";
	
	cum = 0;
	for(auto w = 0u; w < nweek; w++){
		data3 << details.getdate(w*7) << "\t" << w*7;
		auto sum = 0u;
		for(auto l = 0u; l < nloc; l++){
			data3 << "\t" << locnum[l][w];
			sum +=locnum[l][w];
		}
		cum += sum;
		data3 << "\t" << sum << "\t" << cum;
		data3 << "\n";
	}
	
	
	/*
	auto filec = details.output_directory+"/HD_tot_england_wales_plot.txt";
	ofstream datac(filec.c_str());
	for(auto w = 0u; w < nweek; w++){
		auto sum = 0.0;
		for(auto c = 0u; c < area.nrow; c++){
			if(areaused[c] == 1) sum += HDnum[c][w];
		}
		datac << w*7 << " " << sum << "\n";
	}
	
	auto filed = details.output_directory+"/CD_tot_england_wales_plot.txt";
	ofstream datad(filed.c_str());
	for(auto w = 0u; w < nweek; w++){
		auto sum = 0.0;
		for(auto c = 0u; c < area.nrow; c++){
			if(areaused[c] == 1) sum += CDnum[c][w];
		}
		datad << w*7 << " " << sum << "\n";
	}
	
		auto filee = details.output_directory+"/HD_tot_england_wales.txt";
	ofstream datae(filee.c_str());
		datae << "date\tall\n";
	for(auto w = 0u; w < nweek; w++){
		auto sum = 0.0;
		for(auto c = 0u; c < area.nrow; c++){
			if(areaused[c] == 1) sum += HDnum[c][w];
		}
		datae << details.getdate(w*7) << " " << sum << "\n";
	}
	
	auto filef = details.output_directory+"/CD_tot_england_wales.txt";
	ofstream dataf(filef.c_str());
	dataf << "date\tall\n";
	for(auto w = 0u; w < nweek; w++){
		auto sum = 0.0;
		for(auto c = 0u; c < area.nrow; c++){
			if(areaused[c] == 1) sum += CDnum[c][w];
		}
		dataf << details.getdate(w*7) << " " << sum << "\n";
	}
	*/
}

void Data::generatedeathdata_scotland()
{
	Table area = load_table("areadata_Scotland_National.txt");
		
	Table tab = load_table("deaths scotland_11-26.txt");
	
	auto nweek = 36u;
	
	vector< vector <unsigned int> > num;
	vector <string> linedate;
	vector <unsigned int> HDnum;
	vector <unsigned int> IDnum;
	
	vector <unsigned int> carehome(nweek), hospital(nweek), other(nweek), home(nweek);
	
		
	num.resize(nweek); linedate.resize(nweek); HDnum.resize(nweek); IDnum.resize(nweek);
	for(auto w = 0u; w < nweek; w++){
		HDnum[w] = 0; IDnum[w] = 0;
		carehome[w] = 0;
		hospital[w] = 0;
		other[w] = 0;
		home[w] = 0;
		
		num[w].resize(area.nrow);
		for(auto r = 0u; r < area.nrow; r++) num[w][r] = 0;
	}

	auto tmin = 10000u, tmax = 0u;
	
	auto tt = 75u;
	
	auto sum = 0.0, sumI = 0.0, sumH = 0.0;
	for(auto row = 0u; row < tab.nrow; row++){
		auto reg = tab.ele[row][0];
		auto date = tab.ele[row][1];
		
		auto n = stoi(tab.ele[row][4]);
		auto sex = tab.ele[row][5];
		auto age = tab.ele[row][6];
		auto cause = tab.ele[row][7];
		auto loc = tab.ele[row][8];
		if(loc != "All") cout << loc << " loc\n";
		if(sex == "All" && age == "All" && cause == "COVID-19 related" && date.substr(0,4) == "w/c "){
			date = date.substr(4,date.length()-4) ;
			//cout << date << " dat\n";
			//cout << details.start+tt << " " <<  details.gettime(date,"Error") << " yy\n";
			if(details.start+tt > details.gettime(date,"Error")) emsg("Problem");
			auto t = details.gettime(date,"Error") - (details.start+tt);
			if(t < tmin) tmin = t;
			if(t > tmax) tmax = t;
			
			auto w = t/7; if(w > nweek) emsg("out");
				
			auto r = 0u; while(r < area.nrow && area.ele[r][0] != reg) r++;
			if(1 == 0 && r < area.nrow){
				if(loc == "All"){				
					linedate[w] = date;
					num[w][r] += n;
					sum += n;
				}
			}
			else{
				if(reg == "S92000003"){
					if(loc != "All"){
						linedate[w] = date;
						cout << loc << " loc\n";
						if(loc == "Hospital") hospital[w] += n;
						else{
							if(loc == "Care Home") carehome[w] += n;
							else{
								if(loc == "Home / Non-institution") home[w] += n;
								else{
									if(loc == "Other institution") other[w] += n;
									else{
										emsg("P");
									}
								}
							}
						}
						
						if(loc == "Hospital"){ HDnum[w] += n; sumH += n;}
						else{ IDnum[w] += n; sumI += n;}
					}
				}
			}
		}
	}
	cout << tmin << " " << tmax <<" max\n";
	//emsg("P");
	cout << "Total in regions: " << sum << endl;
	cout << "Total from hospitals: " << sumH << endl;
	cout << "Total from community: " << sumI << endl;

	ofstream compare(details.output_directory+"/comp.txt");
	compare << "Week\tHospital\tCare Home\tHome\tOther\n";
	for(auto w = 0u; w < nweek; w++){
		compare << w << "\t" <<  hospital[w] << "\t" <<  carehome[w] << "\t"  <<  home[w] << "\t"  <<  other[w] << "\n";
	}		
	
	string file = details.output_directory+"/D_reg_NRS.txt";
	ofstream deathout(file.c_str());
	deathout << "date"; for(auto r = 0u; r < area.nrow; r++) deathout << "\t" << area.ele[r][0];
	deathout << endl;
	for(auto w = 0u; w < nweek; w++){
		deathout << linedate[w];
		for(auto r = 0u; r < area.nrow; r++){
			//sum = 0; for(auto ww = 0u; ww < w; ww++) sum +=;
			deathout << "\t" <<  num[w][r];
		}
		deathout << endl;
	}
	
	file = details.output_directory+"/deaths_NRS.txt";
	ofstream allout(file.c_str());
	allout << "date\tday\tall" << endl;

	for(auto w = 0u; w < nweek; w++){
		allout << linedate[w] << "\t" << 75 + w*7 << " " << HDnum[w]+IDnum[w]  << endl;
	}

	
	file = details.output_directory+"/HD_NRS.txt";
	ofstream HDout(file.c_str());
	HDout << "date\tall" << endl;
	for(auto w = 0u; w < nweek; w++){
		HDout << linedate[w] << "\t" << HDnum[w] << endl;
	}

	file = details.output_directory+"/CD_NRS.txt";
	ofstream IDout(file.c_str());
	IDout << "date\tall" << endl;
	for(auto w = 0u; w < nweek; w++){
		IDout << linedate[w] << "\t" << IDnum[w] << endl;
	}
	
	
	vector < vector <unsigned int> > agedist;
	const int nagecat = 7;
	const string agecat[nagecat] = {"0 years","1-14 years","15-44 years","45-64 years","65-74 years","75-84 years","85 years and over"};	
	agedist.resize(nagecat);
	for(auto a = 0u; a < nagecat; a++){
		agedist[a].resize(nweek);
		for(auto w = 0u; w < nweek; w++) agedist[a][w] = UNSET;
	}
		
	for(auto row = 0u; row < tab.nrow; row++){
		auto reg = tab.ele[row][0];
		auto date = tab.ele[row][1];
		
		auto n = stoi(tab.ele[row][4]);
		auto sex = tab.ele[row][5];
		auto age = tab.ele[row][6];
		auto cause = tab.ele[row][7];
		auto loc = tab.ele[row][8];
		
		if(sex == "All" && age != "All" && loc == "All" && cause == "COVID-19 related" && date.substr(0,4) == "w/c "){
			date = date.substr(4,date.length()-4) ;
			//cout << date << " dat\n";
		
			if(details.start+tt > details.gettime(date,"Error")) emsg("Problem");
			auto t = details.gettime(date,"Error") - (details.start+tt);
			if(t < tmin) tmin = t;
			if(t > tmax) tmax = t;
			
			auto w = t/7; if(w > nweek) emsg("out");
			linedate[w] = date;
				
			auto a = 0u; for(a = 0; a < nagecat; a++) if(agecat[a] == age) break;
			if(a == nagecat){ cout << age << " age\n"; emsg("Problem2");}
			
			if(agedist[a][w] != UNSET) emsg("P");
			agedist[a][w] = n;
		}
	}
	
	file = details.output_directory+"/D_age.txt";
	ofstream Dageout(file.c_str());
	Dageout << "date";
	for(auto a = 0u; a < nagecat; a++) Dageout << "\t" << agecat[a];
	Dageout<< endl;
	for(auto w = 0u; w < nweek; w++){
		Dageout << linedate[w];
		for(auto a = 0u; a < nagecat; a++) Dageout << "\t" << agedist[a][w];
		Dageout << endl;
	}
	
	file = details.output_directory+"/D_age_plot.txt";
	ofstream Dageplot(file.c_str());
	for(auto w = 0u; w < nweek; w++){
		Dageplot << tt + w*7;
		for(auto a = 0u; a < nagecat; a++) Dageplot << "\t" << agedist[a][w];
		
		vector <double> sumst(nagecat);
		auto sum = 0.0;
		for(auto a = 0u; a < nagecat; a++){
			sumst[a] = sum;
			sum += agedist[a][w];
		}
				
		for(auto a = 0u; a < nagecat; a++) Dageplot << "\t" << sumst[a]/sum;
		Dageplot << endl;
	}
	
	file = details.output_directory+"/age_dist.txt";
	ofstream Dagedist(file.c_str());
	for(auto a = 0u; a < nagecat; a++){
		auto sum = 0.0; for(auto w = 0u; w < nweek; w++) sum += agedist[a][w];
		Dagedist << agecat[a] << "\t" << sum << "\n";
	}
	
	file = details.output_directory+"/age_dist_short.txt";
	ofstream Dagedistshort(file.c_str());
	for(auto a = 0u; a < nagecat; a++){
		auto sum = 0.0; for(auto w = 0u; w < 3; w++) sum += agedist[a][w];
		Dagedistshort << agecat[a] << "\t" << sum << "\n";
	}
	
	for(auto a = 0u; a < nagecat; a++){
		auto sum = 0.0; for(auto w = 0u; w < 3; w++) sum += agedist[a][w];
		Dagedistshort << agecat[a] << "\t" << sum << "\n";
	}
	
	for(auto w = 0u; w < nweek; w++){
		cout << w << " " << linedate[w] << " dat\n";
	}

	for(auto a = 0u; a < nagecat; a++){
		auto sum = 0.0; for(auto w = 0u; w < 20; w++) sum += agedist[a][w];
		cout << agecat[a] << "\t" << sum << "\n";
	}
}


/// Converts the matrix sent by stephen
void Data::convert_Mdata()
{
	ifstream regin(data_directory+"/contact matrix regional names.txt");
	ofstream regdat(data_directory+"/regiondata.txt");
	ofstream areadat(data_directory+"/areadata.txt");
		
	regdat << "name\tcode" << endl;
	areadat << "area\teasting\tnorthing\tregion\tdensity\tage0-14\tage15-44\tage45-64\t65+\tMale\tFemale" << endl;
	do{
		string line, code, name;
		getline(regin,line);
		if(line.length() > 2){
			code = line.substr(1,9);
			int j = line.length()-10; while(line.substr(j,1) != "\"") j++;
			name = line.substr(13,j);
			strip(name);
	
			regdat << name << "\t" << code << endl;
			areadat << code << "\t" << ran() << "\t" << ran() << "\t" << code << "\t1\t50000\t100000\t100000\t100000\t50\t50" << endl;
		}
	}while(!regin.eof());
	
	const unsigned int si=171; 
	double M[si][si];
	
	ifstream ip(data_directory+"/contactmatrix.txt");
	unsigned int i, j;
	for(i = 0; i < si; i++){
		for(j = 0; j < si; j++){
			ip >> M[i][j];
		}
	}
		
	ofstream op(data_directory+"/Mdata.txt");
	op << "'region1'	'region2'	'contact'" << endl;
	for(i = 0; i < si; i++){
		for(j = i; j < si; j++){
			op << i << "\t" << j << "\t" << M[i][j] << endl;
		}
	}
	
	/*
	Table tabarea = load_table("areadata.txt");
	
	Table tabcheck = load_table("datazones.txt");
	
	vector <unsigned int> conv(tabarea.nrow);
	
	for(auto r = 0u; r < tabarea.nrow; r++){
		auto st = tabcheck.ele[r][0];
		auto code = st.substr(st.size()-9,9);
		
		unsigned int a;
		for(a = 0; a < tabarea.nrow; a++) if(tabarea.ele[a][1] == code) break;
		if(a == tabarea.nrow) cout << "cannot find\n";
		conv[r] = a;
		
		cout << tabarea.ele[a][1] << " " << tabcheck.ele[a][0] << "  g:" << code << " samz\n";
	}

	cout << tabarea.nrow << " " << tabcheck.nrow << " jj\n";

	unsigned int i, j;
	double val;
	
	ifstream inp(data_directory+"/M.txt");
	ofstream op(data_directory+"/Mdata.txt");
	op << "DZ1 DS2 contact" << endl;
	op << setprecision(5);
	do{
		inp >> i >> j >> val;
		
		cout << i << " " << j << " " << val << " k\n";
		auto a1 = conv[i-1]; 
		auto a2 = conv[j-1];
		if(a1 > a2){ auto a = a1; a1 =a2; a2= a;}
		op << a1 << "\t" << a2 << "\t" << val << endl;
	}while(!inp.eof());
	*/
}


/// Plots raw data
void Data::plotrawdata()
{
	Table tab = load_table("apple_mobility.txt");
	
	auto row_walk = 0u;
	while(row_walk < tab.nrow && !(tab.ele[row_walk][1] == "Scotland" && tab.ele[row_walk][2] == "walking")) row_walk++;
	
	auto row_drive = 0u;
	while(row_drive < tab.nrow && !(tab.ele[row_drive][1] == "Scotland" && tab.ele[row_drive][2] == "driving")) row_drive++;
	
	auto row_transit = 0u;
	while(row_transit < tab.nrow && !(tab.ele[row_transit][1] == "Scotland" && tab.ele[row_transit][2] == "transit")) row_transit++;
	
	cout << row_walk << " " << row_drive << " " << row_transit << "row\n";
	
	vector <double> val_walk, val_driving, val_transit;
	for(auto col = 9u; col < tab.ncol-4; col++){
		auto sum = 0.0; for(unsigned int c = col-3; c <= col+3; c++) sum += get_double(tab.ele[row_walk][c],"error");
		val_walk.push_back(sum);
		
		sum = 0.0; for(unsigned int  c = col-3; c <= col+3; c++) sum += get_double(tab.ele[row_drive][c],"error");
		val_driving.push_back(sum);
		
		sum = 0.0; for(unsigned int c = col-3; c <= col+3; c++) sum += get_double(tab.ele[row_transit][c],"error");
		val_transit.push_back(sum);
	}
	
	ofstream data(details.output_directory+"/apple_mobility_plot.txt");
	for(auto i = 0u; i < val_walk.size(); i++){
		data << i + 15 << " " <<  100*val_walk[i]/val_walk[0] << " " <<  100*val_driving[i]/val_driving[0] << " " <<  100*val_transit[i]/val_transit[0] << "\n";
	}
	
	/*
	string fileraw = details.output_directory+"/rawdeaths.txt";
	
	ofstream raw(fileraw.c_str());
	for(auto &tdata : transdata){                                // Loads transition data for inference
		cout << tdata.fromstr << " " << tdata.tostr << " fro to\n";
		auto sum = 0.0;
		for(auto col = 0u; col < tdata.num[0].size(); col++){
			//raw << tdata.start+col*tdata.units << " " << sum << "\n";
			auto sum2 = 0.0;
			for(auto row = 0u; row < tdata.num.size(); row++){
				sum += tdata.num[row][col];
				sum2 += tdata.num[row][col];
			}
			raw << tdata.start+(col+0.5)*tdata.units << " " << sum2 << "\n";
		}
	}
	
	Table tab = load_table("areadata.txt");
	auto col = find_column(tab,"density");
	
	string file = details.output_directory+"/death_vs_density.txt";
	ofstream outp(file.c_str());
	
	for(auto &tdata : transdata){                                // Loads transition data for inference
		string file2 = tdata.file; 
		Table tab2 = load_table(file2);
		
		for(auto r = 0u; r < nregion; r++){
			auto col2 = find_column(tab2,region[r].code);
			auto sum = 0.0;
			for(auto row = 0u; row < tab2.nrow; row++){
				sum += get_double(tab2.ele[row][col2].c_str());
			}
			
			outp << get_double(tab.ele[r][col].c_str()) << " " << sum << "\n";
		}
	}
	*/
}


/// Generates a initpop file
void Data::generate_initpop()
{
	ofstream fout(details.output_directory+"/initpop.csv");
	
	fout << "area"; 
	for(auto d = 0u; d < ndemocat; d++){
		if(democat[d].value.size() > 1){
			fout << "," << democat[d].name;
		}
	}
	fout << ",compartment,population" << endl;
	
	for(auto a=0u; a < narea; a++){
		for(auto dp = 0u; dp < ndemocatpos; dp++){
			fout << area[a].code;
			for(auto d = 0u; d < ndemocat; d++){
				if(democat[d].value.size() > 1){
					fout << "," << democat[d].value[democatpos[dp][d]];
				}
			}
			fout << ",S," << area[a].pop_init[0][dp] << endl;
		}
	}
}


// Simulates files for time-varying covariates
void Data::generate_tvcovar()
{
	ofstream fout(details.output_directory+"/tvcovar.csv");
	fout << "Date,Temperature" << endl;
	auto mu = 20.0, sd = 2.0, lam = 0.1;
	auto T = mu;
	for(auto t = 0u; t < details.period; t++){
		T += -lam*(T-mu) + normal_sample(0,sd);
		fout << details.getdate(t) << "," << T << endl; 
	}
	
	ofstream fout2(details.output_directory+"/areatvcovar.csv");
	fout2 << "Date";
	for(auto c = 0u; c < narea; c++){
		fout2 << "," << area[c].code;
	}
	fout2 << endl;
	
	vector <double> Tvec(narea);
	for(auto c = 0u; c < narea; c++) Tvec[c] = mu + normal_sample(0,3*sd);
			
	for(auto t = 0u; t < details.period; t++){
		fout2 << details.getdate(t);
		for(auto c = 0u; c < narea; c++){
			Tvec[c] += -lam*(Tvec[c]-mu) + normal_sample(0,sd);	

			fout2 << "," << Tvec[c]; 
		}
		fout2 << endl;
	}
}


/// Simulates a file that can be used to test algorithm
void Data::generate_level_effect()
{
	auto file = details.output_directory+"/level.csv";

	string levels = "level1,level2";
	
	auto levs = split(levels,',');  
	
	ofstream lout(file.c_str());
	lout << "Date";
	for(auto c = 0u; c < narea; c++) lout << "," << area[c].code; 
	lout << endl;
	
	auto T =  details.end- details.start;
	vector< vector <string> > map;
	map.resize(T);
	for(auto t = 0u; t < T; t++) map[t].resize(narea);

	for(auto c = 0u; c < narea; c++){
		auto t = 0u;
		do{
			auto val = (unsigned int)(levs.size()*ran());
			auto tnext = t + (unsigned int)(30*ran());
			if(tnext > T) tnext = T;
			while(t < tnext){ map[t][c] = levs[val]; t++;}
		}while(t < T);
	}
	
	for(auto t = 0u; t < T; t++){
		lout << details.getdate(t); 
		for(auto c = 0u; c < narea; c++) lout << "," << map[t][c];
		lout << endl;
	}
}


/// Loads up a GEOJSON file and reduces its size (so it can be put on the GitHub repository)
void Data::reducesize_geojson(const string file)
{
	auto filefull = data_directory+"/"+file;

	ifstream boundfile(filefull);
	if(!boundfile) emsg("Cannot open the file '"+file+"'.");
	string s;
	do{
		string st;
		getline( boundfile,st); s += st; 
		if( boundfile.eof()) break;
	}while(true);
	
	json jso = json::parse(s);
	
	auto &h = jso["features"];
	
	vector <unsigned int> rem(h.size());
	
	auto loop= 0u;
	for(auto &it : h.items()){                          // Deletes unused areas
		json val = it.value();
		//cout << val << " val\n";
		auto h2 = val["type"];
		if(h2 == "Feature"){
			auto c = UNSET; 
			auto h3 = val["properties"];
			for(auto &it3 : h3.items()){     // Checks if one of the properties contains the name of the area 
				if(it3.value().is_string()){
					string name = it3.value();
					strip(name);		
					c = 0; while(c < narea && name != area[c].code) c++;
					if(c < narea) break;
				}
			}
		
			if(c == UNSET || c == narea) rem[loop] = 1;
		}
		loop++;
	}

	for(int loop = h.size()-1; loop >= 0; loop--){
		if(rem[loop] == 1) h.erase(h.begin()+loop);
	}

	const double acc = 1000;
	for(auto &it : h.items()){                                  // Removes significant figures
		json &val = it.value();
		auto h2 = val["type"];
		if(h2 == "Feature"){
			auto &h3 = val["geometry"];	
			auto h4 = h3["type"];		
			if(h4 == "Polygon"){
				auto &h5 = h3["coordinates"];
				for(auto &it2 : h5.items()){
					json &val2 = it2.value();
					for(auto &it3 : val2.items()){
						json &val3 = it3.value();
						double sx = val3[0]; val3[0] = int(sx*acc)/acc;
						double sy = val3[1]; val3[1] = int(sy*acc)/acc;
					}
				}
			}
				
			if(h4 == "MultiPolygon"){
				auto &h5 = h3["coordinates"];
				for(auto &it : h5.items()){
					json &val = it.value();
					for(auto &it2 : val.items()){
						json &val2 = it2.value();
						for(auto &it3 : val2.items()){
							json &val3 = it3.value();
							double sx = val3[0]; val3[0] = int(sx*acc)/acc;
					  	double sy = val3[1]; val3[1] = int(sy*acc)/acc;
						}
					}
				}
			}
		}
	}
	
	double dif = 0.004;
	auto nump = 0u;
	for(auto &it : h.items()){                                  // Removes fine details
		json &val = it.value();
		auto &h2 = val["type"];
		if(h2 == "Feature"){
			auto &h3 = val["geometry"];	
			auto &h4 = h3["type"];		
			double xst=0, yst=0;
			
			if(h4 == "Polygon"){
				auto &h5 = h3["coordinates"];
				for(auto &it2 : h5.items()){
					json &val2 = it2.value();
					auto loop = 0u; 
					vector <unsigned int> rem(val2.size());  for(auto i = 0u; i <val2.size(); i++) rem[i] = 0;
					for(auto &it3 : val2.items()){
						json &val3 = it3.value();
						if(loop == 0){
							xst = val3[0]; yst = val3[1]; 
						}
						else{
							double x = val3[0];
							double y = val3[1];
							double dx = x-xst;
							double dy = y-yst;
							if(dx*dx+dy*dy < dif*dif) rem[loop] = 1;
							else{ xst += dx; yst += dy;}
						}
						loop++;
					}
					
					for(int loop = val2.size()-1; loop >= 0; loop--){
						if(rem[loop] == 1) val2.erase(val2.begin()+loop);
					}
					nump += val2.size();
				}
			}
				
			if(h4 == "MultiPolygon"){
				auto &h5 = h3["coordinates"];
				for(auto &it : h5.items()){
					json &val = it.value();
					for(auto &it2 : val.items()){
						json &val2 = it2.value();
						auto loop = 0u; 
						vector <unsigned int> rem(val2.size());  for(auto i = 0u; i <val2.size(); i++) rem[i] = 0;
						for(auto &it3 : val2.items()){
							json &val3 = it3.value();
							if(loop == 0){
								xst = val3[0]; yst = val3[1]; 
							}
							else{
								double x = val3[0];
								double y = val3[1];
								double dx = x-xst;
								double dy = y-yst;
								if(dx*dx+dy*dy < dif*dif) rem[loop] = 1;
								else{ xst += dx; yst += dy;}
							}
							loop++;
						}
						
						for(int loop = val2.size()-1; loop >= 0; loop--){
							if(rem[loop] == 1) val2.erase(val2.begin()+loop);
						}
						nump += val2.size();
					}					
				}
			}
		}
	}
	
	cout << nump <<" nump\n";
	
	auto fi =  details.output_directory+"/Scotland_areas.geojson";
	ofstream fiout(fi);
	fiout << jso.dump() << endl;
	
	//cout << jso.dump() << "dum\n";    
}


/// Generates simulated time variation in age mixing matrix 
void Data::generate_age_mixing_perturb()
{
	auto file =  details.output_directory+"/amp.csv";

	//vector <string> date = {"start","2020-02-26","2020-03-11","2020-03-23","2020-03-30","2020-04-08","2020-04-22","2020-05-06","2020-05-20","2020-06-03","2020-06-17","2020-07-01","2020-07-15","2020-07-29","2020-08-12","2020-08-26","2020-09-09","2020-09-23","2020-10-07","2020-10-21","2020-11-04","2020-11-18","2020-12-02","2020-12-16","2020-12-30","2021-01-06","2021-01-13","2021-01-27","2021-02-10","2021-02-24","2021-03-10","2021-03-24","2021-04-07","2021-04-21","2021-05-05","2021-05-19","end"};

	//vector <string> date = {"start","2020-02-26","2020-03-11","2020-03-23","2020-04-08","2020-05-06","2020-06-03","2020-07-01","2020-07-29","2020-08-26","2020-09-23","2020-10-21","2020-11-18","2020-12-16","2021-01-06","2021-01-27","2021-02-24","2021-03-24","2021-04-21","2021-05-19","end"};

	vector <string> date = {"start","14","28","42","56","70","84","98","112","126","140","154","168","182","196","end"};

	//vector <string> cat = {"0-4","5-9","10-14","15-19","20-24","25-29","30-34","35-39","40-44","45-49","50-54","55-59","60-64","65-69","70-74","75-79","80+","CH"};
	
	//vector <string> cat = {"0-9","10-19","20-29","30-39","40-49","50-59","60-69","70-79","80+","CH"};
	
	//vector <string> cat = {"0-19","20-39","40-59","60-79"};
	vector <string> cat = {"0-19","20-39","40-59","60+"};
	
	ofstream op(file);
	op << "bp";
	for(auto j = 0u; j < cat.size(); j++) op << ",Par" << cat[j];
	for(auto j = 0u; j < cat.size(); j++) op << ",Value" << cat[j];
	op << endl;

	vector <double> val(cat.size());
	for(auto j = 0u; j < cat.size(); j++) val[j] = normal_sample(0,0.2);

	vector <double> frac(cat.size());
	auto sum = 0.0;
	for(auto j = 0u; j < cat.size(); j++){
		frac[j] = area[0].pop_dp[j];
		sum += frac[j];
	}			
	for(auto j = 0u; j < cat.size(); j++) frac[j] /= sum;
			
	auto step = 1;
	
	for(auto i = 0u; i < date.size(); i += step){
		if(i < date.size()-step) op << date[i]; else op << "end";
		
		for(auto j = 0u; j <  cat.size(); j++) op << ",Mixing factor " << cat[j] << " " << date[i];
		
		auto sum = 0.0;
		vector <double> expval(cat.size());
		for(auto j = 0u; j < cat.size(); j++){
			expval[j] = exp(val[j]);
			sum += expval[j]*frac[j];
		}
		
		for(auto j = 0u; j < cat.size(); j++){
			op << ",";
			if(j < cat.size()-1) op << (expval[j]/sum);
			else op << "*";
		}
		op << endl;
	
		for(auto j = 0u; j < cat.size(); j++) val[j] += normal_sample(0,0.2);
	}
		
	cout << "age_mixing_perturb = [" << endl;
	for(auto j = 0u; j < cat.size(); j++){
		cout << "{ agecat='age" << cat[j] << "', bp='[bp:amp.csv]', param='[Par" << cat[j] << ":amp.csv]', value='[Value" << cat[j] << ":amp.csv]', prior='MDir(0.3)'" <<  ", smooth='0.1', smooth_type='log_smooth'}," << endl;
	}
	cout << "]" << endl;	
}

void Data::generate_geo_mixing_matrix()
{
	auto narea = 10u;
	
	vector <unsigned int> pop = {2000,1000,500,3000,2000,1000,4000,2000,1000,3000};
	
	vector < vector <double> > mat;
	mat.resize(narea);
	for(auto j = 0u; j < narea; j++){
		mat[j].resize(narea);
	}

	for(auto j = 0u; j < narea; j++){
		for(auto i = j+1; i < narea; i++){
			mat[j][i] = 0.02+ran()*0.02;
			mat[i][j] = mat[j][i]*pop[i]/pop[j];
		}
	}		
	
	for(auto j = 0u; j < narea; j++){
		auto sum = 0.0;
		for(auto i = 0u; i < narea; i++){
			if(i != j) sum += mat[j][i];
		}
		mat[j][j] = 1-sum;
	}
	
	ofstream fout(details.output_directory+"/mixing2.csv");

	auto tot = 0.0;
	for(auto j = 0u; j < narea; j++){
		for(auto i = 0u; i < narea; i++){
			tot += mat[j][i]/pop[j];
		}
	}
	//auto fac = narea/tot;
	
	for(auto j = 0u; j < narea; j++){
		auto sum = 0.0; for(auto i = 0u; i < narea; i++) sum += mat[j][i]/pop[j];
		auto fac = 1.0/sum;
		
		for(auto i = 0u; i < narea; i++){
			if(i > 0) fout << ",";
			fout << fac*mat[j][i]/pop[j];
		}
		fout << endl;
	}
}

