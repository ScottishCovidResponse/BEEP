// Stores details of the simulation / inference

#include <iostream>
#include <math.h> 
#include <cstring>
#include <sstream>

using namespace std;

#include "details.hh"

Details::Details(const Inputs &inputs)
{
	mode = inputs.mode();
	
	siminf = inputs.get_siminf();

	output_directory = inputs.find_string("outputdir","Ouput");                       // Output directory

	string timeformat = inputs.find_string("timeformat","number");                    // Time format
	if(timeformat == "number"){ time_format = TIME_FORMAT_NUM; time_format_str = "time";}
	else{ 
		if(timeformat == "year-month-day"){ time_format = TIME_FORMAT_YMD; time_format_str = "date";}
		else emsgroot("Do not recognise time format '"+time_format_str+"'.");
	}
	
	efoi_factor = inputs.find_integer("efoi_factor",100000);                          // External force of infection factor
	 
	string startstr = inputs.find_string("start","UNSET");                            // Beginning and ending times
	if(startstr == "UNSET") emsgroot("The 'start' time must be set"); 
	start = gettime(startstr);

	string endstr = inputs.find_string("end","UNSET");
	if(endstr == "UNSET") emsgroot("The 'end' time must be set"); 
	end = gettime(endstr);
	
	period = end - start;
	
	inputs.find_timeplot(timeplot);
	for(auto &tp : timeplot) tp.time = gettime(tp.time_str);

	switch(modeltype){
		case IND_MODEL: division_per_time = 1; break;
		case POP_MODEL: division_per_time = 2; break;
	}
	
	ndivision = division_per_time*period;
	division_time.resize(ndivision+1);
	for(auto s = 0u; s <= ndivision; s++) division_time[s] = double(s*period)/ndivision;
	
	
	if(mode == PMCMC_INF){                                             // This defines when particle filtering is performed 
		obs_combine_type = COMBINE_SECTION;
		pmcmc_obs_period = 14*division_per_time;
		sec_define.resize(ndivision);
		for(auto s = 0u; s < ndivision; s++){
			if(s > 0 && s%pmcmc_obs_period == 0) sec_define[s] = true;
			else sec_define[s] = false;
		}
		//for(auto s = 0u; s < ndivision; s++) cout << s << " " <<  sec_define[s] << "  sec define\n"; emsg("O");
	}
	else obs_combine_type = COMBINE_VALUE;
}

/// Gets the time from a string
unsigned int Details::gettime(string st) const 
{
	if(st == "start") return start;
	if(st == "end") return end;
	
	for(auto &tp : timeplot){
		if(st == tp.name) return tp.time;
	}
	
	auto t = 0u;
	const char *buf = st.c_str();
	
	switch(time_format){
	case TIME_FORMAT_NUM:
		t = atoi(buf);
		if(isnan(t)) emsg("The time '"+st+"' is not a number");
		break;

	case TIME_FORMAT_YMD:
		struct tm result;		
		memset(&result, 0, sizeof(result));
		if(strptime(buf,"%Y-%m-%d",&result) != NULL){
			time_t tt = mktime(&result);
			t = tt/(60*60*24);
		}
		else emsg("'"+st+"' is not regonised as Year-Month-Day format.");
		break;
		
	default:
		emsg("The time format is not recognised.");
		break;
	}

	return t;	
}

/// Returns a date from a time
string Details::getdate(unsigned int t) const
{
	t += start;

	stringstream ss; 
	switch(time_format){
	case TIME_FORMAT_NUM:
		ss << t;
		break;
		
	case TIME_FORMAT_YMD:
		time_t tt = (t + 0.5)*(60*60*24);	
		
		struct tm *timeinfo;
		timeinfo = localtime(&tt);
		
		char buffer[80];
		strftime(buffer,80,"%Y-%m-%d",timeinfo);
		ss << buffer;
		break;
	}
	
	return ss.str();
}
