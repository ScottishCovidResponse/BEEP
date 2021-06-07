// Stores details of the simulation / inference

#include <iostream>
#include <math.h> 
#include <cstring>
#include <sstream>

using namespace std;

#include "details.hh"

Details::Details(Inputs &inputs)
{
	toml_file = inputs.inputfilename;
	
	mode = inputs.mode();
	
	stochastic = true;
	auto dynamics = inputs.find_string("dynamics","stochastic");  
	if(dynamics != "stochastic"){
		if(dynamics != "deterministic") emsg("'dynamics' must be set to 'stochastic' or 'deterministic'");
		stochastic = false;
	}
	
	if(stochastic == false && (mode == ABC_MBP || mode == MC3_INF || mode == MCMC_MBP || mode == PAIS_INF)){
		emsgroot("When running 'dynamics=\"determinstic\"', 'mode' can only use inference methods which do not rely on MPBs: 'abc', 'abcsmc' or 'pmcmc'"); 
	}
	
	siminf = inputs.get_siminf();

	output_directory = inputs.find_string("outputdir","Ouput");       // Output directory

	string timeformat = inputs.find_string("time_format","number");    // Time format
	if(timeformat == "number"){ time_format = TIME_FORMAT_NUM; time_format_str = "time";}
	else{ 
		if(timeformat == "year-month-day"){ time_format = TIME_FORMAT_YMD; time_format_str = "date";}
		else{
			if(timeformat == "day/month/year"){ time_format = TIME_FORMAT_DMY_SLASH; time_format_str = "date";}
			else{
				if(timeformat == "day.month.year"){ time_format = TIME_FORMAT_DMY_DOT; time_format_str = "date";}
				else{
					emsgroot("The value of 'time_format="+time_format_str+"' is not recognised (it must be 'number', 'year-month-day', 'day/month/year' or 'day.month.year').");
				}
			}
		}
	}
	
	efoi_factor = inputs.find_positive_integer("efoi_factor",100000);   // External force of infection factor
	 
	string startstr = inputs.find_string("start","UNSET");              // Beginning and ending times
	if(startstr == "UNSET") emsgroot("A value for 'start' must be set"); 
	start = gettime(startstr,"In 'start'");

	string endstr = inputs.find_string("end","UNSET");
	if(endstr == "UNSET") emsgroot("A value for 'end' must be set"); 
	end = gettime(endstr,"In 'end'");
	
	period = end - start;
	
	pred_start = UNSET;
	pred_end = UNSET;

	if(mode == PREDICTION){ 
		string predstr = inputs.find_string("prediction_start","UNSET");
		if(predstr != "UNSET"){
			pred_start = gettime(predstr,"In 'prediction_start'") - start;
			if(pred_start > period) emsg("'prediction_start' cannot be after 'end'");
		}
	
		predstr = inputs.find_string("prediction_end","UNSET");
		if(predstr != "UNSET"){
			pred_end = gettime(predstr,"In 'prediction_end'") - start;
		}
	}
	if(pred_end != UNSET && pred_end > period) period = pred_end;
	
	inputs.find_timeplot(timeplot);
	for(auto &tp : timeplot) tp.time = gettime(tp.time_str,"In 'time_labels'");

	division_per_time = inputs.find_positive_integer("steps_per_unit_time",2);
	
	ndivision = division_per_time*period;
	division_time.resize(ndivision);
	for(auto s = 0u; s < ndivision; s++){
		division_time[s] = double((s+0.5)*period)/ndivision;
		if(time_format == TIME_FORMAT_NUM) division_time[s] += start;
	}
	
	trans_combine = inputs.find_double("trans_combine",UNSET);         // Determines if trans data is combined
	
	if(mode == PMCMC_INF){                                             // This defines when particle filtering is performed 
		obs_section = true;
		pmcmc_obs_period = 14*division_per_time;
		sec_define.resize(ndivision);
		for(auto s = 0u; s < ndivision; s++){
			if(s > 0 && s%pmcmc_obs_period == 0) sec_define[s] = true;
			else sec_define[s] = false;
		}
	}
	else obs_section = false;
	
	inputs.find_mcmc_update(mcmc_update);
}


/// Gets the time from a string
unsigned int Details::gettime(const string st, const string em) const 
{
	if(st == "start") return start;
	if(st == "end") return end;
	
	for(auto &tp : timeplot){
		if(st == tp.name) return tp.time;
	}
	
	auto t = 0u;
	const char *buf = st.c_str();
	struct tm result;		
		
	switch(time_format){
	case TIME_FORMAT_NUM:
		t = get_int(st,em);
		break;

	case TIME_FORMAT_YMD:
		memset(&result, 0, sizeof(result));
		if(strptime(buf,"%Y-%m-%d",&result) != NULL){
			time_t tt = mktime(&result);
			t = tt/(60*60*24);
		}
		else emsg(em+" the expression '"+st+"' is not regonised as 'year-month-day' format.");
		break;
		
	case TIME_FORMAT_DMY_SLASH:
		memset(&result, 0, sizeof(result));
		if(strptime(buf,"%d/%m/%y",&result) != NULL){
			time_t tt = mktime(&result);
			t = tt/(60*60*24);
		}
		else emsg(em+" the expression '"+st+"' is not regonised as 'day/month/year' format.");
		break;
		
	case TIME_FORMAT_DMY_DOT:
		memset(&result, 0, sizeof(result));
		if(strptime(buf,"%d.%m.%y",&result) != NULL){
			time_t tt = mktime(&result);
			t = tt/(60*60*24);
		}
		else emsg(em+" the expression '"+st+"' is not regonised as 'day.month.year' format.");
		break;
	
	default:
		emsg("The time format is not recognised.");
		break;
	}

	return t;	
}


/// Returns a date from a time
string Details::getdate(const unsigned int t) const
{
	auto tshift = t + start;

	time_t tt = (tshift + 0.5)*(60*60*24);	
	struct tm *timeinfo;
	timeinfo = localtime(&tt);
	char buffer[80];
		
	stringstream ss; 
	switch(time_format){
	case TIME_FORMAT_NUM:
		ss << tshift;
		break;
		
	case TIME_FORMAT_YMD:
		strftime(buffer,80,"%Y-%m-%d",timeinfo);
		ss << buffer;
		break;
		
	case TIME_FORMAT_DMY_SLASH:
		strftime(buffer,80,"%d/%m/%y",timeinfo);
		ss << buffer;
		break;
		
	case TIME_FORMAT_DMY_DOT:
		strftime(buffer,80,"%d.%m.%y",timeinfo);
		ss << buffer;
		break;
		
	default:
		emsg("The time format is not recognised.");
		break;
	}
	
	return ss.str();
}
