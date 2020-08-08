
#include <vector>
#include <string>

using namespace std;

#include "combinetrace.hh"

/// Finishing this is on the TO DO list
void combine_trace(DATA & /* data */, Inputs & /*inputs*/)
{
	vector <string> paramname;
	vector < vector < vector <double> > > vals;
	
		/*
		if (cmdlineparams.count("dirs") == 0) emsg("When using the 'combinetrace' mode, you must set the 'dirs' property");
		vector <string> dirs;
		string output, distfile="";
		unsigned int burnin=UNSET;
		dirs = split(cmdlineparams["dirs"],',');
		
		if(cmdlineparams.count("output") == 0) emsg("When using the 'combinetrace' mode, you must set the 'output' property");
	
		if(cmdlineparams.count("distribution") == 1) distfile = cmdlineparams["distribution"];
		if(cmdlineparams.count("burnin") == 1){
			burnin = atoi(cmdlineparams["burnin"].c_str());
			if(std::isnan(burnin)) emsg("The 'burnin' property must be an integer."); 
		}
	
		output = cmdlineparams["output"];
		data.combinetrace(dirs,output,distfile,burnin);
		*/
	
}

/*
// This function combines results from different trace files to generate overall statistics
void DATA::combinetrace(vector <string> inputdirs, string outputfile, string distfile, unsigned int burnin)
{
	unsigned int inp, i, row, th;
	double v;

	TABLE tab;
	
	vals.resize(inputdirs.size());
	for(inp = 0; inp < inputdirs.size(); inp++){
		tab = loadtable("trace.txt",inputdirs[inp]);
	
		for(i = 1; i < tab.ncol; i++){
			if(tab.heading[i] == "zero") break;
			if(inp == 0) paramname.push_back(tab.heading[i]);
			else{
				if(i-1 >= paramname.size()) emsg("The columns in the input trace files do not match up.");
				if(paramname[i-1] != tab.heading[i]) emsg("The headings in the input trace files do not match up.");
			}
		}
		vals[inp].resize(paramname.size());
		
		for(th = 0; th < paramname.size(); th++){
			for(row = 0; row < tab.nrow; row++){
				v = atof(tab.ele[row][th+1].c_str()); if(std::isnan(v)) emsg("In file '"+inputdirs[inp]+"/trace.txt' the quantity '"+tab.ele[row][th]+"' is not a number.");
				vals[inp][th].push_back(v);
			}
		}
	}
	
	Output output(details,data,model);
	output.combinedtrace(paramname,vals,outputfile,distfile,burnin);
}
*/
