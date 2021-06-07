// Outputs various graphs and statistics

#include <iostream>
#include <vector>
#include <sys/stat.h>
#include <sstream>
#include <algorithm>
#include <cmath>
#include <iomanip> 

using namespace std;

#include "output.hh"
#include "mpi.hh"
#include "details.hh"
#include "model.hh"
#include "state.hh"

Output::Output(const Details &details, const Data &data, const Model &model, Inputs &inputs, const ObservationModel &obsmodel, Mpi &mpi) :  inputs(inputs), details(details), data(data), model(model), obsmodel(obsmodel), mpi(mpi)
{
	if(mpi.core == 0){
		ensure_directories();                                             // Creates directories for outputting
		
		print_model_specification();                                      // Summarises the data and model
		
		inputs.find_outputprop(prop);                                     // Load output properties
		inputs.find_nrun(nrun);
		
		readme(); 	                                                    	// Generates a readme file for outputs
	}
	mpi.barrier();
}
	

/// Generates the output files along with the final pdf report that gives graphical outputs
void Output::generate_graphs(vector <ParamSample> &psamp, const vector <Sample> &opsamp) const
{ 
	for(auto &psa : psamp) psa.paramval = model.dirichlet_correct(psa.paramval);
	
	cout << endl;
	if(suppress_output == false){
		if(details.siminf == SIMULATE) cout << "Outputs in directory '" << details.output_directory << "':" << endl;
		else cout << post_lab << " outputs in directory '" << details.output_directory << "':" << endl;
	}
	
	vector <OutputPlot> op;                                           // Stores the plots for the final pdf

	spline_plots(opsamp,op);                                          // Generates all the different types of plot

	graph_plots(opsamp,op);
	 
	susceptibility_distributions(psamp,op);

	mean_distributions(psamp,op);

	branch_prob_distributions(psamp,op);

	level_effect_distributions(psamp,op);
	
	area_effect_distributions(psamp,op);
	
	derived_parameter_distributions(opsamp,op);

	add_democat_change(op);
	
	set_labelon(op);                                                  // Determines if labels put in description
	
	auto grfile = "Graphs";
	generate_gnuplot_file(op,grfile);                                 // Generates the gnuplot file
	
	generate_pdf(grfile,"Output graphs");                             // Runs gnuplot and generates pdf 
	
	generate_pdf_description(op,grfile);                              // Generates a text descrition of the pdf file

	if(details.siminf == INFERENCE){
		vector <OutputPlot> op_diag;                                    // Stores the plots for the diagnostic pdf

		add_generation_plots(op_diag);
		
		posterior_parameter_distributions(psamp,op_diag);
		
		add_trace_plots(op_diag);
		
		set_labelon(op_diag);                                           // Determines if labels put in description
		
		if(op_diag.size() > 0){
			auto grfile = "Diagnostic_Graphs";
			generate_gnuplot_file(op_diag,grfile);                        // Generates the gnuplot file
		
			generate_pdf(grfile,"Output diagnostic graphs");              // Runs gnuplot and generates pdf 
		
			generate_pdf_description(op_diag,grfile);                     // Generates a text descrition of the pdf file
		}
	
		spatial_R_map(psamp);                                           // Outputs a file giving spatial variation in R
	
		posterior_parameter_estimates(psamp);                           // Outputs a file giving posterior estimates
	}
}


/// Generate the final pdf report
void Output::generate_pdf(const string file, const string desc) const
{
	auto opdir = details.output_directory;
	stringstream ss; ss << "gnuplot '" << opdir << "/gnuplot.txt'" << endl;
	system(ss.str().c_str());
		
	stringstream ss2; ss2 << "ps2pdf " << opdir << "/" << file << ".ps " << opdir << "/" << file << ".pdf" << endl;
	system(ss2.str().c_str());

	stringstream ss3; ss3 << "rm -f " << opdir << "/" << file << ".ps" << endl; 
	system(ss3.str().c_str());
		
	stringstream ss4; ss4 << "rm -f " << opdir << "/gnuplot.txt" << endl;
	system(ss4.str().c_str());

	stringstream ss5; ss5 << "rm -f " << opdir << "/line.csv" << endl;
	system(ss5.str().c_str());

	cout << desc << " are placed into '" << opdir << "/" << file << ".pdf'" << endl;
}


/// Generates a text description of the pdf file
void Output::generate_pdf_description(const vector <OutputPlot> &op, const string grfile) const
{
	auto file = details.output_directory+"/"+grfile+"_description.txt";
	ofstream desc(file);
	if(!desc) emsg("Cannot open the file '"+file+"'");
	
	desc << "This file provides a description of the source files used to generate the graphs in '";
	desc << details.output_directory << "/" << file << ".pdf'" << endl << endl;
	
	for(auto i = 0u; i < op.size(); i++){
		const auto &oppl = op[i];
		desc << "PAGE " << i+1 << ": \"" << oppl.title << "\"" << endl;
		
		for(auto j = 0u; j < oppl.line.size(); j++){
			const auto &line = oppl.line[j];
			if(line.name == "") desc << "  This ";
			else desc << "  '" << line.name << "' ";
			
			switch(oppl.type){
			case OP_SPLINE: case OP_GRAPH: case OP_PARAMDIST: case OP_GENERATION: case OP_LOG_GENERATION:
			case OP_CPU: case OP_TRACE: case OP_ME:
				desc << "uses columns " << line.xcol << " and " << line.ycol << " in '";
				break;
			case OP_GRAPH_MARGINAL:	
				desc << "is generated from ";
				break;
			case OP_MARGINAL:	
				desc << "uses columns 1, 2 and " << line.ycol << " in '";
				break;
			}
			desc << line.file << "'" << endl;
		}
		
		if(oppl.legend_labelson){
			desc << endl << "  ";
			for(auto i =0u; i < oppl.label.size(); i++){
				if(i != 0) desc << ", ";
				desc << reference_label(i) << ": " << oppl.label[i];
			}
			desc << endl;
		}
		desc << endl;
	}
	
	cout << "The source files for these graphs are referenced in '" << file << "'" << endl << endl;
}


/// Saves a file giving time variation in different model quantities
void Output::spline_plots(const vector <Sample> &opsamp, vector <OutputPlot> &op) const
{
	if(opsamp.size() == 0) return;

	vector <SplineOutput> sim_spline_output;
	
	if(plot_sim_param == true){  
		if(details.mode == SIM) sim_spline_output = opsamp[0].spline_output;
		else{
			vector <double> vec; for(auto par : model.param) vec.push_back(par.value);	
			sim_spline_output = model.get_spline_output(vec,vector < vector < vector < vector <double> > > > ());
		}
	}
	
	auto nspline_out = opsamp[0].spline_output.size();
	for(auto sp = 0u; sp < nspline_out; sp++){
		const auto &splout = opsamp[0].spline_output[sp];
		auto name = splout.name;
		auto file = post_dir+"/spline/"+name+".csv";
		auto filefull = details.output_directory+"/"+file;
		ofstream tout(filefull.c_str());
		if(!tout) emsg("Cannot open the file '"+filefull+"'");
	
		if(suppress_output == false) cout << "'" << file << "' gives the time variation in "+name+"." << endl;

		bool plot_sim = false;
		if(plot_sim_param == true && sim_spline_output[sp].splineval.size() > 0) plot_sim = true;
	
		tout << "# Gives the time variation in "+name+"." << endl;	
		tout << "Time,mean,95% CI min,95% CI max,Simulated" << endl;

		auto min = LARGE, max = -LARGE;
		for(auto st = 0u; st < details.ndivision; st++){
			vector <double> vec; 
			for(const auto &opsa : opsamp){
				vec.push_back(opsa.spline_output[sp].splineval[st]);
			}
			auto stat = get_statistic(vec);
			
			tout << details.division_time[st] << "," << stat.mean << "," << stat.CImin << "," << stat.CImax;
			
			double num = stat.CImax; if(num > max) max = num; num = stat.CImin; if(num < min) min = num;
			
			if(plot_sim == true) tout << "," << sim_spline_output[sp].splineval[st];
			tout << endl; 
		}
	
		if(name != "UNSET"){                                                // Generates the final plot
			OutputPlot oppl(OP_SPLINE,splout.desc,"Time",name,min,max);
			string lab;
			if(post_lab != "Simulation"){
				oppl.addline(post_lab,filefull,1,2,RED_SOLID);
				oppl.addline("95% CI min",filefull,1,3,RED_DASHED);
				oppl.addline("95% CI max",filefull,1,4,RED_DASHED);
			  lab = "True";
			}
			
			if(plot_sim == true) oppl.addline(lab,filefull,1,5,BLACK_DASHED);
			op.push_back(oppl);
		}
	}
}


/// Plots a graph giving time variation in transition rate or population or marginal distributions
void Output::graph_plots(const vector <Sample> &opsamp, vector <OutputPlot> &op) const
{
	if(opsamp.size() == 0) return;
	
	vector <GraphMultiPlot> graph_plot;
	
	for(auto gr_num = 0u; gr_num < data.graph.size(); gr_num++){
		auto gr = data.graph[gr_num];
		
		const auto &dt = data.datatable[gr.datatable];
		
		auto file = details.output_directory+"/"+post_dir+"/state/"+gr.file;

		ofstream stateout(file.c_str());
		if(!stateout) emsg("Cannot open the file '"+file+"'");
		
		switch(gr.type){
			case GRAPH_TIMESERIES: stateout << "Time"; break;
			case GRAPH_MARGINAL: stateout << "Ref. Num.,Demographic state"; break;
		} 
		stateout << ",Mean,95% CI min,95% CI max" << endl; 
		
		auto imax = opsamp[0].graph_state[gr_num].size();
		auto min = LARGE, max = -LARGE;
		for(auto i = 0u; i < imax; i++){
			vector <double> vec; for(const auto &opsa : opsamp) vec.push_back(opsa.graph_state[gr_num][i]);		
			auto stat = get_statistic(vec);
			
			switch(gr.type){
				case GRAPH_TIMESERIES: stateout << details.division_time[i*graph_step+graph_step/2]; break;				
				case GRAPH_MARGINAL: stateout << i << "," << dt.demolist[i]; break;
			}
			
			stateout << "," << stat.mean << "," << stat.CImin << "," << stat.CImax << endl; 
			
			if(stat.CImax > max) max = stat.CImax;
			if(stat.CImin < min) min = stat.CImin;
		}
		
		auto thresh_flag = false;
			
		string file_data, file_data_thresh;
		if(dt.optype == DATA){
			file_data = details.output_directory+"/"+post_dir+"/state/"+gr.file.substr(0,gr.file.length()-4)+"-data.csv";
	
			ofstream dataout(file_data.c_str());
			if(!dataout) emsg("Cannot open the file '"+file_data+"'");
			
			if(gr.type == GRAPH_MARGINAL){
				dataout << "Ref. Num,Demographic state,Data,Threshold" << endl;
			}
			else{
				dataout << "Time,Data" << endl;
			}
			
			if(suppress_output == false) cout << file_data << " file" << endl;
		
			file_data_thresh = details.output_directory+"/"+post_dir+"/state/"+gr.file.substr(0,gr.file.length()-4)+"-data_thresh.csv";
		
			ofstream datathreshout;
			if(gr.type != GRAPH_MARGINAL){
				for(const auto &p : gr.point){ if(data.obs[p.obs].value == THRESH) thresh_flag = true;}
				
				if(thresh_flag == true){
					datathreshout.open(file_data_thresh.c_str());
					if(!datathreshout) emsg("Cannot open the file '"+file_data_thresh+"'");
					datathreshout << "Time,Data" << endl;
				}
			}
			
			for(auto i = 0u; i < gr.point.size(); i++){
				const auto &p = gr.point[i];
				auto ob = p.obs;
				
				auto val = data.obs[ob].value;
				
				auto shift = 0.0; if(details.time_format == TIME_FORMAT_NUM) shift = details.start;
				
				switch(dt.type){
					case MARGINAL:
						if(val == THRESH)	dataout << p.xi << "," << dt.demolist[i] << "," << dt.threshold << endl;
						else{
							if(val == UNKNOWN) dataout << p.xi << "," << dt.demolist[i] << "," << endl;
							else dataout << p.xi << "," << dt.demolist[i] << "," << val << endl;
						}
						break;
						
					case TRANS:
						if(val == THRESH){
							dataout << endl;
							if(details.stochastic == true){
								datathreshout << p.xi + shift << "," << dt.threshold << endl;
								datathreshout << p.xf + shift << "," << dt.threshold << endl;
							}
							else{
								datathreshout << 0.5*(p.xi+p.xf) + shift << "," << val << endl;
							}
						}
						else{
							if(val == UNKNOWN){
								dataout << endl;
								if(thresh_flag == true) datathreshout << endl;
							}
							else{
								val /= p.xf - p.xi;
								if(details.stochastic == true){
									dataout << p.xi + shift << "," << val << endl;
									dataout << p.xf + shift << "," << val << endl;
								}
								else{
									dataout << 0.5*(p.xi+p.xf) + shift << "," << val << endl;
								}
								if(thresh_flag == true) datathreshout << endl;
							}
						}
						break;
					
					case POP: case POPFRAC:
						if(val == THRESH){
							dataout << endl;
							datathreshout << 0.5*(p.xi+p.xf) + shift << "," << dt.threshold << endl;
						}
						else{
							if(val == UNKNOWN){
								dataout << endl;
								if(thresh_flag == true) datathreshout << endl;
							}
							else{
								dataout << 0.5*(p.xi+p.xf) + shift << "," << val << endl;
								if(thresh_flag == true) datathreshout << endl;
							}
						}
						break;
				}
				
				if(val != THRESH && val != UNKNOWN){
					if(val > max) max = val; if(val < min) min = val;
				}
			}
		}

		bool new_plot = true;
		
		auto plot_name = dt.plot_name;
		
		if(plot_name == "") plot_name = gr.desc;
		else{
			auto i = 0u; while(i < graph_plot.size() && graph_plot[i].plot_name != plot_name) i++;
			if(i < graph_plot.size() && graph_plot[i].type == GRAPH_TIMESERIES && gr.type == GRAPH_TIMESERIES){
				GraphMultiPlot &gp = graph_plot[i];
				if(min < gp.min) gp.min = min;
				if(max > gp.max) gp.max = max;
				gp.name.push_back(gr.name);
				gp.file.push_back(file);
				gp.line_colour.push_back(dt.line_colour);
				if(gp.line_colour.size() > 8) emsgroot("When outputting graphs cannot have more than 8 lines in one plot");
				
				if(file_data != ""){
					gp.file_data.push_back(file_data);
					if(thresh_flag == true) gp.file_data_thresh.push_back(file_data_thresh);
				}
				new_plot = false;
			}
		}
		
		if(new_plot == true){
			GraphMultiPlot gp;
			gp.plot_name = plot_name;
			gp.name.push_back(gr.name);
			if(!(gr.type == GRAPH_MARGINAL && post_lab == "Simulation")) gp.file.push_back(file);
			gp.type = gr.type;
			gp.line_colour.push_back(dt.line_colour);
			gp.min = min;
			gp.max = max;
			if(file_data != ""){
				gp.file_data.push_back(file_data);
				if(thresh_flag == true) gp.file_data_thresh.push_back(file_data_thresh);
			}
			graph_plot.push_back(gp);
		}
	}
	
	
	for(const auto &gp : graph_plot){	
		auto min = gp.min, max = gp.max;
		if(min == max) max++;
		vector <LineType> lt, lt2;
		get_line_colours(gp.line_colour,lt,lt2);	
				
		switch(gp.type){                                                // Generates the final plot
			case GRAPH_TIMESERIES:
				{
					auto ylabel = gp.name[0]; if(gp.file.size() > 1) ylabel = "Value";
					
					OutputPlot oppl(OP_GRAPH,gp.plot_name,"Time",ylabel,min,max);
					
					for(auto i = 0u; i < gp.file.size(); i++){
						auto legend = post_lab; if(gp.file.size() > 1) legend = gp.name[i];
						auto file = gp.file[i];
						if(post_lab == "Simulation"){
							oppl.addline(legend,file,1,2,lt[i]);
						}
						else{
							oppl.addline(legend,file,1,2,lt[i]);
							oppl.addline("95% CI min",file,1,3,lt2[i]);
							oppl.addline("95% CI max",file,1,4,lt2[i]);
						}
					}
					
					for(auto file_data : gp.file_data) oppl.addline("Data",file_data,1,2,BLACK_SOLID);
					
					for(auto file_data_thresh : gp.file_data_thresh) oppl.addline("Threshold",file_data_thresh,1,2,BLACK_DOTTED);
					
					op.push_back(oppl);
				}
				break;
				
			case GRAPH_MARGINAL:
				{
					if(gp.file.size() != 1) emsgEC("Output",10);
					OutputPlot oppl(OP_GRAPH_MARGINAL,gp.plot_name,"Category",gp.name[0],min,max);
					oppl.addline(post_lab,gp.file[0],1,3,lt[0]);
					
					for(auto file_data : gp.file_data) oppl.addline("Data",file_data,1,3,BLACK_SOLID);
					op.push_back(oppl);
				}
				break;
		}
	}
}


/// Works out the colours for each of the lines
void Output::get_line_colours(vector <LineColour> line_colour, vector <LineType> &lt, vector <LineType> &lt2) const
{
	for(auto i = 0u; i < line_colour.size(); i++){
		if(line_colour[i] == UNSET_COLOUR){
			unsigned int j;
			for(j = 1; j < BLACK; j++){
				auto ii = 0u; while(ii < line_colour.size() && line_colour[ii] != j) ii++;
				if(ii == line_colour.size()){
					line_colour[i] = (LineColour)j; 
					break;
				}
			}
			if(j == BLACK) emsgEC("Output",4); 
		}
	}
	
	for(auto i = 0u; i < line_colour.size(); i++){
		switch(line_colour[i]){
			case RED: lt.push_back(RED_SOLID); lt2.push_back(RED_DASHED); break;
			case GREEN: lt.push_back(GREEN_SOLID); lt2.push_back(GREEN_DASHED); break;
			case BLUE: lt.push_back(BLUE_SOLID); lt2.push_back(BLUE_DASHED); break;
			case YELLOW: lt.push_back(YELLOW_SOLID); lt2.push_back(YELLOW_DASHED); break;
			case CYAN: lt.push_back(CYAN_SOLID); lt2.push_back(CYAN_DASHED); break;
			case MAGENTA: lt.push_back(MAGENTA_SOLID); lt2.push_back(MAGENTA_DASHED); break;
			case BLACK: lt.push_back(BLACK_SOLID); lt2.push_back(BLACK_DASHED); break;
			default: emsgEC("Output",5); break;
		}
	}
}
	
	
/// Outputs a bar chart of relative susceptibility split into demographic groups   
void Output::susceptibility_distributions(const vector <ParamSample> &psamp, vector <OutputPlot> &op) const 
{
	for(auto c = 0u; c < data.ndemocat; c++){
		auto amax = data.democat[c].value.size();
		if(data.democat[c].sus_vari == true && amax > 1){
			auto type = data.democat[c].name;
			
			auto file = type+".csv";
			auto filefull = details.output_directory+"/"+post_dir+"/susceptibility/"+file;
			ofstream distout(filefull.c_str());
			if(!distout) emsg("Cannot open the file '"+filefull+"'");
			
			if(suppress_output == false) cout << "'" << file << "' gives stratified susceptibility for " << type << "." << endl;
			
			distout << "# Stratified susceptibility for " << type << "." << endl;
		
			vector <string> label;
			distout << "Ref.Num.,Ref.Lab," << type << ",Mean,CI min,CImax";
			if(plot_sim_param == true) distout << ",Simulated";
			distout << endl;					
			for(auto a = 0u; a < amax; a++){
				distout << a << "," << reference_label(a) << "," << data.democat[c].value[a];
			
				auto th = data.democat[c].sus_param[a];
				vector <double> vec; for(const auto &psa : psamp) vec.push_back(psa.paramval[th]);
				auto stat = get_statistic(vec);
		
				distout << "," << stat.mean << "," << stat.CImin << "," << stat.CImax;
				if(plot_sim_param == true) distout << "," << model.param[th].value;
				distout << endl;
			}
			
			OutputPlot oppl(OP_MARGINAL,"Susceptibility with "+type,data.democat[c].name,"Relative susceptibility",0,0);
			oppl.label = label;
			if(post_lab != "Simulation") oppl.addline(post_lab,filefull,1,4,RED_SOLID);
			if(plot_sim_param == true) oppl.addline("True",filefull,1,7,BLACK_SOLID);
			op.push_back(oppl);
		}
	}
}


/// Gets information about a dependency from a string 
vector <DemoDep> Output::get_demodep(const string dep) const 
{
	vector <DemoDep> demodep;
	
	vector <unsigned int> dep_order;
	auto demo_ref = inputs.find_demo_ref(dep,data.democat,dep_order);

	for(auto k = 0u; k < demo_ref.size(); k++){
		string name;
		for(auto c : dep_order){
			if(name != "") name += " ";
			name += data.democat[c].value[demo_ref[k][c]];
		}		
		
		unsigned int dp;
		for(dp = 0u; dp < data.ndemocatpos; dp++){
			auto i = 0u; while(i < data.ndemocat && (data.democatpos[dp][i] == demo_ref[k][i] || demo_ref[k][i] == UNSET)) i++;
			if(i == data.ndemocat) break;
		}
		if(dp == data.ndemocatpos) emsgEC("Output",1);
		
		DemoDep dedp; dedp.name = name; dedp.dp = dp;
		demodep.push_back(dedp);
	}
	
	return demodep;
}

	
/// Outputs a bar chart giving compartmental means split by demographic group
void Output::mean_distributions(const vector <ParamSample> &psamp, vector <OutputPlot> &op) const 
{
	for(const auto &co : model.comp){
		auto dep = co.mean_dep;
		if(co.num == 0 && dep != ""){ 
			auto dir = details.output_directory+"/"+post_dir+"/compartment_mean";
			ensure_directory(dir);
			auto file = co.name+"_"+dep+".csv";
			auto filefull = dir+"/"+file;
		
			ofstream distout(filefull.c_str());
			if(!distout) emsg("Cannot open the file '"+filefull+"'");
	
			if(suppress_output == false) cout << "'" << file << "' gives stratified compartment '" << co.name << "' mean residency time." << endl;
	
			distout << "# Stratified mean compartmental residency time for compartment '" << co.name << "'." << endl;

			distout << "Ref.Num.,Ref.Lab,Demographic state";
			distout << ",Mean,95% CI Min,95% CI Max";
			if(plot_sim_param == true) distout << ",Simulated";
			distout << endl;
		
			vector <string> label;
			auto max = -LARGE, min = LARGE;
			auto demodep = get_demodep(dep);
			for(auto i = 0u; i < demodep.size(); i++){
				distout << i << "," << reference_label(i) << "," << demodep[i].name;
					
				auto th = co.param_mean[demodep[i].dp];
			
				vector <double> vec; for(const auto &psa : psamp) vec.push_back(psa.paramval[th]);
				auto stat = get_statistic(vec);

				distout << "," << stat.mean << "," << stat.CImin << "," << stat.CImax;	
				if(stat.CImax > max) max = stat.CImax;
				if(stat.CImin < min) min = stat.CImin;
				
				if(plot_sim_param == true){
					auto num = model.param[th].value; if(num > max) max = num; if(num < min) min = num;
					distout << "," << num;
				}
				distout << endl;
			}
			distout << endl;

			if(min != max){
				OutputPlot oppl(OP_MARGINAL,"Compartment "+co.name+" mean residency time","Demographic classification","Probability",0,0);
				oppl.label = label;
				if(post_lab != "Simulate") oppl.addline(post_lab,filefull,1,4,RED_SOLID);
				if(plot_sim_param == true) oppl.addline("True",filefull,1,7,BLACK_SOLID);
				op.push_back(oppl);	
			}
		}
	}
}


/// Outputs a bar chart giving branching probabilities split by demographic groups
void Output::branch_prob_distributions(const vector <ParamSample> &psamp, vector <OutputPlot> &op) const 
{
	for(const auto tr : model.trans){
		auto dep = tr.prob_dep;
		if(dep != ""){	
			auto dir = details.output_directory+"/"+post_dir+"/branch_prob";
			auto file = tr.name_file+"_"+dep+".csv";
			ensure_directory(dir);
			auto filefull = dir+"/"+file;
		
			ofstream distout(filefull.c_str());
			if(!distout) emsg("Cannot open the file '"+filefull+"'");
	
			if(suppress_output == false) cout << "'" << file << "' gives stratified branching probability for transition '" << tr.name << "'." << endl;
	
			distout << "# Stratified branching probability for transition '" << tr.name << "'." << endl;

			distout << "Ref.Num.,Ref.Lab,Demographic state";
			distout << ",Mean,95% CI Min,95% CI Max";
			if(plot_sim_param == true) distout << ",Simulated";
			distout << endl;
		
			vector <string> label;
			auto max = -LARGE, min = LARGE;
			auto demodep = get_demodep(dep);
			for(auto i = 0u; i < demodep.size(); i++){
				distout << i << "," << reference_label(i) << "," << demodep[i].name;
					
				auto th = tr.param_prob[demodep[i].dp];
				
				vector <double> vec; for(const auto &psa : psamp) vec.push_back(psa.paramval[th]);
				auto stat = get_statistic(vec);

				distout << "," << stat.mean << "," << stat.CImin << "," << stat.CImax;	
				if(stat.CImax > max) max = stat.CImax;
				if(stat.CImin < min) min = stat.CImin;
				
				if(plot_sim_param == true){
					auto num = model.param[th].value; if(num > max) max = num; if(num < min) min = num;
					distout << "," << num;
				}
				distout << endl;
			}
			distout << endl;

			if(min != max){
				OutputPlot oppl(OP_MARGINAL,"Branching probability for transition '"+tr.name+"'","Demographic classification","Probability",0,0);
				oppl.label = label;
				if(post_lab != "Simulation") oppl.addline(post_lab,filefull,1,4,RED_SOLID);
				if(plot_sim_param == true) oppl.addline("True",filefull,1,7,BLACK_SOLID);
				op.push_back(oppl);	
			}
		}
	}
}


/// Outputs a bar chart giving level effects split by demographic groups
void Output::level_effect_distributions(const vector <ParamSample> &psamp, vector <OutputPlot> &op) const 
{
	if(data.level_effect.on == false) return;

	auto file = post_dir+"/level_effect.csv";
	auto filefull = details.output_directory+"/"+file;
		
	ofstream distout(filefull.c_str());
	if(!distout) emsg("Cannot open the file '"+filefull+"'");
	
	if(suppress_output == false) cout << "'" << file << "' gives level effects." << endl;

	distout << "# Level effects." << endl;

	distout << "Ref.Num.,Ref.Lab,Parameter,Mean,95% CI Min,95% CI Max";
	if(plot_sim_param == true) distout << ",Simulated";
	distout << endl;
		
	auto max = -LARGE, min = LARGE;
		
	vector <string> label;
	const auto &list = model.level_effect_param_list;
	for(auto i = 0u; i < list.size(); i++){
		auto th = list[i];
		vector <double> vec; for(const auto &psa : psamp) vec.push_back(psa.paramval[th]);
		auto stat = get_statistic(vec);
		
		distout << i << "," << reference_label(i) << "," << model.param[th].name << ",";
		distout << stat.mean << "," << stat.CImin << "," << stat.CImax;	
		if(plot_sim_param == true){
			auto num = model.param[th].value; if(num > max) max = num; if(num < min) min = num;
			distout << "," << num;
		}
		distout << endl;
			
		if(stat.CImax > max) max = stat.CImax;
		if(stat.CImin < min) min = stat.CImin;
	}
	
	if(min != max){
		OutputPlot oppl(OP_MARGINAL,"Relative transmission for different levels","Level Effect","Relative transmission rate",0,0);
		oppl.label = label;
		if(post_lab != "Simulation") oppl.addline(post_lab,filefull,1,4,RED_SOLID);
		if(plot_sim_param == true) oppl.addline("True",filefull,1,7,BLACK_SOLID);
		op.push_back(oppl);	
	}
}


/// Outputs a bar chart giving relative transmission split by area
void Output::area_effect_distributions(const vector <ParamSample> &psamp, vector <OutputPlot> &op) const 
{
	if(data.area_effect.on == false) return;

	auto file = post_dir+"/area_effect.csv";
	auto filefull = details.output_directory+"/"+file;
		
	ofstream distout(filefull.c_str());
	if(!distout) emsg("Cannot open the file '"+filefull+"'");
	
	if(suppress_output == false) cout << "'" << file << "' gives area effects." << endl;

	distout << "# Area effects." << endl;

	distout << "Ref.Num.,Ref.Lab,Parameter,Mean,95% CI Min,95% CI Max";
	if(plot_sim_param == true) distout << ",Simulated";
	distout << endl;
		
	auto max = -LARGE, min = LARGE;
		
	vector <string> label;
	const auto &list = model.area_effect_param_list;
	for(auto i = 0u; i < list.size(); i++){
		auto th = list[i];
		vector <double> vec; for(const auto &psa : psamp) vec.push_back(psa.paramval[th]);
		auto stat = get_statistic(vec);
		
		label.push_back(model.param[th].name);
		distout << i << "," << reference_label(i) << ",";
		distout << model.param[th].name << "," << stat.mean << "," << stat.CImin << "," << stat.CImax;	
		if(plot_sim_param == true){
			auto num = model.param[th].value; if(num > max) max = num; if(num < min) min = num;
			distout << "," << num;
		}
		distout << endl;
			
		if(stat.CImax > max) max = stat.CImax;
		if(stat.CImin < min) min = stat.CImin;
	}
	
	if(min != max){
		OutputPlot oppl(OP_MARGINAL,"Relative transmission for different areas","Area Effect","Relative transmission rate",0,0);
		oppl.label = label;
		if(post_lab != "Simulation") oppl.addline(post_lab,filefull,1,4,RED_SOLID);
		if(plot_sim_param == true) oppl.addline("True",filefull,1,7,BLACK_SOLID);
		op.push_back(oppl);	
	}
}


/// Outputs the probability distributions for derived quantities (generation time, external infections etc..)
void Output::derived_parameter_distributions(const vector <Sample> &opsamp, vector <OutputPlot> &op) const
{
	if(opsamp.size() == 0) return;
	
	const auto &derived_param = opsamp[0].derived_param;
	for(auto dp = 0u; dp < derived_param.size(); dp++){
		vector <double> vec; for(const auto &opsa : opsamp) vec.push_back(opsa.derived_param[dp].value);
		
		auto paramdist = get_distribution(vec,UNSET,UNSET);

		if(paramdist.variation == true){
			const auto &derpar = derived_param[dp];
			
			string file = derpar.file;
			string filefull = details.output_directory+"/"+post_dir+"/parameter/"+file;
			ofstream distout(filefull.c_str());
			if(!distout) emsg("Cannot open the file '"+filefull+"'");
			
			if(suppress_output == false){
				cout << "'" << file << "' gives the probability distributions for the "+ derpar.desc+"." << endl;
			}
			
			distout << "# " << post_lab << " probability distributions for the generation time." << endl;

			
			distout << "Value,Probability" << endl;
			
			for(auto b = 0u; b < paramdist.value.size(); b++){
				distout << paramdist.value[b] << "," << paramdist.prob[b] << endl;
			}
			
			OutputPlot oppl(OP_PARAMDIST,derpar.desc,derpar.name,"Probability",0,0);
			oppl.addline(post_lab,filefull,1,2,GREEN_SOLID);
			op.push_back(oppl);
		}
	}
}


/// Outputs a file containing posterior distributions for model parameters
void Output::posterior_parameter_distributions(const vector <ParamSample> &psamp, vector <OutputPlot> &op) const
{
	if(model.param.size() == 0) return;
	
	auto file = "model.csv";
	auto filefull = details.output_directory+"/"+post_dir+"/parameter/"+file;
	ofstream distout(filefull.c_str());
	if(!distout) emsg("Cannot open the file '"+filefull+"'");
	
	if(suppress_output == false) cout << "'" << file << "' gives the probability distributions for parameters." << endl;
	
  distout << "# " << post_lab << " probability distributions for model parameters." << endl;
	
	vector <unsigned int> pvary;
	
	vector <Distribution> paramdist(model.param.size());
	for(auto p = 0u; p < model.param.size(); p++){
		vector <double> vec; for(const auto &psa : psamp) vec.push_back(psa.paramval[p]);
		double priormin = UNSET, priormax = UNSET;
		if(model.param[p].priortype == UNIFORM_PRIOR){ priormin = model.param[p].val1; priormax = model.param[p].val2;}
		paramdist[p] = get_distribution(vec,priormin,priormax);
		
		if(paramdist[p].variation == true) pvary.push_back(p);
	}

	for(auto i = 0u; i < pvary.size(); i++){
		auto p = pvary[i];
		distout << model.param[p].name << " (value)" << "," << model.param[p].name << " (probability)";
		if(i < pvary.size()-1) distout << ",";
	}	
	distout << endl; 
	
	if(pvary.size() == 0){
		distout << "Parameters did not vary" << endl;
		return;
	}
	
	for(auto b = 0u; b < paramdist[0].value.size(); b++){
		for(auto i = 0u; i < pvary.size(); i++){
			auto p = pvary[i];
			distout << paramdist[p].value[b] << "," << paramdist[p].prob[b];
			if(i < pvary.size()-1) distout << ",";
		}
		distout << endl;
	}
	
	for(auto i = 0u; i < pvary.size(); i++){
		auto p = pvary[i];
	
		OutputPlot oppl(OP_PARAMDIST,post_lab+" distribution for "+model.param[p].name,model.param[p].name,"Probability",0,0);
		oppl.addline(post_lab,filefull,2*i+1,2*i+2,GREEN_SOLID);
		if(plot_sim_param == true) oppl.addline("True",to_string(model.param[p].value),1,2,BLACK_DASHED);
		op.push_back(oppl);
	}
}


/// Adds the plots which show how quantities change as a function of generation number (ABC-SMC / ABC-MBP)
void Output::add_generation_plots(vector <OutputPlot> &op) const
{
	if(details.mode != ABC_SMC && details.mode != ABC_MBP && details.mode != PAIS_INF) return;
	
	auto file = details.output_directory+"/Diagnostics/Generation.csv";

	{
		OutputPlot oppl(OP_LOG_GENERATION,"EF cut-off as a function of generation","Generation","EF cut-off",UNSET,UNSET);
		oppl.addline("EF cut-off",file,1,3,BLACK_SOLID);
		op.push_back(oppl);
	}
	
	{
		string title, label;
		switch(details.mode){
			case ABC_MBP: case ABC_SMC: title = "Model evidence as a function of EF cut-off"; label = "EF cut-off"; break;
			case PAIS_INF: title = "Model evidence as a function of invT"; label = "invT"; break;
			default: emsgEC("Output",2); break;
		}

		OutputPlot oppl(OP_ME,title,label,"log(Model evidence)",UNSET,UNSET);
		oppl.addline("",file,3,5,BLACK_SOLID);
		oppl.addline("6",file,3,5,BLACK_SOLID);
		op.push_back(oppl);
	}
	
	for(auto p = 0u; p < model.param.size(); p++){                           // These are graphs for the parameters
		if(model.param[p].priortype != FIXED_PRIOR){
			double min = UNSET, max = UNSET; 
			if(model.param[p].priortype == UNIFORM_PRIOR){
				min = model.param[p].val1; max = model.param[p].val2;
			}
			
			OutputPlot oppl(OP_GENERATION,post_lab+" estimate for "+model.param[p].name+" as a function of generation","Generation",model.param[p].name,min,max);
			oppl.addline(post_lab,file,1,3*p+7,BLUE_SOLID);
			oppl.addline("95% CI min",file,1,3*p+8,BLUE_DASHED);
			oppl.addline("95% CI max",file,1,3*p+9,BLUE_DASHED);
			if(plot_sim_param == true) oppl.addline("True",to_string(model.param[p].value),1,2,BLACK_DASHED);
			op.push_back(oppl);
		}
	}
	
	auto num = 0u; for(auto i = 0u; i < data.datatable.size(); i++){ if(data.datatable[i].weight != 0) num++;}			
	if(num > 1){
		auto fileDT = details.output_directory+"/Diagnostics/EF_datatable.csv";
		OutputPlot opplDT(OP_LOG_GENERATION,"Contribution to EF from different data sources","Generation","EF contribution",UNSET,UNSET);
		for(auto i = 0u; i < data.datatable.size(); i++){
			if(data.datatable[i].weight != 0) opplDT.addline(rus(data.datatable[i].file),fileDT,1,i+2,get_linestyle(i));
		}
	
		op.push_back(opplDT);
	}

	{
		OutputPlot oppl(OP_CPU,"CPU time as a function of EF","EF cut-off","CPU time (minutes)",UNSET,UNSET);
		oppl.addline("",file,3,2,BLACK_SOLID);
		op.push_back(oppl);
	}
}


/// Detemines the order of line usage when plotting serveral lines
LineType Output::get_linestyle(const unsigned int i) const
{
	switch(i%16){
		case 0: return RED_SOLID;
		case 1: return GREEN_SOLID;
		case 2: return BLUE_SOLID;
		case 3: return BLACK_SOLID;
		case 4: return RED_DASHED;
		case 5: return GREEN_DASHED;
		case 6: return BLUE_DASHED;
		case 7: return BLACK_DASHED;
		case 8: return RED_DOTTED;
		case 9: return GREEN_DOTTED;
		case 10: return BLUE_DOTTED;
		case 11: return BLACK_DOTTED;
		case 12: return RED_DOTDASH;
		case 13: return GREEN_DOTDASH;
		case 14: return BLUE_DOTDASH;
		case 15: return BLACK_DOTDASH;
	}
	return BLACK_SOLID;
}


/// Adds the trace plots (PMCMC / MC3)
void Output::add_trace_plots(vector <OutputPlot> &op) const
{
	if(details.mode != MC3_INF && details.mode != MCMC_MBP && details.mode != PMCMC_INF) return;
	
	for(auto p = 0u; p < model.param.size(); p++){                           // These are graphs for the parameters
		if(model.param[p].priortype != FIXED_PRIOR){
			OutputPlot oppl(OP_TRACE,"","Sample",model.param[p].name,UNSET,UNSET);
			
			for(auto ru = 0u; ru < nrun; ru++){
				if(nrun == 1){
					auto file = details.output_directory+"/Diagnostics/Trace.csv";
					oppl.addline("",file,1,2+p,get_linestyle(ru));	
				}
				else{
					auto file = details.output_directory+"/Diagnostics/Trace_Run"+to_string(ru+1)+".csv";
					oppl.addline("Run "+to_string(ru+1),file,1,2+p,get_linestyle(ru));	
				}
			}
			if(plot_sim_param == true) oppl.addline("True",to_string(model.param[p].value),1,2,BLACK_DASHED);
			op.push_back(oppl);
		}
	}
}


/// Adds graphs giving any demographic change 
void Output::add_democat_change(vector <OutputPlot> &op) const
{
	for(const auto &dcc : data.democat_change){
		auto d = dcc.d;
		auto desc = "The change in "+data.democat[d].name+" with time";
	
		auto file =	dcc.name;
		if(dcc.democats_filt != ""){ file += "_"+dcc.democats_filt; desc += " for "+dcc.democats_filt;}
		if(dcc.geo_filt != ""){ file += "_"+dcc.geo_filt; desc += " for "+dcc.geo_filt;}
		file += ".csv";
		auto filefull = details.output_directory+"/"+post_dir+"/democat_change/"+file;
		
		ofstream lout(filefull.c_str());
		if(!lout) emsg("Cannot open the file '"+filefull+"'");
	
		lout << "# " << desc << endl;
		lout << "Time";
		for(auto v = 0u; v < data.democat[d].value.size(); v++) lout << "," << data.democat[d].value[v]; 
		lout << endl;
		
		for(auto sett = 0u; sett < details.ndivision; sett++){
			lout << details.division_time[sett]; 
			for(auto v = 0u; v < data.democat[d].value.size(); v++) lout << "," << 100*dcc.frac[sett][v];
			lout << endl;
		}
	
		OutputPlot oppl(OP_GRAPH,desc,"Time","Percent",0,100);
		for(auto v = 0u; v < data.democat[d].value.size(); v++){
			oppl.addline(data.democat[d].value[v],filefull,1,2+v,get_linestyle(v));
		}
		op.push_back(oppl);
	}
}


// For marginals works out if labels go into description
void Output::set_labelon(vector <OutputPlot> &op) const
{
	for(auto i = 0u; i < op.size(); i++){
		auto &oppl = op[i];	
		if(oppl.label.size() != 0){            
			auto max = 0u;
			for(const auto &st : oppl.label){ if(st.length() > max) max = st.length();}
			if(max*oppl.label.size() > 100) oppl.legend_labelson = true;
		}
	}
}


/// Adds a time to pred_timeplot
void Output::add_pred_timeplot(const string name, unsigned int time, vector <TimePlot> &pred_timeplot) const 
{
	TimePlot tp;
	tp.name = name;
	tp.time = time;
	auto i = 0u; while(i < pred_timeplot.size() && pred_timeplot[i].time != tp.time) i++;
	if(i == pred_timeplot.size()) pred_timeplot.push_back(tp);
}


/// Gets pred_timeplot
vector <TimePlot> Output::get_pred_timeplot() const 
{
	vector <TimePlot> pred_timeplot;
	if(details.mode == PREDICTION){		
		switch(details.mode){
			case PREDICTION: add_pred_timeplot("Prediction Start",model.modelmod.pred_start,pred_timeplot); break;
			default: emsgEC("Output",3); break;
		}
	
		for(const auto &mod : data.modification){
			add_pred_timeplot("",mod.start,pred_timeplot);
			add_pred_timeplot("",mod.end,pred_timeplot);
		}
	}
	return pred_timeplot;
}


/// Generates a file for plotting output with gnuplot
void Output::generate_gnuplot_file(const vector <OutputPlot> &op, const string grfile) const
{
	auto opdir = details.output_directory; 
	
	auto linefile = opdir+"/line.csv";
	ofstream line(linefile);       // This provides a way of drawing vertical lines
	if(!line) emsg("Cannot open the file '"+linefile+"'");
	line << "1,0" << endl << "1,1" << endl;
	
	auto file = opdir+"/gnuplot.txt";
	ofstream gnuplot(file);
	if(!gnuplot) emsg("Cannot open the file '"+file+"'");

	gnuplot << "set terminal postscript enhanced color font ',20' dl 3" << endl;
	gnuplot << "set output '" << opdir << "/" << grfile << ".ps'" << endl;

	gnuplot << "set datafile separator ','" << endl;
	
	gnuplot << "set style line 1 lt 1 lc rgb '#ff2222' lw 5" << endl; // RED_SOLID
	gnuplot << "set style line 2 lt 2 lc rgb '#ffaaaa' lw 5" << endl; // RED_DASHED
	gnuplot << "set style line 3 lt 1 lc rgb '#22ff22' lw 5" << endl; // GREEN_SOLID
	gnuplot << "set style line 4 lt 2 lc rgb '#aaffaa' lw 5" << endl; // GREEN_DASHED
	gnuplot << "set style line 5 lt 1 lc rgb '#2222ff' lw 5" << endl; // BLUE_SOLID
	gnuplot << "set style line 6 lt 2 lc rgb '#aaaaff' lw 5" << endl; // BLUE_DASHED
	gnuplot << "set style line 7 lt 1 lc rgb '#222222' lw 5" << endl; // BLACK_SOLID
	gnuplot << "set style line 8 lt 2 lc rgb '#222222' lw 5" << endl; // BLACK_DASHED
	gnuplot << "set style line 9 lt 3 lc rgb '#ff2222' lw 5" << endl; // RED_DOTTED
	gnuplot << "set style line 10 lt 5 lc rgb '#ffaaaa' lw 5" << endl; // RED_DOTDASH
	gnuplot << "set style line 11 lt 4 lc rgb '#22ff22' lw 5" << endl; // GREEN_DOTTED
	gnuplot << "set style line 12 lt 5 lc rgb '#aaffaa' lw 5" << endl; // GREEN_DOTDASH
	gnuplot << "set style line 13 lt 4 lc rgb '#2222ff' lw 5" << endl; // BLUE_DOTTED
	gnuplot << "set style line 14 lt 5 lc rgb '#aaaaff' lw 5" << endl; // BLUE_DOTDASH
	gnuplot << "set style line 15 lt 4 lc rgb '#222222' lw 5" << endl; // BLACK_DOTTED
	gnuplot << "set style line 16 lt 5 lc rgb '#222222' lw 5" << endl; // BLACK_DOTDASH
	gnuplot << "set style line 17 lt 1 lc rgb '#ffff22' lw 5" << endl; // YELLOW_SOLID
	gnuplot << "set style line 18 lt 2 lc rgb '#ffff22' lw 5" << endl; // YELLOW_DASHED
	gnuplot << "set style line 19 lt 1 lc rgb '#22ffff' lw 5" << endl; // CYAN_SOLID
	gnuplot << "set style line 20 lt 2 lc rgb '#22ffff' lw 5" << endl; // CYAN_DASHED
	gnuplot << "set style line 21 lt 1 lc rgb '#ff22ff' lw 5" << endl; // MAGENTA_SOLID
	gnuplot << "set style line 22 lt 2 lc rgb '#ff22ff' lw 5" << endl; // MAGENTA_DASHED

	gnuplot << "set style line 80 lt 1 lc rgb '#222222' lw 2" << endl; // BLACK_THIN
	gnuplot << "set style line 81 lt 2 lc rgb '#ff2222' lw 2" << endl; // RED_THIN
	gnuplot << "set style line 82 lt 2 lc rgb '#22ff22' lw 2" << endl; // GREEN_THIN
	gnuplot << "set style line 83 lt 2 lc rgb '#2222ff' lw 2" << endl; // BLUE_THIN
	
	auto pred_timeplot = get_pred_timeplot();

	auto multiplot = UNSET;
	for(auto i = 0u; i < op.size(); i++){
		const auto &oppl = op[i];
			
		if(oppl.title == "") gnuplot << "set title" << endl;
		else gnuplot << "set title '" << label(oppl.title) << "'" << endl;
		gnuplot << "set autoscale" << endl;
		
		auto keyfontsize = 22u, labelfontsize = 22u, ticfontsize = 20u; 
		if(oppl.type == OP_TRACE){ keyfontsize = 16u; labelfontsize = 18; ticfontsize = 16u;}
		gnuplot << "set tics font ', " << ticfontsize << "'" << endl;
		gnuplot << "set key font ', " << keyfontsize << "'" << endl;
		gnuplot << "set xlabel '" << label(oppl.xaxis) << "' font '," << labelfontsize << "'" << endl;
		gnuplot << "set ylabel '" << label(oppl.yaxis) << "' font '," << labelfontsize << "'" << endl;
	
		if(oppl.type == OP_PARAMDIST)	gnuplot << "unset ytics" << endl;
		else gnuplot << "set ytics" << endl;
			
		if(multiplot == 2 || oppl.type != OP_TRACE){ multiplot = UNSET; gnuplot << "unset multiplot" << endl;}
				
		switch(oppl.type){
			case OP_SPLINE: case OP_GRAPH:
				{
					auto ymax = 1.2*oppl.max;
					gnuplot << "set yrange [0:" << ymax << "]" << endl;
		
					auto shift = details.start; if(details.time_format == TIME_FORMAT_NUM) shift = 0;
					
					for(const auto &tp : details.timeplot){
						gnuplot << "set arrow from " << tp.time - shift << ",0 to " << tp.time - shift << ","
										<< ymax << " lw 3 lc rgb '#000066' lt 1 nohead" << endl;
				
						gnuplot << "set label '" << label(tp.name) << "' at " << tp.time - shift << "," << ymax << " left rotate by 270 offset 1,-1" << endl;
					}
					
					for(const auto &tp : pred_timeplot){
						gnuplot << "set arrow from " << tp.time << ",0 to " << tp.time << ","
										<< ymax << " lw 3 lc rgb '#0000aa' ";

						if(tp.name != "") gnuplot << "lt 1"; else gnuplot << "lt 4"; 
						gnuplot << " nohead" << endl;
				
						if(tp.name != ""){
							gnuplot << "set label '" << label(tp.name) << "' at " << tp.time << "," << ymax << " tc rgb '#0000aa'  left rotate by 270 offset 1,-1" << endl;
						}
					}
					
					gnuplot << "plot ";
					for(auto j = 0u; j < oppl.line.size(); j++){
						const auto &line = oppl.line[j];
						if(j != 0) gnuplot << ", ";
						gnuplot << "'" << line.file << "' using " << line.xcol << ":" << line.ycol << " with lines ls " << line.style << " ";
						if(line.name == "" || line.name == "95% CI min" || line.name == "95% CI max") gnuplot << "notitle"; 
						else gnuplot << "title '" << label(line.name) << "'";
					}
					
					gnuplot << endl;
					gnuplot << "unset label"<< endl;
					gnuplot << "unset arrow"<< endl;
							
					if(false){
						for(auto j = 0u; j < details.timeplot.size(); j++){
							const auto &tp = details.timeplot[j];
							gnuplot << ", '" << opdir << "/line.csv' using ($1*" 
										<< tp.time - details.start << "):($2*"
										<< ymax << ") with lines ls " << 80+j%4 << " title '" << label(tp.name) << "'";
						}	
					}
				}
				break;
				
			case OP_GRAPH_MARGINAL:	
				gnuplot << "set autoscale" << endl;
				gnuplot << "set yrange [0:]" << endl;
				gnuplot << "set boxwidth 0.5" << endl;
				gnuplot << "set style fill solid" << endl;
				{
					auto col = oppl.line[0].ycol;
					gnuplot << "plot '" << oppl.line[0].file << "' using ($1+0.5):" << col << ":xtic(2) with boxes  notitle, ";
					if(details.mode != SIM){
						gnuplot << "'" << oppl.line[0].file << "' using ($1+0.5):" << col << ":" << col+1 << ":" << col+2;
						gnuplot << " with errorbars lt 1  lc rgb '#000000' lw 3 ps 2 title '" << label(oppl.line[0].name) << "'";
					}
					if(oppl.line.size() > 1){
						gnuplot << ",'" << oppl.line[1].file << "' using ($1+0.5):" << oppl.line[1].ycol 
											<< " with point lt 1 lc rgb '#0000ff' lw 5 pt 2 ps 2 title '" <<  label(oppl.line[1].name) << "'";
					}
				}
				break;
				
			case OP_MARGINAL:	
				gnuplot << "set autoscale" << endl;
				gnuplot << "set yrange [0:]" << endl;
				gnuplot << "set boxwidth 0.5" << endl;
				gnuplot << "set style fill solid" << endl;
				{
					auto lab = 3u; if(oppl.legend_labelson == true) lab = 2;
					auto col = oppl.line[0].ycol;
					gnuplot << "plot '" << oppl.line[0].file << "' using ($1+0.5):" << col << ":xtic(" << lab <<") with boxes notitle, ";
					if(details.mode != SIM){
						gnuplot << "'" << oppl.line[0].file << "' using ($1+0.5):" << col << ":" << col+1 << ":" << col+2;
						gnuplot << " with errorbars lt 1  lc rgb '#000000' lw 3 ps 2 title '" << label(oppl.line[0].name) << "', ";
					}
					if(oppl.line.size() > 1){
						gnuplot << "'" << oppl.line[1].file << "' using ($1+0.5):" << oppl.line[1].ycol 
											<< " with point lt 1 lc rgb '#0000ff' lw 5 pt 2 ps 2 title '" <<  label(oppl.line[1].name) << "'";
					}
				}
				break;
				
			case OP_PARAMDIST:
				gnuplot << "set yrange [0:1]" << endl;

				gnuplot << "plot ";
				for(auto j = 0u; j < oppl.line.size(); j++){
					const auto &line = oppl.line[j];
					if(j != 0) gnuplot << ", ";
					
					if(line.name == "True"){
						gnuplot << "'" << linefile  << "' using " << "($1*" << line.file << "):2";
					}
					else{
						gnuplot << "'" << line.file << "' using " << line.xcol << ":" << line.ycol;
					}
					gnuplot << " with lines ls " << line.style;
					
					if(line.name == "") gnuplot << " notitle"; else gnuplot << " title '" << label(line.name) << "'";
				}
				break;
				
			case OP_GENERATION:
				if(oppl.min != UNSET) gnuplot << "set yrange [" << oppl.min << ":" << oppl.max << "]" << endl;
	
				gnuplot << "plot ";
				for(auto j = 0u; j < oppl.line.size(); j++){
					const auto &line = oppl.line[j];
					if(j != 0) gnuplot << ", ";
					
					if(line.name == "True"){
						gnuplot << line.file;
					}
					else{
						gnuplot << "'" << line.file << "' using " << line.xcol << ":" << line.ycol;
					}
					
					gnuplot << " with lines ls " << line.style << " ";
					if(line.name == "" || line.name == "95% CI min" || line.name == "95% CI max") gnuplot << "notitle"; 
					else gnuplot << "title '" << label(line.name) << "'";
				}
				break;
				
			case OP_LOG_GENERATION:
				gnuplot << "set logscale y" << endl;
					
				if(oppl.min != UNSET) gnuplot << "set yrange [" << oppl.min << ":" << oppl.max << "]" << endl;
	
				gnuplot << "plot ";
				for(auto j = 0u; j < oppl.line.size(); j++){
					const auto &line = oppl.line[j];
					if(j != 0) gnuplot << ", ";
					gnuplot << "'" << line.file << "' using " << line.xcol << ":" << line.ycol;
					gnuplot << " with lines ls " << line.style << " ";
					if(line.name == "EF cut-off") gnuplot << "notitle"; 
					else gnuplot << "title '" << label(line.name) << "'";
				}
				gnuplot << endl;
				
				gnuplot << "unset logscale y" << endl;
				break;
				
			case OP_CPU:
				gnuplot << "set logscale y" << endl;
				gnuplot << "set logscale x" << endl;
					
				gnuplot << "plot ";
				for(auto j = 0u; j < oppl.line.size(); j++){
					const auto &line = oppl.line[j];
					if(j != 0) gnuplot << ", ";
					gnuplot << "'" << line.file << "' using " << line.xcol << ":" << line.ycol;
					gnuplot << " with lines ls " << line.style << " ";
					if(line.name == "EF cut-off") gnuplot << "notitle"; 
					else gnuplot << "title '" << label(line.name) << "'";
				}
				gnuplot << endl;
				
				gnuplot << "unset logscale y" << endl;
				gnuplot << "unset logscale x" << endl;
				break;
				
			case OP_TRACE:
				if(multiplot == UNSET){ multiplot = 0; gnuplot << "set multiplot layout 2,1 rowsfirst" << endl;}
				multiplot++;
				
			  gnuplot << "set size ratio 0.4" << endl;
				gnuplot << "plot ";
				for(auto j = 0u; j < oppl.line.size(); j++){
					const auto &line = oppl.line[j];
					if(j != 0) gnuplot << ", ";
					if(line.name == "True"){
						gnuplot << line.file;
					}
					else{
						gnuplot << "'" << line.file << "' using " << line.xcol << ":" << line.ycol;
					}
					gnuplot << " with lines ls " << line.style << " ";
					gnuplot << "title '" << label(line.name) << "'";
				}
				gnuplot << endl;
			  gnuplot << "set size noratio" << endl;
				break;
				
			case OP_ME:
				gnuplot << "set logscale x" << endl;
					
				gnuplot << "plot ";
				for(auto j = 0u; j < oppl.line.size(); j++){
					const auto &line = oppl.line[j];
					if(j != 0) gnuplot << ", ";
					if(line.name == ""){  // Ordinary line
						gnuplot << "'" << line.file << "' using " << line.xcol << ":" << line.ycol;
						gnuplot << " with lines ls " << line.style << " ";
					
					}
					else{                 // Error bar
						gnuplot << "'" << line.file << "' using " << line.xcol << ":" << line.ycol << ":" << line.name;
						gnuplot << " with yerrorbars ls " << line.style << " ";
					}
					gnuplot << "notitle";
				}
				gnuplot << endl;
				
				gnuplot << "unset logscale x" << endl;
				break;
		}
		gnuplot << endl << endl;
	}
}


/// Outputs a file containing posterior estimates for model parameters
void Output::posterior_parameter_estimates(const vector <ParamSample> &psamp) const
{
	auto file = "Parameter_estimates.csv";
	auto filefull = details.output_directory+"/"+file;

	ofstream paramout(filefull.c_str());
	if(!paramout) emsg("Cannot open the file '"+filefull+"'");
	
	if(suppress_output == false) cout << "'" << file << "' gives the model parameters." << endl;
	
	if(details.siminf == SIMULATE) paramout << "# Values for model parameters." << endl;
	else paramout << "# Posterior distributions for model parameters." << endl;
	paramout << endl;
	
	if(details.siminf == SIMULATE) paramout << "Name,Value" << endl;
	else{
		paramout << "Name,Mean,95% CI min,95% CI max,Standard deviation,Prior,";
		paramout << "Effective sample size,Gelman-Rubin statistic" << endl;
	}

	auto GR = get_Gelman_Rubin_statistic(psamp);
	
	auto ESS = get_effective_sample_size(psamp);
	
	vector <double> paramav(model.param.size());
	for(auto p = 0u; p < model.param.size(); p++){
		vector <double> vec; for(const auto psa : psamp) vec.push_back(psa.paramval[p]);
		auto stat = get_statistic(vec);
		paramav[p] = stat.mean;
		if(model.param[p].name != "zero" && model.param[p].name != "one"){
			paramout << model.param[p].name  << "," << stat.mean;
			if(details.siminf == INFERENCE){ 
				paramout << "," << stat.CImin << ","<< stat.CImax << "," << stat.sd << ",";
				paramout << replace(model.print_prior(p),",","|") << ",";
				if(ESS[p] == UNSET) paramout << "---"; else paramout << ESS[p]; paramout << ",";
				if(GR[p] == UNSET) paramout << "---"; else paramout << GR[p];
			}
			paramout << endl; 
		}
	}
	
	if(details.siminf == SIMULATE) cout << "Parameter values are given in '" << filefull << "'" << endl << endl;
	else cout << post_lab+" parameter estimates are given in '" << filefull << "'" << endl << endl;
}


/// Outputs a spatial distribution for R as a function of time 
void	Output::spatial_R_map(const vector <ParamSample> &psamp) const
{
	for(auto st = 0u; st < data.nstrain; st++){
		vector < vector <double> > Rmap;
		
		Rmap.resize(details.period);
		for(auto t = 0u; t < details.period; t++){
			Rmap[t].resize(data.narea); for(auto c = 0u; c < data.narea; c++) Rmap[t][c] = 0;
		}			
		
		for(const auto &info : model.Rspline_info){
			auto sp = info.spline_ref;
			
			for(const auto &psa : psamp){
				const auto &paramv_dir = psa.paramval;
				
				auto Rt = model.create_disc_spline(sp,paramv_dir);
		
				auto Rfac = paramv_dir[data.strain[st].Rfactor_param];
				
				auto areafactor = model.create_areafactor(paramv_dir); 
			
				for(auto t = 0u; t < details.period; t++){
					auto sett = t*details.division_per_time;
					for(auto c = 0u; c < data.narea; c++){
						Rmap[t][c] += Rfac*Rt[sett]*areafactor[t][c]/psamp.size();
					}
				}
			}
		}
		
		string file = "R_map"; if(data.nstrain > 1) file += "_"+data.strain[st].name; file += ".csv";			
		auto filefull = details.output_directory+"/"+post_dir+"/"+file;
		ofstream Rmapout(filefull.c_str());
		if(!Rmapout) emsg("Cannot open the file '"+filefull+"'");
	
		Rmapout << details.time_format_str; for(auto c = 0u; c < data.narea; c++) Rmapout << "," << data.area[c].code; Rmapout << endl;
		for(auto t = 0u; t < details.period; t++){
			Rmapout << details.getdate(t); for(auto c = 0u; c < data.narea; c++) Rmapout << "," << Rmap[t][c]; Rmapout << endl;
		}
	}
}


/// Removes underscore
string Output::rus(string st) const
{
	for(auto i = 0u; i < st.length(); i++){
		if(st.substr(i,1) == "_") st = st.substr(0,i)+" "+st.substr(i+1,st.length()-(i+1));
	}
	return st;
}


/// Simulates a file that can be used to test algorithm
void Output::simulate_level_file() const
{
	auto file = details.output_directory+"/level.csv";

	string levels = "level1,level2";
	
	auto levs = split(levels,',');  
	
	ofstream lout(file.c_str());
	lout << "Date";
	auto narea = data.narea;
	for(auto c = 0u; c < narea; c++) lout << "," << data.area[c].code; 
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


/// Generates simulated data with files provided in the input TOML file
void Output::simulated_data(const vector <double> &obs_value, const string dir) const
{	
	cout << endl << "Simulated data in directory '" << dir <<"':" << endl;
	
	for(auto &dt : data.datatable){
		if(dt.optype == DATA){
			auto file = dt.file;
			auto filefull =  dir+"/"+file;
			cout << "  Generating data file '" << file << "'" << endl;
			
			string sep;
			if(stringhasending(file,".txt")) sep = "\t";
			else{
				if(stringhasending(file,".csv")) sep = ",";
				else emsg("File '"+file+"' is not it '.txt' or '.csv' format");
			}
			
			ofstream dataout(filefull);
			if(!dataout) emsg("Cannot open the file '"+filefull+"'");

			if(dt.type != POPFRAC && details.stochastic == true) dataout << std::fixed << std::setprecision(0);
	
			if(!dataout) emsg("Cannot open the file '"+filefull+"'");
			
			switch(dt.type){
				case POP: case POPFRAC: case TRANS:
					{
						dataout << details.time_format_str;
						for(auto i : dt.graph_ref) dataout << sep << data.graph[i].colname;
						dataout << endl;

						if(dt.graph_ref.size() == 0) emsg("The file '"+file+"' contains no graphs");
						
						auto nrow = data.graph[dt.graph_ref[0]].point.size();
						for(auto row =0u; row < nrow; row++){
							dataout << details.getdate(dt.start + row*dt.timestep);
							for(auto i : dt.graph_ref){
								dataout << sep;
								
								auto val = obs_value[data.graph[i].point[row].obs];
								if(dt.threshold != UNSET && val <= dt.threshold){
									dataout << data.threshold_str;
								}
								else{
									if(val < 0) val = 0;
									dataout	<< val;
								}
							}
							dataout << endl;
						}
					}
					break;
					
				case MARGINAL:
					{
						dataout << "Demographic state" << sep << "Value" << endl;
						auto nrow = dt.demolist.size();
						for(auto row =0u; row < nrow; row++){
							auto i = dt.graph_ref[0];
							dataout << dt.demolist[row] << sep << obs_value[data.graph[i].point[row].obs] << endl;
						}
					}
					break;
			}
		}
	}
}


bool WeightedPoint_ord (WeightedPoint p1,WeightedPoint p2) { return (p1.val < p2.val); }

/// Gets mean and credible interval for a series of weighted samples
Statistics Output::get_statistic_with_weight(vector <WeightedPoint> vec) const  
{
	Statistics stat;
	
	if(vec.size() == 0) emsgEC("Output",4);
	
	auto sum = 0.0, sum2 = 0.0, sumw = 0.0; 
	for(auto v : vec){ sum += v.val*v.w; sum2 += v.val*v.val*v.w; sumw += v.w;}
	
	stat.mean = sum/sumw; 
	auto var = sum2/sumw - (sum/sumw)*(sum/sumw); if(var < 0) var = 0;
	stat.sd = sqrt(var);
	
	sort(vec.begin(),vec.end(),WeightedPoint_ord);
	
	auto n = vec.size();
	if(n >= 2){
		double sumtail = 0;
		auto i = 0; while(i < int(n) && sumtail < 0.025*sumw){ sumtail += vec[i].w; i++;}
		i--;
		double valw = vec[i].w;
		sumtail -= valw;
		double f = (0.025*sumw-sumtail)/valw;

		stat.CImin = vec[i].val*(1-f) + vec[i+1].val*f;
		
		sumtail = 0;
		i = n-1; while(i >= 0 && sumtail < 0.025*sumw){ sumtail += vec[i].w; i--;}
		i++;
		valw = vec[i].w;
		sumtail -= valw;
		f = (0.025*sumw-sumtail)/valw;

		stat.CImax = vec[i].val*(1-f) + vec[i-1].val*f;
	}
	else{
		stat.CImin = vec[0].val;
		stat.CImax = vec[0].val;
	}
	
	return stat;
}
		
		
/// Calculates diagnostic statistics
Statistics Output::get_statistic(const vector <double> &vec) const                       
{
	Statistics stat;
	
	auto n = vec.size();

	auto sum = 0.0, sum2 = 0.0; 	
	for(auto v : vec){ sum += v; sum2 += v*v;}
	
	stat.mean = sum/n; 	
	auto var = sum2/n - (sum/n)*(sum/n); if(var < 0) var = 0;
	stat.sd = sqrt(var);
	
	if(n == 0){
		stat.CImin = UNSET; stat.CImax = UNSET; 
	}
	else{ 
		auto vec2 = vec;
		sort(vec2.begin(),vec2.end());

		if(n >= 2){
			auto i = (unsigned int)((n-1)*0.025); auto f = (n-1)*0.025 - i;
			stat.CImin = vec2[i]*(1-f) + vec2[i+1]*f;
				
			i = (unsigned int)((n-1)*0.975); f = (n-1)*0.975 - i;
			stat.CImax = vec2[i]*(1-f) + vec2[i+1]*f;
		}
		else{
			stat.CImin = vec2[0];
			stat.CImax = vec2[0];
		}
	}
	
	return stat;
}


/// Calculate the effective sample size for a 
vector <unsigned int> Output::get_effective_sample_size(const vector <ParamSample> &psamp) const
{
	auto nparam = model.param.size();
	vector <unsigned int> ESS(nparam);
	for(auto th = 0u; th < nparam; th++) ESS[th] = UNSET;
	
	if(details.mode == PMCMC_INF || details.mode == MC3_INF || details.mode == MCMC_MBP){
		vector <ParamSample> psamp_order;
	
		for(auto ru = 0u; ru < nrun; ru++){
			for(const auto &ps : psamp){ if(ps.run == ru) psamp_order.push_back(ps);}
		}
	
		auto nsamp = psamp_order.size();
		vector <double> store(nsamp);
		for(auto th = 0u; th < nparam; th++){
			if(model.param[th].priortype != FIXED_PRIOR){
				auto av = 0.0, av2 = 0.0;
				for(auto s = 0u; s < nsamp; s++){
					auto val = psamp_order[s].paramval[th];
					av += val; av2 += val*val;
				}
				auto mean = av/nsamp, sd = sqrt(av2/nsamp - (av/nsamp)*(av/nsamp));
				
				for(auto s = 0u; s < nsamp; s++) store[s] = (psamp_order[s].paramval[th]-mean)/sd;
				
				auto sum = 1.0;
				for(auto d = 0u; d < nsamp/2; d++){
					auto a = 0.0; for(auto s = 0u; s < nsamp-d; s++) a += store[s]*store[s+d]; 
					auto cor = a/(nsamp-d); if(cor < 0) break;
					sum += 0.5*cor;			
				}
				
				ESS[th] = (unsigned int)(nsamp/sum);
			}	
		}
	}
	
	return ESS;
}
		
			
/// Calculate the Gelman Rubin statistic for each of the model parameters based on weighted samples
vector <double>	Output::get_Gelman_Rubin_statistic(const vector <ParamSample> &psamp, const vector <double> &w, const unsigned int nrun) const 
{
	vector <ParamSample> psamp_new;

	auto NN = psamp.size();
	double wsum[NN];
	for(auto ru = 0u; ru < nrun; ru++){ // This performs a bootstrap step that account for different particle weights 
		auto sum = 0.0;
		for(auto i = 0u; i < NN; i++){ 
			if(psamp[i].run == ru) sum += w[i];
			wsum[i] = sum;
		}
		
		for(auto j = 0u; j < NN/nrun; j++){
			auto z = ran()*sum;
			auto p = 0u; while(p < NN && z > wsum[p]) p++; if(p == NN) emsgEC("Output",5);
			psamp_new.push_back(psamp[p]);
		}
	}

	return get_Gelman_Rubin_statistic(psamp_new);
}


/// Calculate the Gelman Rubin statistic for each of the model parameters
vector <double>	Output::get_Gelman_Rubin_statistic(const vector <ParamSample> &psamp) const 
{
	vector <double> GR;

	vector < vector < vector <double> > > param_GR;
	param_GR.resize(nrun);
	for(const auto &ps : psamp){
		param_GR[ps.run].push_back(ps.paramval);
	}

	if(false){
		cout << nrun << " nrun" << endl;
		for(auto ru = 0u; ru < nrun; ru++){
			cout << ru << " " << param_GR[ru].size() << "  size" << endl;	
		}
	}
	
	return get_Gelman_Rubin_statistic(param_GR);
}


/// Calculate the Gelman Rubin statistic for each of the model parameters
vector <double> Output::get_Gelman_Rubin_statistic(const vector < vector < vector <double> > > &param_GR) const
{
	auto nparam = model.param.size();
	vector <double> GR(nparam);
	for(auto th = 0u; th < nparam; th++) GR[th] = UNSET;
		
	auto nrun = param_GR.size(); if(nrun < 1) emsgEC("Output",6);
	
	if(nrun > 1){
		auto N = LARGE;
		for(auto ru = 0u; ru < nrun; ru++){
			if(param_GR[ru].size() < N) N = param_GR[ru].size();
		}		
		if(N == 0) emsgEC("Output",7);
		
		for(auto th = 0u; th < nparam; th++){
			if(model.param[th].priortype != FIXED_PRIOR){
				double mu[nrun], vari[nrun];
				
				auto muav = 0.0;
				for(auto ru = 0u; ru < nrun; ru++){ 
					auto valav = 0.0; for(auto i = 0u; i < N; i++) valav += param_GR[ru][i][th]/N;
					auto varr = 0.0; for(auto i = 0u; i < N; i++) varr += (param_GR[ru][i][th]-valav)*(param_GR[ru][i][th]-valav)/(N-1);
					mu[ru] = valav;
					vari[ru] = varr;
					muav += mu[ru]/nrun;
				}
				auto W = 0.0; for(auto ru = 0u; ru < nrun; ru++) W += vari[ru]/nrun;
				auto B = 0.0; for(auto ru = 0u; ru < nrun; ru++) B += (mu[ru]-muav)*(mu[ru]-muav)*N/(nrun-1);
				GR[th] = sqrt(((1-1.0/N)*W + B/N)/W);
			}
		}
	}
	
	return GR;
}
		
		
/// Gets the probability distributions for a given set of samples
Distribution Output::get_distribution(const vector <double> &vec, const double priormin, const double priormax) const
{
	Distribution dist;
	
	auto min = LARGE, max = -LARGE;	
	for(auto v : vec){
		if(v < min) min = v;
		if(v > max) max = v;
	}
	
	dist.variation = true;
	if((max-min)*(max-min) < VTINY){
		dist.variation = false;
		auto d = 0.01*min; if(d < TINY) d = TINY;
		min -= d; max += d;
	}
	
	auto d = max-min;
	min -= 0.2*d; max += 0.2*d;
	if(priormin != UNSET && min < priormin) min = priormin;
	if(priormax != UNSET && max > priormax) max = priormax;
	if(priormin != UNSET && priormax != UNSET && d > 0.6*(priormax-priormin)){ min = priormin; max = priormax;}
	
	vector <double> bin(prop.nbin);
	for(auto b = 0u; b < prop.nbin; b++) bin[b] = 0;
		
	auto bin_max = 0.0;
	switch(prop.type){
	case BINNING:
		for(auto v : vec){
			auto b = (unsigned int)(prop.nbin*(v-min)/(max-min)); if(b >= prop.nbin) emsgEC("Output",8);
			bin[b]++;
		}

		for(auto val : bin){ if(val > bin_max) bin_max = val;}
		bin_max *= 1.3;

		for(auto b = 0u; b < prop.nbin; b++){ 
			dist.value.push_back(to_string(min+(b)*(max-min)/prop.nbin)); 
			dist.prob.push_back(to_string(bin[b]/bin_max));
			
			dist.value.push_back(to_string(min+(b+1)*(max-min)/prop.nbin)); 
			dist.prob.push_back(to_string(bin[b]/bin_max));
		}
		break;
		
	case KDE:
		auto dd = (max-min)/prop.nbin;
		
		for(auto v : vec){
			auto bf = (v-min)/dd;
			auto b = (unsigned int) bf; 
			auto bmin = b-prop.h, bmax = b+prop.h+1;

			for(auto bb = bmin; bb <= bmax; bb++){			
				if(bb >= 0 && bb < prop.nbin){
					auto dist = ((bb+0.5)-bf)/prop.h;
					if(dist >= -1 && dist <= 1){
						if(dist > 0) bin[bb] += 1-dist;
						else bin[bb] += 1+dist;
					}
				}
			}
		}

		for(auto val : bin){ if(val > bin_max) bin_max = val;}
		bin_max *= 1.3;
		
		for(auto b = 0u; b < prop.nbin; b++){ 
			dist.value.push_back(to_string(min+(b+0.5)*(max-min)/prop.nbin)); 
			dist.prob.push_back(to_string(bin[b]/bin_max));
		}
		break;
	}
	
	return dist;
}


/// Outputs how the error function and posterior distributions change with generation number (ABC-SMC and ABC-MBP)
void Output::generation_plot(const string file, const vector <Generation> &generation) const
{
	auto nparam = model.param.size(); 

	string filefull = details.output_directory+"/"+file;
	ofstream genout(filefull.c_str());
	if(!genout) emsg("Cannot open the file '"+filefull+"'");

	genout << "# This shows how the EF and infered posterior distribiutions vary with generation number" << endl;
	genout << "Generation,CPU Time,EF cut-off,InvT,Model evidence,Model evidence SD"; 
	for(const auto &param : model.param) genout << "," << param.name << " (mean)," << param.name << " (95% CI min)," << param.name << " (95% CI max)";
	genout << endl;
	
	long timestart = generation[0].time;
	for(auto g = 1u; g < generation.size(); g++){
		const Generation &gen=generation[g];
		
		genout << g << "," << mpi.ncore*(gen.time - timestart)/(60.0*CLOCKS_PER_SEC) << ",";

		switch(details.mode){
			case PAIS_INF:
				genout << gen.EFmax<< "," << gen.invT;
				break;
			default:
				genout << gen.EFcut << ",UNSET" ;
				break;
		}
		auto stat = get_statistic(gen.model_evidence);
	
		genout << "," << stat.mean << "," << stat.sd;
	
		auto psamp = gen.param_samp;
		auto npsamp = psamp.size();

		for(auto &psa : psamp) psa.paramval =  model.dirichlet_correct(psa.paramval);
		
		for(auto th = 0u; th < nparam; th++){
			vector <WeightedPoint> vec;
			for(auto i = 0u; i < npsamp; i++){
				WeightedPoint pw;
				pw.val = psamp[i].paramval[th];
				switch(details.mode){
					case ABC_SMC: pw.w = gen.w[i]; break;
					case ABC_MBP: case PAIS_INF: case MC3_INF: case MCMC_MBP: pw.w = 1; break;
					default: emsgEC("Output",9); break;
				}
				vec.push_back(pw);
			}

			Statistics stat = get_statistic_with_weight(vec);
			
			genout << "," << stat.mean << "," << stat.CImin << "," << stat.CImax;
		}
		genout << endl;
	}	
}


/// Outputs the error function for different datatables
void Output::EF_datatable_plot(const string file, const vector <Generation> &generation) const
{
	string filefull = details.output_directory+"/"+file;
	ofstream ldtyout(filefull.c_str());
	if(!ldtyout) emsg("Cannot open the file '"+filefull+"'");

	ldtyout << "# This shows the contribution to the error function as a function of generation" << endl;
	
	ldtyout << "Generation";
	for(auto &dt : data.datatable) ldtyout << "," << dt.file;
	ldtyout << endl;
		 
	for(auto g = 1u; g < generation.size(); g++){
		const Generation &gen = generation[g];
		
		ldtyout << g;
		for(auto dt = 0u; dt < data.datatable.size(); dt++){
			auto val = 0.0; for(const auto &sa : gen.EF_datatable) val += sa[dt];
			val /= gen.EF_datatable.size();
			
			ldtyout << "," << val;
		}
		ldtyout << endl;
	}
}


/// Outputs results for the generations (ABC-SMC and ABC-MBP approaches)
void Output::generation_results(const vector <Generation> &generation) const 
{
	timer[TIME_RESULTS].start();
		
	if(mpi.core == 0){
		generation_plot("Diagnostics/Generation.csv",generation);
		EF_datatable_plot("Diagnostics/EF_datatable.csv",generation);
	}
	
	timer[TIME_RESULTS].stop();
}


/// Creates all the directories needed for outputting
void Output::ensure_directories()
{
	post_dir = "Posterior"; post_lab = "Posterior";
	if(details.mode == SIM){ post_dir = "Simulation_Results"; post_lab = "Simulation";}
	if(details.mode == MULTISIM){ post_dir = "Multisim_Results"; post_lab = "Multisim";}
		
	ensure_directory(details.output_directory);
	if(details.siminf != SIMULATE) ensure_directory(details.output_directory+"/Diagnostics");
	if(details.mode == MC3_INF) ensure_directory(details.output_directory+"/Diagnostics/Other Chains");
	auto dir = details.output_directory+"/"+post_dir;
	ensure_directory(dir);
	ensure_directory(dir+"/parameter");
	ensure_directory(dir+"/state");
	ensure_directory(dir+"/susceptibility");
	ensure_directory(dir+"/spline");

	if(data.democat_change.size() > 0) ensure_directory(dir+"/democat_change");
	
	if(details.siminf == INFERENCE) ensure_directory(dir+"/samples");
	
	if(details.siminf == SIMULATE) ensure_directory(details.output_directory+"/Simulated_data");
}

		
/// Create a directory if it doesn't already exist
void Output::ensure_directory(const string &path) const
{
	struct stat st;
	if (stat(path.c_str(), &st) == -1){  	// Directory not found
		int ret = mkdir(path.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
		if(ret == -1) emsg("Error creating the directory '"+path+"'");
	}
}


/// Summarises the data and model
void Output::print_model_specification() const 
{
	auto file = details.output_directory+"/Model_specification.txt";
	ofstream modspec(file); 
	if(!modspec) emsg("Cannot open the file '"+file+"'");
	modspec << data.print();                             
	modspec << model.print();
}
		

/// Outputs the model evidence
void Output::final_model_evidence(const vector <double> &ME_list, const double invT_final, const double cutoff_final) const
{
	cout << endl;
	if(invT_final != UNSET) cout << "For invT = " << invT_final << " ";
	else cout << "For EF cut-off " << cutoff_final << " ";

	auto stat = get_statistic(ME_list);
	cout << "the log of the model evidence is " << stat.mean;
	if(nrun > 1) cout << "   (SD = " << stat.sd << ")";
	cout << endl;  
	cout << "When comparing models a higher value implies a better fit with the data." << endl; 
}


/// Initialises trace plot for parameters
void Output::trace_plot_inititialise(string name, ofstream &trace) const 
{
	auto file = details.output_directory+"/Diagnostics/"+name+".csv";

	trace.open(file);		
	if(!trace) emsg("Cannot open the file '"+file+"'");
	trace << "state";
	for(const auto &par : model.param) trace << "," << par.name; 
	trace << ",Li";
	trace << ",Pri"; 
	trace << endl;
}


/// Outputs trace plot
void  Output::trace_plot(const unsigned int samp, const double Li, const vector <double> &paramval, ofstream &trace) const 
{
	trace << samp; 
	for(const auto &pval : paramval) trace << "," << pval; 

	trace << "," << Li; 
	trace << "," << model.prior(paramval);
	trace << endl;
}


/// Generates a readme file for the output directory
void Output::readme() const
{
	auto file = details.output_directory+"/README.md";
	ofstream read(file);
	if(!read) emsg("Cannot open the file '"+file+"'");
	
	read << "# Outputs" << endl << endl;
	read << "This directory contains the outputs from the BEEPmbp analysis."<< endl;
	read << "We give a brief description below of the various files and folders contained within it:" << endl << endl;
			
	read << "## Parameter_estimates.csv" << endl << endl;
	
	read << "If inference is performed, this file gives the posterior means and 95% credible intervals for each of the model parameters." << endl;
	
	read << endl;
	
	read << "## Graphs.pdf" << endl << endl;
	
	read << "This visualises outputs from BEEPmbp." << endl << endl;
	
	read << "## Graphs_description.txt" << endl << endl;
	
	read << "Provides information about the source data for the graphs." << endl << endl;
			
	read << "## Model_specification.txt" << endl << endl;
	
	read << "This gives a summary of the model used to perform the analysis." << endl;
	
	read << endl;
	
	read << "## " << post_lab << endl << endl;
	
	read << "This folder contains files giving posterior probability distributions for the model parameters." << endl << endl;
	read << "These are arranged in a number of different ways:" << endl << endl;
	
	read << "**parameter** - Contains files for the model parameters." << endl << endl;
	read << "**state** - Contains files which compare the system state with that observed in the actual data files." << endl << endl;
	read << "**spline** - Contains files giving time variation in splines used within the model." << endl << endl;
	read << "**susceptibility** - Contains files giving the variation in susceptibility for different demographic classes within the model." << endl << endl;
	read << "**sample** - Contains files giving raw posterior samples." << endl << endl;
	read << "**Rmap.csv** - For spatial models this gives the variation in R across different regions." << endl << endl;
	read << endl;
	
	read << "## Simulated_data" << endl << endl;
	
	read << "This gives simulated data files corresponding to the specifications provided in the input TOML file." << endl << endl;
	
	read << "## Diagnostics" << endl << endl;
	
	read << "This provides diagnostic information to inform how well the algorithm is performing:" << endl << endl;
	
	read << "**Generation.csv** - Shows how model parameters move from the prior distribution to an approximation of the posterior distribution as a function of the generation number ('generations' are used within the ABC-MBP and ABC-SMC procedures to filter particles such that they provide a better and better approximation to the posterior)." << endl << endl;
	
	read << "**MCMC_proposals.txt** - This shows the performance of any MCMC proposals used." << endl;
}


/// This takes a vector of particle sample from each mpi core and combines them
void Output::generate_graphs(vector <Particle> &particle_store) const 
{
	timer[TIME_RESULTS].start();

	if(mpi.core == 0) cout << endl << "Generating outputs..." << endl;
			
	State state(details,data,model,obsmodel);
	
	auto dir = details.output_directory+"/"+post_dir+"/samples/";
	vector <ParamSample> psamp;
	vector <Sample> opsamp;
	mpi.gather_samples(psamp,opsamp,particle_store,state,dir);

	if(mpi.core == 0) generate_graphs(psamp,opsamp);

	timer[TIME_RESULTS].stop();
}


/// Prints percentage complete
void Output::print_percentage(const double s, const unsigned int nsamp, unsigned int &percentage) const
{
	if(mpi.core == 0){
		unsigned int step = 5;
		
		unsigned int per = step*int((100/step)*s/nsamp);
		if(per != percentage || percentage == UNSET){ cout << " " << per << "%" << endl; percentage = per;}
	}
}


/// Generates a new output plot
OutputPlot::OutputPlot(OutputPlotType type_, string title_, string xaxis_, string yaxis_, double min_, double max_)
{
	type = type_; title = title_; xaxis = xaxis_; yaxis = yaxis_; min = min_; max = max_; legend_labelson = false;
}


/// Generates a reference label
string Output::reference_label(const unsigned int i) const
{
	string alphabet = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz";
	if(i < 52) return alphabet.substr(i,1);
	else return alphabet.substr(i/52,1)+alphabet.substr(i%52,1);
}
	
	
/// Makes the underscore display correctly in gnuplot
string Output::label(string lab) const
{
	auto i=0u;
	while(i < lab.length()){
		if(lab.substr(i,1) == "_"){
			lab = lab.substr(0,i+1)+"{"+lab.substr(i+1);
			while(i < lab.length() && lab.substr(i,1) != " "&& lab.substr(i,1) != "-") i++;
			lab = lab.substr(0,i)+"}"+lab.substr(i);
		}
		i++;
	}
	lab = replace(lab,"->"," -> ");
	return lab;
}


/// Adds a line to an output plot
void OutputPlot::addline(const string name, const string file, const unsigned int xcol, const unsigned int ycol, const LineType style)
{
	OutputLine opl; opl.name = name; opl.file = file; opl.xcol = xcol; opl.ycol = ycol; opl.style = style;
	line.push_back(opl);
};
