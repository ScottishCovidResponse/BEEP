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
		inputs.find_plot_param_values(plot_param_values);
		inputs.find_nrun(nrun);
		inputs.find_stateuncer(stateuncer);
		
		inputs.find_area_plot(area_plot);                                 // Finds information about plotting areas 
		
		boundaries = load_boundaries();
		
		readme(); 	                                                    	// Generates a readme file for outputs
	}
	mpi.barrier();
}
	

/// Generates the output files along with the final pdf report that gives graphical outputs
void Output::generate_graphs(vector <ParamSample> &psamp, const vector <Sample> &opsamp, const double invT) const
{ 
	if(mpi.core == 0){
		cout << endl << "Generating outputs in directory '" << details.output_directory;
		cout << "'..." << endl;
	}	
	
	for(auto &psa : psamp) psa.paramval = model.dirichlet_correct(psa.paramval);
	
	if(suppress_output == false && details.siminf != DATAVIEW){
		if(details.siminf == SIMULATE) cout << "Outputs in directory '" << details.output_directory << "':" << endl;
		else cout << details.analysis_type << " outputs in directory '" << details.output_directory << "':" << endl;
	}

	vector <OutputPlot> op;                                             // Stores the plots for the final pdf

	compartmental_model(op);                                            // Generates all the different types of plot
		
	if(details.mode == DATAONLY){
		posterior_parameter_estimates(psamp,op);
		
		spatial_mixing_map(op);
	
		age_mixing_matrix(opsamp,op);
		
		datatable_maps(opsamp,op);
		
		graph_plots(opsamp,op,invT); 
		
		add_democat_change(op);
	}
	else{
		posterior_parameter_estimates(psamp,op);
		
		posterior_parameter_distributions(psamp,op);
			
		spatial_R_map(opsamp,op);  	
		
		spline_plots(opsamp,op);

		spatial_mixing_map(op);
		
		age_mixing_matrix(opsamp,op);
			
		datatable_maps(opsamp,op);
			
		graph_plots(opsamp,op,invT); 
		
		add_democat_change(op);
		
		susceptibility_distributions(psamp,op);

		mean_distributions(psamp,op);

		branch_prob_distributions(psamp,op);

		age_mixing_perturb_distributions(psamp,op);
				
		level_effect_distributions(psamp,op);
		
		area_effect_distributions(psamp,op);
		
		derived_parameter_distributions(opsamp,op);

		covar_data(op);
		
		level_data(op);
			
		set_labelon(op);                                                  // Determines if labels put in description
		
		set_graph_source(op);                                             // Sets source files
		
		if(false){
			auto grfile = "Graphs";
			generate_gnuplot_file(op,grfile);                                 // Generates the gnuplot file
		
			generate_pdf(grfile,"Output graphs");                             // Runs gnuplot and generates pdf 
		
			generate_pdf_description(op,grfile);                              // Generates a text descrition of the pdf file
		}
	}
	
	if(details.siminf == INFERENCE){
		vector <OutputPlot> op_diag;                                    // Stores the plots for the diagnostic pdf

		add_generation_plots(op_diag);
		
		add_trace_plots(op_diag);
		
		set_labelon(op_diag);                                           // Determines if labels put in description
		
		set_graph_source(op);                                           // Sets source files
		
		if(false && op_diag.size() > 0){
			auto grfile = "Diagnostic_Graphs";
			generate_gnuplot_file(op_diag,grfile);                        // Generates the gnuplot file
		
			generate_pdf(grfile,"Output diagnostic graphs");              // Runs gnuplot and generates pdf 
		
			generate_pdf_description(op_diag,grfile);                     // Generates a text descrition of the pdf file
		}
	
		for(auto opu : op_diag) op.push_back(opu);
	}
	
	generate_visualisation(op,"vis_files/BEEP_output.js");         // Generates a file used by the visualisation software
	
	if(1 == 0) EF_dist(psamp);                                        // This is used for the figure in the paper
}


/// Generates the distribution in error function
void Output::EF_dist(const vector <ParamSample> &psamp) const
{
	const auto nbin = 100u;
	
	double EFmin = LARGE, EFmax = -LARGE;
	for(const auto &ps : psamp){
		if(ps.EF > EFmax) EFmax = ps.EF;
		if(ps.EF < EFmin) EFmin = ps.EF;
	}		
	EFmin = 0;
	
	vector <double> bin(nbin);
	for(auto b = 0u; b < nbin; b++) bin[b] = 0;
	
	for(const auto &ps : psamp){
		auto b = int(0.999999*nbin*(ps.EF-EFmin)/(EFmax-EFmin));
		bin[b]++;
	}

	cout << "Outputing distribution in error function" << endl;

	auto file = details.output_directory+"/EF_dist.txt";
	ofstream dist(file);
	for(auto b = 0u; b < nbin; b++){
		dist << EFmin+(EFmax-EFmin)*double(b)/nbin << "," << bin[b] << endl;
	}
}


/// Generate the final pdf report
void Output::generate_pdf(const string file, const string desc) const
{
	auto opdir = details.output_directory;
	stringstream ss; ss << "gnuplot '" << opdir << "/gnuplot.txt'" << endl;
	system(ss.str().c_str());
		
	stringstream ss2;
	ss2 << "ps2pdf '" << opdir << "/" << file << ".ps' '" << opdir << "/" << file << ".pdf'" << endl;
	system(ss2.str().c_str());

	stringstream ss3; ss3 << "rm -f '" << opdir << "/" << file << ".ps'" << endl; 
	system(ss3.str().c_str());
		
	stringstream ss4; ss4 << "rm -f '" << opdir << "/gnuplot.txt'" << endl;
	system(ss4.str().c_str());

	stringstream ss5; ss5 << "rm -f '" << opdir << "/line.csv'" << endl;
	system(ss5.str().c_str());

	cout << desc << " are placed into '" << opdir << "/" << file << ".pdf'" << endl;
}


/// Sets the source files for the plots
void Output::set_graph_source(vector <OutputPlot> &op) const
{
	for(auto i = 0u; i < op.size(); i++){   // Sets the source files
		const auto &oppl = op[i];
	
		stringstream ss;
		for(auto j = 0u; j < oppl.line.size(); j++){
			const auto &line = oppl.line[j];
			if(line.name == "") ss << "This ";
			else ss << "'" << line.name << "' ";
			
			switch(oppl.type){
			case OP_SPLINE: case OP_GRAPH: case OP_PARAMDIST: case OP_GENERATION: case OP_LOG_GENERATION:
			case OP_CPU: case OP_TRACE: case OP_ME:
				ss << "uses columns " << line.xcol << " and " << line.ycol << " in '";
				break;
			case OP_GRAPH_MARGINAL:	
				ss << "is generated from ";
				break;
			case OP_MARGINAL:	
				ss << "uses columns 1, 2 and " << line.ycol << " in '";
				break;
			case OP_ANIM_MAP: case OP_AGE_MATRIX: case OP_AGE_MATRIX_POST: case OP_AGE_MATRIX_DIF: 
			case OP_PARAM_TABLE: case OP_PRIOR_TABLE: 
			case OP_COMP_MODEL: case OP_FOI_MODEL: case OP_DESC:	
			case OP_MIXING_WITHIN_ANIM_MAP: case OP_MIXING_WITHIN_MAP:
			case OP_MIXING_BETWEEN_ANIM_MAP: case OP_MIXING_BETWEEN_MAP: 
			case OP_MIXING_POP_ANIM_MAP: case OP_MIXING_POP_MAP:
			case OP_AREA_COVAR: case OP_TV_COVAR: case OP_AREA_TV_COVAR: 
			case OP_LEVEL_EFFECT:
				break;
			}
			auto file = line.file;
			data.chop_dir(file,details.output_directory);
			ss << file << "'" << endl;
		}
	
		switch(oppl.type){
			case OP_ANIM_MAP: case OP_AGE_MATRIX: case OP_AGE_MATRIX_POST: case OP_AGE_MATRIX_DIF: 
			case OP_PARAM_TABLE: case OP_PRIOR_TABLE: 
			case OP_COMP_MODEL: case OP_FOI_MODEL: case OP_DESC:
			case OP_MIXING_WITHIN_ANIM_MAP: case OP_MIXING_WITHIN_MAP:
			case OP_MIXING_BETWEEN_ANIM_MAP: case OP_MIXING_BETWEEN_MAP: 
			case OP_MIXING_POP_ANIM_MAP: case OP_MIXING_POP_MAP:
			case OP_AREA_COVAR: case OP_TV_COVAR: case OP_AREA_TV_COVAR: 
			case OP_LEVEL_EFFECT:
				ss << oppl.title << endl; 
				break;
			default: break;
		}

		op[i].source = ss.str();
	}
}
	
	
/// Generates a text description of the pdf file
void Output::generate_pdf_description(vector <OutputPlot> &op, const string grfile) const
{
	auto file = details.output_directory+"/"+grfile+"_description.txt";
	ofstream desc(file);
	if(!desc) emsg("Cannot open the file '"+file+"'");
	
	desc << "This file provides a description of the source files used to generate the graphs in '";
	desc << details.output_directory << "/" << file << ".pdf'" << endl << endl;
	
	auto page = 1u;
	for(auto i = 0u; i < op.size(); i++){
		const auto &oppl = op[i];
		
		if(oppl.type != OP_ANIM_MAP && oppl.type != OP_AGE_MATRIX && oppl.type != OP_AGE_MATRIX_POST
			&& oppl.type != OP_AGE_MATRIX_DIF
			&& oppl.type != OP_COMP_MODEL && oppl.type != OP_FOI_MODEL && oppl.type != OP_DESC
			&& oppl.type != OP_MIXING_WITHIN_ANIM_MAP && oppl.type != OP_MIXING_WITHIN_MAP
			&& oppl.type != OP_MIXING_BETWEEN_ANIM_MAP && oppl.type != OP_MIXING_BETWEEN_MAP
			&& oppl.type != OP_MIXING_POP_ANIM_MAP && oppl.type != OP_MIXING_POP_MAP
			&& oppl.type != OP_AREA_COVAR && oppl.type != OP_TV_COVAR && oppl.type != OP_AREA_TV_COVAR
			&& oppl.type != OP_LEVEL_EFFECT
			){

			desc << "PAGE " << page << ": \"" << oppl.fulldesc << "\"" << endl;
		
			desc << op[i].source;
			
			if(oppl.legend_labelson){
				desc << endl << "  ";
				for(auto i =0u; i < oppl.label.size(); i++){
					if(i != 0) desc << ", ";
					desc << reference_label(i) << ": " << oppl.label[i];
				}
				desc << endl;
			}
			desc << endl;
			page++;
		}
	}
	
	cout << "The source files for these graphs are referenced in '" << file << "'" << endl << endl;
}


/// Generates a text description of the pdf file
void Output::generate_visualisation(const vector <OutputPlot> &op, const string grfile) const
{
	auto opdir = details.output_directory;
	
	ensure_directory(opdir+"/vis_files");              // Copies the output 
	stringstream ss; ss << "cp ./vis_files/* '" << opdir << "/vis_files/'" << endl;
	system(ss.str().c_str());
	stringstream ss2; ss2 << "cp ./visBEEP.html '" << opdir << "/'" << endl;
	system(ss2.str().c_str());
	
	stringstream vis;
	vis << fixed;

	vis << "{";

	vis << "\"time_labels\" : [";
	for(auto i = 0u; i < details.timeplot.size(); i++){
		auto &tp = details.timeplot[i];
		if(i != 0) vis << ",";
		vis << "{\"name\" : \"" << tp.name << "\",\"time\" : "  << tp.time-details.start << "}";
	}
	vis << "]";
	
	auto pred_timeplot = get_pred_timeplot();
		
	vis << ",\"pred_timeplot\" : [";
	for(auto i = 0u; i < pred_timeplot.size(); i++){
		auto &tp = pred_timeplot[i];
		if(i != 0) vis << ",";
		vis << "{\"name\" : \"" << tp.name << "\",\"time\" : " << tp.time << "}";
	}
	vis << "]";

	vis << ",\"plots\" : [";
	for(auto i = 0u; i < op.size(); i++){
		auto opi = op[i];
		
		vis << "{";	
		string ty;
		switch(opi.type){
			case OP_GRAPH: ty = "OP_GRAPH"; break;
			case OP_GRAPH_MARGINAL: ty = "OP_GRAPH_MARGINAL"; break;
			case OP_MARGINAL: ty = "OP_MARGINAL"; break;
			case OP_PARAMDIST: ty = "OP_PARAMDIST"; break;
			case OP_SPLINE: ty = "OP_SPLINE"; break;
			case OP_GENERATION: ty = "OP_GENERATION"; break;
			case OP_LOG_GENERATION: ty = "OP_LOG_GENERATION"; break;
			case OP_CPU: ty = "OP_CPU"; break;
			case OP_TRACE: ty = "OP_TRACE"; break;
			case OP_ME: ty = "OP_ME"; break;
			case OP_ANIM_MAP: ty = "OP_ANIM_MAP"; break;
			case OP_MIXING_WITHIN_ANIM_MAP: ty = "OP_MIXING_WITHIN_ANIM_MAP"; break;
			case OP_MIXING_WITHIN_MAP: ty = "OP_MIXING_WITHIN_MAP"; break;
			case OP_MIXING_BETWEEN_ANIM_MAP: ty = "OP_MIXING_BETWEEN_ANIM_MAP"; break;
			case OP_MIXING_BETWEEN_MAP: ty = "OP_MIXING_BETWEEN_MAP"; break;
			case OP_MIXING_POP_ANIM_MAP: ty = "OP_MIXING_POP_ANIM_MAP"; break;
			case OP_MIXING_POP_MAP: ty = "OP_MIXING_POP_MAP"; break;
			case OP_AGE_MATRIX: ty = "OP_AGE_MATRIX"; break;
			case OP_AGE_MATRIX_POST: ty = "OP_AGE_MATRIX_POST"; break;
			case OP_AGE_MATRIX_DIF: ty = "OP_AGE_MATRIX_DIF"; break;
			case OP_PARAM_TABLE: ty = "OP_PARAM_TABLE"; break;
			case OP_PRIOR_TABLE: ty = "OP_PRIOR_TABLE"; break;
			case OP_COMP_MODEL: ty = "OP_COMP_MODEL"; break;
			case OP_DESC: ty = "OP_DESC"; break;
			case OP_FOI_MODEL: ty = "OP_FOI_MODEL"; break;
			case OP_AREA_COVAR: ty = "OP_AREA_COVAR"; break;
			case OP_TV_COVAR: ty = "OP_TV_COVAR"; break;
			case OP_AREA_TV_COVAR: ty = "OP_AREA_TV_COVAR"; break; 
			case OP_LEVEL_EFFECT: ty = "OP_LEVEL_EFFECT"; break;
		}
		if(false) cout << ty << "type" << endl;

		vis << "\"type\":\"" << ty << "\"";
		vis << ",\"fulldesc\":\"" << opi.fulldesc << "\"";
		vis << ",\"tab\":\"" << opi.tab << "\"";
		vis << ",\"tab2\":\"" << opi.tab2 << "\"";
		vis << ",\"tab3\":\"" << opi.tab3 << "\"";
		vis << ",\"tab4\":\"" << opi.tab4 << "\"";
		
		auto title = opi.title; rem_pagebreak(title);
		vis << ",\"title\":\"" << title << "\"";
		vis << ",\"xaxis\":\"" << opi.xaxis << "\"";
		vis << ",\"yaxis\":\"" << opi.yaxis << "\"";
		vis << ",\"min\":" << opi.min;
		vis << ",\"max\":" << opi.max;
		if(opi.xmin != UNSET){
			vis << ",\"xmin\":" << opi.xmin;
			vis << ",\"xmax\":" << opi.xmax;
		}
		
		string source = opi.source; rem_pagebreak(source);
		vis << ",\"source\":\"" << source << "\"";	
		
		switch(opi.type){
			case OP_ANIM_MAP: case OP_AREA_TV_COVAR: case OP_LEVEL_EFFECT:
				{
					vis << ",\"array\":" << data.get_array_JSON(opi.title,details.output_directory);
					
					auto tab = data.get_table(opi.title,details.output_directory);
					vector <string> dates = data.get_table_column_str(0,tab);
					vis << ",\"dates\":[";
					for(auto i = 0u; i < dates.size(); i++){
						if(i > 0) vis << ",";
						vis << "\"" << dates[i] << "\"";
					}
					vis << "]";
					
					vis << ",\"areas\":[";
					for(auto c = 0u; c < data.narea; c++){
						if(c > 0) vis << ",";
						vis << "\"" << data.area[c].code << "\"";
					}
					vis << "]";		

					if(opi.yaxis == "Greyscale"){
						vis << ",\"popscale\":[";
						for(auto c = 0u; c < data.narea; c++){
							if(c > 0) vis << ",";
							vis << "\"" << data.area[c].total_pop << "\"";
						}
						vis << "]";	
					}
					
					if(opi.type == OP_LEVEL_EFFECT){
						vis << ",\"level_param\":[";
						for(auto i = 0u; i < model.level_effect_param_list.size(); i++){
							if(i > 0) vis << ",";
							vis << "\"" << model.param[model.level_effect_param_list[i]].name << "\"";
						}
						vis << "]";	
					}
				}
				break;
			
			case OP_MIXING_WITHIN_ANIM_MAP: case OP_MIXING_WITHIN_MAP:
			case OP_MIXING_BETWEEN_ANIM_MAP: case OP_MIXING_BETWEEN_MAP: 
			case OP_MIXING_POP_ANIM_MAP: case OP_MIXING_POP_MAP:
				{
					vis << ",\"areas\":[";
					for(auto c = 0u; c < data.narea; c++){
						if(c > 0) vis << ",";
						vis << "\"" << data.area[c].code << "\"";
					}
					vis << "]";
					
					vis << ",\"pop\":[";
					for(auto c = 0u; c < data.narea; c++){
						if(c > 0) vis << ",";
						vis << data.area[c].total_pop;
					}
					vis << "]";
					
					const auto &M = data.genQ.M;
					vis << ",\"N\":" << M.N;
					
					vis << ",\"diag\":[";
					for(auto c = 0u; c < M.N; c++){
						if(c > 0) vis << ",";
						vis << M.diag[c];
					}
					vis << "]";
					
					vis << ",\"to\":[";
					for(auto i = 0u; i < M.to.size(); i++){
						if(i > 0) vis << ",";
						vis << "[";
						for(auto j = 0u; j < M.to[i].size(); j++){
							if(j > 0) vis << ",";
							vis << M.to[i][j];
						}
						vis << "]";
					}
					vis << "]";
					
					vis << ",\"val\":[";
					for(auto i = 0u; i < M.val.size(); i++){
						if(i > 0) vis << ",";
						vis << "[";
						for(auto j = 0u; j < M.val[i].size(); j++){
							if(j > 0) vis << ",";
							vis << M.val[i][j];
						}
						vis << "]";
					}
					vis << "]";
				}
				
				if(opi.type == OP_MIXING_WITHIN_ANIM_MAP || opi.type == OP_MIXING_BETWEEN_ANIM_MAP || opi.type == OP_MIXING_POP_ANIM_MAP){
					vis << ",\"dates\":[";
					for(auto i = 0u; i < details.period; i++){
						if(i > 0) vis << ",";
						vis << "\"" << details.getdate(i) << "\"";
					}
					vis << "]";	
				}
				break;
				
			case OP_AREA_COVAR:
				{
					vis << ",\"areas\":[";
					for(auto c = 0u; c < data.narea; c++){
						if(c > 0) vis << ",";
						vis << "\"" << data.area[c].code << "\"";
					}
					vis << "]";
						
					auto tab = data.get_table(data.areas_file,data.data_directory);
					auto covar = data.get_table_column(opi.yaxis,tab);
				
					vis << ",\"array\":[[";
					for(auto c = 0u; c < data.narea; c++){
						if(c > 0) vis << ",";
						vis << covar[c];
					}
					vis << "]]";
				}
				break;
				
			case OP_AGE_MATRIX: case OP_AGE_MATRIX_POST: case OP_AGE_MATRIX_DIF:
				{
					vis << ",\"ages\":[";
					for(auto a = 0u; a < data.democat[0].value.size(); a++){
						if(a > 0) vis << ",";
						auto val = data.democat[0].value[a];
						if(val.length() > 3){
							if(val.substr(0,3) == "age") val = val.substr(3);
						}
						vis << "\"" << val << "\"";
					}
					vis << "]";
					
					auto N = data.genQ.N[0];
					for(auto i = 0u; i < N.ele.size(); i++){
						for(auto j = 0u; j < N.ele[i].size(); j++) N.ele[j][i] *= data.agedist[j]/N.norm_factor;
					}
					
					Matrix Nplot;
					switch(opi.type){
						case OP_AGE_MATRIX:
							Nplot = N;
							break;
							
						case OP_AGE_MATRIX_POST:
							Nplot = data.get_matrix(opi.title,details.output_directory);
							break;
							
						case OP_AGE_MATRIX_DIF:
							Nplot = data.get_matrix(opi.title,details.output_directory);
							for(auto i = 0u; i < N.ele.size(); i++){
								for(auto j = 0u; j < N.ele[i].size(); j++){
									Nplot.ele[j][i] /= N.ele[j][i];
								}
							}
							break;
							
						default: emsgEC("Output",45); break;
					}
					
					vis << ",\"ele\":[";
					for(auto i = 0u; i < Nplot.ele.size(); i++){
						if(i > 0) vis << ",";
						vis << "[";
						for(auto j = 0u; j < Nplot.ele[i].size(); j++){
							if(j > 0) vis << ",";
							vis << Nplot.ele[j][i];
						}
						vis << "]";
					}
					vis << "]";
				}
				break;
			
			case OP_PARAM_TABLE:
				vis << ",\"param_table\":" << data.get_table_JSON(opi.title,details.output_directory);
				break;
				
			case OP_PRIOR_TABLE:
				{
					vector <unsigned int> cols; 
					cols.push_back(0);
					if(details.mode == ML_INF) cols.push_back(4); else cols.push_back(6);
					vis << ",\"param_table\":" << data.get_table_cols_JSON(opi.title,details.output_directory,cols);
				}
				break;
				
			case OP_COMP_MODEL:
				vis << ",\"model\":" << model.comparmtental_model_JSON();
				break;
			
			case OP_DESC:
				break;
				
			case OP_FOI_MODEL:
				vis << ",\"foi\":" << model.foi_model_JSON();
				break;
			
			default:
				vis << ",\"line\":[";		
				Table tab;
				string file_loaded = "";
				for(auto j = 0u; j < opi.line.size(); j++){						
					vis << "{";
					auto li = opi.line[j];
					auto file = li.file;
					if(file_loaded != file && file.length() > 4){
						auto file_end = file.substr(file.length()-4,4);
						if(file_end == ".csv" || file_end == ".txt"){
							tab = data.get_table(file,details.output_directory);
							file_loaded = file;
						}
					}
				
					vis << "\"name\":\"" << li.name << "\"";
			
					if(opi.type == OP_MARGINAL || opi.type == OP_GRAPH_MARGINAL){
						auto lab = 3u;
						if(opi.type == OP_GRAPH_MARGINAL) lab = 2;
						
						auto vec = data.get_table_column_str(lab-1,tab);
						vis << ",\"label\":[";
						for(auto k = 0u; k < vec.size(); k++){
							auto val = vec[k];
							if(val.length() > 3){
								if(val.substr(0,3) == "age") val = val.substr(3);
							}
							
							vis << "\"" << val << "\""; if(k < vec.size()-1) vis << ",";
						}				
						vis << "]";
					}
					
					if(plot_param_values == true && (opi.type == OP_GENERATION || opi.type == OP_PARAMDIST || opi.type == OP_TRACE) && j != 0 && j == opi.line.size()-1 && li.name == "Value"){
						if(get_double(file,"This parameter value ") != UNSET){
							vis << ",\"true\":" << file;
						}
					}
					else{
						auto xcol = li.xcol;
						if(xcol != UNSET){
							auto vec = data.get_table_column(xcol-1,tab);
							vis << ",\"xcol\":[";
							for(auto k = 0u; k < vec.size(); k++){
								vis << vec[k]; if(k < vec.size()-1) vis << ",";
							}				
							vis << "]";
						}
						
						auto ycol = li.ycol;
						if(ycol != UNSET){
							auto vec = data.get_table_column(ycol-1,tab);
						 
							vis << ",\"ycol\":[";
							for(auto k = 0u; k < vec.size(); k++){
								vis << vec[k]; if(k < vec.size()-1) vis << ",";
							}				
							vis << "]";
		
							auto EBcol = li.EB;
							if(EBcol != UNSET){
								auto vecmin = data.get_table_column(EBcol-1,tab);
								vis << ",\"errbarmin\":[";
								for(auto k = 0u; k < vecmin.size(); k++){
									vis << vecmin[k]; if(k < vecmin.size()-1) vis << ",";
								}				
								vis << "]";
									
								auto vecmax = data.get_table_column(EBcol,tab);
								vis << ",\"errbarmax\":[";
								for(auto k = 0u; k < vecmax.size(); k++){
									vis << vecmax[k]; if(k < vecmax.size()-1) vis << ",";
								}				
								vis << "]";
							}
							
							
							if((opi.type == OP_MARGINAL || opi.type == OP_GRAPH_MARGINAL) && j == 0 && details.siminf == INFERENCE){
								if(details.mode != SIM){
									auto vecmin = data.get_table_column(ycol,tab);
									vis << ",\"errbarmin\":[";
									for(auto k = 0u; k < vecmin.size(); k++){
										vis << vecmin[k]; if(k < vecmin.size()-1) vis << ",";
									}				
									vis << "]";
									
									auto vecmax = data.get_table_column(ycol+1,tab);
									vis << ",\"errbarmax\":[";
									for(auto k = 0u; k < vecmax.size(); k++){
										vis << vecmax[k]; if(k < vecmax.size()-1) vis << ",";
									}				
									vis << "]";
								}
							}
							
							if(opi.type == OP_ME && j == 0 && opi.line.size() > 1){
								auto sd = data.get_table_column(ycol,tab);
								vis << ",\"errbarmin\":[";
								for(auto k = 0u; k < sd.size(); k++){
									vis << vec[k]-sd[k]; if(k < sd.size()-1) vis << ",";
								}				
								vis << "]";
									
								vis << ",\"errbarmax\":[";
								for(auto k = 0u; k < sd.size(); k++){
									vis << vec[k]+sd[k]; if(k < sd.size()-1) vis << ",";
								}				
								vis << "]";
							}
						}
					}
					
					vis << ",\"style\":\"";
					switch(li.style){
						case NOLINE: vis << "NOLINE"; break;
						case RED_SOLID: vis << "RED_SOLID"; break;
						case RED_DASHED: vis << "RED_DASHED"; break;
						case GREEN_SOLID: vis << "GREEN_SOLID"; break;
						case GREEN_DASHED: vis << "GREEN_DASHED"; break;
						case BLUE_SOLID: vis << "BLUE_SOLID"; break;
						case BLUE_DASHED: vis << "BLUE_DASHED"; break;
						case BLACK_SOLID: vis << "BLACK_SOLID"; break;
						case BLACK_DASHED: vis << "BLACK_DASHED"; break;
						case RED_DOTTED: vis << "RED_DOTTED"; break;
						case RED_DOTDASH: vis << "RED_DOTDASH"; break;
						case GREEN_DOTTED: vis << "GREEN_DOTTED"; break;
						case GREEN_DOTDASH: vis << "GREEN_DOTDASH"; break;
						case BLUE_DOTTED: vis << "BLUE_DOTTED"; break;
						case BLUE_DOTDASH: vis << "BLUE_DOTDASH"; break;
						case BLACK_DOTTED: vis << "BLACK_DOTTED"; break;
						case BLACK_DOTDASH: vis << "BLACK_DOTDASH"; break;
						case YELLOW_SOLID: vis << "YELLOW_SOLID"; break;
						case YELLOW_DASHED: vis << "YELLOW_DASHED"; break;
						case CYAN_SOLID: vis << "CYAN_SOLID"; break;
						case CYAN_DASHED: vis << "CYAN_DASHED"; break;
						case MAGENTA_SOLID: vis << "MAGENTA_SOLID"; break;
						case MAGENTA_DASHED: vis << "MAGENTA_DASHED"; break;
						case RED_THIN: vis << "RED_THIN"; break;
						case GREEN_THIN: vis << "GREEN_THIN"; break;
						case BLUE_THIN: vis << "BLUE_THIN"; break;
						case BLACK_THIN: vis << "BLACK_THIN"; break;
					}
					vis << "\"";
					vis << "}";
					if(j < opi.line.size()-1) vis << ",";
				}
				vis << "]";
				break;
		}
		
		if(opi.spline_param_JSON != "") vis << ",\"spline_param\":" << opi.spline_param_JSON;
		
		vis << "}";
		
		if(i < op.size()-1) vis << ",";
	}
	vis << "]";

	vis << ",\"boundaries\" : " << boundaries;
	
	vis << "}";
	
	auto str = vis.str();
	
	for(auto i = 1u; i < str.length(); i++){                    // Ensures that ' character has \ before it
		if(str.substr(i,1) == "'" && str.substr(i-1,1) != "\\"){
			str.insert(i,"\\"); i++;
		}
	}
	
	auto file = details.output_directory+"/"+grfile;

	ofstream visout(file);
	if(!visout) emsg("Cannot open the file '"+file+"'");
		
	visout << "var jsonstr = '" << str << "';" << endl;
	
	cout << "  Open 'visBEEP.html' to visualise results." << endl;
}


/// Loads up the boundaries from the specified file and converts to a string for output
string Output::load_boundaries() const
{
	vector < vector < vector <Coord> > > bound;
	bound.resize(data.narea);
	
	if(area_plot.boundfile != "") data.load_boundaries(area_plot.boundfile,bound);
	else data.create_boundaries(area_plot.xcol,area_plot.ycol,bound);
	
	vector <unsigned int> list;
	for(auto c = 0u; c < data.narea; c++){ if(bound[c].size() == 0) list.push_back(c);}
	
	if(list.size() == data.narea){ warning("Could not find area boundary data so maps cannot be plotted."); return "";}
  
	if(list.size() > 0){
		stringstream ss; ss << "Could not find boundaries for the following areas: ";
		for(auto i = 0u; i < list.size(); i++){
			if(i != 0) ss << ", ";
			ss << data.area[list[i]].code;			
		}
		warning(ss.str()); 
		return "";
	}
	
	if(area_plot.project == EQUI_PROJ){                            // This performs an equirectangular projection
		double latav = 0.0, num = 0.0;
		for(auto c = 0u; c < data.narea; c++){
			for(auto i = 0u; i < bound[c].size(); i++){
				for(auto j = 0u; j < bound[c][i].size(); j++){
					latav += bound[c][i][j].y; num++;
				}
			}
		}
		latav /= num;
		
		if(latav > 180 && latav < -180) emsgroot("The latitute in 'file' is out of range");
		auto fac = cos(2*M_PI*latav/360.0);
		for(auto c = 0u; c < data.narea; c++){
			for(auto i = 0u; i < bound[c].size(); i++){
				for(auto j = 0u; j < bound[c][i].size(); j++){
					bound[c][i][j].x *= fac;
				}
			}
		}
	}
	
	data.rescale_boundary(bound);	
	
	if(area_plot.boundfile == "") data.make_circle_boundary(area_plot.xcol,bound);	
	
	auto d = 0.001;                                                     // Reduces the number of points
	
	auto npo = 0u;
	for(auto c = 0u; c < data.narea; c++){    
		auto xmin = LARGE, xmax = -LARGE, ymin = LARGE, ymax = -LARGE;
		for(auto i = 0u; i < bound[c].size(); i++){
			for(auto j = 0u; j < bound[c][i].size(); j++){
				if(bound[c][i][j].x < xmin) xmin = bound[c][i][j].x;
				if(bound[c][i][j].x > xmax) xmax = bound[c][i][j].x;
				if(bound[c][i][j].y < ymin) ymin = bound[c][i][j].y;
				if(bound[c][i][j].y > ymax) ymax = bound[c][i][j].y;
			}
		}
		auto dx = xmax-xmin, dy = ymax-ymin;
		auto dd = dx*0.01; if(dy > dx) dd = dy*0.01;
	
		for(auto i = 0u; i < bound[c].size(); i++){
			auto j = 1u;
			vector <Coord> boundnew;
			double x = LARGE, y = LARGE;
			for(j = 0; j < bound[c][i].size(); j++){
				auto dx = bound[c][i][j].x - x;
				auto dy = bound[c][i][j].y - y;
				auto r = dx*dx + dy*dy;
				if(r > d*d || r > dd*dd){
					x += dx; y += dy;
					boundnew.push_back(bound[c][i][j]);
				}
			}
			bound[c][i] = boundnew;
			
			npo += bound[c][i].size();
		}
	}
	
	stringstream ss;
	ss << "[";
	for(auto c = 0u; c < data.narea; c++){
		if(c != 0) ss << ",";
		ss << "[";
		for(auto i = 0u; i < bound[c].size(); i++){
			if(i != 0) ss << ",";
			ss << "[";
			for(auto j = 0u; j < bound[c][i].size(); j++){
				if(j != 0) ss << ",";
				ss << "[" << bound[c][i][j].x << "," << bound[c][i][j].y << "]";
			}
			ss << "]";
		}
		ss << "]";
	}
	ss << "]";

	return ss.str();
}


/// Outputs a spatial distribution for R as a function of time 
void Output::spatial_R_map(const vector <Sample> &opsamp, vector <OutputPlot> &op) const
{
	if(data.narea == 1) return; 
	
	auto nmap_out = opsamp[0].Rmap.size();
	for(auto m = 0u; m < nmap_out; m++){
		const auto &map_out = opsamp[0].Rmap[m];
		auto file = post_dir+"/"+map_out.file;
		auto filefull = details.output_directory+"/"+file;
		
		ofstream Rmapout(filefull.c_str());
		if(!Rmapout) emsg("Cannot open the file '"+filefull+"'");
		Rmapout << fixed;
		
		Rmapout << details.time_format_str; for(auto c = 0u; c < data.narea; c++) Rmapout << "," << data.area[c].code; Rmapout << endl;
		for(auto t = 0u; t < details.period; t++){
			Rmapout << details.getdate(t); 
			for(auto c = 0u; c < data.narea; c++){
				// Note credible intervals could be calculated here
				auto av = 0.0; for(const auto &opsa : opsamp) av += opsa.Rmap[m].map[t][c];
				av /= opsamp.size();
				
				Rmapout << "," << av;
			}
			Rmapout << endl;
		}
		
		OutputPlot oppl(OP_ANIM_MAP,file,map_out.fulldesc,map_out.tab,map_out.tab2,map_out.tab3,map_out.tab4,"","Colour",UNSET,UNSET);
		op.push_back(oppl);
	}
}


/// Saves a file giving time variation in different model quantities
void Output::spline_plots(const vector <Sample> &opsamp, vector <OutputPlot> &op) const
{
	if(opsamp.size() == 0) return;
	
	vector <SplineOutput> sim_spline_output;
	
	if(plot_param_values == true){  
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
		tout << fixed;
		
		if(suppress_output == false) cout << "'" << file << "' gives the time variation in "+name+"." << endl;

		bool plot_sim = false;
		if(plot_param_values == true && sim_spline_output[sp].splineval.size() > 0){
			if(sim_spline_output[sp].splineval[0] != UNSET) plot_sim = true;	
		}
		
		tout << "# Gives the time variation in "+name+"." << endl;	
		tout << "Time,mean,95% CI min,95% CI max";
		if(plot_sim == true) tout << ",Simulated";
		tout << endl;
		
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
			//spline_plots
			string line_desc;
			if(details.siminf == SIMULATE){
				if(details.mode == MULTISIM){
					line_desc = "  The red dashed lines denote 95% variation across simulations.";				
				}
				else{
					line_desc = "  The red lines represent the posterior mean, with the dashed lines denoting 95% credible intervals.";		
				}
				if(plot_sim == true){
					line_desc += " The black dashed line shows the true value used to simulate the data.";			
				}
			}
				
			OutputPlot oppl(OP_SPLINE,splout.desc,splout.fulldesc+line_desc,splout.tab,splout.tab2,splout.tab3,splout.tab4,"Time",name,min,max);
			oppl.spline_param_JSON = splout.spline_param_JSON;
			
			string lab;
			if(details.analysis_type != "Simulation"){
				oppl.addline(details.analysis_type,filefull,1,2,RED_SOLID);
				oppl.addline("95% CI min",filefull,1,3,RED_DASHED);
				oppl.addline("95% CI max",filefull,1,4,RED_DASHED);
			  lab = "True";
			}
			
			if(plot_sim == true) oppl.addline(lab,filefull,1,5,BLACK_DASHED);
			op.push_back(oppl);
		}
	}
}


/// Generates maps based on datatables
void Output::datatable_maps(const vector <Sample> &opsamp, vector <OutputPlot> &op) const
{
	if(data.narea == 1) return;
	
	for(auto &dt : data.datatable){
		if(dt.type != MARGINAL && dt.geo_dep != ""){		
			string fulldesc = "Geographical map: This shows a map giving the ";
			if(details.siminf == SIMULATE) fulldesc += "underlying ";
			else fulldesc += "posterior mean estimate for the ";
			fulldesc += "spatial and temporal variation in "+data.observation_description(dt.type,dt.observation)+" (note, this uses a logarithmic greyscale).  Clicking on the play button animates over time (note, the left and right arrow keys can be used to step single time units).";
			
			auto file = post_dir+"/state/map_"+dt.file;
			auto filefull = details.output_directory+"/"+file;
			ofstream resultout(filefull);
			if(!resultout) emsg("Cannot open the file '"+file+"'");
	
			resultout << fixed;
			
			resultout << details.time_format_str;
			for(auto c = 0u; c < data.narea; c++) resultout << "," << data.area[c].code;
			resultout << endl;
			
			for(auto t = 0u; t < details.period; t++){ 
				resultout << details.getdate(t);
				for(auto gr_num : dt.graph_ref){
					vector <double> vec; for(const auto &opsa : opsamp) vec.push_back(opsa.graph_state[gr_num][t]);		
					auto stat = get_statistic(vec);
					resultout << "," << stat.mean;
				}			
				resultout << endl;
			}	
		
			OutputPlot oppl(OP_ANIM_MAP,file,fulldesc,details.analysis_type,"Data Table",dt.file,"Map","","Greyscale",UNSET,UNSET);
			if(details.siminf != DATAVIEW) op.push_back(oppl);
			
			if(details.siminf == INFERENCE || details.siminf == DATAVIEW){
				auto file_data = post_dir+"/state/map_data_"+dt.file;
			
				vector <string> cols;
				cols.push_back(details.time_format_str); for(auto c = 0u; c < data.narea; c++) cols.push_back(data.area[c].code);
							
				data.generate_file_from_table(file_data,dt.file,cols);
				
				fulldesc = "Geographical map: This is a visualisation of the file *"+dt.file+"* which shows the ";
				fulldesc += "spatial and temporal variation in "+data.observation_description(dt.type,dt.observation,dt.timestep)+" (note, this uses a logarithmic greyscale).  Clicking on the play button animates over time (note, the left and right arrow keys can be used to step single time units).";
			
				OutputPlot oppl(OP_ANIM_MAP,file_data,fulldesc,"Data","Data Table",dt.file,"Map","","Greyscale",UNSET,UNSET);
				op.push_back(oppl);
			}
		}
	}
}


/// Plots a graph giving time variation in transition rate or population or marginal distributions
void Output::graph_plots(const vector <Sample> &opsamp, vector <OutputPlot> &op, const double invT) const
{
	if(opsamp.size() == 0) return;
	
	auto ngraph = data.graph.size();
	
	auto EBflag = false;                          // Sets if error bars are used for the obs model
	if(details.siminf == INFERENCE && invT != UNSET) EBflag = true;
	
	vector <GraphMultiPlot> graph_plot;
	
	for(auto gr_num = 0u; gr_num < ngraph; gr_num++){
		auto gr = data.graph[gr_num];
		const auto &dt = data.datatable[gr.datatable];

		auto file = details.output_directory+"/"+post_dir+"/state/"+gr.file;

		ofstream stateout(file.c_str());
		if(!stateout) emsg("Cannot open the file '"+file+"'");
		stateout << fixed;
			
		switch(gr.type){
			case GRAPH_TIMESERIES: stateout << "Time"; break;
			case GRAPH_MARGINAL: stateout << "Ref. Num.,Demographic state"; break;
		} 
		stateout << ",Mean";
		switch(stateuncer){
			case CI: 
				stateout << ",95% CI min,95% CI max" << endl; 
				break;
				
			case CURVES:
				for(auto sa = 0u; sa < opsamp.size(); sa++) stateout << ",Sample " << sa;
				stateout << endl;
				break;
		}
		
		auto imax = opsamp[0].graph_state[gr_num].size();
		auto min = LARGE, max = -LARGE;
		for(auto i = 0u; i < imax; i++){
			vector <double> vec; for(const auto &opsa : opsamp) vec.push_back(opsa.graph_state[gr_num][i]);		
			auto stat = get_statistic(vec);
		
			switch(gr.type){
				case GRAPH_TIMESERIES:
					stateout << details.division_time[i*details.graph_step+details.graph_step/2]; 
					break;				
				case GRAPH_MARGINAL: stateout << i << "," << dt.demolist[i]; break;
			}
			
			stateout << "," << stat.mean;
			switch(stateuncer){
				case CI: 
					stateout << "," << stat.CImin << "," << stat.CImax << endl;
					if(stat.CImax > max) max = stat.CImax;
					if(stat.CImin < min) min = stat.CImin;	
					break;
					
				case CURVES:
					for(const auto &opsa : opsamp){
						auto val = opsa.graph_state[gr_num][i]; 
						stateout << "," << val; 
						if(val > max) max = val; 
						if(val < min) min = val;	
					}
					stateout << endl;
					break;
			}
		}
		
		auto thresh_flag = false;
			
		string file_data, file_data_thresh;
		if(dt.optype == DATA){
			file_data = details.output_directory+"/"+post_dir+"/state/"+gr.file.substr(0,gr.file.length()-4)+"-data.csv";
	
			ofstream dataout(file_data.c_str());
			if(!dataout) emsg("Cannot open the file '"+file_data+"'");
			dataout << fixed;
			
			if(gr.type == GRAPH_MARGINAL){
				dataout << "Ref. Num,Demographic state,Data";
			}
			else{
				dataout << "Time,Data";
			}
			
			if(EBflag == true) dataout << ",EB min,EB max";
	
			dataout << endl;
			
			if(suppress_output == false) cout << file_data << " file" << endl;
		
			file_data_thresh = details.output_directory+"/"+post_dir+"/state/"+gr.file.substr(0,gr.file.length()-4)+"-data_thresh.csv";
		
			ofstream datathreshout;
			if(gr.type != GRAPH_MARGINAL){
				for(const auto &p : gr.point){ if(data.obs[p.obs].value == THRESH) thresh_flag = true;}
				
				if(thresh_flag == true){
					datathreshout.open(file_data_thresh.c_str());
					if(!datathreshout) emsg("Cannot open the file '"+file_data_thresh+"'");
					datathreshout << fixed;
			
					datathreshout << "Time,Data" << endl;
				}
			}
			
			for(auto i = 0u; i < gr.point.size(); i++){
				const auto &p = gr.point[i];
				const auto &ob = data.obs[p.obs];
				
				auto val = ob.value;
				
				auto EBmin = 0.0, EBmax = 0.0;
				if(EBflag == true){
					auto sd = 0.0;
					
					switch(dt.obsmodel){
						case NORMAL_OBSMODEL: sd = ob.sd; break;
						case NORMAL_PERCENT_OBSMODEL: sd = ob.sd; break;
						case POISSON_OBSMODEL: sd = sqrt(val); break;	
						case NEGBINO_OBSMODEL: sd = sqrt(val + val*val/ob.shape); break; 
					}
					
					EBmin = val - sd/sqrt(invT); EBmax = val + sd/sqrt(invT); 
				}
				
				if(dt.type == TRANS && val != UNKNOWN && val != THRESH){
					val /= p.xf - p.xi;
					EBmin /= p.xf - p.xi; EBmax /= p.xf - p.xi;
				}
				
				auto shift = 0.0; if(details.time_format == TIME_FORMAT_NUM) shift = details.start;
				
				string EBstr = "";
				if(EBflag == true) EBstr = ","+to_string(EBmin)+","+to_string(EBmax);
			
				switch(dt.type){
					case MARGINAL:
						if(val == THRESH)	dataout << p.xi << "," << dt.demolist[i] << "," << dt.threshold << EBstr << endl;
						else{
							if(val == UNKNOWN) dataout << p.xi << "," << dt.demolist[i] << "," << "0" << EBstr << endl;
							else dataout << p.xi << "," << dt.demolist[i] << "," << val << EBstr << endl;
						}
						break;
						
					case TRANS:
						if(val == THRESH){
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
								if(thresh_flag == true) datathreshout << endl;
							}
							else{
								if(details.stochastic == true){
									dataout << p.xi + shift << "," << val << EBstr << endl;;
									dataout << p.xf + shift << "," << val << EBstr << endl;;
								}
								else{
									dataout << 0.5*(p.xi+p.xf) + shift << "," << val << EBstr << endl;;
								}
								if(thresh_flag == true) datathreshout << endl;
							}
						}
						break;
					
					case POP: case POPFRAC:
						if(val == THRESH){
							datathreshout << 0.5*(p.xi+p.xf) + shift << "," << dt.threshold << endl;
						}
						else{
							if(val == UNKNOWN){
								if(thresh_flag == true) datathreshout << endl;
							}
							else{
								dataout << 0.5*(p.xi+p.xf) + shift << "," << val << EBstr << endl;
								if(thresh_flag == true) datathreshout << endl;
							}
						}
						break;
				}
								
				if(val != THRESH && val != UNKNOWN){
					if(val > max) max = val; 
					if(val < min) min = val;
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
					gp.file_data_EB.push_back(EBflag);
					if(thresh_flag == true)	gp.file_data_thresh.push_back(file_data_thresh);
				}
				new_plot = false;
			}
		}
		
		if(new_plot == true){
			GraphMultiPlot gp;
			gp.fulldesc = gr.fulldesc; gp.fulldesc_data = gr.fulldesc_data;
			gp.tab = gr.tab; gp.tab2 = gr.tab2; gp.tab3 = gr.tab3; gp.tab4 = gr.tab4;
			gp.plot_name = plot_name;
			gp.name.push_back(gr.name);
			if(!(gr.type == GRAPH_MARGINAL && details.analysis_type == "Simulation")) gp.file.push_back(file);
			gp.type = gr.type;
			gp.line_colour.push_back(dt.line_colour);
			gp.min = min;
			gp.max = max;
			if(file_data != ""){
				gp.file_data.push_back(file_data);
				gp.file_data_EB.push_back(EBflag);
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
				
		auto time_label = "Time";
		
		switch(gp.type){                                                // Generates the final plot
			case GRAPH_TIMESERIES:
				{
					auto ylabel = gp.name[0]; if(gp.file.size() > 1) ylabel = "Value";
		
					OutputPlot oppl(OP_GRAPH,gp.plot_name,gp.fulldesc,gp.tab,gp.tab2,gp.tab3,gp.tab4,time_label,ylabel,min,max);
				
					for(auto i = 0u; i < gp.file.size(); i++){
						auto legend = details.analysis_type; if(gp.file.size() > 1) legend = gp.name[i];
						auto file = gp.file[i];
						if(details.analysis_type == "Simulation"){
							oppl.addline(legend,file,1,2,lt[i]);
						}
						else{
							switch(stateuncer){
								case CI: 
									oppl.addline("95% CI min",file,1,3,lt2[i]);
									oppl.addline("95% CI max",file,1,4,lt2[i]);
									break;
							
								case CURVES:
									for(auto sa = 0u; sa < opsamp.size(); sa++){
										oppl.addline("",file,1,3+sa,lt2[i]);
									}
									break;
							}
							
							oppl.addline(legend,file,1,2,lt[i]);
						}
					}
					
					for(auto i = 0u; i < gp.file_data.size(); i++){
						if(gp.file_data_EB[i] == true){
							//oppl.addline("Data",gp.file_data[i],1,2,BLACK_SOLID,3);
							oppl.addline("Data",gp.file_data[i],1,2,BLACK_SOLID);
							oppl.addline("Obs Model",gp.file_data[i],1,3,BLACK_DASHED);
							oppl.addline("Obs Model",gp.file_data[i],1,4,BLACK_DASHED);
						}
						else oppl.addline("Data",gp.file_data[i],1,2,BLACK_SOLID);
					}
					
					for(auto file_data_thresh : gp.file_data_thresh) oppl.addline("Threshold",file_data_thresh,1,2,BLACK_DOTTED);
					
					if(details.siminf != DATAVIEW) op.push_back(oppl);
					
					if((details.siminf == INFERENCE || details.siminf == DATAVIEW) && gp.file_data.size() > 0){
						OutputPlot oppl(OP_GRAPH,gp.plot_name,gp.fulldesc_data,"Data",gp.tab2,gp.tab3,gp.tab4,time_label,ylabel,min,max);
						oppl.xmin = 0; oppl.xmax = details.period;
						
						for(auto i = 0u; i < gp.file_data.size(); i++){
							oppl.addline("Data",gp.file_data[i],1,2,BLACK_SOLID);
							
							if(gp.file_data_EB[i] == true){
								oppl.addline("Obs Model",gp.file_data[i],1,3,BLACK_DASHED);
								oppl.addline("Obs Model",gp.file_data[i],1,4,BLACK_DASHED);
							}
						}
					
						for(auto file_data_thresh : gp.file_data_thresh){
							oppl.addline("Threshold",file_data_thresh,1,2,BLACK_DOTTED);
						}
						
						op.push_back(oppl);
					}
				}
				break;
				
			case GRAPH_MARGINAL:
				{
					OutputPlot oppl(OP_GRAPH_MARGINAL,gp.plot_name,gp.fulldesc,gp.tab,gp.tab2,gp.tab3,gp.tab4,"Category",gp.name[0],min,max);
					
					if(gp.file.size() > 0){
						if(gp.file.size() != 1) emsgEC("Output",10);
						oppl.addline(details.analysis_type,gp.file[0],1,3,lt[0]);
					}
					
					for(auto file_data : gp.file_data) oppl.addline("Data",file_data,1,3,BLACK_SOLID);
					if(details.siminf != DATAVIEW) op.push_back(oppl);
					
					if((details.siminf == INFERENCE || details.siminf == DATAVIEW) && gp.file_data.size() > 0){
						OutputPlot oppl(OP_GRAPH_MARGINAL,gp.plot_name,gp.fulldesc_data,"Data",gp.tab2,gp.tab3,gp.tab4,"Category",gp.name[0],min,max);
		
						for(auto file_data : gp.file_data) oppl.addline("Data",file_data,1,3,BLACK_SOLID);
						op.push_back(oppl);
					}
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
			case RED:
				lt.push_back(RED_SOLID);
				if(stateuncer == CURVES) lt2.push_back(BLUE_THIN);
				else lt2.push_back(RED_DASHED);
			
				break;
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


/// Outputs a file containing posterior estimates for model parameters
void Output::posterior_parameter_estimates(const vector <ParamSample> &psamp, vector <OutputPlot> &op) const
{
	auto file = "Parameter_estimates.csv";
	auto filefull = details.output_directory+"/"+file;

	ofstream paramout(filefull.c_str());
	if(!paramout) emsg("Cannot open the file '"+filefull+"'");
	
	paramout << fixed << setprecision(3);
	
	if(suppress_output == false) cout << "'" << file << "' gives the model parameters." << endl;
	
	switch(details.siminf){
	case SIMULATE: paramout << "Name,Value" << endl; break;
	case INFERENCE: 
		paramout << "Name,Mean,95% CI,SD";
		if(details.mode != ML_INF) paramout << ",ESS,GRS";
		paramout << ",Prior" << endl;
		break;
	case DATAVIEW: paramout << "Name" << endl; break;
	}
	
	auto GR = get_Gelman_Rubin_statistic(psamp);
	
	auto ESS = get_effective_sample_size(psamp);
	
	vector <double> paramav(model.param.size());
	for(auto p = 0u; p < model.param.size(); p++){
		vector <double> vec; for(const auto psa : psamp) vec.push_back(psa.paramval[p]);
		auto stat = get_statistic(vec);
		paramav[p] = stat.mean;
		if(model.param[p].name != "zero" && model.param[p].name != "one"){
			switch(details.siminf){
			case SIMULATE:
				paramout << model.param[p].name  << "," << stat.mean;
				break;
				
			case INFERENCE:
				paramout << replace(model.param[p].name,"→","->")  << "," << stat.mean;
				paramout << "," << stat.CImin << " --- "<< stat.CImax << "," << stat.sd;
				if(details.mode != ML_INF){
					paramout << ",";
					if(ESS[p] == UNSET) paramout << "---"; else paramout << ESS[p];
					paramout << ",";
					if(GR[p] == UNSET) paramout << "---"; else paramout << GR[p];
				}
				paramout << "," << replace(model.print_prior(p),",","|");
				break;
				
			case DATAVIEW:
				paramout << replace(model.param[p].name,"→","->");
				break;
			}
			paramout << endl; 
			
			if(false){ // Outputs estimates to terminal 
				cout << model.param[p].name << "  " << stat.mean << "," <<  stat.CImin << "," << stat.CImax << "endl\n"; 
			}
		}
	}
	
	switch(details.siminf){
	case SIMULATE:
		cout << "  Parameter values given in '" << file << "'" << endl; 
		break;
	case INFERENCE:
		cout << "  " << details.analysis_type+" parameter estimates given in '" << file << "'" << endl;
		break;
	case DATAVIEW:
		break;
	}
	
	string fulldesc;
	switch(details.siminf){
	case SIMULATE:
		{
			fulldesc = "Simulation parameter values: This table provides the names and values for each model parameter.";
			fulldesc += "  Click on the parameter names (in blue) to link to where they appear in the model.";
		
			OutputPlot oppl(OP_PARAM_TABLE,file,fulldesc,"Parameters","","","","","",UNSET,UNSET);
			op.push_back(oppl);
		}
		break;
	
	case INFERENCE:
		{
			fulldesc = "Parameter priors: The prior captures any knowledge regarding parameter values before the data is considered. This table shows the priors places on each of the model parameters."; 
			fulldesc += "  Click on the parameter names (in blue) to link to where they appear in the model.";
			
			OutputPlot oppl_pr(OP_PRIOR_TABLE,file,fulldesc,"Prior","","","","","",UNSET,UNSET);
			op.push_back(oppl_pr);
		
			fulldesc = "Posterior parameter estimates: This table provides posterior estimates for each of the model parameters. The first column gives the parameter name followed by the posterior mean, 95% credible interval and standard deviation.";

			if(details.mode != ML_INF){
				fulldesc += " The effective sample size (ESS) provides an estimate of the number of independent posterior samples generated (exceeding 200 implies convergence) and the Gelman-Rubin statistic (GRS) measures how well independent runs map out the same posterior distribution (less than around 1.05 usually implies convergence). Note, ESS and GRS are not available for all inference methods.";
			}
			fulldesc += " Finally a column shows the prior (as defined in the TOML file) for comparison.";
		
			fulldesc += "  Click on the parameter names (in blue) to link to where they appear in the model.";
			
			OutputPlot oppl(OP_PARAM_TABLE,file,fulldesc,details.analysis_type,"Parameters","Table","","","",UNSET,UNSET);
			op.push_back(oppl);
		}
		break;
		
	case DATAVIEW:
		{
			fulldesc = "Model parameters: This table provides a list of each model parameter.";
			fulldesc += "  Click on the parameter names (in blue) to link to where they appear in the model.";
		
			OutputPlot oppl(OP_PARAM_TABLE,file,fulldesc,"Parameters","","","","","",UNSET,UNSET);
			op.push_back(oppl);
		}
		break;
	}
}
	
	
/// Outputs a bar chart of relative susceptibility split into demographic groups   
void Output::susceptibility_distributions(const vector <ParamSample> &psamp, vector <OutputPlot> &op) const 
{
	for(auto c = 0u; c < data.ndemocat; c++){
		auto amax = data.democat[c].value.size();
		if(data.democat[c].sus_vari == true && amax > 1){
			auto plot_sim_param = plot_param_values;
			for(auto a = 0u; a < amax; a++){ if(model.param[data.democat[c].sus_param[a]].value == UNSET) plot_sim_param = false;}
	
			auto type = data.democat[c].name;
			
			auto file = type+".csv";
			auto filefull = details.output_directory+"/"+post_dir+"/susceptibility/"+file;
			ofstream distout(filefull.c_str());
			if(!distout) emsg("Cannot open the file '"+filefull+"'");
			distout << fixed;
			
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
			
			string fulldesc = "Relative susceptibility: ";
			if(details.siminf == SIMULATE){
				fulldesc += "The red bars show the relative susceptibility across different *"+type+"* (note the average susceptibility, weighted by population size, is fixed to 1). ";
			}
			else{
				fulldesc += "The red bars show the posterior mean for the relative susceptibility across different *"+type+"* (note the average susceptibility, weighted by population size, is fixed to 1). The error bars show the 95% credible intervals.";
			
				if(plot_sim_param == true) fulldesc += "  The horizontal black lines show the true values used to simulate the data. ";
			}
			
			OutputPlot oppl(OP_MARGINAL,"Susceptibility with "+type,fulldesc,details.analysis_type,"Transmission","Susceptibility",type,type,"Relative susceptibility",UNSET,UNSET);
			oppl.label = label;
			if(details.analysis_type != "Simulation") oppl.addline(details.analysis_type,filefull,1,4,RED_SOLID);
			if(plot_sim_param == true) oppl.addline("Value",filefull,1,7,BLACK_SOLID);
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
	vector <string> previous_plots;
	
	for(const auto &co : model.comp){
		auto dep = co.mean_dep;
		if(co.num == 0 && dep != "" && find_in(previous_plots,co.name) == UNSET){ 
			auto demodep = get_demodep(dep);
			
			auto plot_sim_param = plot_param_values;
			if(details.analysis_type == "Simulation") plot_sim_param = true;
			
			for(auto i = 0u; i < demodep.size(); i++){
				if(model.param[co.param_mean[demodep[i].dp]].value == UNSET) plot_sim_param = false;
			}
		
			auto dir = details.output_directory+"/"+post_dir+"/compartment_mean";
			ensure_directory(dir);
			auto file = co.name+"_"+dep+".csv";
			auto filefull = dir+"/"+file;
		
			ofstream distout(filefull.c_str());
			if(!distout) emsg("Cannot open the file '"+filefull+"'");
	
			if(suppress_output == false) cout << "'" << file << "' gives stratified compartment '" << co.name << "' mean residency time." << endl;
	
			distout << "# Stratified mean compartmental residency time for compartment '" << co.name << "'." << endl;

			distout << "Ref.Num.,Ref.Lab,Demographic state";
			distout << ",Mean,95% CI min,95% CImax";
			if(plot_sim_param == true) distout << ",Simulated";
			distout << endl;
		
			vector <string> label;
			auto max = -LARGE, min = LARGE;
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
			
			if(min != max){
				string fulldesc = "Residency time: ";
				if(details.siminf == SIMULATE){
					fulldesc += "The red bars show the mean residency time in the *"+co.name+"* compartment for different *"+dep+"*. ";
				}
				else{
					fulldesc += "The red bars show the posterior mean for the mean residency time in the *"+co.name+"* compartment for different *"+dep+"*. The error bars show the 95% credible intervals.";
				
					if(plot_sim_param == true) fulldesc += "  The horizontal black lines show the true values used to simulate the data. ";
				}
				
				// This prevents the same plot multiple times
				previous_plots.push_back(co.name);
				
				OutputPlot oppl(OP_MARGINAL,"Compartment "+co.name+" mean residency time",fulldesc,details.analysis_type,"Compartmental model","Residency time",co.name,dep,"Residency time",UNSET,UNSET);
				oppl.label = label; 
				if(details.analysis_type != "Simulation") oppl.addline(details.analysis_type,filefull,1,4,RED_SOLID);
				if(plot_sim_param == true) oppl.addline("Value",filefull,1,7,BLACK_SOLID);
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
			auto demodep = get_demodep(dep);
			
			auto plot_sim_param = plot_param_values;
			if(details.analysis_type == "Simulation") plot_sim_param = true;
			
			for(auto i = 0u; i < demodep.size(); i++){
				if(model.param[tr.param_prob[demodep[i].dp]].value == UNSET) plot_sim_param = false;
			}
		
			auto dir = details.output_directory+"/"+post_dir+"/branch_prob";
			auto file = tr.name_file+"_"+dep+".csv";
			ensure_directory(dir);
			auto filefull = dir+"/"+file;
		
			ofstream distout(filefull.c_str());
			if(!distout) emsg("Cannot open the file '"+filefull+"'");
			distout << fixed;
			
			if(suppress_output == false) cout << "'" << file << "' gives stratified branching probability for transition '" << tr.name << "'." << endl;
	
			distout << "# Stratified branching probability for transition '" << tr.name << "'." << endl;

			distout << "Ref.Num.,Ref.Lab,Demographic state";
			distout << ",Mean,95% CI Min,95% CI Max";
			if(plot_sim_param == true) distout << ",Simulated";
			distout << endl;
		
			vector <string> label;
			auto max = -LARGE, min = LARGE;
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
	
			if(min != max){
				string fulldesc = "Branching probability: ";
				if(details.siminf == SIMULATE){
					fulldesc += "The red bars show branching probabilities for the *"+tr.name+"* transition for different *"+dep+"*. ";
				}
				else{
					fulldesc += "The red bars show the posterior mean branching probabilities for the *"+tr.name+"* transition for different *"+dep+"*. The error bars show the 95% credible intervals.";
				
					if(plot_sim_param == true) fulldesc += "  The horizontal black lines show the true values used to simulate the data. ";
				}
				
				OutputPlot oppl(OP_MARGINAL,"Branching probability for transition *"+tr.name+"*",fulldesc,details.analysis_type,"Compartmental model","Branching",tr.name,dep,"Probability",UNSET,UNSET);

				oppl.label = label;
				if(details.analysis_type != "Simulation") oppl.addline(details.analysis_type,filefull,1,4,RED_SOLID);
				
				if(plot_sim_param == true) oppl.addline("Value",filefull,1,7,BLACK_SOLID);
				op.push_back(oppl);	
			}
		}
	}
}


/// Plots marginal distribitions showing how the age mixing matrix is modified at different break points
void Output::age_mixing_perturb_distributions(const vector <ParamSample> &psamp, vector <OutputPlot> &op) const 
{
	const auto &matmod = data.genQ.matmod; 
	
	if(matmod.size() == 0) return;
	if(matmod.size() != data.nage) emsgEC("Output",32);
	
	auto plot_sim_param = plot_param_values;
	if(details.analysis_type == "Simulation") plot_sim_param = true;
		
	auto nbp = matmod[0].bp.size();
	for(auto bp = 0u; bp < nbp; bp++){
		auto d = matmod[0].bp[bp];

		auto dir = details.output_directory+"/"+post_dir+"/age_mixing";
		auto file = "Day "+to_string(d)+".csv";
		ensure_directory(dir);
		auto filefull = dir+"/"+file;
	
		ofstream distout(filefull.c_str());
		if(!distout) emsg("Cannot open the file '"+filefull+"'");
		distout << fixed;
			
		if(suppress_output == false) cout << "'" << file << "' gives stratified age mixing modification on day " << d << "." << endl;

		distout << "# Stratified age mixing modification on day " << d << "." << endl;

		distout << "Ref.Num.,Ref.Lab,Age";
		distout << ",Mean,95% CI Min,95% CI Max";
		if(plot_sim_param == true) distout << ",Simulated";
		distout << endl;
	
		vector <string> label;
		auto max = -LARGE, min = LARGE;
		for(auto i = 0u; i < data.nage; i++){
			distout << i << "," << reference_label(i) << "," << data.democat[0].value[i];
				
			auto &spl = model.spline[matmod[i].spline_ref];
			
			auto th = spl.p[bp].param;
			
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

			string fulldesc = "Age mixing modification: ";
			/*
			if(details.siminf == SIMULATE){
				fulldesc += "The red bars show branching probabilities for the *"+tr.name+"* transition for different *"+dep+"*. ";
			}
			else{
				fulldesc += "The red bars show the posterior mean branching probabilities for the *"+tr.name+"* transition for different *"+dep+"*. The error bars show the 95% credible intervals.";
			
				if(plot_sim_param == true) fulldesc += "  The horizontal black lines show the true values used to simulate the data. ";
			}
			*/
			
		OutputPlot oppl(OP_MARGINAL,"Age mixing modification",fulldesc,details.analysis_type,"Age Mixing","Break points","Day "+to_string(d),"Age","Value",UNSET,UNSET);

		oppl.label = label;
		if(details.analysis_type != "Simulation") oppl.addline(details.analysis_type,filefull,1,4,RED_SOLID);
		
		if(plot_sim_param == true) oppl.addline("Value",filefull,1,7,BLACK_SOLID);

		op.push_back(oppl);
	}
}
	
	
/// Outputs a bar chart giving level effects split by demographic groups
void Output::level_effect_distributions(const vector <ParamSample> &psamp, vector <OutputPlot> &op) const 
{
	const auto &list = model.level_effect_param_list;
	
	auto plot_sim_param = plot_param_values;
	for(auto i = 0u; i < list.size(); i++){
		if(model.param[list[i]].value == UNSET) plot_sim_param = false;
	}
			
	if(data.level_effect.on == false) return;

	auto file = post_dir+"/level_effect.csv";
	auto filefull = details.output_directory+"/"+file;
		
	ofstream distout(filefull.c_str());
	if(!distout) emsg("Cannot open the file '"+filefull+"'");
	distout << fixed;
			
	if(suppress_output == false) cout << "'" << file << "' gives level effects." << endl;

	distout << "# Level effects." << endl;

	distout << "Ref.Num.,Ref.Lab,Parameter,Mean,95% CI Min,95% CI Max";
	if(plot_sim_param == true) distout << ",Simulated";
	distout << endl;
		
	auto max = -LARGE, min = LARGE;
		
	vector <string> label;
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
		auto fulldesc = "Levels: This shows the relative rates of disease transmission for different levels.";
		
		OutputPlot oppl(OP_MARGINAL,"Relative transmission for different levels",fulldesc,details.analysis_type,"Transmission","Levels","","Level Effect","Relative transmission rate",0,0);
		oppl.label = label;
		if(details.analysis_type != "Simulation") oppl.addline(details.analysis_type,filefull,1,4,RED_SOLID);
		if(plot_sim_param == true) oppl.addline("Value",filefull,1,7,BLACK_SOLID);
		op.push_back(oppl);	
	}
}


/// Outputs a bar chart giving relative transmission split by area
void Output::area_effect_distributions(const vector <ParamSample> &psamp, vector <OutputPlot> &op) const 
{
	const auto &list = model.area_effect_param_list;
	
	auto plot_sim_param = plot_param_values;
	for(auto i = 0u; i < list.size(); i++){
		if(model.param[list[i]].value == UNSET) plot_sim_param = false;
	}
	
	if(data.area_effect.on == false) return;

	auto file = post_dir+"/area_effect.csv";
	auto filefull = details.output_directory+"/"+file;
		
	ofstream distout(filefull.c_str());
	if(!distout) emsg("Cannot open the file '"+filefull+"'");
	distout << fixed;
			
	if(suppress_output == false) cout << "'" << file << "' gives area effects." << endl;

	distout << "# Area effects." << endl;

	distout << "Ref.Num.,Ref.Lab,Parameter,Mean,95% CI Min,95% CI Max";
	if(plot_sim_param == true) distout << ",Simulated";
	distout << endl;
		
	auto max = -LARGE, min = LARGE;
		
	vector <string> label;
	for(auto i = 0u; i < list.size(); i++){
		auto th = list[i];
		vector <double> vec; for(const auto &psa : psamp) vec.push_back(psa.paramval[th]);
		auto stat = get_statistic(vec);
		
		label.push_back(model.param[th].name);
		distout << i << "," << reference_label(i) << ",";
		auto name = model.param[th].name;
		if(name.length() > 12){
			if(name.substr(0,12) == "Area effect ") name = name.substr(12);
		}	
		
		distout << name << "," << stat.mean << "," << stat.CImin << "," << stat.CImax;	
		if(plot_sim_param == true){
			auto num = model.param[th].value; if(num > max) max = num; if(num < min) min = num;
			distout << "," << num;
		}
		distout << endl;
			
		if(stat.CImax > max) max = stat.CImax;
		if(stat.CImin < min) min = stat.CImin;
	}
	
	if(min != max){
		string fulldesc = "Area effects: ";
		if(details.siminf == SIMULATE){
			fulldesc += "This histogram shows the factor &α_a& change in &R& for different areas (note, this factor is defined to have a poplation average of one). ";
		}
		else{
			fulldesc += "This histogram shows the posterior distribution for the factor &α_a& change in &R& for different areas (note, this factor is defined to have a poplation average of one). The error bars give 95% credible intervals.";
		
			if(plot_sim_param == true) fulldesc += "  The horizontal black lines show the true values used to simulate the data. ";
		}
		
		OutputPlot oppl(OP_MARGINAL,"Relative transmission for different areas",fulldesc,details.analysis_type,"Transmission","Area Effect","","Area Effect","Relative transmission rate",0,0);
		oppl.label = label;
		if(details.analysis_type != "Simulation") oppl.addline(details.analysis_type,filefull,1,4,RED_SOLID);
		if(plot_sim_param == true) oppl.addline("Value",filefull,1,7,BLACK_SOLID);
		op.push_back(oppl);	
	}
}


/// Outputs the probability distributions for derived quantities (generation time, external infections etc..)
void Output::derived_parameter_distributions(const vector <Sample> &opsamp, vector <OutputPlot> &op) const
{
	if(details.siminf == SIMULATE) return;
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
			distout << fixed;
			
			if(suppress_output == false){
				cout << "'" << file << "' gives the probability distributions for the "+ derpar.desc+"." << endl;
			}
			
			distout << "# " << details.analysis_type << " probability distributions for the "+ derpar.desc+"." << endl;

			
			distout << "Value,Probability" << endl;
			
			for(auto b = 0u; b < paramdist.value.size(); b++){
				distout << paramdist.value[b] << "," << paramdist.prob[b] << endl;
			}
			
			auto fulldesc = "Posterior distribution: The posertior probability distribution for the "+ toLower(derpar.desc)+".";
			
			OutputPlot oppl(OP_PARAMDIST,derpar.desc,fulldesc,details.analysis_type,"Parameters","Derived Dist.",derpar.name,derpar.name,"Probability",0,0);
			oppl.addline(details.analysis_type,filefull,1,2,GREEN_SOLID);
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
	distout << fixed;
			
	if(suppress_output == false) cout << "'" << file << "' gives the probability distributions for parameters." << endl;
	
	distout << "# " << details.analysis_type << " probability distributions for model parameters." << endl;
	
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
		
		auto plot_sim_param = plot_param_values;
		if(model.param[p].value == UNSET) plot_sim_param = false;
		
		auto title = details.analysis_type+" distribution for "+model.param[p].name;
		auto fulldesc = "Posterior distribution: This shows the "+toLower(details.analysis_type)+" distribution for the model parameter *"+model.param[p].name+"*.";
		fulldesc += "  The vertical dashed line give the true value used to simulate the data.";
	
		OutputPlot oppl(OP_PARAMDIST,title,fulldesc,details.analysis_type,"Parameters","Distributions",model.param[p].name,model.param[p].name,"Probability",0,0);
		oppl.addline(details.analysis_type,filefull,2*i+1,2*i+2,GREEN_SOLID);
		if(plot_sim_param == true) oppl.addline("Value",to_string(model.param[p].value),1,2,BLACK_DASHED);
		op.push_back(oppl);
	}
}


/// Adds the plots which show how quantities change as a function of generation number (ABC-SMC / ABC-MBP)
void Output::add_generation_plots(vector <OutputPlot> &op) const
{
	string title, fulldesc, label, tab3, tab4;
	bool EF_flag = false;
	unsigned int col;
	
	switch(details.mode){
		case ABC_MBP: case ABC_SMC: case ABC_DA: case ABC_CONT:
			EF_flag = true; label = "log(EF cut-off)"; tab3 = "Error function"; tab4 = "Cut-off"; col = 4;
			break;
		
		case PAIS_INF: 
			label = "log(invT)"; tab3 = "Inverse temperature"; tab4 = "Value"; col = 6;
			break;
		
		case ML_INF:
			EF_flag = true; label = "log(Average EF)"; tab3 = "Error function"; tab4 = "Mean"; col = 4;
			break;
			
		default: return;
	}
		
	auto file = details.output_directory+"/Diagnostics/Generation.csv";
	
	if(label != ""){
		title = label+" as a function of generation";
		
		if(EF_flag == true){
			if(details.mode == ML_INF){
				fulldesc = "Error function reduction: The error function is defined as -2 times the log of the posterior (i.e. the sum of the logs of the likelihood and prior). Minimisation of this quantitiy gives the maximum &a& &posteriori& (MAP) estimate. This graph shows how the average error function reduces as a function of generation number as CMA-ES proceeds. Convergence is indicated by the curve plateauing to a minimum value.";
			}
			else{
				fulldesc = "Error function reduction: This graph shows how the cut-off in the error function reduces as a function of generation number. A smaller error function implies a closer match between the states generated by the model and the actual data.";
			}
		}
		else{
			fulldesc = "Inverse temperature increase: This graph shows how the increases as a function of generation number. A larger inverse temperature implies a closer match between the states generated by the model and the actual data.";
		}
		
		OutputPlot oppl(OP_LOG_GENERATION,title,fulldesc,"Posterior","Diagnostics",tab3,tab4,"Generation",label,UNSET,UNSET);
		oppl.addline(tab3,file,1,col,BLACK_SOLID);
		op.push_back(oppl);
	}
	
	for(auto p = 0u; p < model.param.size(); p++){                           // These are graphs for the parameters
		if(model.param[p].priortype != FIXED_PRIOR){
			auto plot_sim_param = plot_param_values;
			if(model.param[p].value == UNSET) plot_sim_param = false;
	
			double min = UNSET, max = UNSET; 
			if(model.param[p].priortype == UNIFORM_PRIOR){
				min = model.param[p].val1; max = model.param[p].val2;
			}
			
			auto name = replace(model.param[p].name,":","=");
			auto title = details.analysis_type+" estimate for "+name+" as a function of generation";
			string fulldesc;

			if(details.mode == ML_INF){
				fulldesc = "Parameter convergence: This shows the convergence of the posterior estimate for the model parameter *"+name+"* onto a maximim &a& &posterior& (MAP) estimate as a function of the number of generations. Note, here the dashed lines DO NOT represent credible intervals. Under CMAES the MAP estimate is first generated and then credible intervals are calculated using the local curvature of the posterior near to the estimate (using an approximation to the inverse Hessian matrix)."; 
			}
			else{
				fulldesc= "Parameter convergence: This shows the convergence of the posterior estimate for the model parameter *"+name+"* as a function of the number of generations.";
			}
			
			if(plot_sim_param == true) fulldesc += "  The horizontal dashed line shows the true parameter value used to generate the data.";
			
			OutputPlot oppl(OP_GENERATION,title,fulldesc,"Posterior","Diagnostics","Convergence",model.param[p].name,"Generation",model.param[p].name,min,max);
			oppl.addline(details.analysis_type,file,1,3*p+9,BLUE_SOLID);
			oppl.addline("95% CI min",file,1,3*p+10,BLUE_DASHED);
			oppl.addline("95% CI max",file,1,3*p+11,BLUE_DASHED);
			if(plot_sim_param == true) oppl.addline("Value",to_string(model.param[p].value),1,2,BLACK_DASHED);
			op.push_back(oppl);
		}
	}
				
	auto fileDT = details.output_directory+"/Diagnostics/EF_datatable.csv";
	string ylab;
	if(EF_flag == true){
		title = "Contribution to EF from different data sources";
		fulldesc = "Stratified error function: This graph shows the error function (EF) contributions from each of the different data sources as a function of generation number. The lines with the highest EF indicate which data sources the model finds hardest to generate from the model.";
		ylab = "log(EF contribution)";
	}
	else{
		title = "Contribution to observation model from different data sources";
		fulldesc = "Stratified observation model: This graph shows the observation model contributions from each of the different data sources as a function of generation number. The lines with the highest value indicate which data sources the model finds hardest to generate from the model.";
		ylab = "log(Obs. model contribution)";
	}
	
	if(false){
		OutputPlot opplDT(OP_LOG_GENERATION,title,fulldesc,"Posterior","Diagnostics",tab3,"By data source","Generation",ylab,UNSET,UNSET);
		auto c = 2u;
		for(auto i = 0u; i < data.datatable.size(); i++){
			if(data.datatable[i].file != ""){
				opplDT.addline(rus(data.datatable[i].file),fileDT,1,c,get_linestyle(i));
				c++;
			}
		}

		op.push_back(opplDT);
	}
	
	if(details.mode != ML_INF){
		auto fulldesc = "CPU time: This graph shows the CPU time required to achieve a given cut-off in the error function. This graph can be used to compare the efficiency of the different inference algorithms.";
		
		OutputPlot oppl(OP_CPU,"CPU time as a function of EF",fulldesc,"Posterior","Diagnostics","Error function","CPU time","EF cut-off","CPU time (minutes)",UNSET,UNSET);
		oppl.addline("",file,3,2,BLACK_SOLID);
		op.push_back(oppl);
	}
	
	if(label != "" && details.mode != ML_INF){
		title = "Model evidence as a function of "+label;
		fulldesc = "Model Evidence: This graph shows the model evidence (ME) as a function of ";
		if(EF_flag == true) fulldesc += "the error function cut-off";
		else fulldesc += "the inverse temperature";
		fulldesc += " (on a log-log plot).";
		if(nrun > 1){
			fulldesc += " The error bars show the standard deviation across multiple runs of the algorithm.";
		}
		
		fulldesc += "  The ME is a measure for comparing how well the data agrees with the model. When performing model comparison care must be taken to ensure it is done so at a consistent value of ";
		if(EF_flag == true) fulldesc += "error function cut-off.";
		else fulldesc += "inverse temperature.";

		fulldesc += "  The ratio in ME gives the Bayes factors between models. A Bayes factor above 3.2 is considered substantial support for one model over another, and over 10 provides strong evidence.";
		
		
		OutputPlot oppl(OP_ME,title,fulldesc,"Model Evidence","","","",label,"log(Model evidence)",UNSET,UNSET);
		oppl.addline("",file,col,7,BLACK_SOLID);
		if(nrun > 1) oppl.addline("Error bar",file,col,7,BLACK_SOLID);
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
			auto fulldesc = "Trace plot: This shows a trace plot for the parameter *"+model.param[p].name+"*.";
			
			OutputPlot oppl(OP_TRACE,"",fulldesc,"Posterior","Diagnostics","Trace",model.param[p].name,"Sample",model.param[p].name,UNSET,UNSET);
			
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
			if(plot_param_values == true && model.param[p].value != UNSET){
				oppl.addline("Value",to_string(model.param[p].value),1,2,BLACK_DASHED);
			}
			
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
		if(dcc.democats_filt != ""){ 
			file += "_"+dcc.democats_filt; desc += " for "+dcc.democats_filt;
			
		}
		if(dcc.geo_filt != ""){ file += "_"+dcc.geo_filt; desc += " for "+dcc.geo_filt;}
		file += ".csv";
		auto filefull = details.output_directory+"/"+post_dir+"/democat_change/"+file;
		
		ofstream lout(filefull.c_str());
		if(!lout) emsg("Cannot open the file '"+filefull+"'");
		lout << fixed;
			
		lout << "# " << desc << endl;
		lout << "Time";
		for(auto v = 0u; v < data.democat[d].value.size(); v++) lout << "," << data.democat[d].value[v]; 
		lout << endl;
		
		for(auto sett = 0u; sett < details.ndivision; sett++){
			lout << details.division_time[sett]; 
			for(auto v = 0u; v < data.democat[d].value.size(); v++) lout << "," << 100*dcc.frac[sett][v];
			lout << endl;
		}
	
		string fulldesc = "Democat change: This graph shows the imposed time variation in the *"+data.democat[d].name+"* demographic category";
		if(dcc.democats_filt != "") fulldesc += " for "+replace(dcc.democats_filt,":","=");
		if(dcc.geo_filt != "") fulldesc += " for "+replace(dcc.geo_filt,":","=");
		fulldesc += ".";
		
		OutputPlot oppl(OP_GRAPH,desc,fulldesc,"Data","Democat Change",data.democat[d].name,dcc.file,"Time","Percent",0,100);
		for(auto v = 0u; v < data.democat[d].value.size(); v++){
			oppl.addline(data.democat[d].value[v],filefull,1,2+v,get_linestyle(v));
		}
		op.push_back(oppl);
	}
}


/// Outputs the compartmental model
void Output::compartmental_model(vector <OutputPlot> &op) const
{
	string fulldesc = "Compartmental model: ";
	fulldesc += "The force of infection &λ_t& gives the probability per unit time of an individual becoming infected. The occupancy time within a compartment (assuming an onwards transition) can either be exponentially distributed (denoted Exp(&m&) where &m& is a mean time) or Erlang distributed (denoted Γ(&m&,&k&) where &m& is a mean time and &k& is an integer shape parameter). The probability an individual goes down one branch or another is determined by &b&. The quantities &λ_t&, &m& and &b& may depend on the demographic classification &d& of an individual.  Clicking on the distributions or branching probabilities shows how they are set by the underlying model parameters."; 
	OutputPlot oppl(OP_COMP_MODEL,"The input TOML file",fulldesc,"Model","Compartments","","","","",UNSET,UNSET);
	op.push_back(oppl);
	
	fulldesc = "The force of infection: This gives the probability per unit time of a susceptible individual becoming infected. This quantity may vary with time, depend on the area the individual resides in, their demographic group, or the strain of the virus.  Click on the terms highlighted in blue to see how model parameters inform them (details of the dependency are given in the manual)."; 
	OutputPlot oppl2(OP_FOI_MODEL,"The input TOML file",fulldesc,"Model","Force of infection","","","","",UNSET,UNSET);
	op.push_back(oppl2);
	
	if(details.description != "UNSET"){
		OutputPlot oppl(OP_DESC,"The input TOML file","Analysis description: "+details.description,"Model","Description","","","","",UNSET,UNSET);
		op.push_back(oppl);
	}
}


/// Outputs covariate data
void Output::covar_data(vector <OutputPlot> &op) const
{
	for(auto i = 0u; i < data.covar.size(); i++){
		const auto &cv = data.covar[i];
		
		string fulldesc, desc, xaxis;
		OutputPlotType type=OP_AREA_COVAR;
		
		if(cv.func == LOG_TRANS) xaxis = "Log transformed";
		switch(cv.type){
			case AREA_COVAR:
				fulldesc = "Covariate *"+cv.name+"*: This is a visualistation of the *"+cv.name+"* area covariate. ";	
				desc = "This uses the *"+cv.name+"* column in the file *"+data.areas_file+"*.";
				type = OP_AREA_COVAR;
				break;
				
			case TV_COVAR:
				fulldesc = "Time-varying covariate *"+cv.name+"*: This is a visualistation of the *"+cv.name+"* time-varying covariate in the file *"+cv.file+"*. ";
				desc = "This uses the *"+cv.name+"* column in the file *"+cv.file+"*.";
				xaxis = "Time";
				type = OP_TV_COVAR;
				break;
				
			case AREA_TV_COVAR:
				{
					auto file_data = post_dir+"/covar/map_data_"+cv.file;
			
					vector <string> cols;
					cols.push_back(details.time_format_str); for(auto c = 0u; c < data.narea; c++) cols.push_back(data.area[c].code);
							
					data.generate_file_from_table(file_data,cv.file,cols);

					fulldesc = "Spatial time-varying covariate: This is a visualistation of *"+cv.name+"*, a spatial and time-varying covariate defined in the file *"+cv.file+"*. ";
					desc = file_data;
					type = OP_AREA_TV_COVAR;
				}
				break;
		}
		
		if(cv.func == LOG_TRANS) fulldesc += "This quantity is log-transformed before being used. ";
		fulldesc += " The parameter *"+model.param[model.covariate_param[i]].name+"* is a fixed effect used to incorporate this covariate into the model. ";
		if(details.siminf == INFERENCE) fulldesc += "This parameter is estimated during inference."; 
		
		OutputPlot oppl(type,desc,fulldesc,"Data","Covariates",cv.name,"",xaxis,cv.name,UNSET,UNSET);
	
		if(cv.type == TV_COVAR){
			auto file = post_dir+"/covar/data_"+cv.name+".csv";
			data.make_table_with_time(cv.file,file,cv.name);
			oppl.addline("",file,1,2,RED_SOLID);
		}
		
		op.push_back(oppl);
	}
}


/// Outputs level effect data
void Output::level_data(vector <OutputPlot> &op) const
{
	const auto &le = data.level_effect;
	if(le.on != true) return;
	
	auto file_data = post_dir+"/covar/map_level_effect.csv";
	auto filefull = details.output_directory+"/"+file_data;
	
	ofstream fout(filefull);
	if(!fout) emsg("Cannot open the file '"+filefull+"'");
	
	fout << details.time_format_str;
	for(auto c = 0u; c < data.narea; c++) fout << "," << data.area[c].code;
	fout << endl;
	
	for(auto t = 0u; t < le.param_map.size(); t++){
		fout << details.getdate(t); for(auto c = 0u; c < data.narea; c++) fout << "," << le.param_map[t][c];
		fout << endl;
	}
				
	auto fulldesc = "Level effects: This is a visualistation of the file *"+le.file+"*, which shows the restriction level each area is in as a function of time. ";
	auto desc = file_data;
	
	OutputPlot oppl(OP_LEVEL_EFFECT,desc,fulldesc,"Data","Level effects","","","","",UNSET,UNSET);
	op.push_back(oppl);
}


/// Outputs a map giving the spatial mixing between individuals
void Output::spatial_mixing_map(vector <OutputPlot> &op) const
{
	if(data.narea == 1) return;
	
	{
		auto fulldesc = "Spatial mixing within areas: This shows the average percentage of daily contacts (which could potentially transmit infection) made by an individual within a given area with other individuals in that area (with all other contacts made elsewhere). This visualises the *"+data.genQ.M_name+"* input file.";
		OutputPlot oppl(OP_MIXING_WITHIN_MAP,"This is a visualistation of the '"+data.genQ.M_name+"' input file",fulldesc,"Data","Spatial Mixing","Within areas","","","",UNSET,UNSET);
		op.push_back(oppl);
	}
	
	string mod_name;	
	if(model.spline[model.geo_spline_ref].name != "UNSET") mod_name = "&m_t&";
	
	if(mod_name != ""){
		auto fulldesc = "Spatial mixing within areas: This shows the average percentage of daily contacts (which could potentially transmit infection) made by an individual within a given area with other individuals in that area (with all other contacts made elsewhere).  Clicking on the play button animates over time (note, the left and right arrow keys can be used to step single time units). This visualises the *"+data.genQ.M_name+"* input file modulated by the spline *"+mod_name+"*.";
		OutputPlot oppl(OP_MIXING_WITHIN_ANIM_MAP,"This is a visualistation of the '"+data.genQ.M_name+"' input file",fulldesc,details.analysis_type,"Spatial Mixing","Within areas","","","",UNSET,UNSET);
		op.push_back(oppl);
	}
	
	{
		auto fulldesc = "Spatial mixing between areas: This shows the average percentage of contacts made by an individual within a given area with individuals in other areas (thicker blue arrows indicate a higher percentage and the black circles represent the geographical centers of each area). Hovering over an area shows results for that area only. This visualises the *"+data.genQ.M_name+"* input file.";
		OutputPlot oppl(OP_MIXING_BETWEEN_MAP,"This is a visualistation of the '"+data.genQ.M_name+"' input file",fulldesc,"Data","Spatial Mixing","Between areas","","","",UNSET,UNSET);
		op.push_back(oppl);
	}
	
	if(mod_name != ""){
		auto fulldesc = "Spatial mixing between areas: This shows the average percentage of contacts made by an individual within a given area with individuals in other areas (thicker blue arrows indicate a higher percentage and the black circles represent the geographical centers of each area).  Clicking on the play button animates over time (note, the left and right arrow keys can be used to step single time units). Hovering over an area shows results for that area only. This visualises the *"+data.genQ.M_name+"* input file modulated by the spline *"+mod_name+"*.";
		OutputPlot oppl(OP_MIXING_BETWEEN_ANIM_MAP,"This is a visualistation of the '"+data.genQ.M_name+"' input file",fulldesc,details.analysis_type,"Spatial Mixing","Between areas","","","",UNSET,UNSET);
		op.push_back(oppl);
	}
	
	/*
	if(mod_name != ""){
		auto fulldesc = "Spatial population flow between areas: This shows an animation giving the effective population movement between areas (thicker dark blue lines indicate a higher percentage). This visualises the *"+data.genQ.M_name+"* input file modulated by the spline *"+mod_name+"*. Hovering over an area shows only results for that area.";
		OutputPlot oppl(OP_MIXING_POP_ANIM_MAP,"This is a visualistation of the '"+data.genQ.M_name+"' input file",fulldesc,details.analysis_type,"Spatial Mixing","Population flow","","",UNSET,UNSET);
		op.push_back(oppl);
	}
	else{
		auto fulldesc = "Spatial population flow between areas: This shows an animation giving the effective population movement between areas (thicker dark blue lines indicate a higher percentage).This visualises the *"+data.genQ.M_name+"* input file modulated by the spline *"+mod_name+"*. Hovering over an area shows only results for that area.";
		OutputPlot oppl(OP_MIXING_POP_MAP,"This is a visualistation of the '"+data.genQ.M_name+"' input file",fulldesc,details.analysis_type,"Spatial Mixing","Population flow","","",UNSET,UNSET);
		op.push_back(oppl);
	}
	*/
}
	
	
/// Outputs a matrix giving the mixing between different age groups
void Output::age_mixing_matrix(const vector <Sample> &opsamp, vector <OutputPlot> &op) const
{
	if(data.nage == 1) return;
	
	auto fulldesc = "Age Mixing matrix: This shows the relative contact rate with which different age groups mix within the population (darker red colours indicate more frequent contacts). This visualises the *"+data.genQ.N_name[0]+"* input file.";
	OutputPlot oppl(OP_AGE_MATRIX,"This is a visualistation of the '"+data.genQ.N_name[0]+"' input file",fulldesc,"Data","Age Mixing","","","","",UNSET,UNSET);
	op.push_back(oppl);
	
		
	auto file = post_dir+"/age_mixing/matrix.csv";
	auto filefull = details.output_directory+"/"+file;
		
	ofstream Nout(filefull.c_str());
	if(!Nout) emsg("Cannot open the file '"+filefull+"'");
	Nout << fixed;
			
	auto Nsize = data.nage;
	
	Nout << "Age category";
	for(auto c = 0u; c < Nsize; c++) Nout << "," << data.democat[0].value[c] << " Mean";
	for(auto c = 0u; c < Nsize; c++) Nout << "," << data.democat[0].value[c] << " CI min";
	for(auto c = 0u; c < Nsize; c++) Nout << "," << data.democat[0].value[c] << " CI max";
	Nout << endl;

	auto N = data.genQ.N[0];
	
	for(auto j = 0u; j < Nsize; j++){
		vector <Statistics> stat(Nsize);
		for(auto i = 0u; i < Nsize; i++){
			vector <double> vec;
			for(const auto &opsa : opsamp){
				vec.push_back(opsa.Nsample[0][j][i]*data.agedist[j]/N.norm_factor);	
			}				
			stat[i] = get_statistic(vec);
		}
		
		Nout << data.democat[0].value[j];
		for(auto i = 0u; i < Nsize; i++) Nout << "," << stat[i].mean;
		for(auto i = 0u; i < Nsize; i++) Nout << "," << stat[i].CImin;
		for(auto i = 0u; i < Nsize; i++) Nout << "," << stat[i].CImax;
		Nout << endl;
	}
				 	
	if(details.siminf == INFERENCE){
		auto fulldesc2 = "Modified age mixing matrix: This shows the inferred relative contact rate with which different age groups mix within the population (darker red colours indicate more frequent contacts).";
		OutputPlot oppl2(OP_AGE_MATRIX_POST,file,fulldesc2,"Posterior","Age Mixing","Matrix","","","",UNSET,UNSET);
		op.push_back(oppl2);
		
		auto fulldesc3 = "Modified age mixing matrix: This shows the inferred relative contact rate with which different age groups mix within the population (darker red colours indicate more frequent contacts).";
		OutputPlot oppl3(OP_AGE_MATRIX_DIF,file,fulldesc3,"Posterior","Age Mixing","Difference","","","",UNSET,UNSET);
		op.push_back(oppl3);
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
	gnuplot << "set style line 23 lt 1 lc rgb '#ff2222' lw 0.25" << endl; // RED_THIN
	gnuplot << "set style line 24 lt 1 lc rgb '#aaffaa' lw 0.25" << endl; // GREEN_THIN
	gnuplot << "set style line 25 lt 1 lc rgb '#aaaaff' lw 0.25" << endl; // BLUE_THIN
	gnuplot << "set style line 26 lt 1 lc rgb '#222222' lw 0.25" << endl; // BLACK_THIN
	
	
	gnuplot << "set style line 80 lt 1 lc rgb '#222222' lw 2" << endl; // BLACK_THIN
	gnuplot << "set style line 81 lt 2 lc rgb '#ff2222' lw 2" << endl; // RED_THIN
	gnuplot << "set style line 82 lt 2 lc rgb '#22ff22' lw 2" << endl; // GREEN_THIN
	gnuplot << "set style line 83 lt 2 lc rgb '#2222ff' lw 2" << endl; // BLUE_THIN
	
	auto pred_timeplot = get_pred_timeplot();

	auto multiplot = UNSET;
	for(auto i = 0u; i < op.size(); i++){
		const auto &oppl = op[i];
			
		if(oppl.tab != "Data" && oppl.type != OP_ANIM_MAP && oppl.type != OP_AGE_MATRIX
			&& oppl.type != OP_AGE_MATRIX_POST && oppl.type != OP_AGE_MATRIX_DIF
			&& oppl.type != OP_COMP_MODEL && oppl.type != OP_FOI_MODEL && oppl.type != OP_DESC
			&& oppl.type != OP_MIXING_WITHIN_ANIM_MAP && oppl.type != OP_MIXING_WITHIN_MAP
			&& oppl.type != OP_MIXING_BETWEEN_ANIM_MAP && oppl.type != OP_MIXING_BETWEEN_MAP
			&& oppl.type != OP_MIXING_POP_ANIM_MAP && oppl.type != OP_MIXING_POP_MAP
			&& oppl.type != OP_AREA_COVAR && oppl.type != OP_TV_COVAR && oppl.type != OP_AREA_TV_COVAR
			&& oppl.type != OP_LEVEL_EFFECT){
			
			if(oppl.title == "") gnuplot << "set title" << endl;
			else gnuplot << "set title '" << label(oppl.title) << "'" << endl;
			gnuplot << "set autoscale" << endl;
			
			auto keyfontsize = 22u, labelfontsize = 22u, ticfontsize = 20u; 
			if(oppl.type == OP_TRACE){ keyfontsize = 16u; labelfontsize = 18; ticfontsize = 16u;}
			gnuplot << "set tics font ', " << ticfontsize << "'" << endl;
			gnuplot << "set key font ', " << keyfontsize << "'" << endl;
			gnuplot << "set xlabel '" << label(replace(oppl.xaxis,"age:age","")) << "' font '," << labelfontsize << "'" << endl;
			gnuplot << "set ylabel '" << label(replace(oppl.yaxis,"age:age","")) << "' font '," << labelfontsize << "'" << endl;
		
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
											<< ymax << " lw 3 lc rgb '#0000aa' lt 1 nohead" << endl;
					
							gnuplot << "set label  '" << label(tp.name) << "' at " << tp.time - shift << "," << ymax << " left rotate by 270 offset 0.7,-0.2 font ',16' textcolor rgb '#0000aa' " << endl;
						}
						
						for(const auto &tp : pred_timeplot){
							gnuplot << "set arrow from " << tp.time << ",0 to " << tp.time << ","
											<< ymax << " lw 3 lc rgb '#00aa00' ";

							if(tp.name != "") gnuplot << "lt 1"; else gnuplot << "lt 4"; 
							gnuplot << " nohead" << endl;
					
							if(tp.name != ""){
								gnuplot << "set label '" << label(tp.name) << "' at " << tp.time << "," << ymax << " tc rgb '#00aa00'  left rotate by 270 offset 1,-1 font ',16' " << endl;
							}
						}
						
						gnuplot << "plot ";
						for(auto j = 0u; j < oppl.line.size(); j++){
							const auto &line = oppl.line[j];
							if(j != 0) gnuplot << ", ";
							gnuplot << "'" << line.file << "' using " << line.xcol << ":" << line.ycol << " with lines ls " << line.style << " ";
							
							if(line.name == "" || line.name == "95% CI min" || line.name == "95% CI max" || oppl.line.size() == 1) gnuplot << "notitle"; 
							else gnuplot << "title '" << label(line.name) << "'";
						}
						gnuplot << endl;

						/*
						gnuplot << "plot ";
						for(auto j = 0u; j < oppl.line.size(); j++){
							const auto &line = oppl.line[j];
							if(j != 0) gnuplot << ", "; 
							
							if(j == oppl.line.size()-1){
								gnuplot << "'" << line.file << "' using " << line.xcol << ":" << line.ycol << " lc rgb '#000000' lw 3 ps 2 ";
							}
							else{
								gnuplot << "'" << line.file << "' using " << line.xcol << ":" << line.ycol << " with lines ls " << line.style << " ";
							}
							gnuplot << "notitle"; 
						}
						*/
						
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
						if(oppl.line.size() == 0) emsg("Marginal graph '"+oppl.title+"' does not contain any lines");
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
						
						if(line.name == "Value"){
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
						
						if(line.name == "Value"){
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
						if(line.name == "Value"){
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
					gnuplot << "plot ";
					for(auto j = 0u; j < oppl.line.size(); j++){
						const auto &line = oppl.line[j];
						if(j != 0) gnuplot << ", ";
						if(line.name == ""){  // Ordinary line
							gnuplot << "'" << line.file << "' using " << line.xcol << ":" << line.ycol;
							gnuplot << " with lines ls " << line.style << " ";
						
						}
						else{                 // Error bar
							gnuplot << "'" << line.file << "' using " << line.xcol << ":" << line.ycol << ":8";
							gnuplot << " with yerrorbars ls " << line.style << " ";
						}
						gnuplot << "notitle";
					}
					gnuplot << endl;
					break;
					
				case OP_ANIM_MAP: case OP_AGE_MATRIX: case OP_AGE_MATRIX_POST: case OP_AGE_MATRIX_DIF:
				case OP_PARAM_TABLE: case OP_PRIOR_TABLE: 
				case OP_COMP_MODEL: case OP_FOI_MODEL: case OP_DESC:
				case OP_MIXING_WITHIN_ANIM_MAP: case OP_MIXING_WITHIN_MAP:
				case OP_MIXING_BETWEEN_ANIM_MAP: case OP_MIXING_BETWEEN_MAP: 
				case OP_MIXING_POP_ANIM_MAP: case OP_MIXING_POP_MAP:
				case OP_AREA_COVAR: case OP_TV_COVAR: case OP_AREA_TV_COVAR: case OP_LEVEL_EFFECT:
					break;
			}
			gnuplot << endl << endl;
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
			dataout << fixed;
			
			if(dt.type != POPFRAC && details.stochastic == true) dataout << std::fixed << std::setprecision(0);
	
			if(!dataout) emsg("Cannot open the file '"+filefull+"'");
			
			switch(dt.type){
				case POP: case POPFRAC: case TRANS:
					{
						dataout << details.time_format_str;
						for(auto i : dt.graph_ref){
							dataout << sep << data.graph[i].colname;
							if(dt.load_sd == true) dataout << sep << data.graph[i].colname << " SD";
						}
						dataout << endl;

						if(dt.graph_ref.size() == 0) emsg("The file '"+file+"' contains no graphs");
						
						auto nrow = data.graph[dt.graph_ref[0]].point.size();
						for(auto row =0u; row < nrow; row++){
							dataout << details.getdate(dt.start + row*dt.timestep);
							for(auto i : dt.graph_ref){
								dataout << sep;
								
								auto ob = data.graph[i].point[row].obs;
								auto val = obs_value[ob];
								if(dt.threshold != UNSET && val <= dt.threshold){
									dataout << data.threshold_str;
								}
								else{
									if(val < 0) val = 0;
									dataout	<< val;
								}
								
								if(dt.load_sd == true){
									dataout << sep << data.obs[ob].sd;
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
	auto var = sum2/sumw - (sum/sumw)*(sum/sumw); if(var < TINY) var = 0;
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
	auto var = sum2/n - (sum/n)*(sum/n); if(var < VTINY) var = 0;
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


/// Calculates diagnostic statistics using 75% credible interval
Statistics Output::get_statistic_75_percent(const vector <double> &vec) const                       
{
	Statistics stat;
	
	auto n = vec.size();

	auto sum = 0.0, sum2 = 0.0; 	
	for(auto v : vec){ sum += v; sum2 += v*v;}
	
	stat.mean = sum/n; 	
	auto var = sum2/n - (sum/n)*(sum/n); if(var < VTINY) var = 0;
	stat.sd = sqrt(var);
	
	if(n == 0){
		stat.CImin = UNSET; stat.CImax = UNSET; 
	}
	else{ 
		auto vec2 = vec;
		sort(vec2.begin(),vec2.end());

		if(n >= 2){
			auto i = (unsigned int)((n-1)*0.25); auto f = (n-1)*0.25 - i;
			stat.CImin = vec2[i]*(1-f) + vec2[i+1]*f;
				
			i = (unsigned int)((n-1)*0.75); f = (n-1)*0.75 - i;
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
				for(auto d = 1u; d < nsamp/2; d++){
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


/// Returns the correntation between two sets of samples		
vector <double> Output::get_correlation(const vector <ParamSample> &before, const vector <ParamSample> &after) const
{
	auto nparam = model.param.size();
	vector <double> cor(nparam);
	
	if(mpi.core == 0){
		auto N = before.size();
		
		if(N != after.size()) emsgEC("Output",13);
		for(auto i = 0u; i < N; i++){
			if(before[i].run != after[i].run) emsgEC("Output",14);
		}
		
		for(auto th = 0u; th < nparam; th++){
			if(model.param[th].priortype == FIXED_PRIOR){
				cor[th] = UNSET;
			}
			else{
				auto av_bef = 0.0, av_bef2 = 0.0, av_aft = 0.0, av_aft2 = 0.0, av_bef_aft = 0.0; 
				for(auto i = 0u; i < N; i++){
					auto val_bef = before[i].paramval[th];
					auto val_aft = after[i].paramval[th];
					av_bef += val_bef; av_bef2 += val_bef*val_bef;
					av_aft += val_aft; av_aft2 += val_aft*val_aft;
					av_bef_aft += val_bef*val_aft;
				}
				av_bef /= N; av_bef2 /= N; av_aft /= N; av_aft2 /= N; av_bef_aft /= N;
				
				cor[th] = (av_bef_aft - av_bef*av_aft + TINY)/(sqrt((av_bef2 - av_bef*av_bef + TINY)*(av_aft2 - av_aft*av_aft + TINY)));
			}
		}
	}
	
	mpi.bcast(cor);
	
	return cor;
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
	//vector <double> GR;

	vector < vector < vector <double> > > param_GR;
	param_GR.resize(nrun);
	for(const auto &ps : psamp){
		auto ru = ps.run; if(ru == UNSET) ru = 0;
		param_GR[ru].push_back(ps.paramval);
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


/// Sets the time, in minutes, a generation is finished by
void Output::set_generation_time(Generation &gen) const 
{
	timer[TIME_ALG].stop();
	auto alg_time = mpi.sum(timer[TIME_ALG].val);
	auto wait_time = mpi.sum(timer[TIME_WAIT].val);
	gen.time = (alg_time - wait_time)/(60.0*CLOCKS_PER_SEC);
	timer[TIME_ALG].start();
}
	
	
/// Outputs how the error function and posterior distributions change with generation number (ABC-SMC and ABC-MBP)
void Output::generation_plot(const string file, const vector <Generation> &generation) const
{
	auto nparam = model.param.size(); 

	string filefull = details.output_directory+"/"+file;
	ofstream genout(filefull.c_str());
	if(!genout) emsg("Cannot open the file '"+filefull+"'");
	genout << fixed;
	
	genout << "# This shows how the EF and infered posterior distribiutions vary with generation number" << endl;
	genout << "Generation,CPU Time,EF cut-off,log(EF cut-off),InvT,log(InvT),Model evidence,Model evidence SD"; 
	for(const auto &param : model.param) genout << "," << param.name << " (mean)," << param.name << " (95% CI min)," << param.name << " (95% CI max)";
	genout << endl;
	
	auto gmin = 1u;
	if(details.mode == ML_INF) gmin = 0;
	
	auto step = int((generation.size()-gmin)/100);
	if(step == 0) step = 1;
	
	for(auto g = gmin; g < generation.size(); g += step){
		const Generation &gen = generation[g];
		
		genout << g << "," << gen.time << ",";

		switch(details.mode){
			case PAIS_INF:
				genout << gen.EFmax << "," << log(gen.EFmax) << "," << gen.invT << "," << log(gen.invT);
				break;
			default:
				auto val = gen.EFcut; if(val < 0.001) val = 0.001;
				genout << gen.EFcut << "," << log(val) << ",UNSET,UNSET";
				break;
		}
		auto stat = get_statistic(gen.model_evidence);
	
		genout << "," << stat.mean << "," << stat.sd;
	
		auto psamp = gen.param_samp;
		auto npsamp = psamp.size();

		for(auto &psa : psamp) psa.paramval = model.dirichlet_correct(psa.paramval);
		
		for(auto th = 0u; th < nparam; th++){
			vector <WeightedPoint> vec;
			for(auto i = 0u; i < npsamp; i++){
				WeightedPoint pw;
				pw.val = psamp[i].paramval[th];
				switch(details.mode){
					case ABC_SMC: pw.w = gen.w[i]; break;
					case ABC_MBP: case ABC_DA: case ABC_CONT: case PAIS_INF: case MC3_INF: case MCMC_MBP:
					case ML_INF:
						pw.w = 1;
						break;
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

	ldtyout << "# This shows the log of the contribution to the error function as a function of generation" << endl;
	
	ldtyout << "Generation";
	for(auto &dt : data.datatable){
		if(dt.file != "") ldtyout << "," << dt.file;
	}
	
	ldtyout << endl;
		 
	for(auto g = 1u; g < generation.size(); g++){
		const Generation &gen = generation[g];
		
		ldtyout << g;
		for(auto dt = 0u; dt < data.datatable.size(); dt++){
			if(data.datatable[dt].file != ""){
				auto val = 0.0; for(const auto &sa : gen.EF_datatable) val += sa[dt];
				val /= gen.EF_datatable.size();
			
				ldtyout << "," << log(val);
			}
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
	
	mpi.barrier();
	
	timer[TIME_RESULTS].stop();
}


/// Creates all the directories needed for outputting
void Output::ensure_directories()
{
	post_dir = "Posterior";
	if(details.mode == SIM) post_dir = "Simulation_Results";
	if(details.mode == MULTISIM) post_dir = "Multisim_Results";
		
	ensure_directory(details.output_directory);
	if(details.siminf != SIMULATE) ensure_directory(details.output_directory+"/Diagnostics");
	if(details.mode == MC3_INF) ensure_directory(details.output_directory+"/Diagnostics/Other Chains");
	auto dir = details.output_directory+"/"+post_dir;
	ensure_directory(dir);
	ensure_directory(dir+"/parameter");
	ensure_directory(dir+"/state");
	ensure_directory(dir+"/susceptibility");
	ensure_directory(dir+"/spline");
	ensure_directory(dir+"/covar");
	
	if(data.democat_change.size() > 0) ensure_directory(dir+"/democat_change");
	
	if(details.siminf == INFERENCE) ensure_directory(dir+"/samples");
	
	if(data.nage > 1) ensure_directory(dir+"/age_mixing");
			
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
	trace << fixed;
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
	read << "This directory contains the outputs from the BEEP analysis."<< endl;
	read << "We give a brief description below of the various files and folders contained within it:" << endl << endl;
			
	read << "## Parameter_estimates.csv" << endl << endl;
	
	read << "If inference is performed, this file gives the posterior means and 95% credible intervals for each of the model parameters." << endl;
	
	read << endl;
	
	read << "## Graphs.pdf" << endl << endl;
	
	read << "This visualises outputs from BEEP." << endl << endl;
	
	read << "## Graphs_description.txt" << endl << endl;
	
	read << "Provides information about the source data for the graphs." << endl << endl;
			
	read << "## Model_specification.txt" << endl << endl;
	
	read << "This gives a summary of the model used to perform the analysis." << endl;
	
	read << endl;
	
	read << "## " << details.analysis_type << endl << endl;
	
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
void Output::generate_graphs(vector <Particle> &particle_store, const double invT) const 
{
	timer[TIME_RESULTS].start();
				
	State state(details,data,model,obsmodel);

	auto dir = details.output_directory+"/"+post_dir+"/samples/";
	vector <ParamSample> psamp;
	vector <Sample> opsamp;
	
	mpi.barrier(); if(mpi.core == 0 && mpi.ncore > 1) cout << "Gathering samples..." << endl << flush;
		
	mpi.gather_samples(psamp,opsamp,particle_store,state,dir);
	mpi.barrier();
	
	if(mpi.core == 0) generate_graphs(psamp,opsamp,invT);
	
	mpi.barrier();
	
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
	lab = replace(lab,"→"," → ");
	return lab;
}


/// Generates a new output plot
OutputPlot::OutputPlot(OutputPlotType type_, string title_, string fulldesc_, string tab_, string tab2_, string tab3_, string tab4_, string xaxis_, string yaxis_, double min_, double max_)
{
	type = type_; title = title_; fulldesc = fulldesc_; tab = tab_; tab2 = tab2_; tab3 = tab3_; tab4 = tab4_;
	xaxis = xaxis_; yaxis = yaxis_; min = min_; max = max_; legend_labelson = false;
	xmin = UNSET; xmax = UNSET;    
}


/// Adds a line to an output plot
void OutputPlot::addline(const string name, const string file, const unsigned int xcol, const unsigned int ycol, const LineType style, const unsigned int EB)
{
	OutputLine opl; opl.name = name; opl.file = file; opl.xcol = xcol; opl.ycol = ycol; opl.style = style; opl.EB = EB;
	
	line.push_back(opl);
};
