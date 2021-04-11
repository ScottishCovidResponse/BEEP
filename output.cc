// Outputs various graphs and statistics

#include <iostream>
#include <vector>
#include <fstream>
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

Output::Output(const Details &details, const Data &data, const Model &model, const Inputs &inputs, const ObservationModel &obsmodel, Mpi &mpi) :  details(details), data(data), model(model), obsmodel(obsmodel), mpi(mpi)
{
	if(mpi.core == 0){
		ensure_directories();                                                       // Creates directories for outputting
		
		print_model_specification();                                                // Summarises the data and model
		
		inputs.find_outputprop(prop);                                               // Load output properties
		
		readme(); 	                                                              	// Generates a readme file for outputs
	}
}
	

/// Generates the output files along with the final pdf report that gives graphical outputs
void Output::generate_graphs(vector <ParamSample> &psamp, const vector <Sample> &opsamp) const
{ 
	dirichlet_normalisation(psamp);                                               // Dirichlet distributions are normalised

	cout << endl;
	if(suppress_output == false){
		if(details.siminf == SIMULATE) cout << "Outputs in directory '" << details.output_directory << "':" << endl;
		else cout << "Posterior outputs in directory '" << details.output_directory << "':" << endl;
	}
	
	vector <OutputPlot> op;                                                       // Stores the plots for the final pdf
	
	spline_plots(opsamp,op);                                                      // Generates all the different types of plot
	
	graph_plots(opsamp,op);
	 
	susceptibility_distributions(psamp,op);
	
	branch_prob_distributions(psamp,op);

	derived_parameter_distributions(opsamp,op);

	posterior_parameter_distributions(psamp,op);

	add_generation_plots(op);
	
	//spatial_R_map(opsamp);
	
	posterior_parameter_estimates(psamp);
	
	generate_gnuplot_file(op);                                                    // Generates the gnuplot file
	
	generate_pdf();                                                               // Runs gnuplot and generates pdf 
	
	generate_pdf_description(op);                                                 // Generates a text descrition of the pdf file
}


/// Generate the final pdf report
void Output::generate_pdf() const
{
	auto opdir = details.output_directory;
	stringstream ss; ss << "gnuplot '" << opdir << "/gnuplot.txt'" << endl;
	system(ss.str().c_str());
		
	stringstream ss2; ss2 << "ps2pdf " << opdir << "/Graphs.ps " << opdir << "/Graphs.pdf" << endl;
	system(ss2.str().c_str());

	stringstream ss3; ss3 << "rm -f " << opdir << "/Graphs.ps" << endl; 
	system(ss3.str().c_str());
		
	stringstream ss4; ss4 << "rm -f " << opdir << "/gnuplot.txt" << endl;
	system(ss4.str().c_str());

	stringstream ss5; ss5 << "rm -f " << opdir << "/line.txt" << endl;
	system(ss5.str().c_str());

	cout << "Output graphs are placed into '" << opdir << "/Graphs.pdf'" << endl << endl;
}


/// Generates a text description of the pdf file
void Output::generate_pdf_description(const vector <OutputPlot> &op) const
{
	auto file = details.output_directory+"/Graphs_description.txt";
	ofstream desc(file);
	if(!desc) emsg("Cannot output the file '"+file+"'");
	
	desc << "This file provides a description of the source files used to generate the graphs in '" << details.output_directory << "/Graphs.pdf'" << endl << endl;
	
	for(auto i = 0u; i < op.size(); i++){
		const auto &oppl = op[i];
		desc << "PAGE " << i+1 << ": \"" << oppl.title << "\"" << endl;
		
		for(auto j = 0u; j < oppl.line.size(); j++){
			const auto &line = oppl.line[j];
			desc << "  '" << line.name << "' ";
			switch(oppl.type){
			case OP_SPLINE: case OP_GRAPH: case OP_PARAMDIST: case OP_GENERATION: case OP_LOG_GENERATION: case OP_CPU:
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
		vector <double> vec;
		for(auto par : model.param) vec.push_back(par.value);			
		sim_spline_output = model.get_spline_output(vec);
	}
		
	auto nspline_out = opsamp[0].spline_output.size();
	for(auto sp = 0u; sp < nspline_out; sp++){
		const auto &splout = opsamp[0].spline_output[sp];
		auto name = splout.name;
		auto file = "Posterior/spline/"+name+".txt";
		auto filefull = details.output_directory+"/"+file;
		ofstream tout(filefull.c_str());
		if(!tout) emsg("Cannot output the file '"+filefull+"'");
	
		if(suppress_output == false) cout << "'" << file << "' gives the time variation in "+name+"." << endl;
	
		tout << "# Gives the time variation in "+name+"." << endl;	
		tout << "# Time from start, mean, minimum of 95% credible interval, maximum of 95% credible interval, estimated sample size" << endl;

		auto min = LARGE, max = -LARGE;
		for(auto st = 0u; st < details.ndivision; st++){
			vector <double> vec; 
			for(const auto &opsa : opsamp) vec.push_back(opsa.spline_output[sp].splineval[st]);
			auto stat = get_statistic(vec);
			
			tout << (st+0.5)*details.period/details.ndivision << " " << stat.mean << " " << stat.CImin << " " << stat.CImax;
			
			double num = atof(stat.CImax.c_str()); if(num > max) max = num; if(num < min) min = num;
			
			if(plot_sim_param == true) tout << " " << sim_spline_output[sp].splineval[st];
			tout << endl; 
		}
		
		if(min != max){   // Generates the final plot
			OutputPlot oppl(OP_SPLINE,splout.desc,"Time",name,min,max);
			string lab;
			if(details.mode != SIM){
				oppl.addline("Posterior",filefull,1,2,RED_SOLID);
				oppl.addline("95% CI minimum",filefull,1,3,RED_DASHED);
				oppl.addline("95% CI maximum",filefull,1,4,RED_DASHED);
			  lab = "True";
			}
			
			if(plot_sim_param == true) oppl.addline(lab,filefull,1,5,BLACK_DASHED);
			op.push_back(oppl);
		}
	}
}


/// Plots a graph giving time variation in transition rate or population or marginal distributions
void Output::graph_plots(const vector <Sample> &opsamp, vector <OutputPlot> &op) const
{
	if(opsamp.size() == 0) return;
	
	for(auto gr_num = 0u; gr_num < data.graph.size(); gr_num++){
		auto gr = data.graph[gr_num];
		
		auto file = details.output_directory+"/Posterior/state/"+gr.filename+".txt";

		ofstream stateout(file.c_str());
		if(!stateout) emsg("Could not open file '"+file+"'");
		
		auto imax = opsamp[0].graph_state[gr_num].size();
		auto min = LARGE, max = -LARGE;
		for(auto i = 0u; i < imax; i++){
			vector <double> vec; for(const auto &opsa : opsamp) vec.push_back(opsa.graph_state[gr_num][i]);		
			auto stat = get_statistic(vec);
			
			switch(gr.type){
				case GRAPH_TIMESERIES:
					stateout << details.division_time[i*graph_step+graph_step/2];
					break;
					
				case GRAPH_MARGINAL:
					stateout << i;
					break;
			}
			
			stateout << " " << stat.mean << " " << stat.CImin << " " << stat.CImax << endl; 
			
			auto num = atof(stat.CImax.c_str()); if(num > max) max = num;
			num = atof(stat.CImin.c_str()); if(num < min) min = num;
		}
		
		string file_data;
		if(data.datatable[gr.datatable].optype == DATA){
			file_data = details.output_directory+"/Posterior/state/"+gr.filename+"-data.txt";

			ofstream dataout(file_data.c_str());
			if(!dataout) emsg("Could not open file '"+file_data+"'");
			
			if(suppress_output == false) cout << file_data << " file" << endl;
			
			for(auto i = 0u; i < gr.point.size(); i++){
				const auto& p = gr.point[i];
				
				auto sum_true = 0.0;
				for(auto ob : p.obs){
					auto val = data.obs[ob].value;
					if(val != UNKNOWN && val != THRESH) sum_true += val;
				}	
			
				if(sum_true >= UNSET) sum_true = 0;
				
				if(gr.type == GRAPH_MARGINAL){
					const auto dt = data.obs[p.obs[0]].datatable;
					dataout << p.xi << " " << data.datatable[dt].demolist[i] << " " << sum_true << endl;
				}
				else{
					sum_true /= p.xf - p.xi;
					dataout << p.xi << " " << sum_true << endl;
					dataout << p.xf << " " << sum_true << endl;
				}
				
				if(sum_true > max) max = sum_true; if(sum_true < min) min = sum_true;
			}
		}

		if(min != max){   // Generates the final plot
			switch(gr.type){
				case GRAPH_TIMESERIES:
					{
						OutputPlot oppl(OP_GRAPH,gr.desc,"Time",gr.name,min,max);
						if(details.mode != SIM){
							oppl.addline("Posterior",file,1,2,RED_SOLID);
							oppl.addline("95% CI minimum",file,1,3,RED_DASHED);
							oppl.addline("95% CI maximum",file,1,4,RED_DASHED);
						}
						if(file_data != "") oppl.addline("Data",file_data,1,2,BLACK_SOLID);
						op.push_back(oppl);
					}
					break;
					
				case GRAPH_MARGINAL:
					{
						OutputPlot oppl(OP_GRAPH_MARGINAL,gr.desc,"Category",gr.name,min,max);
						if(details.mode != SIM){
							oppl.addline("Posterior",file,1,3,RED_SOLID);
						}
						if(file_data != "") oppl.addline("Data",file_data,1,3,BLACK_SOLID);
						op.push_back(oppl);
					}
					break;
			}
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
			
			auto file = type+".txt";
			auto filefull = details.output_directory+"/Posterior/susceptibility/"+file;
			ofstream distout(filefull.c_str());
			if(!distout) emsg("Cannot output the file '"+filefull+"'");
			
			if(suppress_output == false) cout << "'" << file << "' gives stratified susceptibility for " << type << "." << endl;
			
			distout << "# Age-stratified susceptibility for " << type << "." << endl;
		
			distout << "#" << type << "\tMean\tCI min\tCImax";
			if(plot_sim_param == true) distout << "\tSimulated";
			distout << endl;					
			for(auto a = 0u; a < amax; a++){
				distout << a << "\t" << data.democat[c].value[a];
			
				auto th = model.sus_param_list[c][model.sus_ref[c][a]]; 
				vector <double> vec; for(const auto &psa : psamp) vec.push_back(psa.paramval[th]);
				auto stat = get_statistic(vec);
				
				distout << "\t" << stat.mean << "\t" << stat.CImin << "\t" << stat.CImax;
				if(plot_sim_param == true) distout << "\t" << model.param[th].value;
				distout << endl;
			}
			
			OutputPlot oppl(OP_MARGINAL,"Susceptibility with "+type,"Category","Relative susceptibility",0,0);
			if(details.mode != SIM){
				oppl.addline("Posterior",filefull,1,3,RED_SOLID);
			}
			if(plot_sim_param == true) oppl.addline("True",filefull,1,6,BLACK_SOLID);
			op.push_back(oppl);
		}
	}
}


/// Outputs a bar chart giving branching probabilities split by demographic groups
void Output::branch_prob_distributions(const vector <ParamSample> &psamp, vector <OutputPlot> &op) const 
{
	auto file = "Posterior/branching_probability.txt";
	auto filefull = details.output_directory+"/"+file;
	ofstream distout(filefull.c_str());
	if(!distout) emsg("Cannot output the file '"+filefull+"'");
	
	if(suppress_output == false) cout << "'" << file << "' gives stratified branching probability." << endl;
	
  distout << "# Stratified branching probability." << endl;

	distout << "# Demographic state";

	vector <unsigned int> tr_list;
	vector <double> min, max;
	for(auto c = 0u; c < model.comp.size(); c++){
		auto kmax = model.comp[c].trans.size();
		if(kmax > 1){
			for(auto k = 0u; k < kmax; k++){
				tr_list.push_back(model.comp[c].trans[k]);
				min.push_back(LARGE); max.push_back(-LARGE); 
			}
		}
	}
	
	for(auto tr : tr_list){
		distout << "\t" << model.trans_str(tr) << "\tMin\tMax";
		if(plot_sim_param == true) distout << "\tSimulated";
	}
	distout << endl;
	
	for(auto dp = 0u; dp < data.ndemocatpos; dp++){
		distout << dp << "\t";
		for(auto c = 0u; c < data.ndemocat; c++){
			if(c != 0) distout << ",";
			distout << data.democat[c].value[data.democatpos[dp][c]];
		}
			
		for(auto i = 0u; i < tr_list.size(); i++){
			auto tr = tr_list[i];
			auto th = model.trans[tr].probparam[dp];

			vector <double> vec; for(const auto &psa : psamp) vec.push_back(psa.paramval[th]);
			auto stat = get_statistic(vec);

			distout << "\t" << stat.mean << "\t" << stat.CImin << "\t" << stat.CImax;	
			auto num = atof(stat.CImax.c_str()); if(num > max[i]) max[i] = num;
			num = atof(stat.CImin.c_str()); if(num < min[i]) min[i] = num;
			
			if(plot_sim_param == true){
				num = model.param[th].value; if(num > max[i]) max[i] = num; if(num < min[i]) min[i] = num;
				distout << "\t" << num;
			}
		}
		distout << endl;
	}
	
	auto col = 3u;
	for(auto i = 0u; i < tr_list.size(); i++){
		if(min[i] != max[i]){
			auto tr = tr_list[i];
			OutputPlot oppl(OP_MARGINAL,model.trans_str(tr)+" branching probability","Demographic classification","Probability",0,0);
			if(details.mode != SIM) oppl.addline("Posterior",filefull,1,col,RED_SOLID);
			col += 3;
			if(plot_sim_param == true){ oppl.addline("True",filefull,1,col,BLACK_SOLID); col++;}
			op.push_back(oppl);	
		}
	}
	distout << endl;
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
			
			string file = derpar.filename+".txt";
			string filefull = details.output_directory+"/Posterior/parameter/"+file;
			ofstream distout(filefull.c_str());
			if(!distout) emsg("Cannot output the file '"+filefull+"'");
			
			if(suppress_output == false){
				cout << "'" << file << "' gives the probability distributions for the "+ derpar.desc+"." << endl;
			}
			
			distout << "# Posterior probability distributions for the generation time." << endl;

			
			distout << "# Value\t" << "Probability" << endl;
			
			for(auto b = 0u; b < paramdist.value.size(); b++){
				distout << paramdist.value[b] << "\t" << paramdist.prob[b] << endl;
			}
			
			OutputPlot oppl(OP_PARAMDIST,derpar.desc,derpar.name,"Probability",0,0);
			oppl.addline("Posterior",filefull,1,2,GREEN_SOLID);
			op.push_back(oppl);
		}
	}
}


/// Outputs a file containing posterior distributions for model parameters
void Output::posterior_parameter_distributions(const vector <ParamSample> &psamp, vector <OutputPlot> &op) const
{
	auto file = "model.txt";
	auto filefull = details.output_directory+"/Posterior/parameter/"+file;
	ofstream distout(filefull.c_str());
	if(!distout) emsg("Cannot output the file '"+filefull+"'");
	
	if(suppress_output == false) cout << "'" << file << "' gives the probability distributions for parameters." << endl;
	
  distout << "# Posterior probability distributions for model parameters." << endl;
	distout << endl;

	vector <Distribution> paramdist(model.param.size());
	for(auto p = 0u; p < model.param.size(); p++){
		vector <double> vec; for(const auto &psa : psamp) vec.push_back(psa.paramval[p]);
		auto priormin = UNSET, priormax = UNSET;
		if(model.param[p].priortype == UNIFORM_PRIOR){ priormin = model.param[p].val1; priormax = model.param[p].val2;}
		paramdist[p] = get_distribution(vec,priormin,priormax);
	}
	
	for(auto p = 0u; p < model.param.size(); p++){
		distout << model.param[p].name << "\t" << "Probability\t";
	}	
	distout << endl;
	
	if(paramdist.size() > 0){
		for(auto b = 0u; b < paramdist[0].value.size(); b++){
			for(auto p = 0u; p < model.param.size(); p++){
				distout << paramdist[p].value[b] << "\t" << paramdist[p].prob[b] << "\t";
			}
			distout << endl;
		}
	}
	
	for(auto p = 0u; p < model.param.size(); p++){
		if(paramdist[p].variation == true){
			OutputPlot oppl(OP_PARAMDIST,"Posertior distribution for "+model.param[p].name,model.param[p].name,"Probability",0,0);
			oppl.addline("Posterior",filefull,2*p+1,2*p+2,GREEN_SOLID);
			if(plot_sim_param == true) oppl.addline("True",to_string(model.param[p].value),1,2,BLACK_DASHED);
			op.push_back(oppl);
		}
	}
}


/// Adds the plots which show how quantities change as a function of generation number (ABC-SMC / ABC-MBP)
void Output::add_generation_plots(vector <OutputPlot> &op) const
{
	if(details.mode != ABC_SMC && details.mode != ABC_MBP && details.mode != ABC_MBP_GR && details.mode != PAIS_INF) return;
	
	auto file = details.output_directory+"/Diagnostics/Generation.txt";
	for(auto p = 0u; p < model.param.size(); p++){                           // These are graphs for the parameters
		if(model.param[p].priortype != FIXED_PRIOR){
			auto min = 0.0, max = 0.0; 
			switch(model.param[p].priortype){
				case UNIFORM_PRIOR: min = model.param[p].val1; max = model.param[p].val2; break;
				case EXP_PRIOR: min = UNSET; max = UNSET; break;
				default: emsgEC("Output",66); break;
			}
			
			OutputPlot oppl(OP_GENERATION,"Posterior estimate for "+model.param[p].name+" as a function of generation","Generation",model.param[p].name,min,max);
			oppl.addline("Posterior",file,1,3*p+4,BLUE_SOLID);
			oppl.addline("95% CI minimum",file,1,3*p+5,BLUE_DASHED);
			oppl.addline("95% CI maximum",file,1,3*p+6,BLUE_DASHED);
			if(plot_sim_param == true) oppl.addline("True",to_string(model.param[p].value),1,2,BLACK_DASHED);
			op.push_back(oppl);
		}
	}
	
	auto fileDT = details.output_directory+"/Diagnostics/EF_datatable.txt";
	OutputPlot opplDT(OP_LOG_GENERATION,"Contribtion to EF from different data sources","Generation","EF contribution",UNSET,UNSET);
	for(auto i = 0u; i < data.datatable.size(); i++){
		if(data.datatable[i].weight != 0) opplDT.addline(rus(data.datatable[i].file),fileDT,1,i+2,LineType(1+(i+1)%8));
	}
	op.push_back(opplDT);
	
	{
		OutputPlot oppl(OP_LOG_GENERATION,"EF cut-off as a function of generation","Generation","EF cut-off",UNSET,UNSET);
		oppl.addline("EF cut-off",file,1,3,BLACK_SOLID);
		op.push_back(oppl);
	}
	
	{
		OutputPlot oppl(OP_CPU,"CPU time as a function of EF","EF cut-off","CPU time (minutes)",UNSET,UNSET);
		oppl.addline("",file,3,2,BLACK_SOLID);
		op.push_back(oppl);
	}
}


/// Generates a file for plotting output with gnuplot
void Output::generate_gnuplot_file(const vector <OutputPlot> &op) const
{
	auto opdir = details.output_directory; 
	
	auto linefile = opdir+"/line.txt";
	ofstream line(linefile);       // This provides a way of drawing vertical lines
	if(!line) emsg("Cannot output the file '"+linefile+"'");
	line << "1 0" << endl << "1 1" << endl;
	
	auto file = opdir+"/gnuplot.txt";
	ofstream gnuplot(file);
	if(!gnuplot) emsg("Cannot output the file '"+file+"'");

	gnuplot << "set terminal postscript enhanced color font ',20' dl 3" << endl;
	gnuplot << "set output '" << opdir << "/Graphs.ps'" << endl;

	gnuplot << "set style line 1 lt 1 lc rgb '#ff2222' lw 5" << endl; // RED_SOLID
	gnuplot << "set style line 2 lt 3 lc rgb '#ffaaaa' lw 5" << endl; // RED_DASHED
	gnuplot << "set style line 3 lt 1 lc rgb '#22ff22' lw 5" << endl; // GREEN_SOLID
	gnuplot << "set style line 4 lt 3 lc rgb '#aaffaa' lw 5" << endl; // GREEN_DASHED
	gnuplot << "set style line 5 lt 1 lc rgb '#2222ff' lw 5" << endl; // BLUE_SOLID
	gnuplot << "set style line 6 lt 3 lc rgb '#aaaaff' lw 5" << endl; // BLUE_DASHED
	gnuplot << "set style line 7 lt 1 lc rgb '#222222' lw 5" << endl; // BLACK_SOLID
	gnuplot << "set style line 8 lt 3 lc rgb '#222222' lw 5" << endl; // BLACK_DASHED
	
	gnuplot << "set style line 20 lt 1 lc rgb '#222222' lw 2" << endl; // BLACK_THIN
	gnuplot << "set style line 21 lt 2 lc rgb '#ff2222' lw 2" << endl; // RED_THIN
	gnuplot << "set style line 22 lt 2 lc rgb '#22ff22' lw 2" << endl; // GREEN_THIN
	gnuplot << "set style line 23 lt 2 lc rgb '#2222ff' lw 2" << endl; // BLUE_THIN
	
	for(auto i = 0u; i < op.size(); i++){
		const auto &oppl = op[i];
		
		gnuplot << "set title '" << oppl.title << "'" << endl;
		gnuplot << "set autoscale" << endl;
		gnuplot << "set xlabel '" << oppl.xaxis << "' font ',22'" << endl;
		gnuplot << "set ylabel '" << oppl.yaxis << "' font ',22'" << endl;
	
		if(oppl.yaxis == "Probability")	gnuplot << "unset ytics" << endl;
		else gnuplot << "set ytics" << endl;
			
		switch(oppl.type){
			case OP_SPLINE: case OP_GRAPH:
				{
					auto ymax = 1.2*oppl.max;
					gnuplot << "set yrange [0:" << ymax << "]" << endl;
		
					for(auto j = 0u; j < details.timeplot.size(); j++){
						const auto &tp = details.timeplot[j];
						gnuplot << "set arrow from " << tp.time - details.start << ",0 to " << tp.time - details.start << ","
										<< ymax << " lw 2 lc rgb '#000066' nohead" << endl;
				
						gnuplot << "set label '" << tp.name << "' at " << tp.time - details.start << "," << ymax << " left rotate by 270 offset 1,-1" << endl;
					}
					
					gnuplot << "plot ";
					for(auto j = 0u; j < oppl.line.size(); j++){
						const auto &line = oppl.line[j];
						if(j != 0) gnuplot << ", ";
						gnuplot << "'" << line.file << "' using " << line.xcol << ":" << line.ycol << " with lines ls " << line.style << " ";
						if(line.name == "" || line.name == "95% CI minimum" || line.name == "95% CI maximum") gnuplot << "notitle"; 
						else gnuplot << "title '" << line.name << "'";
					}
					
					gnuplot << endl;
					gnuplot << "unset label"<< endl;
					gnuplot << "unset arrow"<< endl;
							
					if(false){
						for(auto j = 0u; j < details.timeplot.size(); j++){
							const auto &tp = details.timeplot[j];
							gnuplot << ", '" << opdir << "/line.txt' using ($1*" 
										<< tp.time - details.start << "):($2*"
										<< ymax << ") with lines ls " << 20+j%4 << " title '" << tp.name << "'";
						}	
					}
				}
				break;
				
			case OP_GRAPH_MARGINAL:	
				gnuplot << "set autoscale" << endl;
				gnuplot << "set yrange [0:]" << endl;
				gnuplot << "set boxwidth 0.5" << endl;
				gnuplot << "set style fill solid" << endl;
				gnuplot << "plot '" <<  oppl.line[1].file  << "' using ($1+0.5):3:xtic(2) with boxes notitle, ";
				gnuplot << "'" << oppl.line[0].file << "' using ($1+0.5):2:3:4 with errorbars lt 1  lc rgb '#000000' lw 3 ps 1 notitle";
				break;
				
			case OP_MARGINAL:	
				gnuplot << "set autoscale" << endl;
				gnuplot << "set yrange [0:]" << endl;
				gnuplot << "set boxwidth 0.5" << endl;
				gnuplot << "set style fill solid" << endl;
				{
					auto col = oppl.line[0].ycol;
					gnuplot << "plot '" << oppl.line[0].file << "' using ($1+0.5):" << col << ":xtic(2) with boxes notitle, ";
					gnuplot << "'" << oppl.line[0].file << "' using ($1+0.5):" << col << ":" << col+1 << ":" << col+2;
					gnuplot << " with errorbars lt 1  lc rgb '#000000' lw 3 ps 1 notitle";
					if(oppl.line.size() > 1){
						gnuplot << ",'" << oppl.line[1].file << "' using ($1+0.5):" << oppl.line[1].ycol 
											<< " with point lt 1 lc rgb '#0000ff' lw 3 pt 2 ps 1 notitle";
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
					
					if(line.name == "") gnuplot << " notitle"; else gnuplot << " title '" << line.name << "'";
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
					if(line.name == "" || line.name == "95% CI minimum" || line.name == "95% CI maximum") gnuplot << "notitle"; 
					else gnuplot << "title '" << line.name << "'";
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
					else gnuplot << "title '" << line.name << "'";
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
					else gnuplot << "title '" << line.name << "'";
				}
				gnuplot << endl;
				
				gnuplot << "unset logscale y" << endl;
				gnuplot << "unset logscale x" << endl;
				break;
		}
		gnuplot << endl << endl;
	}
}


/// Outputs a file containing posterior estimates for model parameters
void Output::posterior_parameter_estimates(vector <ParamSample> &psamp) const
{
	auto file = "Parameter_estimates.txt";
	auto filefull = details.output_directory+"/"+file;
	ofstream paramout(filefull.c_str());
	if(!paramout) emsg("Cannot output the file '"+filefull+"'");
	
	if(suppress_output == false) cout << "'" << file << "' gives the model parameters." << endl;
	paramout << "# Posterior distributions for model parameters." << endl;
	paramout << endl;
	
	paramout << "# Name, mean, 95% credible interval" << endl;
	
	vector <double> paramav(model.param.size());
	for(auto p = 0u; p < model.param.size(); p++){
		vector <double> vec; for(const auto psa : psamp) vec.push_back(psa.paramval[p]);
		auto stat = get_statistic(vec);
		paramav[p] = atof(stat.mean.c_str());
		if(model.param[p].name != "zero" && model.param[p].name != "one"){
			paramout << model.param[p].name  << "\t" << stat.mean << "\t(" << stat.CImin << " - "<< stat.CImax << ") " << endl; 
		}
	}
	
	cout << "Parameter posterior estimates are given in '" << filefull << "'" << endl << endl;
}


/// Outputs a spatial distribution for R
void	Output::spatial_R_map( const vector <Sample> &opsamp) const
{
	cout << opsamp[0].graph_state[0][0] << "\n";
	auto file = "Rmap.txt";
	auto filefull = details.output_directory+"/Posterior/"+file;
	ofstream Rmapout(filefull.c_str());
	if(!Rmapout) emsg("Cannot output the file '"+filefull+"'");
	
	if(suppress_output == false) cout << "'" << file << "' gives the time and spatial variation in R0." << endl;
	Rmapout << "# Gives the time and spatial variation in R." << endl;	
	Rmapout << "# Area, Day: ";
	for(auto st = 0u; st < details.ndivision; st++){ Rmapout << (st+1); if(st < details.ndivision-1) Rmapout << ", ";}
	Rmapout << endl;

// TO DO
/*	
	vector <double> areafactor = model.create_areafactor(paramav);
	auto areaav = 0.0; for(const auto& areafa : areafactor) areaav += areafa/data.narea;

	vector <double> Rav(details.ndivision);
	for(auto sett = 0u; sett < details.ndivision; sett++){
		auto sum = 0.0; for(auto s = opsampmin; s < nopsamp; s++) sum += opsamp[s].R0[sett];
		Rav[sett] = sum/(nopsamp-opsampmin);
	}
	
	for(auto c = 0u; c < data.narea; c++){
		Rmapout << data.area[c].code;
		for(auto st = 0u; st < details.ndivision; st++) Rmapout << "\t" << (areafactor[c]/areaav)*Rav[st];
		Rmapout << endl;
	}
	*/
	
	//if(details.mode != SIM){
	//	if(suppress_output == false) cout << "'trace.txt' gives trace plots for model parameters." << endl;
	//}
}


// For the Dirclet distributed quantities this provides the normalisation
void Output::dirichlet_normalisation(vector <ParamSample> &psamp) const
{
	vector <bool> param_used(model.param.size());
		
	for(auto s = 0u; s < psamp.size(); s++){  
		for(auto c = 0u; c < data.ndemocat; c++){  
			if(data.democat[c].sus_vari == true){
				auto sum = 0.0;	for(auto th : model.sus_param_list[c]) sum += psamp[s].paramval[th];
			
				for(auto j = 0u; j < model.sus_param_list[c].size(); j++){
					auto th = model.sus_param_list[c][j];
					psamp[s].paramval[th] /= (sum*model.sus_param_popfrac[c][j]);
				}
			}
		}
		
		for(auto c = 0u; c < model.comp.size(); c++){
			auto kmax = model.comp[c].trans.size();
			if(kmax > 1){
				for(auto dp = 0u; dp < data.ndemocatpos; dp++){
					for(auto th = 0u; th < model.param.size(); th++) param_used[th] = false; 
					
					auto sum = 0.0;
					for(auto k = 0u; k < kmax; k++){
						auto th = model.trans[model.comp[c].trans[k]].probparam[dp];
						sum += psamp[s].paramval[th];
						param_used[th] = true;
					}
					
					for(auto th = 0u; th < model.param.size(); th++) if(param_used[th] == true) psamp[s].paramval[th] /= sum;
				}
			}
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
void Output::simulated_data(const vector <double>& obs_value, string dir) const
{	
	cout << endl << "Simulated data in directory '" << dir <<"':" << endl;
	
	for(auto& dt : data.datatable){
		if(dt.optype == DATA){
			auto file = dt.file;
			auto filefull =  dir+"/"+file;
			cout << "  Generating data file '" << file << "'" << endl;
			
			string sep;
			if(stringhasending(file, ".txt")) sep = "\t"; else sep = ",";
					
			ofstream dataout(filefull);
			if(!dataout) emsg("Cannot output the file '"+filefull+"'");
			
			dataout << std::fixed << std::setprecision(0);
			if(!dataout) emsg("Cannot output the file '"+filefull+"'");
			
			switch(dt.type){
				case POP: case TRANS:
					{
						switch(details.time_format){
						case TIME_FORMAT_NUM: dataout << "time"; break;
						case TIME_FORMAT_YMD: dataout << "date"; break;
						}
						for(auto i : dt.graph_ref) dataout << sep << data.graph[i].colname;
						dataout << endl;

						if(dt.graph_ref.size() == 0) emsg("'"+file+"' contains no graphs");
						
						auto nrow = data.graph[dt.graph_ref[0]].point.size();
						for(auto row =0u; row < nrow; row++){
							dataout << details.getdate(dt.start - dt.shift + row*dt.units);
							for(auto i : dt.graph_ref){
								int val = obs_value[data.graph[i].point[row].obs[0]];
								if(val < 0) val = 0;
								dataout << sep << val;
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
							dataout << dt.demolist[row] << sep << obs_value[data.graph[i].point[row].obs[0]] << endl;
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
	
	if(vec.size() == 0) emsgEC("Output",45);
	
	double sum = 0, sumw = 0; 
	for(auto v : vec){ sum += v.val*v.w; sumw += v.w;}
	
	stat.mean = to_string(sum/sumw); 
	
	sort(vec.begin(),vec.end(),WeightedPoint_ord);
	
	auto n = vec.size();
	if(n >= 2){
		double sumtail = 0;
		auto i = 0; while(i < int(n) && sumtail < 0.025*sumw){ sumtail += vec[i].w; i++;}
		i--;
		double valw = vec[i].w;
		sumtail -= valw;
		double f = (0.025*sumw-sumtail)/valw;

		stat.CImin = to_string(vec[i].val*(1-f) + vec[i+1].val*f);
		
		sumtail = 0;
		i = n-1; while(i >= 0 && sumtail < 0.025*sumw){ sumtail += vec[i].w; i--;}
		i++;
		valw = vec[i].w;
		sumtail -= valw;
		f = (0.025*sumw-sumtail)/valw;

		stat.CImax = to_string(vec[i].val*(1-f) + vec[i-1].val*f);
	}
	else{
		stat.CImin = to_string(vec[0].val);
		stat.CImax = to_string(vec[0].val);
	}
	
	return stat;
}
		
		
/// Calculates diagnostic statistics
Statistics Output::get_statistic(const vector <double> &vec) const                       
{
	Statistics stat;
	
	auto n = vec.size();
	if(n == 0){
		stat.mean = "---"; stat.CImin = "---"; stat.CImax = "---"; 
	}
	else{
		auto sum = 0.0, sum2 = 0.0; for(auto v : vec){ sum += v; sum2 += v*v;}
		sum /= n; sum2 /= n;
		stat.mean = to_string(sum); 
		
		auto vec2 = vec;
		sort(vec2.begin(),vec2.end());
	
		if(n >= 2){
			auto i = (unsigned int)((n-1)*0.025); auto f = (n-1)*0.025 - i;
			stat.CImin = to_string(vec2[i]*(1-f) + vec2[i+1]*f);
				
			i = (unsigned int)((n-1)*0.975); f = (n-1)*0.975 - i;
			stat.CImax = to_string(vec2[i]*(1-f) + vec2[i+1]*f);
		}
		else{
			stat.CImin = to_string(vec2[0]);
			stat.CImax = to_string(vec2[0]);
		}
	}
	
	return stat;
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
			auto b = (unsigned int)(prop.nbin*(v-min)/(max-min)); if(b >= prop.nbin) emsgEC("Output",3);
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
void Output::generation_plot(string file, const vector <Generation> &generation) const
{
	auto nparam = model.param.size(); 

	string filefull = details.output_directory+"/"+file;
	ofstream genout(filefull.c_str());
	if(!genout) emsg("Cannot output the file '"+filefull+"'");

	genout << "# Generation, Time, EF cut-off"; 
	for(const auto& param : model.param) genout << "," << param.name << " (mean,CImin,CImax) ";
	genout << endl;
	
	long timestart = generation[0].time;
	for(auto g = 1u; g < generation.size(); g++){
		const Generation &gen=generation[g];
		
		genout << g << " " << (gen.time - timestart)/(60.0*CLOCKS_PER_SEC) << " ";

		switch(details.mode){
			case PAIS_INF:
				genout << gen.EFmax;
				break;
			default:
				genout << gen.EFcut;
				break;
		}
		
		auto npsamp = gen.param_samp.size();
		vector <ParamSample> psamp(npsamp);
		for(auto i = 0u; i < npsamp; i++) psamp[i].paramval = gen.param_samp[i];
		
		dirichlet_normalisation(psamp);
		
		for(auto th = 0u; th < nparam; th++){
			vector <WeightedPoint> vec;
			for(auto i = 0u; i < npsamp; i++){
				WeightedPoint pw;
				pw.val = psamp[i].paramval[th];
				switch(details.mode){
					case ABC_SMC: pw.w = gen.w[i]; break;
					case ABC_MBP: case ABC_MBP_GR: case PAIS_INF: case MC3_INF: pw.w = 1; break;
					default: emsgEC("Output",11); break;
				}
				vec.push_back(pw);
			}

			Statistics stat = get_statistic_with_weight(vec);
			
			genout << " " << stat.mean << " " << stat.CImin << " " << stat.CImax;
		}
		genout << endl;
	}	
}


/// Outputs the error function for different datatables
void Output::EF_datatable_plot(string file, const vector <Generation> &generation) const
{
	string filefull = details.output_directory+"/"+file;
	ofstream ldtyout(filefull.c_str());
	if(!ldtyout) emsg("Cannot output the file '"+filefull+"'");

	ldtyout << "# Generation";
	for(auto& dt : data.datatable) ldtyout << "\t" << dt.file;
	ldtyout << endl;
		 
	for(auto g = 0u; g < generation.size(); g++){
		const Generation &gen = generation[g];
		
		ldtyout << g;
		for(auto dt = 0u; dt < data.datatable.size(); dt++){
			auto val = 0.0; for(const auto &sa : gen.EF_datatable) val += sa[dt];
			val /= gen.EF_datatable.size();
			
			ldtyout << "\t" << val;
		}
		ldtyout << endl;
	}
}


/// Outputs results for the generations (ABC-SMC and ABC-MBP approaches)
void Output::generation_results(const vector <Generation> &generation) const 
{
	timer[TIME_RESULTS].start();
		
	if(mpi.core == 0){
		generation_plot("Diagnostics/Generation.txt",generation);
		EF_datatable_plot("Diagnostics/EF_datatable.txt",generation);
	}
	
	timer[TIME_RESULTS].stop();
}


/// Creates all the directories needed for outputting
void Output::ensure_directories() const
{
	ensure_directory(details.output_directory);
	ensure_directory(details.output_directory+"/Diagnostics");
	ensure_directory(details.output_directory+"/Posterior");
	ensure_directory(details.output_directory+"/Posterior/parameter");
	ensure_directory(details.output_directory+"/Posterior/state");
	ensure_directory(details.output_directory+"/Posterior/susceptibility");
	ensure_directory(details.output_directory+"/Posterior/spline");
	ensure_directory(details.output_directory+"/Posterior/samples");
	
	if(details.siminf == SIMULATE) ensure_directory(details.output_directory+"/Simulated_data");
}

		
/// Create a directory if it doesn't already exist
void Output::ensure_directory(const string &path) const
{
	struct stat st;
	if (stat(path.c_str(), &st) == -1){  	// Directory not found
		int ret = mkdir(path.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
		if(ret == -1) emsg("Error creating directory '"+path+"'");
	}
}


/// Summarises the data and model
void Output::print_model_specification() const 
{
	auto file = details.output_directory+"/Model_specification.txt";
	ofstream modspec(file); 
	if(!modspec) emsg("Cannot output the file '"+file+"'");
	modspec << data.print();                             
	modspec << model.print();
}
		
		
/// Initialises trace plot for parameters
void Output::trace_plot_inititialise(string name, ofstream &trace) const 
{
	auto file = details.output_directory+"/Diagnostics/"+name+".txt";

	trace.open(file);		
	if(!trace) emsg("The file '"+file+"' could not be opened");
	trace << "state";
	for(const auto& par : model.param) trace << "\t" << par.name; 
	trace << "\tLi";
	trace << "\tPri"; 
	trace << endl;
}


/// Outputs trace plot
void  Output::trace_plot(unsigned int samp, double Li, const vector <double> &paramval, ofstream &trace) const 
{
	trace << samp; 
	for(const auto& pval : paramval) trace << "\t" << pval; 

	trace << "\t" << Li; 
	trace << "\t" << model.prior(paramval);
	trace << endl;
}


/// Generates a readme file for the output directory
void Output::readme() const
{
	auto file = details.output_directory+"/README.md";
	ofstream read(file);
	if(!read) emsg("Cannot output the file '"+file+"'");
	
	read << "# Outputs" << endl << endl;
	read << "This directory contains the outputs from the BEEPmbp analysis."<< endl;
	read << "We give a brief description below of the various files and folders contained within it:" << endl << endl;
			
	read << "## Parameter_estimates.txt" << endl << endl;
	
	read << "If inference is performed, this file gives the posterior means and 95% credible intervals for each of the model parameters." << endl;
	
	read << endl;
	
	read << "## Graphs.pdf" << endl << endl;
	
	read << "This visualises outputs from BEEPmbp." << endl << endl;
	
	read << "## Graphs_description.txt" << endl << endl;
	
	read << "Provides information about the source data for the graphs." << endl << endl;
			
	read << "## Model_specification.txt" << endl << endl;
	
	read << "This gives a summary of the model used to perform the analysis." << endl;
	
	read << endl;
	
	read << "## Posterior" << endl << endl;
	
	read << "This folder contains files giving posterior probability distributions for the model parameters." << endl << endl;
	read << "These are arranged in a number of different ways:" << endl << endl;
	
	read << "**parameter** - Contains files for the model parameters." << endl << endl;
	read << "**state** - Contains files which compare the system state with that observed in the actual data files." << endl << endl;
	read << "**spline** - Contains files giving time variation in splines used within the model." << endl << endl;
	read << "**susceptibility** - Contains files giving the variation in susceptibility for different demographic classes within the model." << endl << endl;
	read << "**sample** - Contains files giving raw posterior samples." << endl << endl;
	read << "**Rmap.txt** - For spatial models this gives the variation in R across different regions." << endl << endl;
	read << endl;
	
	read << "## Simulated_data" << endl << endl;
	
	read << "This gives simulated data files corresponding to the specifications provided in the input TOML file." << endl << endl;
	
	read << "## Diagnostics" << endl << endl;
	
	read << "This provides diagnostic information to inform how well the algorithm is performing:" << endl << endl;
	
	read << "**Generation.txt** - Shows how model parameters move from the prior distribution to an approximation of the posterior distribution as a function of the generation number ('generations' are used within the ABC-MBP and ABC-SMC procedures to filter particles such that they provide a better and better approximation to the posterior)." << endl << endl;
	
	read << "**MCMC_proposals.txt** - This shows the performance of any MCMC proposals used." << endl;
}


/// This takes a vector of particle sample from each mpi core and combines them
void Output::generate_graphs(vector <Particle> &particle_store) const 
{
	timer[TIME_RESULTS].start();
	
	State state(details,data,model,obsmodel);
	
	auto dir = details.output_directory+"/Posterior/samples/";
	vector <ParamSample> psamp;
	vector <Sample> opsamp;
	mpi.gather_samples(psamp,opsamp,particle_store,state,dir);
	
	if(mpi.core == 0) generate_graphs(psamp,opsamp);
	
	timer[TIME_RESULTS].stop();
}


/// Prints percentage complete
void Output::print_percentage(const double s, const unsigned int nsamp, unsigned int &percentage) const
{
	auto stot = mpi.sum(s);
	if(mpi.core == 0){
		unsigned int step = 10;
		if(details.mode == PMCMC_INF || details.mode == MC3_INF){ step = 1; stot /= mpi.ncore;}
		
		unsigned int per = step*int((100/step)*stot/nsamp);
		if(per != percentage || percentage == UNSET){ cout << " " << per << "%" << endl; percentage = per;}
	}
}


/// Generates a new output plot
OutputPlot::OutputPlot(OutputPlotType type_, string title_, string xaxis_, string yaxis_, double min_, double max_)
{
	type = type_; title = title_; xaxis = xaxis_; yaxis = yaxis_; min = min_; max = max_;
}

	
/// Adds a line to an output plot
void OutputPlot::addline(string name, string file, unsigned int xcol, unsigned int ycol, LineType style)
{
	OutputLine opl; opl.name = name; opl.file = file; opl.xcol = xcol; opl.ycol = ycol; opl.style = style;
	line.push_back(opl);
};
