/// Generates the actual plots 
function generate_plots()
{
	plot = get_plot();
	var vis = visjson.plots[plot];
	
	switch(vis.type){
	case "OP_ANIM_MAP": get_colourbar(vis); break;
	case "OP_MIXING_WITHIN_MAP": case "OP_MIXING_WITHIN_ANIM_MAP": get_colourbar_linear(vis); break;
	case "OP_MIXING_BETWEEN_ANIM_MAP": case "OP_MIXING_BETWEEN_MAP": break;
	case "OP_MIXING_POP_ANIM_MAP": case "OP_MIXING_POP_MAP": break;
	case "OP_AREA_COVAR": case "OP_AREA_TV_COVAR":
		if(vis.xaxis == "Log transformed") get_colourbar(vis);
		else get_colourbar_linear(vis);
		break;
	case "OP_LEVEL_EFFECT": break;
	case "OP_AGE_MATRIX": case "OP_PARAM_TABLE": case "OP_PRIOR_TABLE": break;
	case "OP_COMP_MODEL": case "OP_FOI_MODEL": break; case "OP_DESC": break;
	case "OP_CPU": get_logaxrange(vis); break;
	default: get_axrange(vis); break;
	}

	switch(vis.type){
	case "OP_ANIM_MAP": 
	case "OP_MIXING_WITHIN_ANIM_MAP": case "OP_MIXING_WITHIN_MAP":
	case "OP_MIXING_BETWEEN_ANIM_MAP": case "OP_MIXING_BETWEEN_MAP": 
	case "OP_MIXING_POP_ANIM_MAP": case "OP_MIXING_POP_MAP":
	case "OP_AREA_COVAR": case "OP_AREA_TV_COVAR": 
	case "OP_LEVEL_EFFECT":
		mapframe(vis);
		break;
	
	case "OP_COMP_MODEL":
		compmodel_frame(vis);
		break;
		
	case "OP_FOI_MODEL":
		display_foi();
		break;
		
	case "OP_AGE_MATRIX":
		matrixframe(menux+10,10,vis);
		break;
	
	case "OP_SPLINE": case "OP_GRAPH": case "OP_GENERATION": case "OP_ME":
	case "OP_LOG_GENERATION": case "OP_TRACE": case "OP_TV_COVAR":
		graphframe(menux+10,10,40+w,60,20,40,w,vis.xaxis,vis.yaxis,vis);
		break;
		
	case "OP_CPU":
		loggraphframe(menux+10,10,40+w,60,20,40,w,vis.xaxis,vis.yaxis);
		break;
		
	case "OP_MARGINAL":	case "OP_GRAPH_MARGINAL":	
		histogram_frame(menux+10,10,40+w,60,20,40,w,vis.xaxis,vis.yaxis,vis.line[0].label);
		break;
		
	case "OP_PARAMDIST":
		paramdist_frame(menux+10,10,40,60,20,40,vis.xaxis,vis.yaxis);
		break;
		
	case "OP_PARAM_TABLE": case "OP_PRIOR_TABLE":
		table(menux+30,30,vis.param_table);
		break;
	
	case "OP_DESC":
		break;
		
	default: alertp("This is not currently supported");
	}
	
	plot_results();
		 
	x = menux+resdx+30; 
	if(vis.type == "OP_DESC") x = menux+30;
	wid = width-x-30;
	
	sourcey = height-30;

	if(sourceon == 1){                               // Incorporates the button for the source files
		var source = vis.source.replace(/\*/g,"'");
		var spl = source.split("|");
	
		var shift=[];
		var yy = 20;
		for(j = 0; j < spl.length-1; j++){
			alignparagraph(spl[j],wid,SOURCEFONT);
			shift.push(yy); yy += nlines*13;
		}
		yy -= 13;
		sourcey -= yy;
			
		for(j = 0; j < spl.length; j++){
			addbutton(spl[j],x,sourcey+shift[j],wid,20,-1,SOURCEFILEBUT,-1,-1);
		}
	}
	addbutton("SOURCE FILES",x,sourcey,wid,20,SOURCEBUT,SOURCEBUT,plot,-1);

	var topline;                                            // Adds the key
	switch(vis.type){
	case "OP_ANIM_MAP": case "OP_MIXING_WITHIN_MAP": case "OP_MIXING_WITHIN_ANIM_MAP":
	case "OP_AREA_COVAR": case "OP_AREA_TV_COVAR": 
		topline = sourcey - 40;
		addbutton("",x,topline,wid-10,20,-1,COLOURSCALEBUT,col_map.type,-1);
		break;
		
	case "OP_LEVEL_EFFECT": 
		topline = sourcey - vis.level_param.length*30-40;
		for(var i = 0; i < vis.level_param.length; i++){
			addbutton(vis.level_param[i],x+20,topline+30*i+20,wid-30,20,-1,LEVELKEYBUT,collist[i],-1);
		}
		break;
		
	case "OP_COMP_MODEL":
		topline = sourcey - 20;
		if(selparam_table != null){ 
			topline -= 20;
			var num = selparam_table.length; if(num > sellinemax) num = sellinemax;
			topline -= num*seltabledy;
			selparam_info(x,topline,wid);
		}
		break;
		
	case "OP_FOI_MODEL":
		topline = sourcey - 20;
		if(seleq != null){ 
			topline -= 20;
			var num = seleq.length; if(num > sellinemax) num = sellinemax;
			topline -= num*seltabledy;
			seleq_info(x,topline,wid);
		}
		break;
		
	case "OP_ME": case "OP_AGE_MATRIX":	case "OP_PARAM_TABLE": case "OP_PRIOR_TABLE": case "OP_MARGINAL": 
	case "OP_MIXING_BETWEEN_ANIM_MAP": case "OP_MIXING_BETWEEN_MAP":
	case "OP_MIXING_POP_ANIM_MAP": case "OP_MIXING_POP_MAP":
	case "OP_TV_COVAR":
	case "OP_GRAPH_MARGINAL":
	case "OP_LOG_GENERATION":
		topline = sourcey - 20;
		break;
		
	case "OP_DESC": 
		topline = sourcey - 20;
		break;
		
	default:
		if(vis.tab == "Dataa"){
			topline = sourcey - 20;
		}
		else{
			var list=[];
			for(li = 0; li < vis.line.length; li++){
				if(vis.line[li].name && vis.line[li].name != "" &&
					 vis.line[li].name != "95% CI min"& vis.line[li].name != "95% CI max"){
					list.push(li);
				}
			}
			topline = sourcey - list.length*20-40;
		
			var y = topline; 
			if(list.length > 0){
				for(i = 0; i < list.length; i++){
					li = list[i];
					if(vis.line[li].name){
						addbutton(vis.line[li].name,x+20,y,wid,0,-1,KEYBUT,plot,li); y += 20;
					}
				}		
			}
			topline -= 10;
		}
		break;
	}
	
	if(vis.popscale != null){
		topline -= 50;
		var ww = Math.floor(wid/2)-5;
		addbutton("Per 100K individuals",x+5,topline,ww,15,RADIOBUT,RADIOBUT,"rate",RADIORATE); 
		addbutton("Absolute",x+5,topline+20,ww,15,RADIOBUT,RADIOBUT,"absolute",RADIORATE);
		
		addbutton("No highlight",x+10+ww,topline,ww,15,RADIOBUT,RADIOBUT,"off",RADIOHIGH); 
		addbutton("Lower highlight",x+10+ww,topline+20,ww,15,RADIOBUT,RADIOBUT,"on",RADIOHIGH);
	}
	
	if(vis.spline_param != null){
		topline -= 20;
		addbutton("Show parameters",x+15,topline,17,17,CHECKBUT,CHECKBUT,-1,-1);
		topline -= 10;
	}

	var menuname = tree[page].child[ps].child[pss].child[psss].name;
	if(menuname != ""){
		var i = 0; while(i < param_link.length && param_link[i].name != menuname) i++;
		if(i < param_link.length){
			topline -= 20;
			addbutton("Model link: '"+menuname+"'",x+5,topline,wid,30,MENULINKBUT,MENULINKBUT,menuname,-1);	
			topline -= 10;
		}
	}

	y = 30;
	
	var sp = vis.fulldesc.split(":");                          // Draws the title and description
	addbutton(sp[0],x,y,wid,0,-1,TITLEBUT,-1,-1); y += 24;

	var desc = sp[1];
	if(vis.popscale != null){
		if(rateradio == "rate") desc += "  Note, this shows results per 100,000 individuals."; 
		if(lowerhighlight == "on"){
			desc += " Lower highlighting is turned on (values near to the smallest are coloured blue, ";
			desc += "which can be an effective way to visualise when cases first appear).";
		}
	}
	
	topline -= 5;
	dy = topline-y;
	dy = Math.floor(dy/lineheight)*lineheight;

	var content = desc.substr(1).replace(/\*/g, "'");
	addbutton(content,x,y,wid,dy,-1,SLIDEPARAGRAPHBUT,-1,-1);

	alignparagraph(desc.substr(1),wid);
		
	nlines_disp = Math.floor(dy/lineheight);
	slidefrac = nlines_disp/nlines; 
	if(slidefrac < 1){	
		addbutton("",x+wid,y,10,dy,SLIDEBUT,SLIDEBUT,-1,-1);
	}
}


/// Returns the plot from the menu 
function get_plot()
{
	return tree[page].child[ps].child[pss].child[psss].plot;	
}


/// Plots the results on the "resultcan" canvas
function plot_results()
{
	plot = get_plot();
	var vis = visjson.plots[plot];
	
	switch(vis.type){
	case "OP_ANIM_MAP": drawmap(vis);	break;
	case "OP_MIXING_WITHIN_MAP": case "OP_MIXING_WITHIN_ANIM_MAP": drawmap(vis);	break;
	case "OP_AREA_COVAR": case "OP_AREA_TV_COVAR": drawmap(vis); break;
	case "OP_LEVEL_EFFECT": drawmap(vis); break;
	case "OP_COMP_MODEL": drawcompmodel(vis); break;
	case "OP_FOI_MODEL": case "OP_DESC": break;
	case "OP_MIXING_BETWEEN_ANIM_MAP": case "OP_MIXING_BETWEEN_MAP": drawmixingmap(vis); break;
	case "OP_MIXING_POP_ANIM_MAP": case "OP_MIXING_POP_MAP": drawmixingmap(vis); break;
	case "OP_AGE_MATRIX": drawmatrix(vis); break;
	case "OP_PARAM_TABLE": case "OP_PRIOR_TABLE": break;
	case "OP_SPLINE": case "OP_GRAPH": drawlines(vis); draw_timelabels(); break;
	case "OP_GENERATION": case "OP_LOG_GENERATION": case "OP_TRACE": drawlines(vis); break;
	case "OP_TV_COVAR": drawlines(vis); break;
	case "OP_CPU": drawloglines(vis); break;
	case "OP_ME":	drawME(vis); break;
	case "OP_MARGINAL": drawmarginal(vis); break;
	case "OP_GRAPH_MARGINAL": drawmarginal(vis); break;
	case "OP_PARAMDIST": drawparamdist(vis); break;
	default: alertp("plot not supported"); break;
	}
}


/// Gets the bounds for the axes
function get_ax_bound(vis)
{
	axxmin = large; axxmax = -large;
	axymin = large; axymax = -large;
	for(li = 0; li < vis.line.length; li++){
		var visli =  vis.line[li];
		var xcol = visli.xcol;
		if(xcol){
			for(i = 0; i < xcol.length; i++){
				if(xcol[i] < axxmin) axxmin = xcol[i];
				if(xcol[i] > axxmax) axxmax = xcol[i];
			}
		}
		
		var ycol = visli.ycol;
		if(ycol){
			for(i = 0; i < ycol.length; i++){
				if(ycol[i] < axymin) axymin = ycol[i];
				if(ycol[i] > axymax) axymax = ycol[i];
			}
		}
		
		var ebmin = visli.errbarmin;
		var ebmax = visli.errbarmax;
		if(ebmax){
			for(i = 0; i < ebmin.length; i++){
				if(ebmin[i] < axymin) axymin = ebmin[i];
				if(ebmax[i] > axymax) axymax = ebmax[i];
			}
		}
	}
}


/// Sets up the x and y axes
function get_axrange(vis)
{
	get_ax_bound(vis)
	
	//pr(vis.type);
	switch(vis.type){
		case "OP_MARGINAL":	case "OP_GRAPH_MARGINAL": case "OP_SPLINE": case "OP_GRAPH":
			if(axymin > 0) axymin = 0; 
			break;
	}
	
	d = 0.1*(axymax - axymin);
	
	if(axymax > -d && axymax <= 0) axymax = 0; else axymax += d;
	
	if(axymin < d && axymin >= 0) axymin = 0; else axymin -= d;
	
	if(axymin >= 0 && axymax >= 0){
		if(axymin < 0.5*axymax) axymin = 0;
	}
	
	d = 0.05*(axxmax - axxmin);
	
	if(axxmax > -d && axxmax <= 0) axxmax = 0; else axxmax += d;
	
	if(axxmin < d && axxmin >= 0) axxmin = 0; else axxmin -= d;
	
	if(axxmin >= 0 && axxmax >= 0){
		if(axxmin < 0.3*axxmax) axxmin = 0;
	}
	
	if(axymin == axymax){ axymin = 0; axymax *= 1.2;}
	
	setxtics();
	w = setytics();   
}


/// Sets up log scale axes
function get_logaxrange(vis)
{
	get_ax_bound(vis)
	setlogxtics();
	w = setlogytics();   
}


/// Generate array trans based on scale per 100k individuals
function generate_arraytrans(vis)
{
	vis.arraytrans=[];
	for(t = 0; t < vis.array.length; t++){
		vis.arraytrans[t] = [];
		for(c = 0; c < vis.array[t].length; c++){
			if(vis.popscale == null || rateradio != "rate") vis.arraytrans[t][c] = vis.array[t][c];
			else vis.arraytrans[t][c] = vis.array[t][c]*100000/vis.popscale[c];
		}
	}
}


/// Sets up a colour bar (used for maps)
function get_colourbar(vis)
{	
	generate_arraytrans(vis);

	axymin = large; axymax = -large;
	for(t = 0; t < vis.arraytrans.length; t++){
		for(c = 0; c < vis.arraytrans[t].length; c++){
			val = vis.arraytrans[t][c];
			if(val < axymin) axymin = val;
			if(val > axymax) axymax = val;
		}
	}

	if(axymin == axymax){
		if(axymin == 0) axymax = 1;
		else{ axymin *= 0.9; axymax *= 1.1;} 
	}
	
	var shift = 0;
	if(axymin <= 0.05*axymax) shift = 0.05*axymax - axymin;
	
	var light = 220;
	
	var collist=[];
	
	if(vis.yaxis == "Colour"){
		if(axymin < 1) collist.push({R:0, G:255, B:0, val:Math.log(axymin+shift)});
		else collist.push({R:light, G:255, B:light, val:Math.log(axymin+shift)});
	
		if(axymin < 1 && axymax > 1){
			collist.push({R:light, G:255, B:light, val:0});
			collist.push({R:255, G:light, B:light, val:0});
		}
	
		if(axymax > 1) collist.push({R:255, G:0, B:0, val:Math.log(axymax+shift)});
		else collist.push({R:light, G:255, B:light, val:Math.log(axymax+shift)});
	}
	else{
		collist.push({R:255, G:255, B:255, val:Math.log(axymin+shift)});
		collist.push({R:0, G:0, B:0, val:Math.log(axymax+shift)});
	}
			
	var ncoldiv = 100;

	col_map = {type:"log", div:ncoldiv, min:Math.log(axymin+shift), max:Math.log(axymax+shift), shift:shift, map:[]};
	
	var i=0;
	for(div = 0; div <= ncoldiv; div++){
		val = Math.log(axymin+shift) + ((Math.log(axymax+shift) - Math.log(axymin+shift))*div)/ncoldiv;
		while(i < collist.length-2 && val > collist[i+1].val) i++; 

		frac = (val - collist[i].val)/(collist[i+1].val - collist[i].val);
		
		if(vis.yaxis != "Colour" && lowerhighlight == "on" && div < 1){
			col_map.map.push({ R:0, G:0, B:255, val:val});
		}
		else{
			col_map.map.push({ R:Math.floor(collist[i].R*(1-frac)) + Math.floor(collist[i+1].R*frac), 
									 G:Math.floor(collist[i].G*(1-frac)) + Math.floor(collist[i+1].G*frac),
									 B:Math.floor(collist[i].B*(1-frac)) + Math.floor(collist[i+1].B*frac),
									 val:val});
		}
	}
	
	nticky = 0;
	ticky = [];
	var min = Math.exp(col_map.min)-shift, max = Math.exp(col_map.max)-shift;
	for(range = 5; range >= -5; range--){
		add_tick(10 ** range,min,max,shift);
		add_tick((10 ** range)*2,min,max,shift);
		add_tick((10 ** range)*5,min,max,shift);
	}

	if(nticky <= 2){
		var nearmin = min*0.99+max*0.01;
		add_tick(nearmin.toPrecision(2),min,max,shift);
		var nearmax = min*0.01+max*0.99;
		add_tick(nearmax.toPrecision(2),min,max,shift);
	}
	if(min < 0.00000001){ ticky.push(0); nticky++;}
}


/// Adds ticks to the colour bar
function add_tick(val,min,max,shift)
{
	if(val >= min && val <= max && ((nticky < 10 && val > shift) || nticky < 6)){ ticky.push(val); nticky++;} 
}


function get_colourbar_linear(vis)
{
	generate_arraytrans(vis)
	
	axymin = large; axymax = -large;
	for(t = 0; t < vis.arraytrans.length; t++){
		for(c = 0; c < vis.arraytrans[t].length; c++){
			val = vis.arraytrans[t][c];
			if(val < axymin) axymin = val;
			if(val > axymax) axymax = val;
		}
	}
	
	var light = 220;
	
	var collist=[];

	collist.push({R:255, G:255, B:255, val:axymin});
	collist.push({R:0, G:0, B:0, val:axymax});
			
	var ncoldiv = 100;

	col_map = {type:"linear", div:ncoldiv, min:axymin, max:axymax, shift:0, map:[]};
	
	var i=0;
	for(div = 0; div <= ncoldiv; div++){
		val = axymin + ((axymax - axymin)*div)/ncoldiv;
		while(i < collist.length-2 && val > collist[i+1].val) i++; 

		frac = (val - collist[i].val)/(collist[i+1].val - collist[i].val);
		col_map.map.push({ R:Math.floor(collist[i].R*(1-frac)) + Math.floor(collist[i+1].R*frac), 
									 G:Math.floor(collist[i].G*(1-frac)) + Math.floor(collist[i+1].G*frac),
									 B:Math.floor(collist[i].B*(1-frac)) + Math.floor(collist[i+1].B*frac),
									 val:val});
	}
	
	setytics();
}


/// Sets up ticks along the y axis
function setytics()                                      
{
	i = nticksize-1; while(i >= 0 && Math.floor((axymax-axymin)/ticksize[i]) < 3) i--;
	axticky = ticksize[i];

	mm = Math.floor(axymin/axticky + 0.9999);
	nticky = 0;
	while(mm*axticky < axymax){
		ticky[nticky] = rn(mm*axticky); nticky++;
		mm++;
	}

	w = 0;
	for(i = 0; i < nticky; i++){
		ww = textwidth(ticky[i],TICKFONT);
		if(ww > w) w = ww;
	}
	
	return w;
}


/// Sets up ticks along the x axis
function setxtics()
{
	i = nticksize-1; while(i >= 0 && Math.floor((axxmax-axxmin)/ticksize[i]) < 4) i--;
	axtickx = ticksize[i];

	tx = Math.floor(axxmin/axtickx + 0.999999)*axtickx;
	ntickx = 0;
	while(tx < axxmax){
		tickx[ntickx] = rn(tx); ntickx++;	
		tx += axtickx;
	}
}


/// Sets up logarithmically distributed ticks along the y axis
function setlogytics() 
{
	nticky = 0;
	ticky = [];
	for(range = -5; range < 5; range++){
		val = 10 ** range; if(val >= axymin && val <= axymax){ ticky.push(val); nticky++;}
		val = (10 ** range)*2; if(val >= axymin && val <= axymax){ ticky.push(val); nticky++;}
		val = (10 ** range)*5; if(val >= axymin && val <= axymax){ ticky.push(val); nticky++;}
	}
	
	w = 0;
	for(i = 0; i < nticky; i++){
		ww = textwidth(ticky[i],TICKFONT);
		if(ww > w) w = ww;
	}
	
	return w;
}


/// Sets up logarithmically distributed ticks along the x axis
function setlogxtics()
{
	ntickx = 0;
	tickx = [];
	for(range = -5; range < 5; range++){
		val = 10 ** range; if(val >= axxmin && val <= axxmax){ tickx.push(val); ntickx++;}
		val = (10 ** range)*2; if(val >= axxmin && val <= axxmax){ tickx.push(val); ntickx++;}
		val = (10 ** range)*5; if(val >= axxmin && val <= axxmax){ tickx.push(val); ntickx++;}
	}
}


/// Rounds a number when creating ticks
function rn(n)                                          
{
	var x;
	n = parseFloat(n);

	x = 0;
	do{
		num = n.toFixed(x);
		if((num-n)*(num-n) < 0.000000000000001) break;
		x++;
	}while(1 == 1);
	return num;
}


/// Generates a frame for a map
function mapframe(vis)
{
	graphdx = mapdx; graphdy = mapdy;
	addbutton("",menux,0,mapdx,mapdy,CANVASBUT,CANVASBUT,-1,-1);
	addcanbutton("",0,0,graphdx,graphdy,-1,RESULTBUT2,-1,-1);	
	
	if(vis.type == "OP_ANIM_MAP" || vis.type == "OP_MIXING_WITHIN_ANIM_MAP" ||
	   vis.type == "OP_MIXING_BETWEEN_ANIM_MAP" || vis.type == "OP_AREA_TV_COVAR" || vis.type == "OP_LEVEL_EFFECT"){
		addcanbutton("",15,graphdy-45,30,30,PLAYBUT,PLAYBUT,-1,-1);		
		addcanbutton("",60,graphdy-37,graphdx-60-10-60,15,PLAYLINEBUT,PLAYLINEBUT,-1,-1);		
	}
	
	var dx = 22, dy = 26;
	addcanbutton("",mapdx-dx-dx-15,graphdy-dy-10,dx,dy,ZOOMINBUT,ZOOMINBUT,-1,-1);
	addcanbutton("",mapdx-dx-10,graphdy-dy-10,dx,dy,ZOOMOUTBUT,ZOOMOUTBUT,-1,-1); 
}


/// Sets up a frame for drawing the compartmental model
function compmodel_frame(vis)
{
	graphdx = mapdx; graphdy = mapdy;
	addbutton("",menux,0,mapdx,mapdy,CANVASBUT,CANVASBUT,-1,-1);
	addcanbutton("",0,0,graphdx,graphdy,-1,RESULTBUT2,-1,-1);	

	var comp = model.comp;	
	var shiftx = (mapdx-mapdy)/2;
	
	var h = 25 - comp.length;
	if(h > 20) h = 20; if(h < 15) h = 15;
	
	for(k = 0; k < comp.length; k++){  
		if(comp[k].label != null){
			var x = shiftx + mapdy*comp[k].x;
			var y = mapdy*(0.5 + (comp[k].y+comp[k].h/2));
						
			var texth = Math.floor(0.8*h); 
			var font = texth+"px times";
			var text = comp[k].label;
		
			lab = splitintosub(text,font);
			w = lab[lab.length-1].w+lab[lab.length-1].dw;
			w += 6;
			
			addcanbutton(text,x-w/2,y,w,h,VARIABLEBUT,VARIABLEBUT,font,comp[k].mean);
		}
	}		
	
	var h = 30 - comp.length;
	if(h > 24) h = 24; if(h < 15) h = 15;
	
	
	var trans = model.trans;	
	for(k = 0; k < trans.length; k++){  
		if(trans[k].label != null){
			var x = shiftx + mapdy*trans[k].labelx;
			var y = mapdy*(0.5 + (trans[k].labely));
				
			var texth = Math.floor(0.8*h); 
			var font = texth+"px times";
			var text = trans[k].label;
		
			w = textwidth(text,font);
			orth = h*0.7;
			addcanbutton(text,x-w/2 +orth*trans[k].labelnx ,y-h/2+orth*trans[k].labelny,w,h,VARIABLEBUT,VARIABLEBUT,font,trans[k].prob);
		}
	}		
	
	//var dx = 22, dy = 26;
	//addcanbutton("",mapdx-dx-dx-15,10,dx,dy,ZOOMINBUT,ZOOMINBUT,-1,-1);
	//addcanbutton("",mapdx-dx-10,10,dx,dy,ZOOMOUTBUT,ZOOMOUTBUT,-1,-1); 
}


/// Displays the equation for the force of infection
function display_foi()
{
	y = 50;
	addbutton("",menux,0,mapdx,mapdy,CANVASBUT,CANVASBUT,-1,-1);
	
	var si = 30;
	var w;
	
	do{
		var fl = false;
		
		var font = si+"px times";
	
		w = textwidth(foi.name,font)+3; 
		w += textwidth("=",font)+3; 
		for(i = 0; i < foi.eq.length; i++) w += textwidth(foi.eq[i].text,font)+3; 
		
		if(w < mapdx-100) x = Math.floor(mapdx/2-w/2);
		else{ si--; fl = true;}
	}while(fl == true);
	
	var text = foi.name;
	w = textwidth(text,font)+3; 
	h = 25;
			
	addcanbutton(text,x,y,w,h,-1,EQBUT3,font,-1);
	x += w;
	
	text = "=";
	w = textwidth(text,font)+3; 
	addcanbutton(text,x,y,w,h,-1,EQBUT3,font,-1);
	x += w;
	
	for(i = 0; i < foi.eq.length; i++){
		text = foi.eq[i].text;
		w = textwidth(text,font)+3; 
		if(foi.eq[i].param != null){
			addcanbutton(text,x,y,w,h,EQBUT,EQBUT,font,foi.eq[i].param);
		}
		else{
			addcanbutton(text,x,y,w,h,-1,EQBUT3,font,-1);
		}
		x += w;
	}
		
	si = 15;
	
	do{
		var fl = false;
		
		var font = si+"px arial";
		var wmax = 0;
		for(i = 0; i < foi.desc.length; i++){
			w = textwidth(foi.desc[i],font);
			if(w > wmax) wmax = w;
		}
		
		if(wmax < mapdx-100) x = Math.floor(mapdx/2-wmax/2);
		else{ si--; fl = true;}
	}while(fl == true);
	
	y = 110; 
	var dy = (mapdy - y-30)/foi.desc.length; if(dy > 30) dy = 30; if(dy < 17) dy = 17;
	
	for(i = 0; i < foi.desc.length; i++){
		addcanbutton(foi.desc[i],x,y+dy*i,w,h,-1,EQBUT2,font,-1);
	}
}


/// Generates a frame for displaying the age mixing matrix
function matrixframe(x,y,vis)
{
	w = 0;
	for(i = 0; i < vis.ages.length; i++){
		ww = textwidth(vis.ages[i],TICKFONT);
		if(ww > w) w = ww;
	}
	w = Math.floor(w+20);
	
	x2 = 20; y1 = w; y2 = 20;
	
	graphdy = resdy - (y1+y2);
	graphdx = graphdy; 
	x1 = mapdx/2-graphdx/2;

	addbutton("",x,y,resdx,resdy,CANVASBUT,CANVASBUT,-1,-1);
	addcanbutton("",x1,y2,graphdx,graphdy,-1,RESULTBUT4,-1,-1);		

	d = graphdx/vis.ages.length;

	for(i = 0; i < vis.ages.length; i++){
		addcanbutton(vis.ages[i],x1+Math.floor(i*d),y2 + graphdy,Math.floor(d),w,-1,MATRIXXBUT,-1,-1);
		addcanbutton(vis.ages[i],x1-w,y2 + Math.floor(i*d),w,Math.floor(d),-1,MATRIXYBUT,-1,-1);
	}
}


/// Draws the heading for the table
function table(x,y,tab)
{
	plottype = "table"; 
	 
	var linemax = Math.floor((resdy-60)/tabledy);

	graphdx = resdx-40; 
	graphdy = (1+linemax)*tabledy;
	
	var head = tab.heading;
	var ele = tab.ele;
	
	var xx = 0, xxst=[];
	for(c = 0; c < head.length; c++){
		var wmax = 0;
		var w = textwidth_simp(head[c],TABLEHEADFONT);
		if(w > wmax) wmax = w;
		for(row = 0; row < ele.length; row++){
			var w = textwidth_simp(ele[row][c],TABLEFONT);
			if(w > wmax) wmax = w;
		}		
		wmax += 10;
		xxst[c] = xx;
		xx += wmax;
	}
	
	var gap = 8;
	
	tabx_slidefrac = graphdx/xx; 
	if(tabx_slidefrac < 1){	
		addbutton("",x,y+graphdy+gap,graphdx,10,TABXSLIDEBUT,TABXSLIDEBUT,-1,-1);
	}
	
	var dx = Math.floor(tabx_slide*xx);
	
	for(c = 0; c < head.length; c++){
		var xx = xxst[c] - dx;
		addcanbutton(head[c],xx,0,wmax,tabledy,-1,TABLEHEADBUT,-1,-1);
		
		var jjmin = Math.floor(taby_slide*ele.length*1.001);
		var jjmax = jjmin + linemax; if(jjmax > ele.length) jjmax = ele.length;
		for(j = jjmin; j < jjmax; j++){
			if(c == 0) addcanbutton(ele[j][c],xx,tabledy*(j-jjmin+1),xxst[c+1]-xxst[c],tabledy,LINKBUT,LINKBUT,j,-1);
			if(c == 0) addcanbutton(ele[j][c],xx,tabledy*(j-jjmin+1),xxst[c+1]-xxst[c],tabledy,LINKBUT,LINKBUT,j,-1);
			else addcanbutton(ele[j][c],xx,tabledy*(j-jjmin+1),0,tabledy,-1,TABLEBUT,-1,-1);
		}
	}

	addbutton("",x,y,graphdx,graphdy,CANVASBUT,CANVASBUT,-1,-1);
	
	taby_slidefrac = linemax/ele.length; 
	if(taby_slidefrac < 1){	
		addbutton("",x+graphdx+gap,y+tabledy,10,graphdy-tabledy,TABYSLIDEBUT,TABYSLIDEBUT,-1,-1);
	}
}


/// Plots information about a selected model parameter
function selparam_info(x,y,wid)
{
	var wmax = textwidth_simp("Symbol",NORMTABLEHEADFONT); 
	for(row = 0; row < selparam_table.length; row++){
		lab = splitintosub(selparam_table[row].param,NORMTABLEFONT);
		w = lab[lab.length-1].w+lab[lab.length-1].dw;
		if(w > wmax) wmax = w;
	}
	wmax += 5;
	
	addbutton("Symbol",x,y,0,seltabledy,-1,NORMTABLEHEADBUT,-1,-1);
	addbutton("Parameter name",x+wmax,y,0,seltabledy,-1,NORMTABLEHEADBUT,-1,-1);
	
	taby_slidefrac = sellinemax/selparam_table.length; 
	if(taby_slidefrac < 1){	
		addbutton("",x+wid,y+seltabledy,10,sellinemax*seltabledy,TABYSLIDEBUT,TABYSLIDEBUT,-1,-1);
	}
	
	var rowmin = Math.floor(taby_slide*selparam_table.length*1.001);
	var rowmax = rowmin + sellinemax; 
	if(rowmax > selparam_table.length) rowmax = selparam_table.length;
	for(row = rowmin; row < rowmax; row++){
		var yy = y+seltabledy*(row-rowmin+1);
		addbutton(selparam_table[row].param,x,yy,wmax,seltabledy,-1,NORMTABLEBUT,-1,-1);
		addbutton(selparam_table[row].name,x+wmax,yy,wid-wmax,seltabledy,-1,NORMTABLEBUT,-1,-1);
	}
}


/// Plots information about a selected equation
function seleq_info(x,y,wid)
{	
	addbutton("Parameter name",x,y,0,seltabledy,-1,NORMTABLEHEADBUT,-1,-1);
	
	taby_slidefrac = sellinemax/seleq.length; 
	if(taby_slidefrac < 1){	
		addbutton("",x+wid,y+seltabledy,10,sellinemax*seltabledy,TABYSLIDEBUT,TABYSLIDEBUT,-1,-1);
	}
	
	var rowmin = Math.floor(taby_slide*seleq.length*1.001);
	var rowmax = rowmin + sellinemax; 
	if(rowmax > seleq.length) rowmax = seleq.length;
	for(row = rowmin; row < rowmax; row++){
		var yy = y+seltabledy*(row-rowmin+1);
		addbutton(seleq[row],x,yy,wid,seltabledy,-1,NORMTABLEBUT,-1,-1);
	}
}


/// Generates a frame for drawing a graph
function graphframe(x,y,x1,y1,x2,y2,w,labx,laby,vis)
{	
	graphdx = resdx - (x1+x2); graphdy = resdy - (y1+y2);

	addbutton("",x,y,resdx,resdy,CANVASBUT,CANVASBUT,-1,-1);
	addcanbutton("",x1,y2,graphdx,graphdy,-1,RESULTBUT,-1,-1);		

	addcanbutton(labx,x1,y2+graphdy+30,graphdx,30,-1,XLABELBUT,-1,-1);
	addcanbutton(laby,0,y2,30,graphdy,-1,YLABELBUT,-1,-1);
	
	if(ntickx > 0) addcanbutton("X ticks",x1,y2+graphdy,graphdx,30,-1,XTICKBUT,-1,-1);
	if(nticky > 0) addcanbutton("Y ticks",x1-w,y2,w,graphdy,-1,YTICKTRBUT,-1,-1);
	
	if(check == 1 && vis.spline_param != null){
		for(i = 0; i < vis.spline_param.length; i++){
			var x = x1+Math.floor(graphdx*(vis.spline_param[i].time-axxmin)/(axxmax-axxmin));
			
			addcanbutton(vis.spline_param[i].param,x,y2,15,graphdy,PARAMLINEBUT,PARAMLINEBUT,-1,-1);
		}
	}
}


/// Generates a frame for drawing a graph with logarithmic axes
function loggraphframe(x,y,x1,y1,x2,y2,w,labx,laby)
{	
	graphdx = resdx - (x1+x2); graphdy = resdy - (y1+y2);

	addbutton("",x,y,resdx,resdy,CANVASBUT,CANVASBUT,-1,-1);
	addcanbutton("",x1,y2,graphdx,graphdy,-1,RESULTBUT,-1,-1);		

	addcanbutton(labx,x1,y2+graphdy+30,graphdx,30,-1,XLABELBUT,-1,-1);
	addcanbutton(laby,0,y2,30,graphdy,-1,YLABELBUT,-1,-1);
	
	if(ntickx > 0) addcanbutton("X ticks",x1,y2+graphdy,graphdx,30,-1,LOGXTICKBUT,-1,-1);
	if(nticky > 0) addcanbutton("Y ticks",x1-w,y2,w,graphdy,-1,LOGYTICKTRBUT,-1,-1);
}


/// Generates a frame for drawing a parameter distribution
function paramdist_frame(x,y,x1,y1,x2,y2,labx,laby)
{	
	graphdx = resdx - (x1+x2); graphdy = resdy - (y1+y2);

	addbutton("",x,y,resdx,resdy,CANVASBUT,CANVASBUT,-1,-1);
	addcanbutton("",x1,y2,graphdx,graphdy,-1,RESULTBUT,-1,-1);		

	addcanbutton(labx,x1,y2+graphdy+30,graphdx,30,-1,XLABELBUT,-1,-1);
	addcanbutton(laby,0,y2,30,graphdy,-1,YLABELBUT,-1,-1);
	
	if(ntickx > 0) addcanbutton("X ticks",x1,y2+graphdy,graphdx,30,-1,XTICKBUT,-1,-1);
}


/// Generates a frame for drawing a histogram
function histogram_frame(x,y,x1,y1,x2,y2,w,labx,laby,cats)
{
	dx = Math.floor(graphdx/cats.length);

	var wmax = 0;
	for(i = 0; i < cats.length; i++){
		var w = textwidth_simp(cats[i],CATFONT)
		if(w > wmax) wmax = w;
	}
	var dy = 0;
	if(wmax > dx) dy = wmax-20;
	
	graphdx = resdx - (x1+x2); graphdy = resdy - (y1+y2+dy);
	
	if(wmax > dx){    // Labels are placed vertically
		for(i = 0; i < cats.length; i++){
			xx = Math.floor(x1 + graphdx*i/cats.length);		
			addcanbutton(cats[i],xx,y2+graphdy+10,dx,10,-1,CATVERTBUT,-1,-1);
		}
	}
	else{             // Labels 
		for(i = 0; i < cats.length; i++){
			xx = Math.floor(x1+graphdx*i/cats.length);		
			addcanbutton(cats[i],xx,y2+graphdy+5,dx,10,-1,CATBUT,-1,-1);
		}
	}

	addbutton("",x,y,resdx,resdy,CANVASBUT,CANVASBUT,-1,-1);
	addcanbutton("",x1,y2,graphdx,graphdy,-1,RESULTBUT3,-1,-1);		

	addcanbutton(labx,x1,y2+graphdy+30+dy,graphdx,30,-1,XLABELBUT,-1,-1);
	addcanbutton(laby,0,y2,30,graphdy,-1,YLABELBUT,-1,-1);
	
	if(nticky > 0) addcanbutton("Y ticks",x1-w,y2,w,graphdy,-1,YTICKTRBUT,-1,-1);
}


/// Zooms in or out from a map
function mapzoom(fac,x,y)
{
	zoomx += (1-1.0/fac)*map_ratio*(x/mapdx)/zoomfac;
	zoomy += (1-1.0/fac)*(1-y/mapdy)/zoomfac;
	zoomfac *= fac;
	
	buttoninit();
}


/// Fires when an animation on a map is started
function playanim()
{
	playtime = Math.floor(playtimeinit + (getsec()- playstartt)*playtimemax/6);
	
	if(playtime >= playtimemax){ playtime = playtimemax-1; playing = 0;}
	
	plot_results();
	buttonplot();

	if(playing == 1 && playtime < playtimemax-1) setTimeout(function(){ playanim();}, 10);
}


/// Sets the outline for an area
function set_outline(v)
{
	cv.beginPath();
	for(i = 0; i <= v.length; i++){
		var ii = i%v.length;
		var x = mapdy*zoomfac*(v[ii][0]-zoomx), y = mapdy - mapdy*zoomfac*(v[ii][1]-zoomy);
		if(i == 0) cv.moveTo(x,y);
		else cv.lineTo(x,y);
	}	
}


/// Draws a map onto "resultscan"
function drawmap(vis)
{
	plottype = "map";
			
	cv = resultcv; 
	cv.clearRect(0,0,graphdx,graphdy);
	
	playtimemax = vis.array.length;
	
	var b = visjson.boundaries;
	
	for(c = 0; c < b.length; c++){                                  // Fills the areas
		for(p = 0; p < b[c].length; p++){
			set_outline(b[c][p]);
	
			if(vis.type == "OP_LEVEL_EFFECT"){
				cv.fillStyle = collist[vis.array[playtime][c]];
			}
			else{
				var val;
				if(col_map.type == "log") val = Math.log(vis.arraytrans[playtime][c] + col_map.shift);		
				else val = vis.arraytrans[playtime][c];		
				
				var div = Math.floor(0.999999*col_map.div*(val-col_map.min)/(col_map.max-col_map.min));
				cv.fillStyle = "rgb("+col_map.map[div].R+","+col_map.map[div].G+","+col_map.map[div].B+")";
			}
			cv.fill();
		}
	}
	
	cv.lineWidth = 0.4;
	for(c = 0; c < b.length; c++){                                  // Draws the boundaries
		for(p = 0; p < b[c].length; p++){
			set_outline(b[c][p]);
			cv.strokeStyle = BLACK;
			cv.stroke();
		}
	}
	
	if(areaover != -1){
		var val;
		if(vis.type == "OP_LEVEL_EFFECT") val = vis.level_param[vis.array[playtime][areaover]];
		else val = vis.arraytrans[playtime][areaover];
		highlight_area(vis,areaover,val,vis.areas[areaover]);
	}
	
	if(vis.dates != null) centertext(vis.dates[playtime],mapdx/2,30,MAPDATEFONT,BLACK);
}


/// Draws the compartmental model
function drawcompmodel(vis)
{
	plottype = "compmodel";
			
	cv = resultcv; 
	cv.clearRect(0,0,graphdx,graphdy);
	cv.lineWidth = 0.5;
	
	var shiftx = (mapdx-mapdy)/2;
	
	for(i = 0; i < model.comp.length; i++){
		var co = model.comp[i];
		var x = shiftx + mapdy*(co.x-co.w/2);
		var y = mapdy*(0.5 + (co.y-co.h/2));
		var w = Math.floor(mapdy*co.w);
		var h = Math.floor(mapdy*co.h);
		drawroundrect(x,y,w,h,h/5,co.col,darkcol(co.col));  
		centertext(co.name,x+w/2,y+Math.floor(h*0.75),Math.floor(h*0.8)+"px Times",WHITE);
	}
	
	for(i = 0; i < model.trans.length; i++){
		var tr = model.trans[i];
		var arx=[], ary=[];
		var jmax = tr.arrow.length;
		for(j = 0; j < jmax; j++){
			arx[j] = shiftx + mapdy*tr.arrow[j].x;
			ary[j] = mapdy*(0.5 + tr.arrow[j].y);
		}
		for(k = 0; k < arx.length-1; k++){
			drawline(arx[k],ary[k],arx[k+1],ary[k+1],BLACK,1,0); 
		}
		drawarrow(arx[jmax-1],ary[jmax-1],arx[jmax-2],ary[jmax-2],7,BLACK);   
	}
}


/// Highlights an area on the map
function highlight_area(vis,c,val,name)
{
	var b = visjson.boundaries;
	
	cv.lineWidth = 1;
	for(p = 0; p < b[c].length; p++){
		cv.beginPath();
		for(i = 0; i <= b[c][p].length; i++){
			var ii = i%b[c][p].length;
			var x = mapdy*zoomfac*(b[c][p][ii][0]-zoomx), y = mapdy - mapdy*zoomfac*(b[c][p][ii][1]-zoomy);
			if(i == 0) cv.moveTo(x,y);
			else cv.lineTo(x,y);
		}
	
		cv.strokeStyle = BLACK;
		cv.stroke();
	}
	var ra = visjson.click_bound_range[c];
	var x = mapdy*zoomfac*(ra.xav-zoomx);

	if(val != ""){
		if(vis.type != "OP_LEVEL_EFFECT") val = val.toPrecision(3);
		
		ww = textwidth(val,MAPBOLDFONT)+5;
		y = mapdy - mapdy*zoomfac*(ra.yav-zoomy);
		cv.globalAlpha = 0.9;
		drawroundrect(x-ww/2,y-10,ww,13,3,WHITE,WHITE); 
		cv.globalAlpha = 1;
		centertext(val,x,y+1,MAPBOLDFONT,BLACK);
	}
	
	ww = textwidth(vis.areas[c],MAPFONT)+5;
	y = mapdy - mapdy*zoomfac*(ra.ymin-zoomy);
	cv.globalAlpha = 0.9;
	drawroundrect(x-ww/2,y+2,ww,13,3,WHITE,WHITE); 
	cv.globalAlpha = 1;
	centertext(name,x,y+14,MAPFONT,BLACK);
}
	
	
/// Draws a map showing mixing between regions onto "resultscan"
function drawmixingmap(vis)
{
	plottype = "map";
			
	playtimemax = vis.mat_list.length;
	
	var mat = vis.mat_list[playtime];
	
	cv = resultcv; 
	cv.clearRect(0,0,graphdx,graphdy);
	cv.lineWidth = 0.5;
	
	var b = visjson.boundaries;
	for(c = 0; c < b.length; c++){
		for(p = 0; p < b[c].length; p++){
			cv.beginPath();
			for(i = 0; i <= b[c][p].length; i++){
				var ii = i%b[c][p].length;
				var x = mapdy*zoomfac*(b[c][p][ii][0]-zoomx), y = mapdy - mapdy*zoomfac*(b[c][p][ii][1]-zoomy);
				if(i == 0) cv.moveTo(x,y);
				else cv.lineTo(x,y);
			}
		
			cv.strokeStyle = BLACK;
			cv.stroke();
		}
	}
	
	for(c = 0; c < visjson.click_bound_range.length; c++){
		var ra = visjson.click_bound_range[c];
		var x = mapdy*zoomfac*(ra.xav-zoomx);
		var y = mapdy - mapdy*zoomfac*(ra.yav-zoomy);
		fillcircle(x,y,3,BLACK,BLACK,1);   	
	}

	for(i = 0; i < vis.areas.length; i++){
		if(areaover == -1 || i == areaover){
			var rai = visjson.click_bound_range[i];
			var x1 = mapdy*zoomfac*(rai.xav-zoomx);
			var y1 = mapdy - mapdy*zoomfac*(rai.yav-zoomy);
			
			var max;
			if(i == areaover){
				max = 0; for(j = 0; j < vis.areas.length; j++){ if(mat[i][j] > max) max = mat[i][j];}
			}
			else max = vis.max;
			
			var j;
			for(j = 0; j < vis.areas.length; j++){
				if(mat[i][j] != null){
					var raj = visjson.click_bound_range[j];
					var x2 = mapdy*zoomfac*(raj.xav-zoomx);
					var y2 = mapdy - mapdy*zoomfac*(raj.yav-zoomy);
			
					var w = 7*mat[i][j]/max;
					
					var dx = x2-x1, dy = y2-y1;
					var nx = dy, ny = -dx;
					var r = Math.sqrt(dx*dx+dy*dy);
					nx /= r; ny /= r;
					var end = 3/r; if(end > 0.3) end = 0.3;
					var neck = end + 2*w/r; if(neck > 0.5) neck = 0.5;
					var fac = 2.5;
					var mid = 0;
					
					if(areaover == -1){
						var wtot = 7*(mat[i][j]+ mat[j][i])/max;
						mid = wtot/2;
						if(mid < 3) mid = 3;
					}
					
					cv.fillStyle = BLUE;
					cv.beginPath();	
					cv.moveTo(x1+end*dx+(mid+w/2)*nx,y1+(mid+w/2)*ny+end*dy);
					cv.lineTo(x2-neck*dx+(mid+w/2)*nx,y2+(mid+w/2)*ny-neck*dy);
					cv.lineTo(x2-neck*dx+(mid+fac*w/2)*nx,y2+(mid+fac*w/2)*ny-neck*dy);
					cv.lineTo(x2-end*dx+(mid)*nx,y2+(mid)*ny-end*dy);
					cv.lineTo(x2-neck*dx+(mid-fac*w/2)*nx,y2+(mid-fac*w/2)*ny-neck*dy);
					cv.lineTo(x2-neck*dx+(mid-w/2)*nx,y2+(mid-w/2)*ny-neck*dy);
					cv.lineTo(x1+end*dx+(mid-w/2)*nx,y1+(mid-w/2)*ny+end*dy);
					cv.fill();
					
					if(i == areaover && mat[i][j] > 0.5*max){
						cv.save();
						cv.translate(0.5*(x1+x2) + 0*(2+0.7*w)*nx,0.5*(y1+y2) + 0*(2+0.7*w)*ny);
						cv.rotate(Math.atan(dy/dx));
						cv.textAlign = 'center';
						cv.fillStyle = RED;
						cv.fillText(mat[i][j].toPrecision(2), 0, -(2+0.7*w));
						cv.restore();
					}
				}
			}
		}
	}
	
	if(areaover != -1) highlight_area(vis,areaover,"",vis.areas[areaover])
}
	
	
/// Draws a matrix giving the mixing between different age groups
function drawmatrix(vis)
{
	plottype = "matrix";
	
	cv = resultcv; 
	cv.clearRect(0,0,graphdx,graphdy);
	
	var nage = vis.ages.length;

	var max = 0;		
	for(i = 0; i < nage; i++){
		for(j = 0; j < nage; j++){
			if(vis.ele[i][j] > max) max = vis.ele[i][j];
		}
	}	
	
	var d = graphdx/nage;
	for(i = 0; i < nage; i++){
		for(j = 0; j < nage; j++){
			cv.globalAlpha = vis.ele[i][j]/max;
			fillrect(i*d,j*d,d,d,RED);     
		}
	}
	cv.globalAlpha = 1;
}


/// Draw lines onto a graph
function drawlines(vis)
{
	plottype = "graph";
	
	cv = resultcv; 
	cv.clearRect(0,0,graphdx,graphdy);
	cv.lineWidth = 2;
	
	for(li = 0; li < vis.line.length; li++){
		var visli = vis.line[li];
		var xcol = visli.xcol;
		if(xcol){
			var ycol = visli.ycol;
		
			cv.beginPath(); 
			for(i = 0; i < xcol.length; i++){
				x = Math.floor(graphdx*(xcol[i]-axxmin)/(axxmax-axxmin));
				y = Math.floor(graphdy-graphdy*(ycol[i]-axymin)/(axymax-axymin));
				
				if(i == 0) cv.moveTo(x,y);
				else cv.lineTo(x,y);
			}	
		}
		else{
			if(visli.true){
				y = Math.floor(graphdy-graphdy*(visli.true-axymin)/(axymax-axymin));
				cv.beginPath(); 
				cv.moveTo(0,y);
				cv.lineTo(graphdx,y);
			}
			else alertp("Problem");
		}
		setstyle(visli.style);
		cv.stroke();
	} 	
	
	setdash(0);
}


/// Draws the time labels onto the graphs
function draw_timelabels()
{
	for(i = 0; i < visjson.time_labels.length; i++){
		var x = Math.floor(graphdx*(visjson.time_labels[i].time-axxmin)/(axxmax-axxmin));
		add_verticle_line(visjson.time_labels[i].name,x,BLACK);
	}

	for(i = 0; i < visjson.pred_timeplot.length; i++){
		var x = Math.floor(graphdx*(visjson.pred_timeplot[i].time-axxmin)/(axxmax-axxmin));
		add_verticle_line(visjson.pred_timeplot[i].name,x,DBLUE);
	}
}


/// Draws lines using log axes
function drawloglines(vis)
{
	plottype = "graph";
	
	cv = resultcv; 
	cv.clearRect(0,0,graphdx,graphdy);
	cv.lineWidth = 2;
	
	for(li = 0; li < vis.line.length; li++){
		var visli = vis.line[li];
		var xcol = visli.xcol;
		var ycol = visli.ycol;
	
		cv.beginPath(); 
		for(i = 0; i < xcol.length; i++){
			x = Math.floor(graphdx*(Math.log(xcol[i])-Math.log(axxmin))/(Math.log(axxmax)-Math.log(axxmin)));
			y = Math.floor(graphdy-graphdy*(Math.log(ycol[i])-Math.log(axymin))/(Math.log(axymax)-Math.log(axymin)));
			if(i == 0) cv.moveTo(x,y);
			else cv.lineTo(x,y);
		}
	
		setstyle(visli.style);
		cv.stroke();
	} 	
	
	setdash(0);
}


/// Draws the model evidence plot
function drawME(vis)
{
	plottype = "graph";
	
	cv = resultcv; 
	cv.clearRect(0,0,graphdx,graphdy);
	cv.lineWidth = 2;
	
	var visli = vis.line[0];
	var xcol = visli.xcol;
	var ycol = visli.ycol;
	var sd = visli.name;

	cv.beginPath(); 
	for(i = 0; i < xcol.length; i++){
		x = Math.floor(graphdx*(xcol[i]-axxmin)/(axxmax-axxmin));
		y = Math.floor(graphdy-graphdy*(ycol[i]-axymin)/(axymax-axymin));
		
		if(i == 0) cv.moveTo(x,y);
		else cv.lineTo(x,y);
	}	
	setstyle(visli.style);
	cv.stroke();
	
	if(visli.errbarmin != null){
		for(i = 0; i < xcol.length; i++){
			x = Math.floor(graphdx*(xcol[i]-axxmin)/(axxmax-axxmin));
			draw_errorbar(x,ycol[i],visli.errbarmin[i],visli.errbarmax[i],6);
		}
	}
	
	setdash(0);
}


/// Draws posterior distributions for the parameters
function drawparamdist(vis)
{
	plottype = "graph";
	
	cv = resultcv; 
	cv.clearRect(0,0,graphdx,graphdy);
	cv.lineWidth = 2;
	
	var visli = vis.line[0];
	var xcol = visli.xcol;
	var ycol = visli.ycol;
	
	cv.beginPath(); 
	for(i = 0; i < xcol.length; i++){
		x = Math.floor(graphdx*(xcol[i]-axxmin)/(axxmax-axxmin));
		y = Math.floor(graphdy-graphdy*(ycol[i]-axymin)/(axymax-axymin));
		
		if(i == 0) cv.moveTo(x,y);
		else cv.lineTo(x,y);
	}	
	
	setstyle(visli.style);
	cv.stroke();
	
	cv.lineTo(x,graphdy);
	x = Math.floor(graphdx*(xcol[0]-axxmin)/(axxmax-axxmin));
	cv.lineTo(x,graphdy);
	cv.globalAlpha = 0.3;
	cv.fillStyle = "#22ff22";
	cv.fill();
	cv.globalAlpha = 1;
		
	if(vis.line.length > 1){
		var visli = vis.line[1];
		if(visli.true){
			cv.beginPath(); 
			x = Math.floor(graphdx*(visli.true-axxmin)/(axxmax-axxmin));

			cv.moveTo(x,0);
			cv.lineTo(x,graphdy);
			setstyle(vis.line[1].style);
			cv.stroke();
		}
	}	

	setdash(0);
}


/// Draws a vertical line with a label
function add_verticle_line(text,x,col)
{
	cv.strokeStyle = "#000000"; setdash(0);
	cv.lineWidth = 1;
	cv.beginPath(); 
	cv.moveTo(x,0);
	cv.lineTo(x,graphdy);
	cv.stroke();
	verticaltext(text,x+4,5,VERTFONT,col); 
}


/// Links from a parameter name to where it appears in the model
function follow_param_link(param_name)
{
	var i = 0; while(i < param_link.length && param_link[i].name != param_name) i++;
	if(i == param_link.length){ alertp("Cannot find parameter"); return;}
	
	if(param_link[i].selparam_table != null){
		changepage(0,0,-1,-1);
		paramsel = param_link[i].name; 
		selparam_table = param_link[i].selparam_table;
		selparam_name = param_link[i].selparam_name;
	}
	else{
		changepage(0,1,-1,-1);
		paramsel = param_link[i].name;
		seleq = param_link[i].seleq;
		seleq_name = param_link[i].seleq_name;
	}
}


/// Draws a histogram representing a marginal distribution
function drawmarginal(vis)
{
	plottype = "marginal";
	
	cv = resultcv; 
	cv.clearRect(0,0,graphdx,graphdy);

	var li = vis.line[0];
	var xcol = li.xcol;
	var ycol = li.ycol;
	var imax = xcol.length;
	for(i = 0; i < imax; i++){
		var x1 = Math.floor((i+0.25)*graphdx/imax);
		var x2 = Math.floor((i+0.75)*graphdx/imax);
		var y = Math.floor(graphdy-graphdy*(ycol[i]-axymin)/(axymax-axymin));
		fillrect(x1,y,x2-x1,graphdy-y,RED);  
	}
	
	if(li.errbarmin){	
		var dx = graphdx/imax;
		var tic = 4;
		if(tic > dx/4) tic = dx/4;

		for(i = 0; i < imax; i++){
			draw_errorbar(Math.floor((i+0.5)*graphdx/imax),ycol[i],li.errbarmin[i],li.errbarmax[i],tic);
		}
	}
	
	if(vis.line.length > 1){
		var li = vis.line[1];
		var ycol = li.ycol;
		for(i = 0; i < imax; i++){
			var x1 = Math.floor((i)*graphdx/imax);
			var x2 = Math.floor((i+1)*graphdx/imax);
			var y = Math.floor(graphdy-graphdy*(ycol[i]-axymin)/(axymax-axymin));
			drawline(x1,y,x2,y,BLACK,2,0);
		}
	}
}


/// Draws an error bar
function draw_errorbar(xx,y,min,max,tic)
{
	var yy = Math.floor(graphdy-graphdy*(y-axymin)/(axymax-axymin));
	var yymin = Math.floor(graphdy-graphdy*(min-axymin)/(axymax-axymin));
	var yymax = Math.floor(graphdy-graphdy*(max-axymin)/(axymax-axymin));

	var si = 0.7*tic;
	drawline(xx,yymin,xx,yymax,BLACK,1,0);
	drawline(xx-tic,yymin,xx+tic,yymin,BLACK,1,0);
	drawline(xx-tic,yymax,xx+tic,yymax,BLACK,1,0);
	drawline(xx-si,yy-si,xx+si,yy+si,BLACK,1,0);
	drawline(xx-si,yy+si,xx+si,yy-si,BLACK,1,0);
}


/// Sets the style for a line
function setstyle(style)
{
	switch(style){
		case "RED_SOLID": cv.strokeStyle = "#ff2222"; setdash(0); break;
		case "RED_DASHED": cv.strokeStyle = "#ffaaaa"; setdash(1); break;
		case "GREEN_SOLID": cv.strokeStyle = "#22ff22"; setdash(0); break;
		case "GREEN_DASHED": cv.strokeStyle = "#aaffaa"; setdash(1); break;
		case "BLUE_SOLID": cv.strokeStyle = "#2222ff"; setdash(0); break;
		case "BLUE_DASHED": cv.strokeStyle = "#aaaaff"; setdash(1); break;
		case "BLACK_SOLID": cv.strokeStyle = "#000000"; setdash(0); break;
		case "BLACK_DASHED": cv.strokeStyle = "#000000"; setdash(1); break;
		case "RED_DOTTED": cv.strokeStyle = "#ff2222"; setdash(2); break;
		case "RED_DOTDASH": cv.strokeStyle = "#ffaaaa"; setdash(3); break;
		case "GREEN_DOTTED": cv.strokeStyle = "#22ff22"; setdash(2); break;
		case "GREEN_DOTDASH": cv.strokeStyle = "#aaffaa"; setdash(3); break;
		case "BLUE_DOTTED": cv.strokeStyle = "#2222ff"; setdash(2); break;
		case "BLUE_DOTDASH": cv.strokeStyle = "#aaaaff"; setdash(3); break;
		case "BLACK_DOTTED": cv.strokeStyle = "#000000"; setdash(2); break;
		case "BLACK_DOTDASH": cv.strokeStyle = "#000000"; setdash(3); break;
		case "YELLOW_SOLID": cv.strokeStyle = "#ffff22"; setdash(0); break;
		case "YELLOW_DASHED": cv.strokeStyle = "#ffff22"; setdash(1); break;
		case "CYAN_SOLID": cv.strokeStyle = "#22ffff"; setdash(0); break;
		case "CYAN_DASHED": cv.strokeStyle = "#22ffff"; setdash(1); break;
		case "MAGENTA_SOLID": cv.strokeStyle = "#ff22ff"; setdash(0); break;
		case "MAGENTA_DASHED": cv.strokeStyle = "#ff22ff"; setdash(1); break;
	}
}


/// Initialises maps show they can be clicked (or hovered over)
function init_clickable_map()
{
	var b = visjson.boundaries;
	visjson.inter = [];
	visjson.click_bound_range=[];
		
	for(c = 0; c < b.length; c++){
		var xmin = large, xmax = -large, ymin = large, ymax = -large; 
		var xav = 0, yav = 0, nav = 0; 
		for(p = 0; p < b[c].length; p++){          
			for(k = 0; k < b[c][p].length; k++){
				var x = b[c][p][k][0], y = b[c][p][k][1];
				if(x > xmax) xmax = x; if(x < xmin) xmin = x;
				if(y > ymax) ymax = y; if(y < ymin) ymin = y;
				xav += x; yav += y; nav++;
			}
		}
		visjson.click_bound_range[c] = {xmin:xmin, xmax:xmax, ymin:ymin, ymax:ymax, xav:xav/nav, yav:yav/nav};
		
		visjson.inter[c] = [];
		for(yi = Math.floor(clickline*ymin); yi <= Math.floor(clickline*ymax); yi++){
			visjson.inter[c][yi]=[];
		}
	
		for(p = 0; p < b[c].length; p++){          
			for(k = 0; k < b[c][p].length; k++){
				var po = b[c][p][k];
				xi = b[c][p][k][0];
				xf = b[c][p][(k+1)%b[c][p].length][0];
				bx = xf-xi; 
				
				yi = clickline*b[c][p][k][1];
				yf = clickline*b[c][p][(k+1)%b[c][p].length][1];
				by = yf-yi;
				if(by > 0){
					for(j = Math.floor(yi+1); j < yf; j++) visjson.inter[c][j].push(xi + bx*(j-yi)/by);
				}
				
				if(by < 0){
					for(j = Math.floor(yi); j > yf; j--) visjson.inter[c][j].push(xi + bx*(j-yi)/by);
				}
				//if(by == 0) visjson.inter[c][Math.floor(yi)].push(xi);
			}
		}
	}	
}
