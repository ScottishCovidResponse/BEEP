/// Initialises everything before the app is started
function init()
{
	if(!jsonstr){
		ById("loading").innerHTML = "visBEEPmbp must be loaded from the output directory of an analysis";
		return;
	}
	visjson = JSON.parse(jsonstr);
	
	menu_structure();     // Initialises the menu structure on the left
	
	init_page();          // Initialises HTML elements on the page

	init_canvas();        // Initialises canvases for plotting results
	
	init_compmodel();     // Initialises the compartmental model
	
	init_mouse_key();     // Initialises input operations using the mouse and keys
	
	init_variables();     // Initialises various parameters
	
	init_clickable_map(); // Initialises clickable maps
	
	initdone = 1;
	
	setsize();
}


/// Gets the compartment number from the name
function get_comp(name)
{
	var c = 0; while(c < model.comp.length && model.comp[c].name != name) c++;
	if(c == model.comp.length) alertp("Cannot find name:"+name);
	return c;
}


/// Gets permutations 
function permutator(inputArr) {
  var results = [];

  function permute(arr, memo) {
    var cur, memo = memo || [];

    for (var i = 0; i < arr.length; i++) {
      cur = arr.splice(i, 1);
      if (arr.length === 0) {
        results.push(memo.concat(cur));
      }
      permute(arr.slice(), memo.concat(cur));
      arr.splice(i, 0, cur[0]);
    }

    return results;
  }

  return permute(inputArr);
}

/*
function get_permute(vec,vec2)
{
	for(i = 0; i < vec.length; i++){
		var vecnew = [];
		for(j = 0; j < 
	}
}
*/

/// Initialises how the comparmental model is displayed on a page
function init_compmodel()
{
	var i = 0; while(i < visjson.plots.length && visjson.plots[i].type != "OP_COMP_MODEL") i++;
	if(i == visjson.plots.length) alertp("Cannot find comparmtental model");
	model = visjson.plots[i].model;
	
	i = 0; while(i < visjson.plots.length && visjson.plots[i].type != "OP_FOI_MODEL") i++;
	if(i == visjson.plots.length) alertp("Cannot find force of infection");
	foi = visjson.plots[i].foi;
	
	var comp = model.comp;
	var trans = model.trans;
	
	for(k = 0; k < trans.length; k++){
		trans[k].fromc = get_comp(trans[k].from);
		trans[k].toc = get_comp(trans[k].to);
	}		
	
	var i = 0; while(i < trans.length && !trans[i].infection) i++;
	if(i == trans.length) alertp("Cannot find the infection transition");
	
	var sus = trans[i].fromc;
	
	var wmax = 0;
	for(i = 0; i < comp.length; i++){
		var w = textwidth_simp(comp[i].name,textsize+"px arial");
		if(w > wmax) wmax = w;
	}
	var h = textsize+6;
	wmax += 10;
	if(wmax < 1.5*h) wmax = 1.5*h;
	
	var colcomp=[], colcompy=[];
	
	comp[sus].x = 0;
	comp[sus].w = wmax;
	comp[sus].h = h;
	comp[sus].col = 0;
	comp[sus].colnum = 0;
	colcomp[0]=[];
	colcomp[0].push(sus);
	colcompy[0]=[];
	colcompy[0].push(0);
	
	col = 1;
	while(1 == 1){
		colcomp[col]=[]; colcompy[col]=[];
		
		for(cc = 0; cc < comp.length; cc++){
			if(comp[cc].col == null){
				for(k = 0; k < trans.length; k++){
					if(trans[k].toc == cc){
						c = trans[k].fromc;
						if(comp[c].col != null){
							if(comp[c].col == col-1) break;
						}
					}
				}
				
				if(k < trans.length){
					comp[cc].x = wmax*2*col;
					comp[cc].w = wmax;
					comp[cc].h = h;
					comp[cc].col = col;
					comp[cc].colnum = colcomp[col].length;
					colcomp[col].push(cc);
				}
			}
		}
		jmax = colcomp[col].length;
		
		if(jmax == 0) break;
		if(jmax == 1) colcompy[col][0] = 0;
		else{
			for(j = 0; j < jmax; j++) colcompy[col][j] = ((j/(jmax-1))-0.5)*2*h*(jmax-1);
		}
		col++;
	}

	var perm = [];                             // This works out the best y positions by doing permutations
	for(c = 0; c < col; c++){ 	
		var arr=[]; for(i = 0; i < colcomp[c].length; i++) arr.push(i);
		perm[c] = permutator(arr);
	}
	
	var index = [];
	for(c = 0; c < col; c++) index[c] = 0;
	var index_best = [];
	var arrowmax = large;
	do{
		var arrow = 0;
		for(k = 0; k < trans.length; k++){
			var from = trans[k].fromc;
			var to = trans[k].toc;
			var dx = comp[from].x-comp[to].x;
			var fromcol = comp[from].col, fromcolnum = comp[from].colnum;
			var yfr = colcompy[fromcol][perm[fromcol][index[fromcol]][fromcolnum]];
			var tocol = comp[to].col, tocolnum = comp[to].colnum;
			var yto = colcompy[tocol][perm[tocol][index[tocol]][tocolnum]];
			var dy = yfr-yto;
	
			arrow += Math.sqrt(dx*dx+dy*dy);
		}		
		
		if(arrow < arrowmax){
			for(c = 0; c < col; c++) index_best[c] = index[c];
			arrowmax = arrow;
		}		
	
		j = 0;
		var fl;
		do{
			fl = false
			index[j]++; if(index[j] == perm[j].length){ index[j] = 0; j++; fl = true;}
		}while(j < col && fl == true);
	}while(j < col);

	for(c = 0; c < col; c++){
		for(j = 0; j < colcomp[c].length; j++){
			comp[colcomp[c][j]].y = colcompy[c][perm[c][index_best[c]][j]];
		}
	}
	
	var coltop = [], colbot = [];
	for(c = 0; c < col; c++){
		colbot[c] = colcompy[c][0] - h;
		coltop[c] = colcompy[c][colcomp[c].length-1] + h;
	}
	
	for(k = 0; k < trans.length; k++){
		var from = trans[k].fromc;
		var fromx = comp[from].x, fromy = comp[from].y; 
		var fromw = comp[from].w, fromh = comp[from].h; 
		
		var to = trans[k].toc;
		var tox = comp[to].x, toy = comp[to].y; 
	 	var tow = comp[to].w, toh = comp[to].h; 
	
		trans[k].arrow = [];
		
		if(comp[to].col >= comp[from].col){
			var fromint = roundrect_intersect(fromx,fromy,tox,toy,fromx-fromw/2,fromy-fromh/2,fromw,fromh,fromh/5);
			var toint = roundrect_intersect(tox,toy,fromx,fromy,tox-tow/2,toy-toh/2,tow,toh,toh/5);
			
			trans[k].arrow.push({x:fromint.x,y:fromint.y});
			trans[k].arrow.push({x:toint.x,y:toint.y});		
		}
		else{
			var y;
			if(fromy < 0){	
				y = large;
				for(c = comp[to].col; c <= comp[from].col; c++){
					if(colbot[c] < y) y = colbot[c];
				}
				
				for(c = comp[to].col; c <= comp[from].col; c++) colbot[c] = y - 0.5*h;
			}
			else{
				y = -large;
				for(c = comp[to].col; c <= comp[from].col; c++){
					if(coltop[c] > y) y = coltop[c];
				}
				
				for(c = comp[to].col; c <= comp[from].col; c++) coltop[c] = y + 0.5*h;
			}
			
			var fromint = roundrect_intersect(fromx,fromy,fromx,y,fromx-fromw/2,fromy-fromh/2,fromw,fromh,fromh/5);
			var toint = roundrect_intersect(tox,toy,tox,y,tox-tow/2,toy-toh/2,tow,toh,toh/5);
		
			trans[k].arrow.push({x:fromint.x,y:fromint.y});
			trans[k].arrow.push({x:fromint.x,y:y});
			trans[k].arrow.push({x:toint.x,y:y});
			trans[k].arrow.push({x:toint.x,y:toint.y});
		}
	}
	
	for(k = 0; k < comp.length; k++){   // Adds labels to compartments
		if(comp[k].mean != null){
			if(comp[k].shape == 1){
				if(comp[k].mean.length > 1) comp[k].label = "Exp(&m^{"+comp[k].name+"}_d&)";
				else comp[k].label = "Exp(&m^{"+comp[k].name+"}&)";	
			}
			else{
				if(comp[k].mean.length > 1) comp[k].label = "Γ(&m&^{"+comp[k].name+"}_d,"+comp[k].shape+")";
				else comp[k].label = "Γ(&m&^{"+comp[k].name+"},"+comp[k].shape+")";	
			}
		}		
	}
	
	for(k = 0; k < trans.length; k++){   // Adds labels to transitions
		if(trans[k].prob != null || trans[k].infection != null){
			if(trans[k].infection != null) trans[k].label = foi.name;
			else{
				if(trans[k].prob.length > 1) trans[k].label = "&b&_d^{"+trans[k].from+"→"+trans[k].to+"}";
				else trans[k].label = "&b&^{"+trans[k].from+"→"+trans[k].to+"}";
			}
			
			var i = Math.floor((trans[k].arrow.length-1)/2);
			trans[k].labelx = (trans[k].arrow[i].x+trans[k].arrow[i+1].x)/2;
			trans[k].labely = (trans[k].arrow[i].y+trans[k].arrow[i+1].y)/2;
			var nx = trans[k].arrow[i+1].y - trans[k].arrow[i].y;
			var ny = -(trans[k].arrow[i+1].x - trans[k].arrow[i].x);
			var d = Math.sqrt(nx*nx+ny*ny);
			nx /= d; ny /= d;
			if(ny > 0){ nx *= -1; ny *= -1;}
			trans[k].labelnx = nx;
			trans[k].labelny = ny;
		}		
		
		if(trans[k].infection != null){
			
		}
	}

	var xmin = large, xmax = -large;   // Scales
	
	for(k = 0; k < comp.length; k++){
		var x = comp[k].x+comp[k].w/2;
		if(x > xmax) xmax = x;
		var x = comp[k].x-comp[k].w/2;
		if(x < xmin) xmin = x;
	}
	
	scale = 1.0/(xmax-xmin);
	
	for(k = 0; k < comp.length; k++){
		comp[k].x = (comp[k].x-xmin)*scale;
		comp[k].y = comp[k].y*scale;
		comp[k].w *= scale;
		comp[k].h *= scale;
		comp[k].col = collist[k];
	}
	
	for(k = 0; k < trans.length; k++){
		for(i = 0; i < trans[k].arrow.length; i++){
			trans[k].arrow[i].x = (trans[k].arrow[i].x-xmin)*scale;
			trans[k].arrow[i].y = trans[k].arrow[i].y*scale;
		}
		if(trans[k].label != null){
			trans[k].labelx = (trans[k].labelx-xmin)*scale;
			trans[k].labely = trans[k].labely*scale;
		}
	}
}


/// Intialises the HTML elements on the page
function init_page()
{
	ById("container").style.visibility = "hidden";          // Sets up elements on page
	ById("container").style.height = "0px";
	ById("main").style.visibility = "visible";		
	ById("bod").style.backgroundColor = "#ddddff";
		
	logopic = new Image();                                  // Loads the logo picture
	logopic.src = "vis_files/logo_menu.png";
	logopic.onload = function(){ buttonplot();};	
}


/// Initialises canvas objects for plotting the results 
function init_canvas()
{
	myc = ById("myCanvas");
	maincv = myc.getContext("2d");
	cv = maincv;

	graphcan = document.createElement('canvas');
  graphcv = graphcan.getContext('2d');

	resultcan = document.createElement('canvas');
  resultcv = resultcan.getContext('2d');	
}


/// Initialises global variables 
function init_variables()
{
	nticksize = 0;
	for(sh = -10; sh <= 10; sh++){
		for(i = 0; i < 3; i++){   
			f = tickpo[i];
			if(sh < 0){ for(j = 0; j < -sh; j++) f /= 10;}
			else{ for(j = 0; j < sh; j++) f *= 10;}
			ticksize[nticksize] = f; nticksize++;
		}
	}
	
	zoomfac = 1; zoomx = 0; zoomy = 0;
}


/// Based on the plots on visjson, this constructs the menus on the left hand edge
function menu_structure()
{
	for(i = 0; i < visjson.plots.length; i++){
		var vis = visjson.plots[i];
		tab = vis.tab; tab2 = vis.tab2; tab3 = vis.tab3;
	
		var j = 0; while(j < tree.length && tab != tree[j].name) j++;
		if(j == tree.length){
			tree.push({ name:tab, child:[]});
		}
		
		var j2 = 0; while(j2 < tree[j].child.length && tab2 != tree[j].child[j2].name) j2++;
		if(j2 == tree[j].child.length ){
			tree[j].child.push({ name:tab2, child2:[], frac:0, y:0});
		}

		var j3 = 0; while(j3 < tree[j].child[j2].child2.length && tab3 != tree[j].child[j2].child2[j3].name) j3++;
		if(j3 == tree[j].child[j2].child2.length){
			tree[j].child[j2].child2.push({ name:tab3, plot:i});
		}
		else emsg("Problem");
	}
	
	for(i = 0; i < tree.length; i++){ 
		pagesub[i] = 0; pagesubsub[i]=[]; 
		for(j = 0; j < tree[i].child.length; j++){
			pagesubsub[i][j] = 0;
			tree[i].child[j].frac = submax/tree[i].child[j].child2.length;
		}
	}
}


/// This sets up listeners for the mouse and keys
function init_mouse_key()
{
	a = ById("main");
	
	a.addEventListener('mousemove', function(evt) {
		var mousePos = getMousePos(myc, evt);
		mousemove(mousePos.x,mousePos.y);
	}, false);

	a.addEventListener('mousedown', function(evt) {
		var mousePos = getMousePos(myc, evt);
		mdown(mousePos.x,mousePos.y); 
	}, false);

	a.addEventListener ("mouseout", function(evt) {
		drag = 0; 
		buttoninit();
	}, false);
	
	a.addEventListener('mouseup', function(evt) {
		var mousePos = getMousePos(myc, evt); 
		ctrlkey = evt.ctrlKey;
		if(evt.altKey) location.reload(true);
		mup(mousePos.x,mousePos.y);
	}, false);

	document.onkeydown = checkKey;
}


/// Determines what is done when a key is pressed
function checkKey(e) 
{
  e = e || window.event;

	if (e.keyCode == '37'){  // Left arrow
		playtime--; if(playtime < 0) playtime = 0;
	}
	if (e.keyCode == '39'){  // Right arrow
		playtime++; if(playtime > playtimemax-1) playtime = playtimemax-1;
	}	
	
	var num = tree[page].child[pagesub[page]].child2.length;
	var val = pagesubsub[page][pagesub[page]];
	
	if (e.keyCode == '38'){  // Up arrow
		if(val > 0) changepage(-1,-1,val-1);
	}
	if (e.keyCode == '40'){  // Down arrow
		if(val < num-1) changepage(-1,-1,val+1);
	}
	
	plot_results();
	buttonplot();
}


/// Dynamically sets the size of page objects (activates when the window is resized)
function setsize()                                         
{
	if(initdone == 0) return;
	
	var w = window.innerWidth;
	var h = window.innerHeight;

	height = h-50; if(height < 400) height = 400;
	width = Math.floor(menux+(height*aspect_ratio));
	
	if(width > w){
		height = Math.floor(height*(w-40)/width); if(height < 400) height = 400;
		width = Math.floor(menux+(height*aspect_ratio));
	}
	
	mar = Math.floor((w-width)/2); if(mar < 0) mar = 0;
	
	ById("main").style.marginLeft = mar+"px";
	ById("main").style.width = width+"px";
	
	graphcan.width = width;
  graphcan.height = height;
 
  resultcan.width = width;
  resultcan.height = height;
	
	myc.width = width;
  myc.height = height;

  canw = width; canh = height;

	resdx = Math.floor((width - menux)*0.6);
	resdy = height-20; 
	
	mapdx = Math.floor(height*map_ratio);
	mapdy = height;
	
	buttoninit();
	widthold = width; heightold = height;
}


/// Gets the mouse position
function getMousePos(canvas, evt)   
{
	var rect = canvas.getBoundingClientRect();
	return {
		x: evt.clientX - rect.left,
		y: evt.clientY - rect.top
	};
}


/// Used as a shortcut to access document elements
function ById(a){ return document.getElementById(a);}     // Gets an element in DOM


/// Fires if mouse button is clicked down
function mdown(xx,yy) 
{
   var d = new Date(); timedown = d.getTime();
	
	if(timedown-timeup < 300 && timeclick < 300){ mousedblclick(); timedown = 0; return;}
	
	drag = 0; 
	
	if(plottype == "map" && canbut != -1 && over != -1 && buttype[over] == CANVASBUT){
		drag = 1; mxst = mx; myst = my;
	}
	
	switch(buttype[over]){
	case MENUSLIDEBUT:
		var ob = tree[page].child[pagesub[page]];

		if(my >= sliy1 && my <= sliy2){ mxst = mx; myst = my; menuslidest = ob.y; menuslidesize = butdy[over]; drag = 2;}
		else{
			if(my > sliy2){ ob.y += 0.9*ob.frac; if(ob.y > 1-ob.frac) ob.y = 1-ob.frac; buttoninit();}
			else{
				if(my < sliy1){ ob.y -= 0.9*ob.frac; if(ob.y < 0) ob.y = 0; buttoninit();}
			}
		}
		break;	
		
	case SLIDEBUT:
		if(my >= slidey1 && my <= slidey2){ mxst = mx; myst = my; menuslidest = slidey; menuslidesize = butdy[over]; drag = 3;}
		else{
			if(my > slidey2){ slidey += 0.9*slidefrac; if(slidey > 1-slidefrac) slidey = 1-slidefrac; buttoninit();}
			else{
				if(my < slidey1){ slidey -= 0.9*slidefrac; if(slidey < 0) slidey = 0; buttoninit();}
			}
		}
		break;	
		
	case TABYSLIDEBUT:
		if(my >= taby1 && my <= taby2){ mxst = mx; myst = my; menuslidest = taby_slide; menuslidesize = butdy[over]; drag = 4;}
		else{
			if(my > taby2){ taby_slide += 0.9*taby_slidefrac; if(taby_slide > 1-taby_slidefrac) taby_slide = 1-taby_slidefrac; buttoninit();}
			else{
				if(my < taby1){ taby_slide -= 0.9*taby_slidefrac; if(taby_slide < 0) taby_slide = 0; buttoninit();}
			}
		}
		break;	
		
	case TABXSLIDEBUT:
		if(mx >= tabx1 && mx <= tabx2){ mxst = mx; myst = my; menuslidest = tabx_slide; menuslidesize = butdx[over]; drag = 5;}
		else{
			if(mx > tabx2){ tabx_slide += 0.9*tabx_slidefrac; if(tabx_slide > 1-tabx_slidefrac) tabx_slide = 1-tabx_slidefrac; buttoninit();}
			else{
				if(mx < tabx1){ tabx_slide -= 0.9*tabx_slidefrac; if(tabx_slide < 0) tabx_slide = 0; buttoninit();}
			}
		}
		break;	
	}
}


/// Fires when mouse button is released
function mup(xx,yy) 
{
	timeup = (new Date()).getTime();
	timeclick  = timeup-timedown; 
	if(drag != 0){ drag = 0; if(timeclick > 200) buttoninit();}
	if(timeclick < 500){ mouseclick(xx,yy);}
}


/// Prints to the console
function pr(te)
{ 
	console.log(te);
}          
