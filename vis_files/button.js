/// Initialises the buttons used on pages
function buttoninit()
{
	over = -1; canover = -1;
	nbut = 0; ncanbut = 0;
	
	var menubut = set_menubut();
	
	setsize(menubut);
	
	generate_plots();

	addbutton("",0,0,menux,heightpage,-1,MENUBACKBUT,-1,-1);
	addbutton("",menux,canh,width-menux,heightpage-canh,-1,MENUBACKBUT,-1,-1);

	addbutton("",menux,0,width-menux,canh,-1,CORNERSBUT,-1,-1);
	
	addbutton("",menux,canh,width-menux,25,-1,LOWERBARBUT,-1,-1);
	
  addbutton("",15,20,0,0,-1,LOGOBUT,-1,-1);
	
	add_menubut(menubut);
	
	buttonplot();
	
	if(drag == 0) mousemove(mx,my);
}


/// Sets the menu buttons
function set_menubut()
{
	var menubut=[];
	var y = 80; 
  var xx = 15, xx2 = 28; xx3 = 47;
	var ddy = 34, dby = 8;
	var dysub = 23, dysubsub = 19, dysubsubsub = 19;
	for(var i = 0; i < tree.length; i++){
		menubut.push({type:"tab", name:tree[i].name, x:0, y:y, dx:menux, dy: ddy, i:i}); y+= ddy;
		if(page == i){
			for(j = 0; j < tree[i].child.length; j++){
				var name = tree[i].child[j].name;
				if(name != ""){
					ii = i; 
					menubut.push({type:"sub", name:name, x:xx, y:y, dx:menux-xx, dy: dysub, j:j});
					y += dysub;
			
					if(pagesub[page] == j){
						for(jj = 0; jj < tree[i].child[j].child.length; jj++){
							var name = tree[i].child[j].child[jj].name;
							if(name != ""){
								menubut.push({type:"subsub", name:name, x:xx2, y:y, dx:menux-xx2, dy: dysubsub, j:jj});
								y += dysubsub;
							
								if(pagesubsub[page][ps] == jj){
									var len = tree[i].child[j].child[jj].child.length;
									if(len > submax){
										var yst = y;
										var ob = tree[page].child[ps].child[pss];
										var j2min = Math.floor(len*ob.y*1.001);
										if(j2min > len-submax) j2min = len-submax;
								
										for(j2 = j2min; j2 < j2min+submax; j2++){
											name2 = tree[i].child[j].child[jj].child[j2].name;
											if(name2 != ""){
												menubut.push({type:"subsubsub", name:name2, x:xx3, y:y, dx:menux-xx3, dy: dysubsubsub, j:j2});
												y += dysubsubsub;
											}
										}
										menubut.push({type:"menuslide", name:"", x:xx+4+5, y:yst+4, dx:10, dy: submax*dysubsubsub-8});	
									}
									else{
										for(j2 = 0; j2 < tree[i].child[j].child[jj].child.length; j2++){
											name2 = tree[i].child[j].child[jj].child[j2].name;
											if(name2 != ""){
												menubut.push({type:"subsubsub", name:name2, x:xx3, y:y, dx:menux-xx3, dy: dysubsubsub, j:j2});
												y += dysubsubsub;
											}
										}
									}
								}
							}
						}
					}
				}
			}
		}
		y += dby;
	}
	
	return menubut;
}


/// Adds the menu buttons to the page
function add_menubut(menubut)
{
	for(var i = 0; i < menubut.length; i++){
		var mb = menubut[i];
		switch(mb.type){
			case "tab":
				addbutton(mb.name,mb.x,mb.y,mb.dx,mb.dy,TABBUT,TABBUT,mb.i,-1);
				break;
				
			case "sub":
				addbutton("",mb.x,mb.y,mb.dx,mb.dy,-1,PAGESUBBACKBUT,-1,-1);
				addbutton(mb.name,mb.x,mb.y,mb.dx,mb.dy,PAGESUBBUT,PAGESUBBUT,mb.j,-1);
				break;
				
			case "subsub":
				addbutton(mb.name,mb.x,mb.y,mb.dx,mb.dy,PAGESUBSUBBUT,PAGESUBSUBBUT,mb.j,-1);
				break;
				
			case "subsubsub":
				addbutton(mb.name,mb.x,mb.y,mb.dx,mb.dy,PAGESUBSUBSUBBUT,PAGESUBSUBSUBBUT,mb.j,-1);
				break;
				
			case "menuslide":
				addbutton(mb.name,mb.x,mb.y,mb.dx,mb.dy,MENUSLIDEBUT,MENUSLIDEBUT,-1,-1);
				break;
		}
	}
}

	
/// Fires when a page is changed
function changepage(pagenew, s, ss, sss)
{
	var pagest=page;
	
	playtime = 0; playing = 0;
	sourceon = 0; slidey = 0;
	tabx_slide = 0; taby_slide = 0;
	selparam_table = null; selparam_name = "";
	seleq = null; seleq_name = ""; 
	paramsel = ""; check = 0;
	
	if(pagenew != -1) page = pagenew;
	var p = page;
	ps = pagesub[p]; if(s != -1){ ps = s; pagesub[page] = ps;}
	pss = pagesubsub[p][ps]; if(ss != -1){ pss = ss; pagesubsub[p][ps] = pss;}
	psss = pagesubsubsub[p][ps][pss]; if(sss != -1){ psss = sss; pagesubsubsub[p][ps][pss] = psss;}
	
	var ob = tree[p].child[ps].child[pss];
	var len = ob.child.length;
	if(len > submax){
		var j2min = Math.floor(len*ob.y*1.001);
		if(psss < j2min) ob.y = psss/len;
		if(psss >= j2min+submax-1) ob.y = (psss-submax+1+0.1)/len;
	}
	
	buttoninit();
}


/// Adds a button onto a page
function addbutton(text,x,y,dx,dy,ac,type,val,val2)       
{
	buttext[nbut] = text;
	butx[nbut] = Math.floor(x);
	buty[nbut] = Math.floor(y);
	butdx[nbut] = Math.floor(dx);
	butdy[nbut] = Math.floor(dy);
	butac[nbut] = ac; 	
	buttype[nbut] = type;
	butover[nbut] = -1;
	butval[nbut] = val;
	butval2[nbut] = val2;
	nbut++;
}


/// Returns the font size
function get_fontsize(font)
{
	var pos2 = font.search("px");
	return parseInt(font.substr(pos2-2,2));
}


/// Resizes the font
function resize_font(font,sinew)
{
	var pos2 = font.search("px");
	var si = parseInt(font.substr(pos2-2,2));
	return font.substr(0,pos2-2)+sinew+font.substring(pos2);
}


/// Splits a string into substrings to accout for subscripts
function splitintosub(st,font)                             
{
	var lab=[];

	var si = get_fontsize(font);
	var fsi = Math.floor(0.7*si);
	
	var font2 = resize_font(font,fsi);
	
	var italicfont = "italic "+si+"px Times";
	var italicfont2 = "italic "+fsi+"px Times";
	
	if(font.includes("bold") == true){ italicfont = "bold "+italicfont; italicfont2 = "bold "+italicfont2;}
	 
	var j = 0, w = 0, italic = 0;
	do{
		if(st.substr(j,1) != "_" && st.substr(j,1) != "^" && st.substr(j,1) != "&"){
			jst = j;
			while(j < st.length && st.substr(j,1) != "_" && st.substr(j,1) != "^" && st.substr(j,1) != "&") j++;
			frag = st.substr(jst,j-jst);

			if(italic == 0) fo = font; else fo = italicfont;
			if(frag == "Σ") fo = resize_font(font,Math.floor(1.2*si));
			dw = textwidth_simp(frag,fo)+1;
		
			lab.push({frag:frag, type:"norm", size:"big", w:w, dw:dw, font:fo, dy:0, italic:italic});
			w += dw;
		}
		else{
			if(st.substr(j,1) == "&"){
				italic = 1-italic;
				j++;
			}
			else{
				var type;
				if(st.substr(j,1) == "_"){ type = "sub"; dy = 0.2*si;}
				else{ type = "sup"; dy = -0.4*si;}
				j++; 
				
				if(st.substr(j,1) == "{"){
					j++;
					jst = j;
					while(j < st.length && st.substr(j,1) != "}") j++;
					if(j ==  st.length) alertp("Bracket problem");
					frag = st.substr(jst,j-jst);
					j++;
				}
				else{
					jst = j;
					while(j < st.length && st.substr(j,1) != " " &&  st.substr(j,1) != "=" && st.substr(j,1) != "]" && st.substr(j,1) != ")" && st.substr(j,1) != "}" && st.substr(j,1) != "," && st.substr(j,1) != "." && st.substr(j,1) != "_" && st.substr(j,1) != "^" && st.substr(j,1) != "&") j++;
					frag = st.substr(jst,j-jst);
				}
			
				if(frag != ""){
					dw = textwidth_simp(frag,font2);
					
					if(italic == 0) fo = font2; else fo = italicfont2;
					lab.push({frag:frag, type:type, w:w, dw:dw, font:fo, dy:dy, italic:italic});
					w += dw;
				}
			}
		}
	}while(j < st.length);

	
	var i;
	for(i = 0; i < lab.length-1; i++){
		if((lab[i].type == "sub" && lab[i+1].type == "sup") || (lab[i].type == "sup" && lab[i+1].type == "sub")){
			var dw = lab[i].dw, dw2 = lab[i+1].dw;
			if(dw2 > dw){
				for(j = i+1; j < lab.length; j++) lab[j].w -= dw;
			}
			else{
				lab[i+1].w -= dw;
				for(j = i+2; j < lab.length; j++) lab[j].w -= dw2;
			}
		}
	}
	
	for(i = 1; i < lab.length; i++){
		if(lab[i].type == "sub"){
			if(lab[i-1].type == "norm" && lab[i-1].italic == 1){		
				dw = 2;
				for(j = i; j < lab.length; j++) lab[j].w -= dw;
			}
		}
	}
	
	return lab;
}


/// Plots left-aligned text
function plottext(text,x,y,font,col,width)               
{
	if(text == undefined){ text = "undefined"; pr("undefined text");}
	if(text.length == 0) return;
	if(!width) width = large;
	
	text = text.toString();
	if(text.search("_") == -1 && text.search("^") == -1 && text.search("&") == -1){   // looks for subscripts
		cv.font = font;
		if(width){
			if(cv.measureText(text).width > width-5){
				while(cv.measureText(text+"...").width > width-5){
					text = text.substr(0,text.length-1);
				}
				text += "...";
			}	
		}
		cv.textAlign = 'left';
		cv.fillStyle = col;
		cv.fillText(text, x, y);
	}
	else{	
		if(text == "[" || text == "]" || text == "(" || text == ")"){
			var fac;
			if(text == "[" || text == "]"){ fac = 1.6; y += 3;} else{ fac = 1.3; y+=2;}
			var si = get_fontsize(font);
			cv.font = resize_font(font,Math.floor(si*fac));
			cv.textAlign = 'left';
			cv.fillStyle = col;
			cv.fillText(text, x, y);
		}
		else{
			var lab = splitintosub(text,font);

			cv.textAlign = 'left';
			cv.fillStyle = col;
			var j;
			for(j = 0; j < lab.length; j++){
				cv.font = lab[j].font;
				if(lab[j].w+lab[j].dw > width-5){
					te = lab[j].frag;

					while(lab[j].w+cv.measureText(te+"...").width > width-5 && te.length > 0){
						te = te.substr(0,te.length-1);
					}
					te += "...";

					cv.fillText(te,x+lab[j].w,y+lab[j].dy);
					break;
				}
				cv.fillText(lab[j].frag,x+lab[j].w,y+lab[j].dy);
			}
		}
	}
}


/// Plots left-aligned vertical text
function plotverticaltext(text,x,y,font,col,width)               
{	
	if(text == undefined){ text = "undefined"; pr("undefined text");}
	if(text.length == 0) return;
	if(!width) width = large;
	
	text = text.toString();
	
	var lab = splitintosub(text,font);
	
	cv.save();
	cv.translate(x, y);
	cv.rotate(Math.PI/2);
	cv.textAlign = 'left';
	cv.fillStyle = col;
	
	var j;
	for(j = 0; j < lab.length; j++){
		cv.font = lab[j].font;
		cv.fillText(lab[j].frag,lab[j].w,lab[j].dy);
	}
	cv.restore();
}


/// Plots text at an angle
function plotangletext(text,x,y,th,font,col,width) 
{
	cv.font = font;
	if(width){
		if(cv.measureText(text).width > width-5){
			while(cv.measureText(text+"...").width > width-5) text = text.substr(0,text.length-1);
			text += "...";
		}	
	}
	
	cv.save();
	cv.translate(x, y);
	cv.rotate(-th);
	cv.textAlign = 'left';
	cv.fillStyle = col;
	cv.fillText(text, 0, 0);
	cv.restore();
}


/// Plots centered text at an angle
function centerplotangletext(text,x,y,th,font,col,width) 
{
	text = text.toString();
	if(text.search("_") == -1 && text.search("^") == -1 && text.search("&") == -1){   // looks for subscripts
		cv.font = font;
		if(width){
			if(cv.measureText(text).width > width-5){
				while(cv.measureText(text+"...").width > width-5) text = text.substr(0,text.length-1);
				text += "...";
			}	
		}
		
		cv.save();
		cv.translate(x, y);
		cv.rotate(-th);
		cv.textAlign = 'center';
		cv.fillStyle = col;
		cv.fillText(text, 0, 0);
		cv.restore();
	}
	else{
		lab = splitintosub(text,font);
		w = lab[lab.length-1].w+lab[lab.length-1].dw;
		cv.textAlign = 'left';
		cv.fillStyle = col;
		cv.save();
		cv.translate(x, y);
		cv.rotate(-th);	
		for(j = 0; j < lab.length; j++){
			cv.font = lab[j].font;
			cv.fillText(lab[j].frag, lab[j].w-w/2,lab[j].dy);
		}
		cv.restore();
	}
}


/// Plots vertical text
function verticaltext(text,x,y,font,col)  
{
	text = text.toString();
	cv.font = font;	
	cv.save();
	cv.translate(x, y);
	cv.rotate(Math.PI/2);
	cv.textAlign = 'left';
	cv.fillStyle = col;
	cv.fillText(text, 0, 0);
	cv.restore();
}


/// Plots centered text in an area
function centertext(text,x,y,font,col,width)        
{
	if(text == undefined) text = "undefined";
	
	text = text.toString();
	if(text.search("_") == -1 && text.search("^") == -1 && text.search("&") == -1){   // looks for subscripts
		cv.font = font;

		if(width){
			if(cv.measureText(text).width > width-10){
				while(cv.measureText(text+"...").width > width-10) text = text.substr(0,text.length-1);
				text += "...";
			}	
		}
		
		cv.textAlign = 'center';
		cv.fillStyle = col;
		cv.fillText(text, x, y);
	}
	else{
		var lab = splitintosub(text,font);
		var w = lab[lab.length-1].w+lab[lab.length-1].dw;
		cv.textAlign = 'left';
		cv.fillStyle = col;
		var j;
		for(j = 0; j < lab.length; j++){
			cv.font = lab[j].font;
			cv.fillText(lab[j].frag,x+lab[j].w-w/2,y+lab[j].dy);
		}
	}
}


/// Right-aligns text
function righttext(text,x,y,font,col,width)
{
	cv.font = font;
	if(width){
		if(cv.measureText(text).width > width-5){
			while(cv.measureText(text+"...").width > width-5) text = text.substr(0,text.length-1);
			text += "...";
		}	
	}
	cv.textAlign = 'right';
	cv.fillStyle = col;
	cv.fillText(text, x, y);
}

function textwidth_simp(text,font)  
{
	cv.font = font;
	return cv.measureText(text).width;
}


/// Returns the width of specified text
function textwidth(text,font)
{ 
	text = text.toString();
	if(text.search("_") == -1 && text.search("^") == -1 && text.search("&") == -1){ 
		cv.font = font;
		return cv.measureText(text).width;
	}
	else{
		var lab = splitintosub(text,font);
		var j = lab.length-1;
		return lab[j].w+lab[j].dw;
	}
}


/// Draws a rectangle
function drawrect(x1, y1, x2, y2, col, style)
{
	cv.beginPath();
	cv.lineWidth = style;
	cv.rect(x1,y1,x2,y2);
	cv.strokeStyle = col;
	cv.stroke();
}


/// Draws reload symbol
function reloadsign(x,y,col)
{
	var dx = 20, dy = 20;
	cv.lineWidth = 3; 
	cv.beginPath();
	cv.arc(x+dx/2, y+dy/2+1, 8, -0.1, 1.5*Math.PI);
	cv.strokeStyle = col;
	cv.stroke(); 
	drawarrow(x+dx/2+7,y+3,x-100,y+6,8,col);
	centertext("Reload",x+dx/2,y+dy+9,"bold 10px georgia",col);
}


/// Draws a filled rectangle
function fillrect(x1, y1, x2, y2, col)          
{
	cv.beginPath();
	cv.rect(x1,y1,x2,y2);
	cv.fillStyle = col;
	cv.fill();
}


/// Sets the dash type
function setdash(dash)                                  
{	
	var sm = 3, md = 9, bi = 14;
	
	switch(dash%10){
	case 0: cv.setLineDash([]); break;
	case 1: cv.setLineDash([md, sm]); break;
	case 2: cv.setLineDash([sm, sm]); break;
	case 3: cv.setLineDash([md, sm, sm, sm]); break;
	case 4: cv.setLineDash([bi, md]); break;
	case 5: cv.setLineDash([bi, sm, sm, sm]); break;
	case 6: cv.setLineDash([bi, sm, md, sm]); break;
	case 7: cv.setLineDash([bi, sm, sm, sm, sm, sm]); break;
	case 8: cv.setLineDash([bi, sm, md, sm, sm, sm]); break;
	case 9: cv.setLineDash([bi, sm, md, sm, md, sm]); break;
	}
}


/// Draws a line
function drawline(x1, y1, x2, y2, col, style, dash)      
{
	cv.lineWidth = style;

	if(dash) setdash(dash);

	cv.beginPath();
	cv.moveTo(x1,y1);
	cv.lineTo(x2,y2);
	cv.strokeStyle = col;
	cv.stroke();
	
	if(dash) cv.setLineDash([]);
}


/// Draws a polygon
function drawpolygon(npoint,col,col2, style)  
{
	var i;

	cv.lineWidth = style;

	cv.beginPath();
	cv.moveTo(polypoint[0][0],polypoint[0][1]);
	for(i = 1; i < npoint; i++) cv.lineTo(polypoint[i][0],polypoint[i][1]);
	cv.closePath();
	cv.fillStyle = col;
	cv.fill();
	cv.strokeStyle = col2;
	cv.stroke();
}


/// Draws a rectangle with rounded corners
function drawroundrect(x,y,dx,dy,r,col,col2)            
{
	var th, i, nth = Math.floor(r/2);
	if(nth < 1){
		cv.lineWidth = 1;	
		cv.beginPath();
		cv.beginPath();
		cv.moveTo(x+r,y);
		cv.lineTo(x+dx-r,y);
		cv.lineTo(x+dx,y+r);
		cv.lineTo(x+dx,y+dy-r);
		cv.lineTo(x+dx-r,y+dy);
		cv.lineTo(x+r,y+dy);
		cv.lineTo(x,y+dy-r);
		cv.lineTo(x,y+r);
		cv.closePath();
		cv.fillStyle = col; 
		cv.fill();
		cv.strokeStyle = col2;
		cv.stroke();
		return;
	}
	
	cv.lineWidth = 1;
	cv.beginPath();
 
	cv.moveTo(x+r,y);

	for(i = 0; i <= nth; i++){
		th = i*Math.PI/(2*nth);
		cv.lineTo(x+dx-r+r*Math.sin(th),y+r-r*Math.cos(th));
	}

	for(i = 0; i <= nth; i++){
		th = i*Math.PI/(2*nth);
		cv.lineTo(x+dx-r+r*Math.cos(th),y+dy-r+r*Math.sin(th));
	}

	for(i = 0; i <= nth; i++){
		th = i*Math.PI/(2*nth);
		cv.lineTo(x+r-r*Math.sin(th),y+dy-r+r*Math.cos(th));
	}

	for(i = 0; i <= nth; i++){
		th = i*Math.PI/(2*nth);
		cv.lineTo(x+r-r*Math.cos(th),y+r-r*Math.sin(th));
	}

	cv.closePath();
	cv.fillStyle = col; 
	cv.fill();
	cv.strokeStyle = col2;
	cv.stroke();
}

function roundrect_intersect(x1,y1,x2,y2,x,y,dx,dy,r)            
{
	var th, i, nth = Math.floor(r/2);

	var line = [];
	line.push({x:x+r,y:y});
	for(i = 0; i <= nth; i++){
		th = i*Math.PI/(2*nth);
		line.push({x:x+dx-r+r*Math.sin(th),y:y+r-r*Math.cos(th)});
	}

	for(i = 0; i <= nth; i++){
		th = i*Math.PI/(2*nth);
		line.push({x:x+dx-r+r*Math.cos(th),y:y+dy-r+r*Math.sin(th)});
	}

	for(i = 0; i <= nth; i++){
		th = i*Math.PI/(2*nth);
		line.push({x:x+r-r*Math.sin(th),y:y+dy-r+r*Math.cos(th)});
	}

	for(i = 0; i <= nth; i++){
		th = i*Math.PI/(2*nth);
		line.push({x:x+r-r*Math.cos(th),y:y+r-r*Math.sin(th)});
	}

	var ax = x1, ay = y1, nax = x2-x1, nay = y2-y1;
	
	var jmax = line.length;
	for(j = 0; j < jmax; j++){
		var bx = line[j].x, by = line[j].y, nbx = line[(j+1)%jmax].x-bx, nby = line[(j+1)%jmax].y-by;

		var beta = ((ax-bx)*nay - (ay-by)*nax)/(nbx*nay - nby*nax + tiny);
		if(beta > 0 && beta < 1){
			var alpha = ((bx-ax)*nby - (by-ay)*nbx)/(nax*nby - nay*nbx + tiny);
			if(alpha > 0){
				return {x:ax + alpha*nax, y: ay + alpha*nay};
			}
		}			
	}
	alertp("Cannot find interesection");
}



/// Draws a circle
function circle(x,y,r,col,style)                       
{
	cv.lineWidth = style;
	cv.beginPath();
	cv.arc(x,y,r,0,2*Math.PI);
	cv.strokeStyle = col;
	cv.stroke();
}


/// Draws a filled circle
function fillcircle(x,y,r,col,col2,style)                
{
	cv.lineWidth = style;
	cv.beginPath();
	cv.arc(x,y,r,0,2*Math.PI);
	cv.fillStyle = col;
	cv.fill();
  
	cv.strokeStyle = col2;
	cv.stroke();
}


/// Draws the corners of the main page
function drawcorners(x,y,dx,dy,r,col) 
{
	var th, i, nth = Math.floor(r/3);
	if(nth < 1) nth = 1;

	cv.lineWidth = 2;
	cv.beginPath();
	cv.fillStyle = col; 

	cv.moveTo(x+dx,y);
	for(i = 0; i <= nth; i++){
		th = i*Math.PI/(2*nth);
		cv.lineTo(x+dx-r+r*Math.sin(th),y+r-r*Math.cos(th));
	} 
	cv.closePath();
	cv.fill();

	/*
	cv.moveTo(x+dx,y+dy);
	for(i = 0; i <= nth; i++){
		th = i*Math.PI/(2*nth);
		cv.lineTo(x+dx-r+r*Math.cos(th),y+dy-r+r*Math.sin(th));
	}
	cv.closePath();
	cv.fill();

	/*
	cv.moveTo(x,y+dy);
	for(i = 0; i <= nth; i++){
		th = i*Math.PI/(2*nth);
		cv.lineTo(x+r-r*Math.sin(th),y+dy-r+r*Math.cos(th));
	}
	cv.closePath();
	cv.fill();
*/

	cv.moveTo(x,y);
	for(i = 0; i <= nth; i++){
		th = i*Math.PI/(2*nth);
		cv.lineTo(x+r-r*Math.cos(th),y+r-r*Math.sin(th));
	}
	cv.closePath();
	cv.fill();
	
}

/// Draws the corners of the main page
function drawlowerbar(x,y,dx,dy,r,col) 
{
	fillrect(x,y,dx,dy,"#a4a3e6");
	centertext("© BEEPmbp 2022",x+dx/2,y+18,"bold 16px arial",WHITE);   
	
	var th, i, nth = Math.floor(r/3);
	if(nth < 1) nth = 1;

	cv.lineWidth = 2;
	cv.beginPath();
	cv.fillStyle = col; 

	
	cv.moveTo(x+dx,y+dy);
	for(i = 0; i <= nth; i++){
		th = i*Math.PI/(2*nth);
		cv.lineTo(x+dx-r+r*Math.cos(th),y+dy-r+r*Math.sin(th));
	}
	cv.closePath();
	cv.fill();

	cv.moveTo(x,y+dy);
	for(i = 0; i <= nth; i++){
		th = i*Math.PI/(2*nth);
		cv.lineTo(x+r-r*Math.sin(th),y+dy-r+r*Math.cos(th));
	}
	cv.closePath();
	cv.fill();	
}


/// Aligns lines of a paragraph within a certain width
function alignparagraph(text,dx,fontset) 
{
	var i, ist, len, yy, font, di;

	dx -= 8;
	if(fontset) font = fontset;
	else font = HELPFONT;

	text = ""+text;

	var textsplit = text.split("  ");
	
	nlines = 0;
	for(var spl = 0; spl < textsplit.length; spl++){
		var te = textsplit[spl];
		
		i = 0; len = te.length;
		di = Math.floor(len*dx/textwidth_simp(te,font));

		yy = 0;
		while(i < len){
			ist = i;
			i += di; if(i > len) i = len; //while(i < len && text.substr(i,1) != " ") i++;
			while(i < len && textwidth_simp(te.substr(ist,i-ist),font) < dx) i++;
			 
			if(textwidth_simp(te.substr(ist,i-ist),font) >= dx){
				do{
					i--; while(te.substr(i,1) != " " && i > ist) i--;
				}while(textwidth_simp(te.substr(ist,i-ist),font) > dx);

				if(i == ist){
					i = ist + di;
					while(textwidth_simp(te.substr(ist,i-ist),font) < dx) i++;
					while(textwidth_simp(te.substr(ist,i-ist),font) > dx) i--;
				}
				else i++;
			}
			
			lines[nlines] = te.substr(ist,i-ist);
			while(i < len && te.substr(i,1) == " ") i++;
			linesy[nlines] = yy;
			yy += lineheight;
			nlines++;
		}
		if(spl < textsplit.length-1){
			lines[nlines] = ""; yy += lineheight; nlines++;
		}
	}
	
	hsto = yy;
}


/// Draws an arrow
function drawarrow(x,y,x2,y2,size,col)                  
{
	var nx, ny, px, py, r;

	nx = x2-x; ny = y2-y;
	r = Math.sqrt(nx*nx+ny*ny);
	if(size > r/5) size = r/5;
	nx *= size/r; ny *= size/r; 
	px = 0.5*ny; py = -0.5*nx;

	polypoint[0][0] = x; polypoint[0][1] = y; 
	polypoint[1][0] = x+Math.round(nx*0.8); polypoint[1][1] = y+Math.round(ny*0.8);
	polypoint[2][0] = x+Math.round(nx+px); polypoint[2][1] = y+Math.round(ny+py);
	drawpolygon(3,col,col,NORMLINE);
	polypoint[2][0] = x+Math.round(nx-px); polypoint[2][1] = y+Math.round(ny-py);
	drawpolygon(3,col,col,NORMLINE);
}


/// Plots the x label on a graph
function plotxlabel(text,x,y,font,col) 
{
	centertext(text,x,y,font,col); 
}


/// Plots the y label on a graph
function plotylabel(text,x,y,font,col)
{
	centerplotangletext(text,x,y,Math.PI/2,font,col);
}


/// Plots an individual button
function indbuttonplot(i)                                 
{
	cv.clearRect (butx[i],buty[i],butdx[i],butdy[i]);
	butplot(i,0);
}


/// Converts decimal to hexidecimal
function hex(c) {                                         
	var hex = (Math.floor(c)).toString(16);
	return hex.length == 1 ? "0" + hex : hex;
}


/// Darkens a colour
function darkcol(col)  
{ 
	if(col == BLACK) return GREY;
	var bigint, r, g, b, frac = 0.7;
	bigint = parseInt(col.substring(1), 16);	
	r = (bigint >> 16) & 255;
	g = (bigint >> 8) & 255;
	b = bigint & 255;
	return "#" + hex(frac*r) + hex(frac*g) + hex(frac*b);
}
