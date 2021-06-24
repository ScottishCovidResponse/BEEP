/// Fires when the mouse moves
function mousemove(x,y)
{
	mx = x; my = y;

	switch(drag){  
	case 1:  // move map
		zoomx += (mxst-mx)/(zoomfac*mapdy); mxst = mx;
		zoomy -= (myst-my)/(zoomfac*mapdy); myst = my;
		plot_results(); buttonplot();
		break;
	
	case 2: // menu slider
		var ob = tree[page].child[pagesub[page]];
		ob.y = menuslidest + (my-myst)/menuslidesize;
		if(ob.y > 1-ob.frac) ob.y = 1-ob.frac; if(ob.y < 0) ob.y = 0;
		buttoninit();
		break;
		
	case 3: // description slider
		slidey = menuslidest + (my-myst)/menuslidesize;
		if(slidey > 1-slidefrac) slidey = 1-slidefrac; if(slidey < 0) slidey = 0;
		plot_results(); buttonplot();
		break;
	
	case 4: // table y slider
		taby_slide = menuslidest + (my-myst)/menuslidesize;
		if(taby_slide > 1-taby_slidefrac) taby_slide = 1-taby_slidefrac; if(taby_slide < 0) taby_slide = 0;
		buttoninit();
		break;
	
	case 5: // table x slider
		tabx_slide = menuslidest + (mx-mxst)/menuslidesize;
		if(tabx_slide > 1-tabx_slidefrac) tabx_slide = 1-tabx_slidefrac; if(tabx_slide < 0) tabx_slide = 0;
		buttoninit();
		break;
	}
	
	var overnew = -1;
	canbut = -1;
	for(i = nbut-1; i >= 0; i--){                   // Determines if mouse over a button
		if(buttype[i] == CANVASBUT) canbut = i;
		if(butac[i] >= 0 && mx >= butx[i] && mx <= butx[i]+butdx[i] && my >= buty[i] && my <= buty[i]+butdy[i]){
			overnew = i; break;
		}
	}

	if(canbut != -1){                               // Determines if mouse over a canvas button
		canovernew = -1;
		if(overnew == over && over == canbut){
			xx =  mx-butx[over]; yy = my-buty[over];
			canmx = xx; canmy = yy;
			for(i = ncanbut-1; i >= 0; i--){
				if(canbutac[i] >= 0 && xx >= canbutx[i] && xx < canbutx[i]+canbutdx[i] &&
										           yy >= canbuty[i] && yy < canbuty[i]+canbutdy[i]){
					canovernew = i;
					break;
				}
			}
		}
	}	
	

	if(canbut != -1){                                 // Replots the button, if neccessary
		if(canovernew != canover && over != -1 && buttype[over] == CANVASBUT){
			canover = canovernew;
			butplot(canbut,-1);
		}
	}

	var areaovernew = -1;                            // Works out if over an area on a map
	if(canbut != -1 && over != -1 && buttype[over] == CANVASBUT && plottype == "map"){
		var mapx = zoomx + (mx-menux)/(mapdy*zoomfac);
		var mapy = zoomy + (mapdy - (my))/(mapdy*zoomfac);
		
		var ra = visjson.click_bound_range;
		for(c = 0; c < ra.length; c++){
			if(mapx > ra[c].xmin && mapx < ra[c].xmax && mapy > ra[c].ymin && mapy < ra[c].ymax){		
				var inter = visjson.inter[c][Math.floor(mapy*clickline)];
				if(inter){
					num = 0; for(k = 0; k < inter.length; k++){ if(inter[k] < mapx) num++;}
					if(num%2 == 1) areaovernew = c;
				}
			}			
		}
	}
	
	if(areaovernew != areaover){                            // Replots the map if the area changes
		areaover = areaovernew;
		plot_results();
		buttonplot();
	}
	
	if(overnew != over){                                    // Replots the buttons, if neccesary
		overold = over; over = overnew;
		if(overold >= 0) butplot(overold,-1);
		if(overnew >= 0) butplot(overnew,-1);
	}

	arrownew = 0;                                           // Sets the cursor (i.e. a hand is used for grabbing)
	if(over == canbut && plottype == "map" && canover == -1 && areaover == -1) arrownew = 1;
	if(drag != 0) arrownew = 1;
	
	if(arrownew != arrow){ 
		arrow = arrownew;
		if(arrow > 0) document.getElementById("bod").style.cursor = "grab"; 
		else document.getElementById("bod").style.cursor = "";
	}
}
