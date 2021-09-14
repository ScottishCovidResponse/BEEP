/// Plots all the buttons on the current page
function buttonplot()                                     
{
	var i, ov; 
	allbutplot = 1;

	cv = maincv;
	cv.clearRect(0,0,canw,canh);
	for(i = 0; i < nbut; i++){
		ov = 0; if(i == over) ov = 1;
		butplot(i,ov);
	}
	allbutplot = 0;
}


// Plots individual button
function butplot(i,ov)   
{ 
	var x, y, dx, dy, col, col2, text, off, w, val, val2, xx, yy, dd;

	x = Math.floor(butx[i]); y = Math.floor(buty[i]); dx = Math.floor(butdx[i]); dy = Math.floor(butdy[i]);
	text = buttext[i]; val = butval[i]; val2 = butval2[i];

	switch(buttype[i]){ 
	case TABBUT:
		arrowdown = 0;
	  if(butval[i] == page){
			fillrect(x,y,dx,dy,WHITE); 
			col = "#3366ff";
			col = "#003366";
			font = MENUFONT;
			off = 0;
			if(tree[page].child.length > 1) arrowdown = 1; 
		}
		else{
			fillrect(x,y,dx,dy,"#ddddff"); 
			font= MENUFONT2;
			off = -1;
			
			col = "#003366";
			if(over==i) col = "#66aaff";
		} 
		
		plottext(text,x+28,y+24+off,font,col);
		
		if(arrowdown == 0){
			ddx = 15; ddx2 = 5; ddy = 5;
			drawline(x+ddx,y+dy/2-ddy,x+ddx+ddx2,y+dy/2,col,THICKLINE)
			drawline(x+ddx,y+dy/2+ddy,x+ddx+ddx2,y+dy/2,col,THICKLINE)
		}
		else{
			ddx = 15; ddx2 = 4; ddy = 4;
			drawline(x+ddx-ddx2,y+dy/2-ddy+1,x+ddx,y+dy/2+ddy+1,col,THICKLINE);
			drawline(x+ddx+ddx2,y+dy/2-ddy+1,x+ddx,y+dy/2+ddy+1,col,THICKLINE);
		}
		break;
		
	case PAGESUBBACKBUT:
		fillrect(x, y, dx, dy, "#ededff");
		break;
		
	case PAGESUBBUT:
		arrowdown = 0;
		if(butval[i] == pagesub[page]){
			if(tree[page].child[pagesub[page]].child[0].name != "") arrowdown = 1; 
			
			if(arrowdown != 1) fillrect(x+dx-8,y+dy/2-6,8,12,WHITE); 
			col = "#cc2222"; 
			font = MENUFONTSMALL;
			off = 0;
		}
		else{
			font = MENUFONTSMALL2;
			off = -1;
			col = "#ff2222"; if(over==i) col = "#dd8888";
		} 

		fillrect(x,y,26+textwidth(text,font), dy, "#ededff");
		off += 1;
		plottext(text,x+24,y+20+off-4,font,col);
		if(arrowdown == 0){
			ddx = 11; ddx2 = 4; ddy = 4;
			drawline(x+ddx,y+dy/2-ddy,x+ddx+ddx2,y+dy/2+1,col,THICKLINE);
			drawline(x+ddx,y+dy/2+ddy,x+ddx+ddx2,y+dy/2+1,col,THICKLINE);
		}
		else{
			ddx = 11; ddx2 = 4; ddy = 3;
			drawline(x+ddx-ddx2,y+dy/2-ddy,x+ddx,y+dy/2+ddy+2,col,THICKLINE);
			drawline(x+ddx+ddx2,y+dy/2-ddy,x+ddx,y+dy/2+ddy+2,col,THICKLINE);
		}
		break;

	case PAGESUBSUBBUT:
		fillrect(x, y, dx, dy, "#ddddff");
		if(butval[i] == pagesubsub[page][pagesub[page]]){
			fillrect(x+dx-8,y+dy/2-6,8,12,WHITE); 
			col = "#222222"; 
			font = MENUFONTSMALL;
			off = 0;
		}
		else{
			font = MENUFONTSMALL2;
			off = -1;
			col = "#222222";
			if(over==i) col = "#888888";
		} 
		
		plottext(text,x+1,y+16+off-2,font,col,dx);
		break;
		
	case PAGESUBSUBSUBBUT:
		fillrect(x, y, dx, dy, "#ddddff");
		if(butval[i] == pagesubsubsub[page][ps][pss]){
			//fillrect(x+dx-8,y+dy/2-6,8,12,WHITE); 
			col = DDBLUE; 
			font = MENUFONTSMALL;
			off = 0;
			fillcircle(x-4,y+10,2,DDBLUE,DDBLUE,1);     
		}
		else{
			font = MENUFONTSMALL2;
			off = -1;
			col = DDBLUE;
			if(over==i) col = "#888888";
		} 
		
		plottext(text,x+1,y+16+off-2,font,col,dx);
		
		break;
		
	case MENUSLIDEBUT:
		var ob = tree[page].child[ps].child[pss];
		
		fillrect(x-1,y-1,dx+2,dy+2,WHITE);
		sliy1 = y + ob.y*dy;
		sliy2 = sliy1 +ob.frac*dy;
	
		col = GREY; if((over == i && my >= sliy1 && my <= sliy2) || drag == 2) col = DDGREY;
		drawrect(x,y,dx,dy,DDGREY,NORMLINE);
		
		fillrect(x,sliy1,dx,sliy2-sliy1,col);
		break;
		
	case SLIDEBUT:
		fillrect(x-1,y-1,dx+2,dy+2,WHITE);
		slidey1 = y + slidey*dy;
		slidey2 = slidey1 +slidefrac*dy;
		
		col = GREY; if((over == i && my >= slidey1 && my <= slidey2) || drag == 3) col = DDGREY;
		drawrect(x,y,dx,dy,DDGREY,NORMLINE);
		
		fillrect(x,slidey1,dx,slidey2-slidey1,col);
		break;
		
	case TABYSLIDEBUT:
		fillrect(x-1,y-1,dx+2,dy+2,WHITE);
		taby1 = y + taby_slide*dy;
		taby2 = taby1 +taby_slidefrac*dy;
		
		col = GREY; if((over == i && my >= taby1 && my <= taby2) || drag == 4) col = DDGREY;
		drawrect(x,y,dx,dy,DDGREY,NORMLINE);
		
		fillrect(x,taby1,dx,taby2-taby1,col);
		break;
		
	case TABXSLIDEBUT:
		fillrect(x-1,y-1,dx+2,dy+2,WHITE);
		tabx1 = x + tabx_slide*dx;
		tabx2 = tabx1 +tabx_slidefrac*dx;

		col = GREY; if((over == i && mx >= tabx1 && mx <= tabx2) || drag == 5) col = DDGREY;
		drawrect(x,y,dx,dy,DDGREY,NORMLINE);
		
		fillrect(tabx1,y,tabx2-tabx1,dy,col);
		break;
		
	case TITLEBUT:
		plottext(text,x+4,y+18,"bold 15px arial",BLACK);	
		break; 
	
	case LEVELKEYBUT:
		var ddx = dy;
		fillrect(x,y,ddx,dy,val);
		drawrect(x,y,ddx,dy,BLACK,1);
		plottext(text,x+ddx+7,y+15,KEYFONT,BLACK);	
		break;
		
	case KEYBUT:
		setstyle(visjson.plots[val].line[val2].style);
		cv.lineWidth = 2;
		cv.beginPath();
		cv.moveTo(x,y+12);
		cv.lineTo(x+50,y+12);
		cv.stroke();
		setdash(0);
		plottext(text,x+60,y+18,KEYFONT,BLACK);
		break; 
	
	case COLOURSCALEBUT: 
		for(i = 0; i < col_map.div; i++){
			var xx = Math.floor(x + (dx*i)/col_map.div);
			var xx2 = Math.floor(x + (dx*(i+1))/col_map.div);
			fillrect(xx,y,xx2-xx,dy,"rgb("+col_map.map[i].R+","+col_map.map[i].G+","+col_map.map[i].B+")");  
		}
		drawrect(x,y,dx,dy,DGREY,1);  

		for(i = 0; i < nticky; i++){
			if(val == "log") xx = Math.floor(x + dx*(Math.log(ticky[i]+col_map.shift)-col_map.min)/(col_map.max-col_map.min));
			else xx = Math.floor(x + dx*(ticky[i]-col_map.min)/(col_map.max-col_map.min));
			drawline(xx, y+dy, xx, y+dy+5,BLACK,2,0); 
			centertext(ticky[i],xx,y+dy+14,"11px arial",BLACK); 
		}
		break;
	
	case SOURCEBUT:
		fillrect(x,y,dx,dy,WHITE);
		col = BLACK; if(over == i) col = GREY;
		x += 4;
		if(sourceon == 0){
			ddx = 2; ddx2 = 3; ddy = 3;
			drawline(x+ddx,y+dy/2-ddy,x+ddx+ddx2,y+dy/2+1,col,THICKLINE);
			drawline(x+ddx,y+dy/2+ddy,x+ddx+ddx2,y+dy/2+1,col,THICKLINE);
		}
		else{
			ddx = 3; ddx2 = 3; ddy = 2;
			drawline(x+ddx-ddx2,y+dy/2-ddy,x+ddx,y+dy/2+ddy+2,col,THICKLINE);
			drawline(x+ddx+ddx2,y+dy/2-ddy,x+ddx,y+dy/2+ddy+2,col,THICKLINE);
		}
	
		plottext(text,x+10,y+14,SOURCEFONT,col);	
		break;
		
	case SOURCEFILEBUT:
		fillrect(x,y,dx,dy,WHITE);
		alignparagraph(text,dx,SOURCEFONT);
		for(j = 0; j < nlines; j++){
			cv.textAlign = 'left';
			cv.style = SOURCEFONT;
			cv.fillStyle = BLACK;
			cv.fillText(lines[j],x+2,y+10+j*13);
		}
		break;

	case SLIDEPARAGRAPHBUT:
		alignparagraph(text,dx);
		nlines_disp = Math.floor(dy/lineheight);
		jmin = Math.floor(slidey*nlines*1.00001);
		jmax = jmin+nlines_disp;
		if(jmax > nlines) jmax = nlines;
		for(j = jmin; j < jmax; j++) plottext(lines[j],x+5,y+18+(j-jmin)*lineheight,HELPFONT,BLACK);
		break;
		
	case MENUBACKBUT:
		fillrect(x,y,dx,dy,"#ddddff");
		break;
	
	case CORNERSBUT:
		drawcorners(x,y,dx,dy,25,"#ddddff");
		break;
		
	case LOWERBARBUT:
		drawlowerbar(x,y,dx,dy,25,"#ddddff");
		break;
		
	case LOGOBUT:
		cv.drawImage(logopic,x,y);
		break;
		
	case CANVASBUT:
		canbuttonplot();
		fillrect(x,y,dx,dy,WHITE); 
	
		cv.drawImage(graphcan,0,0,dx,dy,x,y,dx,dy);
		if(dy == height) drawcorners(menux,0,canw-menux,dy,25,"#ddddff");
		canvasdx = dx; canvasdy = dy;
		break;
	
	case NORMTABLEHEADBUT:
		plottext(text,x+4,y+17,NORMTABLEHEADFONT,BLACK);
		break;
		
	case NORMTABLEBUT:
		var col = BLACK; if(text == paramsel) col =	RED;	
		plottext(text,x+4,y+17,NORMTABLEFONT,col,dx);
		break;
		
	case CHECKBUT:
		fillrect(x,y,dx+150,dy,WHITE);
		drawrect(x+2,y+2,dy-4,dy-4,BLACK,NORMLINE);

		if(check == 1){
			if(over == i) fillrect(x+4,y+4,dy-8,dy-8,GREY);
			else fillrect(x+4,y+4,dy-8,dy-8,BLACK);
		}
		else{
			if(over == i) fillrect(x+4,y+4,dy-8,dy-8,LLGREY);
		}
		
		plottext(text,x+dy+4,y+13,"12px arial",BLACK);
		break;
		
	case MENULINKBUT:
		fillrect(x,y,dx,dy,WHITE);
		col = BLUE; if(over == i){ col = LBLUE;}
		plottext(text,x+7,y+20,"12px arial",col);
		break;
		
	case RADIOBUT:
		fillrect(x,y,dx,dy+4,WHITE);
	
		col=DGREY;
		switch(val2){
			case RADIORATE: selval = rateradio; break;
			case RADIOHIGH: selval = lowerhighlight; break;
		}
		
		if(val == selval) col2 = BLACK;
		else{
			if(over == i) col2 = DGREY;
			else col2 = WHITE
		}

		r = 7;
		drawroundrect(x+1,y+2,2*r,2*r,r,WHITE,col);

		if(val == selval) col3 = BLACK; 
		else{ col3 = GREY;  if(over == i) col3 = DGREY;}
    
		plottext(buttext[i],x+20,y+13,"12px arial",col3,dx-20);
		
		if(col2 != WHITE){
			x += 3; y += 5; r -= 3;
			drawroundrect(x+1,y,2*r,2*r,r,col2,col2);
		}
		break;
		
	default: alertp("Error code EC1: "+buttype[i]); break;
  }
}
