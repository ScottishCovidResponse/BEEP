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
			if(tree[page].child[pagesub[page]].child2[0].name != "") arrowdown = 1; 
			
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
		off += 2;
		
		plottext(text,x+1,y+16+off-2,font,col,dx);
		break;
		
	case MENUSLIDEBUT:
		var ob = tree[page].child[pagesub[page]];
		
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
			xx = Math.floor(x + dx*(Math.log(ticky[i]+col_map.shift)-col_map.min)/(col_map.max-col_map.min));
			drawline(xx, y+dy, xx, y+dy+5,BLACK,2,0); 
			centertext(ticky[i],xx,y+dy+14,"11px arial",BLACK); 
		}
		break;
	
	case SOURCEBUT:
		fillrect(x,y,dx,dy,WHITE);
		col = BLACK; if(over == i) col = GREY;
		
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
		drawcorners(menux,0,canw-menux,dy,30,"#ddddff","#9999ff")
		break;
		
	case LOGOBUT:
		cv.drawImage(logopic,x,y);
		break;
		
	case CANVASBUT:
		canbuttonplot();
		fillrect(x,y,dx,dy,WHITE); 
	
		cv.drawImage(graphcan,0,0,dx,dy,x,y,dx,dy);
		if(dy == height) drawcorners(menux,0,canw-menux,dy,30,"#ddddff","#9999ff");
		canvasdx = dx; canvasdy = dy;
		break;
	
	case NORMTABLEHEADBUT:
		plottext(text,x+4,y+17,NORMTABLEHEADFONT,BLACK);
		break;
		
	case NORMTABLEBUT:
		plottext(text,x+4,y+17,NORMTABLEFONT,BLACK,dx);
		break;
		
	default: alertp("Error code EC1: "+buttype[i]); break;
  }
}
