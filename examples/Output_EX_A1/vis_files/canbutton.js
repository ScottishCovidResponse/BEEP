/// Plots all the canvas buttons on the page
function canbuttonplot()                                
{
	var i, ov;

	cv = graphcv;
	cv.clearRect(0,0,graphcan.width,graphcan.height);
	for(i = 0; i < ncanbut; i++){
		ov = 0; if(i == canover) ov = 1;
		canbutplot(i,ov);
	}

	cv = maincv;
}


/// Click on a canvas button
function canfinalaction(i)
{
	val = canbutval[i]; val2 = canbutval2[i]; 
	x = Math.floor(canbutx[i]); y = Math.floor(canbuty[i]); dx = Math.floor(canbutdx[i]); dy = Math.floor(canbutdy[i]);
	text = canbuttext[i]; 

	switch(canbutac[i]){
	case ZOOMINBUT:
	case ZOOMOUTBUT:
		fac = 1.3; if(canbutac[i] == ZOOMOUTBUT) fac = 1.0/fac;
		mapzoom(fac,mapdx/2,mapdy/2);
		break;

	case PLAYBUT: 
		playing = 1-playing;
		if(playing == 1){
			if(playtime == playtimemax-1) playtime = 0;
			playstartt = getsec();
			playtimeinit = playtime;
			playanim();
		}
		break;
		
	case PLAYLINEBUT:
		playtime = Math.floor(playtimemax*(mx-x-menux)/dx); 
		if(playtime == playtimemax) playtime = playtimemax-1;
		playing = 0;
		break;
		
	case VARIABLEBUT:
		if(text == foi.name) changepage(-1,1,-1,-1);
		else{
			selparam_table = val2;
			selparam_name = text;
			taby_slide = 0;
		}
		paramsel = "";
		break;
		
	case EQBUT:
		seleq = val2;
		seleq_name = text;
		paramsel = "";
		taby_slide = 0;
		break;
		
	case LINKBUT:
		follow_param_link(text);
		break;
		
	case PARAMLINEBUT:
		follow_param_link(text);
		break;
		
	default: alertp("Error code EC2 "+canbutac[i]); break;
	}
}


// Returns the time in seconds
function getsec(){ return new Date().getTime() / 1000;}


/// Plots an individual canvas button
function canbutplot(i,ov)                                
{	
	x = Math.floor(canbutx[i]); y = Math.floor(canbuty[i]); 
	dx = Math.floor(canbutdx[i]); dy = Math.floor(canbutdy[i]);
	text = canbuttext[i]; ty = canbuttype[i];
	val = canbutval[i]; val2 = canbutval2[i];
	
	switch(ty){ 
	case RESULTBUT:
		cv.drawImage(resultcan,0,0,dx,dy,x,y,dx,dy);
		drawline(x,y+dy,x,y-10,BLACK,THICKLINE);
		drawarrow(x,y-20,x,y+dy,15,BLACK);
		drawline(x,y+dy,x+dx+10,y+dy,BLACK,THICKLINE);
		drawarrow(x+dx+20,y+dy,x,y+dy,15,BLACK);
		break;
	
	case RESULTBUT2:
		cv.drawImage(resultcan,0,0,dx,dy,x,y,dx,dy);
		break;
		
	case RESULTBUT3:
		cv.drawImage(resultcan,0,0,dx,dy,x,y,dx,dy);
		drawline(x,y+dy,x,y-10,BLACK,THICKLINE);
		drawarrow(x,y-20,x,y+dy,15,BLACK);
		drawline(x,y+dy,x+dx,y+dy,BLACK,THICKLINE);
		break;
		
	case RESULTBUT4:
		cv.drawImage(resultcan,0,0,dx,dy,x,y,dx,dy);
		drawrect(x,y,dx,dy,BLACK,1);
		break;
		
	case PLAYBUT:
		if(ov == 1) col = BLACK; else col = GREY;

		fillcircle(x+dx/2,y+dy/2,dx/2,col,DGREY,NORMLINE);
		if(playing == 0){
			if(playtime < playtimemax-1){
				ddx = 9; ddy = 9;
				polypoint[0][0] = x+ddx; polypoint[0][1] = y+ddy; 
				polypoint[1][0] = x+dx-ddx; polypoint[1][1] = y+dy/2; 
				polypoint[2][0] = x+ddx; polypoint[2][1] = y+dy-ddy; 
				drawpolygon(3,WHITE,WHITE,THICKLINE);
			}
			else{
				cv.lineWidth = 3; 
				cv.beginPath();
				cv.arc(x+dx/2, y+dy/2+1, 8, -0.1, 1.5*Math.PI);
				cv.strokeStyle = WHITE;
				cv.stroke(); 
				drawarrow(x+dx/2+7,y+9,x-100,y+6,8,WHITE);
			}
		}
		else{
			x1 = 8, y1 = 8, dx1 = 5;
			fillrect(x+x1,y+y1,dx1,dy-2*y1,WHITE);
			fillrect(x+dx-x1-dx1,y+y1,dx1,dy-2*y1,WHITE);
		}
		break;
	
	case PLAYLINEBUT:
		pad = 4;
		fillrect(x,y+pad,dx,dy-2*pad,WHITE);
		xx = Math.floor(dx*playtime/(playtimemax-1));
		fillrect(x,y+pad,xx,dy-2*pad,LGREY);
		drawline(x+xx,y+pad,x+xx,y+dy-pad,DGREY,NORMLINE);
		drawrect(x,y+pad,dx,dy-2*pad,BLACK,NORMLINE);
		break;
		
	case ZOOMINBUT:
	case ZOOMOUTBUT:
		cv.globalAlpha = 0.95;
		drawroundrect(x+1,y+1,dx-2,dy-2,4,WHITE,WHITE);
		cv.globalAlpha = 1;
	
		r = 8; r2 = 5; r3 = 23;
		col = BLACK; if(ov == 1) col = DGREY;
		
		xx = x+dx-r; yy = y+r;
		circle(xx,yy,r,col,THICKLINE);
		drawline(xx-r2,yy,xx+r2,yy,col,THICKLINE);
		if(canbuttype[i] == ZOOMINBUT) drawline(xx,yy-r2,xx,yy+r2,col,THICKLINE);
		th = 2.2;
		drawline(xx+r*Math.cos(th),yy+r*Math.sin(th),xx+r3*Math.cos(th),yy+r3*Math.sin(th),col,VTHICKLINE);
		break;
	
	case CATBUT:
		centertext(text,x+dx/2,y+18,CATFONT,col);
		break;
		
	case CATVERTBUT:
		verticaltext(text,x+dx/2-5,y,CATFONT,col)   
		break;
	
	case XTICKBUT:
		col = BLACK; if(ov == 1 && selectbub != XTICKBUT) col = DGREY;
		for(j = 0; j < ntickx; j++){
			xx = Math.floor(x+dx*(tickx[j]-axxmin)/(axxmax-axxmin));
			centertext(tickx[j],xx,y+dy/2+8,TICKFONT,col); 
			drawline(xx,y,xx,y-10,BLACK,NORMLINE); 
			if(val == -2) drawline(xx,y-graphdy,xx,y-graphdy+10,BLACK,NORMLINE); 
		}
		break;
		
	 case YTICKTRBUT:
		col = BLACK; if(ov == 1 && selectbub != YTICKTRBUT) col = DGREY;
		for(j = 0; j < nticky; j++){
			yy = Math.floor(y+dy-dy*(ticky[j]-axymin)/(axymax-axymin));
			righttext(ticky[j],x+dx-5,yy+6,TICKFONT,col); 
			drawline(x+dx,yy,x+dx+10,yy,BLACK,NORMLINE); 
		}
		break;
		
	case LOGXTICKBUT:
		col = BLACK; if(ov == 1 && selectbub != XTICKBUT) col = DGREY;
		for(j = 0; j < ntickx; j++){
			xx = Math.floor(x+dx*(Math.log(tickx[j])-Math.log(axxmin))/(Math.log(axxmax)-Math.log(axxmin)));
			centertext(tickx[j],xx,y+dy/2+8,TICKFONT,col); 
			drawline(xx,y,xx,y-10,BLACK,NORMLINE); 
			if(val == -2) drawline(xx,y-graphdy,xx,y-graphdy+10,BLACK,NORMLINE); 
		}
		break;
		
	 case LOGYTICKTRBUT:
		col = BLACK; if(ov == 1 && selectbub != YTICKTRBUT) col = DGREY;
		for(j = 0; j < nticky; j++){
			yy = Math.floor(y+dy-dy*(Math.log(ticky[j])-Math.log(axymin))/(Math.log(axymax)-Math.log(axymin)));
			righttext(ticky[j],x+dx-5,yy+6,TICKFONT,col); 
			drawline(x+dx,yy,x+dx+10,yy,BLACK,NORMLINE); 
		}
		break;
	
	case MATRIXXBUT:
		verticaltext(text,x+dx/2-3,y+4,TICKFONT,BLACK);
		break;
		
	case MATRIXYBUT:
		righttext(text,x+dx-5,y+dy/2+6,TICKFONT,BLACK);
		break;
		
	case XLABELBUT:
	 	col=BLACK; if(ov == 1){ col = GREY;}
		plotxlabel(text,x+dx/2,y+dy/2+4,LABELFONT,col);
		break;
		
	case YLABELBUT:
		col=BLACK; if(ov == 1){ col = GREY;}
		plotylabel(text,x+dx/2+5,y+dy/2,LABELFONT,col);
		break;
	
	case TABLEHEADBUT:
		plottext(text,x+7,y+20,TABLEHEADFONT,BLACK);
		break;
		
	case TABLEBUT:
		plottext(text,x+7,y+20,TABLEFONT,BLACK);
		break;
		
	case LINKBUT:
		col = BLUE; if(ov == 1){ col = LBLUE;}
		plottext(text,x+7,y+20,TABLEFONT,col);
		break;
		
	case VARIABLEBUT:
		var col = BLACK; if(ov == 1){ col = GREY;}
		if(text == selparam_name) col = RED;
		
		centertext(text,x+dx/2,y+0.8*dy,val,col);
		break;
		
	case EQBUT:
		var col = BLUE; if(ov == 1){ col = LBLUE;}
		if(text == seleq_name) col = RED;
		
		plottext(text,x+5,y+0.8*dy,val,col);
		break;
		
	case EQBUT2:
		plottext(text,x+5,y+0.8*dy,val,DDGREEN);
		break;
		
	case EQBUT3:
		plottext(text,x+5,y+0.8*dy,val,BLACK);
		break;	
		
	case PARAMLINEBUT:
		col = BLUE; if(ov == 1) col = LBLUE;
		drawline(x,y,x,y+dy,col,NORMLINE,2);
		
		plotverticaltext(text,x+5,y+5,VERTFONT,col);  
		break;
		
	default: alertp(ty+"Error code EC3"); break;
	}
}

function addcanbutton(text,x,y,dx,dy,ac,type,val,val2)   // Adds a canvas button
{
	canbuttext[ncanbut] = text;
	canbutx[ncanbut] = Math.floor(x);
	canbuty[ncanbut] = Math.floor(y);
	canbutdx[ncanbut] = Math.floor(dx);
	canbutdy[ncanbut] = Math.floor(dy);
	canbutac[ncanbut] = ac; 	
	canbuttype[ncanbut] = type;
	canbutover[ncanbut] = -1;
	canbutval[ncanbut] = val;
	canbutval2[ncanbut] = val2;
	ncanbut++;
}
