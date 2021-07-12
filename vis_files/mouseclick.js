// Fires when the mouse is clicked 
function mouseclick(x,y) 
{
	mx = x; my = y;
	if(over >= 0){ finalaction(over); over = -1; canover = -1;}
	
	buttoninit(); mousemove(mx,my);
}


/// Fires for a double click
function mousedblclick()
{
	if(buttype[over] == CANVASBUT){
		if(plottype == "map" && canover == -1) mapzoom(1.3,mx-menux,my);
	}
}
	
	
/// What happens when a button is clicked
function finalaction(i)
{
	var x, y, dx, dy, text, val, val2, ac; 

	x = butx[i]; y = buty[i]; dx = butdx[i]; dy = butdy[i], text = buttext[i]; 
	val = butval[i]; val2 = butval2[i]; ac = butac[i];
	
	if(mx < x || mx > x+dx || my < y || my > y+dy) return;
	
	switch(ac){
	case TABBUT:		
		changepage(val,-1,-1,-1);
		break;
    
	case PAGESUBBUT:
		changepage(-1,val,-1,-1);
		break;
	
	case PAGESUBSUBBUT:
		changepage(-1,-1,val,-1);
		break;
	
	case PAGESUBSUBSUBBUT:
		changepage(-1,-1,-1,val);
		break;
	
	case SOURCEBUT:
		slidey = 0;
		sourceon = 1-sourceon; 
		break;
		
	case CANVASBUT:
		if(canover != -1){
			if(canbutac[canover] != -1) canfinalaction(canover);
		}
		else{
			selparam_table = null; selparam_name = "";
			seleq = null; seleq_name = ""; 
		}
		break;
		
	case CHECKBUT:
		check = 1-check;
		break;
		
	case MENULINKBUT:
		follow_param_link(val);
		break;
		
	case RADIOBUT:
		switch(val2){
			case RADIORATE: rateradio = val; break;
			case RADIOHIGH: lowerhighlight = val; break;
		}
		break;
		
	case MENUSLIDEBUT: case SLIDEBUT: case TABXSLIDEBUT: case TABYSLIDEBUT: break;
	
	default: alertp("Error code EC4"); break;
	}
}

function alertp(st)
{
	if(alertdone == 0){ alert(st); alertdone = 1;}
}
