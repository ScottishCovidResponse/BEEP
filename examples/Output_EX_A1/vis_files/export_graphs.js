function putImage()
{
	var canvas = document.createElement('canvas');
	
	var DX, DY;
	for(i = 0; i < nbut; i++){
		if(buttype[i] == CANVASBUT){ DX = butdx[i]; DY = butdy[i];}
	}
	
	var graph_X = 2;  // The number of graphs on the x axis
	
	var num = tree[page].child[ps].child[pss].child.length;
	
  num = 6;
	
	canvas.width = DX*graph_X;
	canvas.height = Math.floor((num+graph_X-1)/graph_X)*DY;
	outcv = canvas.getContext('2d');
	outcv.beginPath();
	outcv.rect(0,0,canvas.width,canvas.height);
	outcv.fillStyle = "#ffffff";
	outcv.fill();
	
	var psss_st = psss;
	for(pagej = psss_st; pagej < psss_st+num; pagej++){
		var i = (pagej-psss_st)%graph_X;
		var j = Math.floor((pagej-psss_st)/graph_X);
		
		changepage(-1,-1,-1,pagej);
	
		outcv.drawImage(graphcan,0,0,DX,DY,i*DX,j*DY,DX,DY);
	}		
		 
	var image = canvas.toDataURL("image/png").replace("image/png", "image/octet-stream");  // here is the most important part because if you dont replace you will get a DOM 18 exception.

	window.location.href=image; 
}  
