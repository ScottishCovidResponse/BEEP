/// This contains code for loading boundary file information

#include <iostream>
#include <vector>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <math.h> 
#include <iomanip>
#include "json.hpp"

using namespace std;
using json = nlohmann::json;

#include "data.hh"	


/// Creates boundary data from coordinates
void Data::create_boundaries(string x, string y, vector < vector < vector <Coord> > > &bound) const 
{
	if(x == ""){
		auto LY= 1u;
		while(LY*(int)(LY*map_ratio) < narea) LY++;
		auto LX = (int)(LY*map_ratio);
		
		for(auto c = 0u; c < narea; c++){
			vector <Coord> poly;
			Coord ll; 
			ll.x = c%LX;
			ll.y = c/LX;
			poly.push_back(ll);
			bound[c].push_back(poly);
		}
	}
	else{
		string file = inputs.find_string("areas","UNSET");
		Table tab = load_table(file);
		auto xcol = find_column(tab,x); if(xcol == UNSET) emsg("Cannot find the column '"+x+"' in '"+file+"'");
		auto ycol = find_column(tab,y); if(ycol == UNSET) emsg("Cannot find the column '"+y+"' in '"+file+"'");
	
		for(auto c = 0u; c < narea; c++){
			vector <Coord> poly;
			Coord ll; 
			ll.x = get_double(tab.ele[c][xcol],"In file '"+file+"' ");
			ll.y = get_double(tab.ele[c][ycol],"In file '"+file+"' ");
			poly.push_back(ll);
			bound[c].push_back(poly);
		}
	}
}


/// Returns the file type
FileType Data::filetype(const string file) const 
{
	auto spl = split(file,'.');
	if(spl.size() != 2) emsg("File '"+file+"' should either be '.kml' or '.geojson'");
	
	auto type = toLower(spl[1]);
	if(type == "kml") return KML;
	if(type == "geojson") return GEOJSON;
	emsg("File '"+file+"' should either be '.kml' or '.geojson'");
}

/// Loads up the boundaries
void Data::load_boundaries(const string file, vector < vector < vector <Coord> > > &bound) const
{
	auto type = filetype(file);
	
	cout << "Loading boundary file '" << file << "'" << endl;

	switch(type){
		case KML: load_KML(file,bound); break;
		case GEOJSON: load_geojson(file,bound); break;
		default: emsg("File '"+file+"' should either be '.kml' or '.geojson'"); break;
	}
}


/// Loads up a KML file and generates boundary data
void Data::load_KML(const string file, vector < vector < vector <Coord> > > &bound) const
{ 
	XMLDocument doc;
		
	auto filefull = data_directory+"/"+file;
  ifstream inFile(filefull);

  stringstream strStream;
  strStream << inFile.rdbuf();
  string str = strStream.str();
  if(str.size() == 0) emsg("Could not open '"+file+"'.");
	
	char *xml=new char[long(str.size())+1];
 
  xml[str.size()]=0;
  memcpy(xml,str.c_str(),str.size());
  doc.Parse(xml);

  XMLElement* root = doc.FirstChildElement();
	
	placemark_KML(root,bound);
}


/// Loads up a placemark
void Data::placemark_KML(XMLNode* root, vector < vector < vector <Coord> > > &bound) const 
{
	for(auto child = root->FirstChild(); child; child = child->NextSibling()){           
		string s = child->Value();
		if(s == "Placemark"){
			auto c = get_area_name_KML(child);
			if(c != UNSET) get_polygon_KML(child,c,bound);
		}
		else placemark_KML(child,bound);
	}
}


/// iteratively searchs the tree for the name of the area
unsigned int Data::get_area_name_KML(XMLNode* child) const 
{	
	for(auto child2 = child->FirstChild(); child2; child2 = child2->NextSibling()){
		string s2 = child2->Value();
		if(s2 != "Polygon"){	
			XMLText* textNode = child2->ToText();
			if(textNode){
				string name = textNode->Value();
				strip(name);		
				auto c = 0u; while(c < narea && name != area[c].code) c++;
				if(c < narea) return c;
			}		 
			auto c = get_area_name_KML(child2);
			if(c != UNSET) return c;
		}
	}
	
	return UNSET;
}


/// iteratively searchs the tree for polygons
void Data::get_polygon_KML(XMLNode* child, const unsigned int c, vector < vector < vector <Coord> > > &bound) const 
{	
	for(auto child2 = child->FirstChild(); child2; child2 = child2->NextSibling()){
		string s2 = child2->Value();
		if(s2 == "Polygon") get_coordinate_KML(child2,c,bound);
		else get_polygon_KML(child2,c,bound);
	}
}


/// Within polygon iteratively search for coordinates
void Data::get_coordinate_KML(XMLNode* child, const unsigned int c, vector < vector < vector <Coord> > > &bound) const 
{	
	for(auto child2 = child->FirstChild(); child2; child2 = child2->NextSibling()){
		string s2 = child2->Value();
		if(s2 == "coordinates"){	
			XMLText* coord = child2	->FirstChild()->ToText();
			if(coord){
				string name = coord->Value();
				auto spl = split(name,' ');
				
				vector <Coord> poly;
				for(auto i = 0u; i < spl.size(); i++){
					auto spl2 = split(spl[i],',');
					Coord ll; 
					ll.x = get_double(spl2[0],"Reading boundary data, "); 
					ll.y = get_double(spl2[1],"Reading boundary data, ");
					poly.push_back(ll);
				}
				bound[c].push_back(poly);
				return;
			}
		}
		else get_coordinate_KML(child2,c,bound);
	}
}


/// Loads up a GEOJSON file and generates boundary data
void Data::load_geojson(const string file, vector < vector < vector <Coord> > > &bound) const
{
	auto filefull = data_directory+"/"+file;

	ifstream boundfile(filefull);
	if(!boundfile) emsg("Cannot open the file '"+file+"'.");
	string s;
	do{
		string st;
		getline( boundfile,st); s += st; 
		if( boundfile.eof()) break;
	}while(true);
	
	json jso = json::parse(s);
	
	auto h = jso["features"];
	
	for(auto &it : h.items()){
		json val = it.value();
		
		auto h2 = val["type"];
		if(h2 == "Feature"){
			auto c = UNSET; 
			auto h3 = val["properties"];
			for(auto &it3 : h3.items()){     // Checks if one of the properties contains the name of the area 
				if(it3.value().is_string()){
					string name = it3.value();
					strip(name);		
					c = 0; while(c < narea && name != area[c].code) c++;
					if(c < narea) break;
				}
			}
		
			if(c != UNSET && c < narea){
				auto h3 = val["geometry"];	
				auto h4 = h3["type"];		
				if(h4 == "Polygon"){
					auto h5 = h3["coordinates"];
					for(auto &it2 : h5.items()){
						json val2 = it2.value();
						vector <Coord> poly;
						for(auto &it3 : val2.items()){
							json val3 = it3.value();
							Coord ll; ll.x = val3[0]; ll.y = val3[1];
							poly.push_back(ll);
						}
						bound[c].push_back(poly);
					}
				}
				
				if(h4 == "MultiPolygon"){
					auto h5 = h3["coordinates"];
					for(auto &it : h5.items()){
						json val = it.value();
						for(auto &it2 : val.items()){
							json val2 = it2.value();
							vector <Coord> poly;
							for(auto &it3 : val2.items()){
								json val3 = it3.value();
								Coord ll; ll.x = val3[0]; ll.y = val3[1];
								poly.push_back(ll);
							}
							bound[c].push_back(poly);
						}
					}
				}
			}
		}
	}
}


/// Rescales and shifts the boundary to fit into a box map_ratio by 1
void Data::rescale_boundary(vector < vector < vector <Coord> > > &bound) const
{
	double xmin = LARGE, xmax = -LARGE;                   
	double ymin = LARGE, ymax = -LARGE;
	for(auto c = 0u; c < narea; c++){
		for(auto i = 0u; i < bound[c].size(); i++){
			for(auto j = 0u; j < bound[c][i].size(); j++){
				auto valx = bound[c][i][j].x; if(valx > xmax) xmax = valx; if(valx < xmin) xmin = valx; 
				auto valy = bound[c][i][j].y; if(valy > ymax) ymax = valy; if(valy < ymin) ymin = valy; 
			}
		}
	}
	
	auto dx = xmax-xmin, dy = ymax-ymin;                                 // Shifts and scales boundaries
	
	double fac;
	if(dx/dy > map_ratio) fac = map_ratio/dx; else fac = 1.0/dy;
	fac *= 0.9;
	
	auto xav = 0.5*(xmin+xmax), yav = 0.5*(ymin+ymax);
	
	auto npo = 0u;
	for(auto c = 0u; c < narea; c++){                      
		for(auto i = 0u; i < bound[c].size(); i++){
			for(auto j = 0u; j < bound[c][i].size(); j++){
				npo++;
				bound[c][i][j].x = map_ratio/2 + fac*(bound[c][i][j].x-xav); 
				bound[c][i][j].y = 0.5 + fac*(bound[c][i][j].y-yav); 
			}
		}	
	}
}


/// For cases when there is no boundary data then uses circles
void Data::make_circle_boundary(const string xcol, vector < vector < vector <Coord> > > &bound) const
{
	const unsigned int N = 50;

	vector <Coord> center(narea);
	for(auto c = 0u; c < narea; c++) center[c] = bound[c][0][0]; 
	
	auto distmax = 0.0;
	for(auto c = 0u; c < narea; c++){
		auto distmin = LARGE;
		for(auto cc = 0u; cc < narea; cc++){
			if(cc != c){
				auto dx = center[cc].x - center[c].x;
				auto dy = center[cc].y - center[c].y;
				auto dist = dx*dx+dy*dy;
				if(dist < distmin) distmin = dist;
			}
		}
		if(distmin > distmax) distmax = distmin;
	}
	auto r = 0.5*sqrt(distmax);
	
	if(xcol != "") r *= 1.45;
	else r *= 0.9;

	for(auto c = 0u; c < narea; c++){
		auto x = center[c].x; 
		auto y = center[c].y;
		
		vector <Coord> poly;
		for(auto i = 0u; i < N; i++){
			Coord ll; 
			ll.x = x+r*cos((2*M_PI*i)/N);
			ll.y = y+r*sin((2*M_PI*i)/N);
			poly.push_back(ll);
		}
		bound[c][0] = poly;
	}
	
	for(auto c = 0u; c < narea; c++){                         // Clips circle in a Voronoi way
		for(auto cc = 0u; cc < narea; cc++){
			if(cc != c){
				auto x = (center[c].x + center[cc].x)/2;
				auto y = (center[c].y + center[cc].y)/2;
				auto dx = center[cc].y - center[c].y;
				auto dy = -(center[cc].x - center[c].x);
			
				vector <Coord> inter;
				vector <unsigned int> direction;
				vector <unsigned int> pos;
				const auto &p = bound[c][0];
				for(auto i = 0u; i < p.size(); i++){
					auto ii = (i+1)%p.size();
					auto abx = p[ii].x-p[i].x, aby = p[ii].y-p[i].y;
					auto acx = x-p[i].x, acy = y-p[i].y;
					
					auto alpha = (acx*dy - acy*dx)/(abx*dy - aby*dx);
					if(alpha >= 0 && alpha < 1){
						pos.push_back(i);
						Coord co; 
						co.x = p[i].x + alpha*abx;
						co.y = p[i].y + alpha*aby;
						inter.push_back(co);
						if(abx*dy - aby*dx > 0) direction.push_back(1);
						else direction.push_back(0);
					}
				}

				switch(pos.size()){
					case 0: case 1: break;
					case 2:
						{
							vector <Coord> poly;
							if(direction[0] == 1 && direction[1] == 0){
								poly.push_back(inter[1]); poly.push_back(inter[0]);
								for(auto i = pos[0]+1; i <= pos[1]; i++) poly.push_back(p[i]);
							}
							else{
								if(direction[0] == 0 && direction[1] == 1){
									poly.push_back(inter[0]); poly.push_back(inter[1]);
									for(auto i = pos[1]+1; i <= pos[0]+p.size(); i++) poly.push_back(p[i%p.size()]);
								}
								else{
									emsgEC("Data",64);
								}
							}
							bound[c][0] = poly;
						}
						break;
					default: emsgEC("Data",65); break;
				}
			}
		}
	}
	
	rescale_boundary(bound);	
}
