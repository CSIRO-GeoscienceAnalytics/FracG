#include "model.h"
#include "graph.h"
#include "GeoRef.h"
#include "geometrie.h"

MODEL::MODEL ()
{
}

//need a specific comparator, not just selecting an appropriate element from the tuple
struct DistanceCompare
{
	bool operator()(std::pair<long double, point_type> const &t1, std::tuple<long double, point_type> const &t2)
	{
		return get<0>(t1) < get<0>(t2); //if they're both on the end point or (more likely) middle sectoin, add them increasing order of distance, so smaller distances first
	}
};

int AddPoint(Graph G, map<point_type, int, geometry::less<point_type>>& M, point_type new_P, ofstream& Mtxt, std::map<int, point_type>& vertex_map)
{
 string scale = "lc";
 map<point_type, int>::iterator it;
 int index = M.size() + 1;
 it = M.find(new_P);
 
 for (auto Ve : make_iterator_range(vertices(G))) 
	{
		if (geometry::equals(new_P, G[Ve].location))
		{
			if (degree(Ve, G) > 2)
				scale = "lc3";
			else
				scale = "lc2";
		}
	}
 
 if (it != M.end())
	 return(M.find(new_P)->second);
 else
 {
	 M[new_P] = index;
	 vertex_map[index] = new_P;
	 Mtxt << " Point(" << index << ") = {" << new_P.x() << ", " << new_P.y() << ",0 , "<< scale << "}; "<< endl;
	 
 }
 return(index);
}

void AddLines(std::map<pair<int, int>, int>& line_map, map<point_type, int, geometry::less<point_type>> M, ofstream& Mtxt)
{
	std::map<pair<int, int>, int>::iterator it;
	for (int i = 1; i < M.size(); i++)
	{
		if (i != line_map.size())
			{
				it = line_map.find(make_pair(i, i +1));
				if (it == line_map.end())
				{
					int index = line_map.size();
					line_map[make_pair(i, i+1)] = line_map.size() + 1;
					Mtxt << "Line(" << index + 1 << ") = {" << i <<", " << i+1 << "};" << endl;
				}
			}
	}
}

void CreateBoundingBox(double bufffer, box& AOI, vector<line_type>& BoundingBox, vector< std::pair<line_type, std::string> >& DomainLabels)
{
	point_type A, B, C , D;
	line_type AB, BC, CD, DA;
	double min_x = geometry::get<geometry::min_corner, 0>(AOI) - bufffer;
	double min_y = geometry::get<geometry::min_corner, 1>(AOI) - bufffer;
	double max_x = geometry::get<geometry::max_corner, 0>(AOI) + bufffer; 
	double max_y = geometry::get<geometry::max_corner, 1>(AOI) + bufffer; 
	
	geometry::set<0>(A, geometry::get<geometry::min_corner, 0>(AOI) - bufffer); 
	geometry::set<1>(A, geometry::get<geometry::min_corner, 1>(AOI) + bufffer);
	geometry::set<0>(B, geometry::get<geometry::max_corner, 0>(AOI) - bufffer); 
	geometry::set<1>(B, geometry::get<geometry::min_corner, 1>(AOI) + bufffer); 
	geometry::set<0>(C, geometry::get<geometry::max_corner, 0>(AOI) - bufffer); 
	geometry::set<1>(C, geometry::get<geometry::max_corner, 1>(AOI) + bufffer); 
	geometry::set<0>(D, geometry::get<geometry::min_corner, 0>(AOI) - bufffer); 
	geometry::set<1>(D, geometry::get<geometry::max_corner, 1>(AOI) + bufffer); 
	
	geometry::append(AB, A);
	geometry::append(AB, B);
	BoundingBox.push_back(AB);
	
	geometry::append(BC, B);
	geometry::append(BC, C);
	BoundingBox.push_back(BC);
	
	geometry::append(CD, C);
	geometry::append(CD, D);
	BoundingBox.push_back(CD);
	
	geometry::append(DA, D);
	geometry::append(AB, A);
	BoundingBox.push_back(DA);
	
	//define lables
	DomainLabels.push_back(make_pair(AB, "bottom"));
	DomainLabels.push_back(make_pair(CD, "top"));
	DomainLabels.push_back(make_pair(BC, "rigth"));
	DomainLabels.push_back(make_pair(DA, "left"));
	
	box aoi(A,C);
	AOI = aoi;
}

void MODEL::WriteGeo(vector<line_type> faults, string filename)
{
	double buffer = 500;
	double minSize = 10;
	float lc  = 1000, lc2 = 500, lc3 = 100;
	
	box AOI;
	Graph G;
	GRAPH g_class;
	GEOMETRIE geom;
	line_type newLine;
	map_vertex_type map;
	vector<point_type> I;
	polygon_type all_line;
	vector<line_type> BoundingBox;
	std::map<int, point_type> vertex_map;
	std::map<pair<int, int>, int> line_map;
	vector <std::string> SideSets(faults.size());
	vector <std::string> DomainSides(4); // bottom - rigth - top - left
	std::map<point_type, int, geometry::less<point_type>> point_map;
	
	vector< std::pair<line_type, std::string> > DomainLabels;

	string output = (string) filename.c_str() + "_4mesh.geo";
	ofstream Mtxt (output);	
	
//create a boundinx box around fault set--------------------------------
	for (auto L : faults) 
		geometry::append(all_line, L);
	
	AOI = boost::geometry::return_envelope<box>(all_line); 
	CreateBoundingBox(buffer, AOI, BoundingBox, DomainLabels);
	
//first we build  new graph and crop the edges to the bounding box------
	g_class.ReadVEC4MODEL(G, map, faults, AOI); //convert the faults into a graph
	g_class.SplitFaults(G, map, buffer); //split the faults in the graph into fault segments, according to the intersections of the faults
	g_class.RemoveSpurs(G, map, minSize); //remove any spurs from the graph network
	
	vector <pair<long double, point_type>> newEdge;
	for (int i; i < BoundingBox.size();i++)
	{
		for (auto Eg : make_iterator_range(edges(G)))
		{
			if (geometry::intersects(BoundingBox.at(i), G[Eg].trace))
			{
				geometry::intersection(BoundingBox.at(i), G[Eg].trace, I);
				newEdge.push_back(make_pair(0, BoundingBox.at(i).front()));
				for(auto x: I)
					newEdge.push_back(make_pair(geometry::distance(BoundingBox.at(i).front(), x), x));
				newEdge.push_back(make_pair(geometry::distance(BoundingBox.at(i).front(), BoundingBox.at(i).back()), BoundingBox.at(i).back()));
			}
		}
		if (newEdge.size() > 0)
		{
			sort(newEdge.begin(), newEdge.end(), DistanceCompare());
			
			pair<long double, point_type> p;
			BOOST_FOREACH(p, newEdge) 
				geometry::append(newLine, p.second);
			BoundingBox.at(i) = newLine;
		}
	geometry::clear(newLine);
	I.clear();
	newEdge.clear();
	}
	
//now we begin to write the geo file------------------------------------
	
//header with length scales---------------------------------------------
	int P = 0;
	if (Mtxt.is_open())  
	{ 
		Mtxt << "lc = "  << lc  << ";\n" 
			 << "lc2 = " << lc2 << ";\n" 
			 << "lc3 = " << lc3 << ";" << endl;

//create the bounding box of the domain---------------------------------
	for (auto L : BoundingBox)
	{
		BOOST_FOREACH(point_type p, L) 
			AddPoint(G, point_map, p, Mtxt, vertex_map);
	}
	AddLines(line_map, point_map, Mtxt);
	Mtxt << "Line(" << line_map.size() + 1 << ") = {" << point_map.size() <<", " << 1 << "};" << endl;
	line_map[make_pair(point_map.size(), 1)] = line_map.size() + 1;
	
	Mtxt <<  "Curve Loop(1) = {";
	for (std::map<pair<int, int>, int>::iterator it = line_map.begin(); it != line_map.end(); ++it)
	{
		auto nx = std::next(it, 1);
		if (nx == line_map.end())
			Mtxt << it->second << "};" << endl;
		else
			Mtxt << it->second << ", " << endl;
	}
	Mtxt << "Plane Surface(1) = {1};" << endl;
		 
		 
		 
		 
		 
	int domain_l = line_map.size();

	
	
	
	
	int i = 0;
	for (std::map<pair<int, int>, int>::iterator it = line_map.begin(); it != line_map.end(); ++it)
	{

	
		line_type Edge;
		
		
		if(geometry::within(Edge, geom.DefineLineBuffer(DomainLabels[0].first, 10)))
			DomainSides[0] += std::to_string(i) + ", ";
			
		if(geometry::within(Edge, geom.DefineLineBuffer(DomainLabels[1].first, 10)))
			DomainSides[1] += std::to_string(i) + ", ";
		
		if(geometry::within(Edge, geom.DefineLineBuffer(DomainLabels[2].first, 10)))
			DomainSides[2] += std::to_string(i) + ", ";
		
		if(geometry::within(Edge, geom.DefineLineBuffer(DomainLabels[3].first, 10)))
			DomainSides[3] += std::to_string(i) + ", ";
	i++;
	}
	for (auto &s : DomainSides)
		s = s.substr(0, s.length() - 2);
	
	
	
	

//add the lineamnts to the domain---------------------------------------
	for (auto Eg : make_iterator_range(edges(G))) 
	{
		vertex_type U = source(Eg, G);
		vertex_type	u = target(Eg, G);
		int p1 = AddPoint(G, point_map, G[U].location, Mtxt, vertex_map);
		int p2 = AddPoint(G, point_map, G[u].location, Mtxt, vertex_map);
		line_map[make_pair(p1, p2)] = line_map.size()+1;
		Mtxt << "Line(" << line_map.size() << ") = {" << AddPoint(G, point_map, G[U].location, Mtxt, vertex_map) <<", " << AddPoint(G, point_map, G[u].location, Mtxt, vertex_map) << "};" << endl;
		SideSets[G[Eg].FaultNb] += std::to_string(line_map.size()) + ", ";
	}
	
	for (auto &s : SideSets)
		s = s.substr(0, s.length() - 2);

// now that we defined the domain, we add the geometries----------------
	Mtxt <<"SetFactory(\"OpenCASCADE\");\n"
		 <<"Curve{";
	for (std::map<std::pair<int,int>,int>::iterator it = line_map.begin(); it != line_map.end(); ++it ) 
	{
		auto nx = std::next(it, 1);
		if (std::distance(line_map.begin(), it) < domain_l)
			continue;
		if (nx == line_map.end())
			Mtxt << it->second << "}";
		else
			Mtxt << it->second << ", ";
	}
	Mtxt << " In Surface{1};\n"
		 << "Physical Surface(\"block\") = {1}; " << endl;
		 
	for (int i = 0; i < SideSets.size(); i++)
		Mtxt << "Physical Curve(\"ss" << i << "\") = {" << SideSets[i] << "};" << endl;

	for (int i = 0; i < DomainSides.size(); i++)
		Mtxt << "Physical Curve(\"" << DomainLabels[i].second << "\") = {" << DomainSides[i] << "};" << endl;

	Mtxt.close();
	}
	else cout << "Ne results written" << endl;
}
