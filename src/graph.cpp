#include "graph.h"
#include "geometrie.h"
#include "GeoRef.h"
#include "stats.h"

#include <random>

GRAPH::GRAPH ()

{
}


// create and open filestream in folder "statistics"
ofstream CreateGraphFileStream(string folder, string name)
{
	string folder_name = folder + "/graph/";
	const char* stats_dir = folder_name.c_str();
		if(!opendir(stats_dir))
		
	mkdir(stats_dir, 0777);
	string tsv_file = stats_dir + name + (string)"_g_statistics.tsv";
	ofstream txtF; 
	txtF.open (tsv_file, ios::out | ios::app); 
	return(txtF);
}


//take a location, and return the corresponding graph vertex, and the list of possible vertices
//remember that locations can be slightly different due to floating point calculations, so this function checks for vertices within a certain distance of the specified point
//this is a helper function for getting vertices from the graph
void GetVertexData(map_vertex_type& map, point_type const &key, Graph &graph, std::vector<vertex_type> **vertex_list, vertex_type **found)
{
	*vertex_list = nullptr;
	*found = nullptr;
	const double threshold = 1; //verticies are considered to be in the same if the key and existing vertex are within this distance of each other
	//check.set<0>();
	//check.set<1>();
	//check the "native" integer key first, then check all other integer keys around it
	const long long xPos = (long long)key.x();
	const long long yPos = (long long)key.y();
	const long int check_threshold = lrint(ceil(threshold));
	std::vector<point_int> key_list;
	point_int base_key(xPos, yPos);
	key_list.push_back(base_key);
	for (long long xoff = -check_threshold; xoff <= check_threshold; xoff++){
		for (long long yoff = -check_threshold; yoff <= check_threshold; yoff++){
			if (xoff == 0 && yoff == 0) continue; //don't add the base again
			point_int check(xPos+xoff, yPos+yoff);
			key_list.push_back(check);
		}
	}
	
	for (auto cit = key_list.begin(); cit != key_list.end(); cit++) //iterating through list of keys to check
	{
		typename map_vertex_type::iterator it = map.find(*cit); //search the map for the current key
		if (it == map.end()) continue; // no vertex found here, check other locations
		std::vector<vertex_type> *possible = &(it->second);
		
		for (auto pit = possible->begin(); pit != possible->end(); pit++) //check the list of possible existing verticies
		{
			if (geometry::distance(graph[*pit].location, key) <= threshold){
				//if we find a match, return it
				if (vertex_list != NULL) *vertex_list = possible;
				*found = &(*pit);
				return;
			}
		}
	}
	//haven't found a vertex
	if (vertex_list == NULL) return;
	typename map_vertex_type::iterator it = map.find(base_key);
	if (it != map.end()) *vertex_list = &(it->second); //we don't have a vertex, but we do have a list where it should go
}

//add new vertex to graph-----------------------------------------------
//or return an existing vertex, if there is one
vertex_type GRAPH::AddNewVertex(map_vertex_type& map, point_type const& key, Graph& graph)
 {
	std::vector<vertex_type> *possible = nullptr;
	vertex_type *found = nullptr;
	GetVertexData(map, key, graph, &possible, &found);
	if (found != NULL)
	{
		return *found; //we have a vertex, return it
	}
	if (possible == NULL)
	{
		//no list of possible vertices, so make one
		point_int insert_key((long long)key.x(), (long long)key.y());
		std::vector<vertex_type> possible_vector;
		map[insert_key] = possible_vector;
		possible = &(map.find(insert_key)->second);
	}
	//we haven't found an existing vertex, so make a new one
	vertex_type new_vertex = add_vertex(FVertex<point_type>(key), graph);
	possible->push_back(new_vertex);
	return new_vertex;
	//~ return it->second; 
 }
 
//get a vertex from the graph, or null if there is no vertex at the specified location
vertex_type GRAPH::GetVertex(map_vertex_type& map, point_type const& key, Graph& graph)
{
	std::vector<vertex_type> *possible = nullptr;
	vertex_type *found = nullptr;
	GetVertexData(map, key, graph, &possible, &found);
	return *found;
}

//add new edge to the graph---------------------------------------------
//Graph G, Source S, Target T, FaultSeg is the segment of the fault that makes up this edge
 void GRAPH::AddNewEdge(Graph& G, vertex_type S, vertex_type T, line_type FaultSeg)
{
	GRAPH::AddNewEdge(G, S, T, FaultSeg, geometry::length(FaultSeg));
}

//add a new edge to the graph, and record the full length of the fault that comprises this fault segment
void GRAPH::AddNewEdge(Graph& G, vertex_type S, vertex_type T, line_type FaultSeg, double FaultLength)
{
	//cout << "Assing edge between " << S << " and " << T << " with number of vertices " << num_vertices(G) << endl;
	edge_type edge;
	bool success = false;
		if (!boost::edge(S, T, G).second && !geometry::equals(G[S].location, G[T].location))
		{
			std::tie(edge, success) = add_edge(S, T,
				{geometry::length(FaultSeg)/1000, FaultSeg}, G);
			G[edge].fault_length = FaultLength;
	}
}

//read vector data into graph------------------------------------------
void GRAPH::ReadVEC(Graph& graph, map_vertex_type& map, std::vector<line_type> &faults)
{
	GEO g;
	GEOMETRIE geom;
	vertex_type VA, VB;

	//create edges and nodes for faults in the graph-------------------------
	int nb = 0;
	BOOST_FOREACH(line_type f, faults)
	{
		VA = AddNewVertex(map, f.front(), graph);
		VB = AddNewVertex(map, f.back(), graph);

		AddNewEdge(graph, VA, VB, f);
		
		if (boost::edge(VA,VB, graph).second)
		{
			edge_type e = boost::edge(VA,VB, graph).first;
				graph[e].FaultNb = nb;
		}
		nb++;
	}
	cout << "	Converted " << faults.size() << " lineaments into " << num_edges(graph) << " edges." << endl;
}

void GRAPH::ReadVEC4raster(Graph& graph, map_vertex_type& map, std::vector<line_type> &faults)
{
	GEO g;
	GEOMETRIE geom;
	geometry::model::multi_linestring<line_type> intersection;
	vector<std::tuple< std::pair<point_type, point_type>, line_type, unsigned int, double >> G;
	std::pair <point_type, point_type> NODES;
	line_type TRACE;
	long double Xmax, Xmin, Ymax, Ymin;
	vertex_type VA, VB;
	unsigned int type;

	const unsigned int FRONT_ENODE = 1 << 0;
	const unsigned int BACK_ENODE  = 1 << 1;
	const double endpoint_threshold=0.5;
	
	polygon_type pl = g.BoundingBox(GeoTransform, 100);

	//check for edge nodes (crop faults by raster area)--------------------

	BOOST_FOREACH(line_type const& Fault, faults)
	{
		if (geometry::disjoint(Fault, pl)) continue;
		
		geometry::intersection(Fault, pl, intersection);
		
		if (intersection.size() <= 0) continue;
		type = 0;
		const double fault_length = geometry::length(Fault);
		for (auto it = intersection.begin(); it != intersection.end(); it++)
		{
			line_type &segment = *it;
			//use the distance from the original fault's end points to determine if this segment's end points are from the original fault segment, or are cut off by the boundaries of the area of interest
			if (geometry::distance(segment.front(), Fault.front()) > endpoint_threshold && geometry::distance(segment.front(), Fault.back()) > endpoint_threshold) type |= FRONT_ENODE;
			if (geometry::distance(segment.back() , Fault.front()) > endpoint_threshold && geometry::distance(segment.back() , Fault.back()) > endpoint_threshold) type |=  BACK_ENODE;
			G.push_back(make_tuple(make_pair(segment.front(),segment.back()), segment, type, fault_length));
		}
		intersection.clear();
	}

	//create edges and nodes for faults in the graph-------------------------
	for (typename vector <std::tuple< std::pair<point_type, point_type>, line_type, unsigned int, double>>::const_iterator it = G.begin(); it != G.end(); ++it)
	{
		NODES = get<0> (*it);
		TRACE = get<1> (*it);
		type  = get<2> (*it);
		const double fault_length = get<3>(*it);
		
		VA = AddNewVertex(map, NODES.first, graph);
		VB = AddNewVertex(map, NODES.second, graph);

		//if( (degree(VA, graph) == 0) && (degree(VB, graph) == 0) )
		AddNewEdge(graph, VA, VB, TRACE, fault_length);
				
		if (type & FRONT_ENODE) graph[VA].Enode = true;
		if (type &  BACK_ENODE) graph[VB].Enode = true;
	}

	cout << " Converted " << faults.size() << " faults into " << num_edges(graph) << " edges" << endl;
	
	faults.clear();
	for (auto Eg : make_iterator_range(edges(graph)))
		faults.push_back(graph[Eg].trace);
}

Graph GRAPH::ReadVEC4MODEL(std::vector<line_type> faults, box bx)
{
	GEO g;
	Graph graph;
	GEOMETRIE geom;
	map_vertex_type map;
	
	geometry::model::multi_linestring<line_type> intersection;
	vector<std::tuple< std::pair<point_type, point_type>, line_type, unsigned int, int >> G;
	std::pair <point_type, point_type> NODES;
	line_type TRACE;
	long double Xmax, Xmin, Ymax, Ymin;
	vertex_type VA, VB;
	unsigned int type;
		
	const unsigned int FRONT_ENODE = 1 << 0;
	const unsigned int BACK_ENODE  = 1 << 1;
	const double endpoint_threshold=0.5;
	
	//check for edge nodes (crop faults by  area)--------------------
	int nb = 0;
	BOOST_FOREACH(line_type const& Fault, faults)
	{
		if (geometry::disjoint(Fault, bx)) continue;
		
		geometry::intersection(Fault, bx, intersection);
		
		if (intersection.size() <= 0) continue;
		type = 0;
		for (auto it = intersection.begin(); it != intersection.end(); it++)
		{
			line_type &segment = *it;
			//use the distance from the original fault's end points to determine if this segment's end points are from the original fault segment, or are cut off by the boundaries of the area of interest
			if (geometry::distance(segment.front(), Fault.front()) > endpoint_threshold && geometry::distance(segment.front(), Fault.back()) > endpoint_threshold) type |= FRONT_ENODE;
			if (geometry::distance(segment.back() , Fault.front()) > endpoint_threshold && geometry::distance(segment.back() , Fault.back()) > endpoint_threshold) type |=  BACK_ENODE;
			G.push_back(make_tuple(make_pair(segment.front(),segment.back()), segment, type, nb));
		}
		intersection.clear();
		nb++;
	}

	//create edges and nodes for faults in the graph-------------------------
	for (typename vector <std::tuple< std::pair<point_type, point_type>, line_type, unsigned int, int >>::const_iterator it = G.begin(); it != G.end(); ++it)
	{
		NODES = get<0> (*it);
		TRACE = get<1> (*it);
		type  = get<2> (*it);
		nb = get<3> (*it);

		VA = AddNewVertex(map, NODES.first, graph);
		VB = AddNewVertex(map, NODES.second, graph);

		AddNewEdge(graph, VA, VB, TRACE);
				
		if (type & FRONT_ENODE) graph[VA].Enode = true;
		if (type &  BACK_ENODE) graph[VB].Enode = true;
	}

	if (num_edges(graph) > 0)
	{
	
		cout << " Converted " << faults.size() << " faults into " << num_edges(graph) << " edges" << endl;
		
		faults.clear();
		for (auto Eg : make_iterator_range(edges(graph)))
			faults.push_back(graph[Eg].trace);
			
		return (graph);
	}
	else
	{
		cout << "No lineaments found in samling window" << endl;
		return(graph);
	}
}

//find and remove spurs from the network
//a "spur" is a fault segment that extends only slightly past a fault intersection
//these are assumed to be an issue of the digitisation/data quality, and that real physical faults either end at the fault intersection, or continue for a large enough distance past the fault intersection
//this function assumes that the fault graph has already been split into the fault segments between intersections
void GRAPH::RemoveSpurs(Graph& G, map_vertex_type& map, double minDist)
{
	 point_int rmP;
	 vertex_type U ,u;
//  cout << "Checking network" << endl;
	restart:
// 	cout << "Map size = " << map.size() << endl;
	for (auto eg : boost::make_iterator_range(edges(G))) 
	{
		if (geometry::length(G[eg].trace) < minDist)
		{
			U = source(eg, G);
			u = target(eg, G);
			
			if (degree(U, G) == 1 || 
				degree(u, G) == 1)
			{
				remove_edge(U, u, G);
				if (degree(U,G) == 0)
				{
// 					cout << "WARNING: Removing edge to vertex " << endl;
					rmP.set<0>((long long)G[U].location.x());
					rmP.set<1>((long long)G[U].location.y());
					remove_vertex(U, G);
					map.erase(rmP);
				}
				if (degree(u,G) == 0)
				{
// 					cout << "WARNING: Removing edge to vertex " << endl;
					rmP.set<0>((long long)G[u].location.x());
					rmP.set<1>((long long)G[u].location.y());
					remove_vertex(u, G);
					map.erase(rmP);
				}
				goto restart;
			}
		}
	}
}

//return whether the given point should be attached to the front, middle, or end of this fault
//this is used to sort the other faults that intersect with this one
AttachPoint LocateAttachPoint(line_type &fault, point_type &point, double distance_along, double threshold){
	double distance_from_fault = geometry::distance(fault, point);
	const double fault_distance_threshold = 0.5;
	if (distance_from_fault <= fault_distance_threshold) return AttachPoint::middle; //if this point is close enough to the fault, then attach it to the middle
	//from here below, the point is *not* positioned directly on the fault, so figure out if it needs to be attached to either the front, middle, or back
	if (distance_along <= threshold) return AttachPoint::front;
	double distance_to_end = geometry::length(fault);
	if (abs(distance_along - distance_to_end) <= threshold) return AttachPoint::back;
	return AttachPoint::middle;
}

//take a graph and split the faults up into fault segments, according to where the faults intersect
void GRAPH::SplitFaults(Graph& graph, map_vertex_type& map, double minDist)
{
	Graph g;
	map_vertex_type m; //the map translates physical coordinates and nodes in the graph.
	GEOMETRIE geom;
	long double Distance;
	line_type fault, fault2;
	edge_iter Eg, Eg_end, Eg2, Eg2_end;
	vertex_type NewV;
	vector<point_type> Intersec;
	vector<point_type> rmVert;

	vector <std::pair<vertex_type, vertex_type >> removeEdge;
	std::pair<vertex_type, vertex_type > removeS, removeT;
	vector <std::tuple<long double, point_type, AttachPoint>> cross; //point type is the start/end point of the fault, the double is the distance along the (other) fault that the intersection occurrs at, AttachPoint is whether the other fault attaches to the front or back of this fault (for cases where the intersection point is outside of the fault but withing the threshold distance) or middle, if the faults intersect
	//the AttachPoint is used in the sorting of the intersections

	typedef vector <std::tuple < edge_iter, vector< std::tuple<long double, point_type, AttachPoint>> >> UpGraph;
	UpGraph GraphBuild;
	
	for (tie(Eg, Eg_end) = edges(graph); Eg != Eg_end; ++Eg)
	{
		fault  = graph[*Eg].trace;
		int N  = graph[*Eg].FaultNb;
		const double fault_length = geometry::length(fault);
		cross.push_back(make_tuple(0, fault.front(), AttachPoint::middle));
		for (tie(Eg2, Eg2_end) = edges(graph); Eg2 != Eg2_end; ++Eg2)
		{
			fault2 = graph[*Eg2].trace;
			Intersec.clear();
			if(Eg == Eg2) continue;
			
			//check for X-node------------------------------------------------------
			if(geometry::intersects(fault, fault2)) 
			{
				geometry::intersection(fault, fault2, Intersec);
				Distance =  geometry::length(geom.GetSegment(fault, fault.front(), Intersec.at(0)));
				cross.push_back(make_tuple(Distance, Intersec.at(0), AttachPoint::middle));
			}
			//check for Y-node------------------------------------------------------
			else
			{
				if (geometry::distance( fault, fault2.front()) < minDist )
				{
					Distance =  geometry::length(geom.GetSegment(fault, fault.front(), fault2.front()));
					AttachPoint attach_to = LocateAttachPoint(fault, fault2.front(), Distance, minDist);
					cross.push_back(make_tuple(Distance, fault2.front(), attach_to));
				}
				if (geometry::distance(fault, fault2.back()) < minDist )
				{
					Distance =  geometry::length(geom.GetSegment(fault, fault.front(), fault2.back()));
					AttachPoint attach_to = LocateAttachPoint(fault, fault2.back(), Distance, minDist);
					cross.push_back(make_tuple(Distance, fault2.back(), attach_to));
				}
			}
		}
		cross.push_back(make_tuple(geometry::length(fault), fault.back(), AttachPoint::middle));
		
		geom.SortDist(cross); //sort the vertices so that we get them in order of how far along the fault they appear (while taking into account that some intersection points appear off the fault line itself)
		vertex_type prev_vertex = AddNewVertex(m, std::get<1>(cross.front()), g);
		for (vector<std::tuple<long double, point_type, AttachPoint>>::const_iterator I = cross.begin(); I != cross.end(); I++)
		{
			if (I == cross.begin()) continue;
			point_type intersect = get<1>(*I);
			//bool is_start = geometry::distance(fault.front(), intersect) <= minDist;
			//bool is_end   = geometry::distance(fault.back(),  intersect) <= minDist;
			//if (is_start || is_end) continue;
			NewV = AddNewVertex(m, intersect, g);
			if (NewV == prev_vertex) continue;
			
			AddNewEdge(g, prev_vertex, NewV, geom.GetSegment(fault, g[prev_vertex].location, g[NewV].location), fault_length);//also remember the length of the original fault
			
			if(boost::edge(prev_vertex, NewV, g).second) //if an edge exists between prev_vertex and NewV
			{
				auto e = boost::edge(prev_vertex, NewV, g).first; //edge descriptor of the above
				graph[e].FaultNb = N; //set the fault number
			}
			prev_vertex = NewV;
		}

		cross.clear();
	}
	
// 	cout << endl;
	graph = g;
	map = m;
}

void ClassifyEdges(Graph &G, int &Inodes, int &Ynodes, int &Xnodes, int &Enodes, int &II, int &IC, int &CC, int &NbB, double &totalLength, int &numK, double &Area)
{
	box envBox;
	map_type COUNT;
	vertex_type U, u;
	polygon_type all_lines;
	typename map_type::const_iterator it;
	std::map<vertex_type, int> components;
	std::map<vertex_type, int>::iterator comp_it;
	numK = boost::connected_components(G, boost::make_assoc_property_map(components)); 
	
		//classify edges
	for (auto Eg : make_iterator_range(edges(G)))
	{
		U = source(Eg, G);
		u = target(Eg, G);
		
		comp_it = components.find(U);
		if (comp_it != components.end())
		{
			if(components.find(U) == components.find(u))
			cout << "hmm" << endl;
			G[U].component  = comp_it->second;
			G[u].component  = comp_it->second;
			G[Eg].component = comp_it->second;
		}

		if (G[U].Enode == false || G[u].Enode == false)
		{
			G[Eg].component = G[U].component;
			it = COUNT.find(G[U].location);
			if (it == COUNT.end())
			{
				COUNT[G[U].location] = U;
				if (!G[U].Enode)
				{
					if (degree(U, G) == 1 && !G[U].Enode)
						Inodes++;
					
					if (degree(U, G) == 3 && !G[U].Enode) 
						Ynodes++;
						
					if (degree(U, G) == 4 && !G[U].Enode)
						Xnodes++;
				}
				else if (G[U].Enode)
					Enodes++;
			}
			
			it = COUNT.find(G[u].location);
			if (it == COUNT.end())
			{
				COUNT[G[u].location] = u;
				if (!G[u].Enode)
				{
					if (degree(u, G) == 1 && !G[u].Enode)
						Inodes++;
					if (degree(u, G) == 3 && !G[u].Enode)
						Ynodes++;
					if (degree(u, G) == 4 && !G[u].Enode)
						Xnodes++;
				}
				else if (G[u].Enode)
					Enodes++;
			}

			if (degree(U, G) == 1 && degree(u, G) == 1) 
			{
				G[Eg].BranchType = "I-I";
				II++;
			}
		
			if ( degree(U, G) == 3 && degree(u, G) == 3) 
			{
				G[Eg].BranchType = "Y-Y";
				CC++;
			}

			if ( degree(U, G) == 4 && degree(u, G) == 4) 
			{
				G[Eg].BranchType = "X-X";
				CC++;
			}

			if ((degree(U, G) == 1 && degree(u, G) == 3) ||
				(degree(U, G) == 3 && degree(u, G) == 1))
				{
					G[Eg].BranchType = "I-Y";
					IC++;
				}

			if ((degree(U, G) == 1 && degree(u, G) == 4) ||
				(degree(U, G) == 4 && degree(u, G) == 1))
				{
					G[Eg].BranchType = "I-X";
					IC++;
				}

			if ((degree(U, G) == 3 && degree(u, G) == 4) ||
				(degree(U, G) == 4 && degree(u, G) == 3))
				{
					G[Eg].BranchType = "Y-X";
					CC++;
				}
					
//deal with edge nodes--------------------------------------------------
			if (G[U].Enode == true || G[u].Enode == true)
				Enodes++;
		
			if ((G[U].Enode == true && degree(u, G) == 1) ||
				(G[u].Enode == true && degree(U, G) == 1) )
					G[Eg].BranchType = "E-I";

			if ((G[U].Enode == true && degree(u, G) == 3) ||
				(G[u].Enode == true && degree(U, G) == 3) )
					G[Eg].BranchType = "E-Y";

			if ((G[U].Enode == true && degree(u, G) == 4) ||
				(G[u].Enode == true && degree(U, G) == 4) )
					G[Eg].BranchType = "E-X";

//total length of fault traces and area to be analysed------------------
			totalLength += G[Eg].length;
			geometry::append(all_lines, G[Eg].trace);
			NbB++;
		}
	}
	envBox = boost::geometry::return_envelope<box>(all_lines); 
	Area = geometry::area(envBox) * 1e-6 ; 
}



//Topolgy analysis of graph---------------------------------------------
void GRAPH::GraphAnalysis(Graph& G, VECTOR lines, int nb)
{
	assert (num_vertices(G) != 0 && num_edges(G) != 0);
	cout<< "\nGRAPH'S ANALYSIS OF NETWORK" << endl;
	/***********************************************************************
	* Graph analysis based on:												*
	* Sanderson, D. J., Peacock, D. C., Nixon, C. W., & Rotevatn, A. (2018)* 
	* Graph theory and the analysis of fracture networks.					*
	* Journal of Structural Geology.										*
	**********************************************************************/
	STATS stat;
	ofstream txtG;
	vector<line_type> faults = lines.data;
	
	int Inodes = 0,  Xnodes = 0, Ynodes = 0, Enodes = 0, NbB = 0, Nc, NbN, Nl, numK, II = 0 , IC = 0, CC = 0;
	float Cl, Cb, NCfreq, B20,  P20, P21 , B22;
	double Area, totalLength = 0, avLenL, avLenB;

	vertex_type U, u;
	vector <int> Fault_in_component;
	vector < double> Fault_length;
	vector <line_type> lineaments;

	//test whether graph is planer
	if (boyer_myrvold_planarity_test(G))
	{
		ClassifyEdges(G, Inodes, Ynodes, Xnodes, Enodes, II, IC, CC, NbB, totalLength, numK, Area);
		
		cout <<"Area: " << Area << endl;
		//Number of connections
		Nc = Ynodes + Xnodes;
		
		//Number of branches NB = (I + 3Y + 4X) / 2 
		NbN = (Inodes + 3*Ynodes + 4*Xnodes) / 2;

		//Number of lines NL = (I + 2Y) / 2 
		Nl = (Inodes + 2*Ynodes) / 2;

		//Average Length of lines L / NL 
		if (Nl != 0)
		{
			//average length of lines
			avLenL = totalLength / Nl;
			
			//Connections per line 2(Y+X) / NL 
			Cl = 2*(Ynodes + Xnodes) / Nl;
		}

		if (Area != 0)
		{
			//Connecting node frequency (NC km-2)
			NCfreq = Nc / Area;

			//Branch frequency (B20)
			B20 = NbN / Area;

			//Line frequency (P20)
			P20 = Nl / Area;

			//2D intesity (P21 = L / A 
			P21 = totalLength / Area;
		}

		if (NbN != 0)
		{
			//Average length of branches L/NB 
			avLenB = totalLength / NbN;
			
			if (avLenB != 0)
				//Dimensionless intesity B22 
				B22 = P21 * avLenB;
			
			//Conections per branch (3Y + 4X)/ NB 
			Cb = (3*Ynodes + 4*Xnodes) / NbN;
		}

//write results---------------------------------------------------------
		txtG = CreateGraphFileStream(lines.folder, lines.name);
		if (txtG.is_open())  
		{ 
			txtG<< "Nodes: " << "\t" 			 << num_vertices(G) << "\n"
				<< "EDGES: " << "\t" 			 << num_edges(G) << "\n"
				<< "Edgenodes: " 				 << "\t" << Enodes << "\n"
				<< "Inodes: " 			 		 << "\t" << Inodes <<  "\n" 
				<< "Ynodes: " 			 		 << "\t" << Ynodes <<  "\n" 
				<< "Xnodes: " 			 		 << "\t" << Xnodes  << "\n" 
				<< "Enodes: " 			 		 << "\t" << Enodes  << "\n"
				<< "II-Branches: "				 << "\t" << II << "\n"
				<< "IC-Branches: "				 << "\t" << IC << "\n"
				<< "CC-Branches: "				 << "\t" << CC << "\n"
				<<" Connections (Y +X): "		 << "\t" << (Ynodes + Xnodes) << "\n"
				<<" Total length:				"<< "\t" << totalLength << "\n"
				<<" Average length: "			 << "\t" << totalLength / NbB<< "\n" 
				<< "Connecting node frequency: " << "\t" << NCfreq << "\n"
				<< "Branch frequency: " << "\t"  << B20 << "\n"
				<< "Line frequency: " << "\t" 	 << P20 << "\n"
				<< "2D Intesnsity: " << "\t" 	 << P21 << "\n"
				<< "Dimensionless intesity: " 	 << "\t" << B22 << "\n"
				<< "Average degree of network: " << "\t" << (float) 2 * num_edges(G)/ num_vertices(G) << "\n"
				<< "Average connections: " 		 << "\t" << (float) 2 * (Xnodes + Ynodes) / num_vertices(G) << "\n"
				<< "Number of components (c): "  << "\t" << numK << "\n"
				<< "Number of faces (f): " << "\t" << num_edges(G) + numK - num_vertices(G) +1 << "\n" 
				<< "Density (d): " << "\t" << num_edges(G)*(totalLength*totalLength) /(4*Area) << "\n"
				<< "KDE edge orientation" << "\n";
				//	stat.KDE_estimation_strikes(Edges, txtG);
			txtG << endl;

//Analyse the different component of the network------------------------
		cout << " Analysing components of graph " << endl;
		
			for (int i =0; i < numK; i++)
			{
				Inodes = 0, Ynodes = 0, Xnodes = 0, Enodes = 0, NbB = 0, totalLength = 0;
				II = 0, IC =0, CC = 0;
				for (auto Eg : make_iterator_range(edges(G)))
				{
					U = source(Eg, G);
					u = target(Eg, G);
					
					if (G[U].component == i && G[u].component == i)
					{
						Fault_in_component.push_back(G[Eg].FaultNb);
						if (G[U].Enode == false || G[u].Enode == false)
						{
							if (G[Eg].BranchType == "I-I")
								II++;
					
							if (G[Eg].BranchType == "I-Y" || G[Eg].BranchType == "I-X")
								IC++;
					
							if (G[Eg].BranchType == "Y-Y" || G[Eg].BranchType == "X-X")
								CC++;
								
							if (degree(U, G) == 1 && G[U].Enode == false)
							Inodes++;

							if (degree(u, G) == 1 && G[u].Enode == false)
							Inodes++;

							if (degree(U, G) == 3 && G[U].Enode == false) 
							Ynodes++;

							if (degree(u, G) == 3 && G[u].Enode == false)
							Ynodes++;

							if (degree(U, G) == 4 && G[U].Enode == false)
							Xnodes++;

							if (degree(u, G) == 4 && G[u].Enode == false)
							Xnodes++;

							if (G[U].Enode == true || G[u].Enode == true)
							Enodes++;
						}
						totalLength += G[Eg].length;
						NbB++;
					}
				}
				
				NbN = (Inodes + 3*Ynodes + 4*Xnodes) / 2;
				if (NbN > nb)
				{
					sort( Fault_in_component.begin(), Fault_in_component.end() );
					Fault_in_component.erase( unique( Fault_in_component.begin(), Fault_in_component.end() ), Fault_in_component.end() );
					
					for (unsigned ii = 0; ii < Fault_in_component.size(); ii++)
					{
						Fault_length.push_back(geometry::length(faults.at(ii)));
						lineaments.push_back(faults.at(ii));
					}
					txtG<< "COMPONENT NO." << "\t" << i << "\n"
						<< "Branches:" << "\t" << NbB << "\n" 
						<< "Branches (calc):" << "\t" << NbN << "\n" 
						<< "Inodes: " << "\t" << Inodes <<  "\n" 
						<< "Ynodes: " << "\t" << Ynodes <<  "\n" 
						<< "Xnodes: " << "\t" << Xnodes  << "\n" 
						<< "Enodes: " << "\t" << Enodes  << "\n"
						<< "II-Branches: " << "\t" << II << "\n"
						<< "IC-Branches: " << "\t" << IC << "\n"
						<< "CC-Branches: " << "\t" << CC << "\n"
						<< "Connections (Y +X): " << "\t" << (Ynodes + Xnodes) << "\n"	
						<< " Total length:" << "\t" << totalLength << "\n"
						<< " Average length:" << "\t" << totalLength / NbB << "\n" 
						<< " Strike of underlying lineaments" << "\n";
					txtG<< "Length distribution of underlying lineaments" << "\n";
					txtG<< " KDE edge orientation" << "\n";
					txtG << endl;
				}
				Fault_in_component.clear();
				lineaments.clear();
				Fault_length.clear();
			}
		}
		else 
			cout << "ERROR: FAILED TO WRITE RESULTS!" << endl;
	}
	else 
	{
		cout << " ERROR: GRAPH IS NOT PLANAR!" << endl;
		exit (EXIT_FAILURE);
	}
	txtG.close();
}
  
//find shortest path between source and target--------------------------
void GRAPH::ShortPath(Graph G, map_vertex_type m, point_type source, point_type target, double radius)
{ 
	GEOMETRIE geom;
	GEO Gref;
	line_type SV , TV;
	Graph shortP;
	
	vertex_type S = AddNewVertex(m, source, G);
	vertex_type T = AddNewVertex(m, target, G);
	BUFFER BS = geom.DefinePointBuffer(source, radius);
	BUFFER BT = geom.DefinePointBuffer(target, radius);
	
	for (auto vd : boost::make_iterator_range(vertices(G))) 
	{
		if (geometry::within(G[vd].location, BS))
		{
			geometry::append(SV, source);
			geometry::append(SV, G[vd].location);
			AddNewEdge(G, S, vd, SV);
		}
		if(geometry::within(G[vd].location, BT))
		{
			geometry::append(TV, target);
			geometry::append(TV, G[vd].location);
			AddNewEdge(G, T, vd, TV);
		}
		geometry::clear(TV);
		geometry::clear(SV);
	}
	
	vector<double> distances( boost::num_vertices(G));
	vector<edge_type> path;
	vector<vertex_type> predecessors(boost::num_vertices(G));
	
	boost::dijkstra_shortest_paths(G, S,
									boost::weight_map(boost::get(&FEdge::length, G))
									.distance_map(boost::make_iterator_property_map(distances.begin(), boost::get(boost::vertex_index,G)))
									.predecessor_map(boost::make_iterator_property_map(predecessors.begin(), boost::get(boost::vertex_index,G)))
									);

	cout << "shortest path from " << S << " to " << T << ": ";
	vertex_type previous = T;
	bool finished = false;
	for(vertex_type current = predecessors[T]; !finished; previous = current, current = predecessors[current])
	{
		std::pair<edge_type, bool> edge_pair = boost::edge(current,previous,G);
		path.push_back( edge_pair.first );
		cout << edge_pair.first << ", ";
		finished = (current == S) || (current == predecessors[current]);
	}
	cout << endl;

	double distance = 0;
	for(std::vector<edge_type>::size_type i = 0; i != path.size(); i++) 
	{
		vertex_type u_tmp = boost::source(path[i], G);
		vertex_type v_tmp = boost::target(path[i], G);
		edge_type   e_tmp = boost::edge(u_tmp, v_tmp, G).first;
		
		add_vertex(G[u_tmp], shortP);
		add_vertex(G[v_tmp], shortP);
		add_edge(u_tmp, v_tmp, G[e_tmp], shortP);
		distance += G[e_tmp].length;
	}
	cout <<"Dijkstra shortest paths: " << distance << " ("<<path.size() << ") \n" << endl;

//	if (distance > 0)
		//Gref.WriteSHP(shortP, "ShortestPath.shp");
}

//create minimum spanning tree (backbone of network)--------------------
void  GRAPH::MinTree (Graph G)
{  
	Graph min_graph;
	GEO Gref;  
	int i = 0;
	std::string line;
	map_vertex_type map;
	line_type fault;
	
	std::vector<edge_type> spanning_tree;
	kruskal_minimum_spanning_tree( G, std::back_inserter(spanning_tree), 
									weight_map( get(&FEdge::length, G)) );

	for (std::vector<edge_type>::iterator ei = spanning_tree.begin();
		ei != spanning_tree.end(); ++ei) 
	{
		vertex_type S = source(*ei, G);
		vertex_type T = target(*ei, G);
		AddNewVertex(map, G[S].location, min_graph); 
		AddNewVertex(map, G[T].location, min_graph); 
		add_edge(S, T, 
				{geometry::length(G[*ei].trace)/1000, 
					G[*ei].trace},
					min_graph);
		i++;  
	}
	//Gref.WriteSHP(min_graph, "MinTree.shp");
	cout <<"minTree Total eges: " << i << endl;
}

VECTOR GRAPH::ComponentExtract(Graph G, VECTOR lines, int component)
{
	VECTOR extr_lines;
	vector <line_type> Extractor;
//here we extract the initial lineaments--------------------------------
	bool extract;
	for (auto Eg : make_iterator_range(edges(G)))
	{
		extract = true;
		if (G[Eg].component == component)
		{
			line_type extr_line = G[Eg].trace;
			for (auto ext_line : Extractor)
			{
				if (geometry::equals(ext_line, extr_line))
				{
					extract = false;
					break;
				}
			}
			if (extract)
				Extractor.push_back(lines.data.at(G[Eg].FaultNb));
		}
	}
	cout << " Found " << Extractor.size() << " lineaments for component " << component << endl;
	
	//build the struct from the extracted data and return it------------
	extr_lines.name = lines.name + "_component_" + to_string(component);
	extr_lines.refWKT = lines.refWKT;
	for (auto i : Extractor)
		extr_lines.data.push_back(i);
	return(extr_lines);
}
