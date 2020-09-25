#include "../include/graph.h"
#include "../include/geometrie.h"
#include "../include/GeoRef.h"
#include "../include/stats.h"
#include "../include/util.h"

#include <random>

const std::string graph_subdir="graph";

GRAPH::GRAPH ()

{
}

namespace fs = boost::filesystem;

// // create and open filestream in folder "statistics"
// ofstream CreateGraphFileStream(string name)
// {
// // 	string folder_name = folder + "/graph/";
// // 	const char* stats_dir = folder_name.c_str();
// // 		if(!opendir(stats_dir))
// 		
// // 	mkdir(stats_dir, 0777);
//     
//     fs::path file(name);
//     system::error_code ec;
//     fs::create_directories(file.parent_path(), ec);
//     
// // 	string tsv_file = stats_dir + name + (string);
// 	ofstream txtF; 
// 	txtF.open (name, ios::out | ios::app); 
// 	return(txtF);
// }


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
 bool GRAPH::AddNewEdge(Graph& G, vertex_type S, vertex_type T, line_type FaultSeg)
{
	return GRAPH::AddNewEdge(G, S, T, FaultSeg, geometry::length(FaultSeg));
}

//add a new edge to the graph, and record the full length of the fault that comprises this fault segment
bool GRAPH::AddNewEdge(Graph& G, vertex_type S, vertex_type T, line_type FaultSeg, double FaultLength)
{
	edge_type edge;
	bool added_new = false;
    if (!boost::edge(S, T, G).second && !geometry::equals(G[S].location, G[T].location))
    {
			std::tie(edge, added_new) = boost::add_edge(S, T,
				{geometry::length(FaultSeg)/1000, FaultSeg}, G);
			G[edge].fault_length = FaultLength;
	}
// 	cout << "Adding edge between " << S << " and " << T << " with number of vertices " << num_vertices(G) << ", successful ? " << added_new << endl;
	return added_new;
}

//read vector data into graph------------------------------------------
graph_map<point_type, vertex_type, Graph> GRAPH::ConvertLinesToGraph(std::vector<line_type> &faults, const char *refWKT, double distance_threshold)
{
    graph_map<point_type, vertex_type, Graph> map(distance_threshold, refWKT);
    Graph &graph = map.get_graph();
	vertex_type VA, VB;

	//create edges and nodes for faults in the graph-------------------------
	int nb = 0;
	BOOST_FOREACH(line_type f, faults)
	{
		VA = map.add_vertex(f.front());
		VB = map.add_vertex(f.back());

		AddNewEdge(graph, VA, VB, f);
//         std::cout <<" Adding edge from " << VA << " to " << VB << ", there are now " << num_edges(graph) << " in the graph" << std::endl;
		
		if (boost::edge(VA,VB, graph).second)
		{
			edge_type e = boost::edge(VA,VB, graph).first;
				graph[e].FaultNb = nb;
		}
		nb++;
	}
	cout << "Converted " << faults.size() << " lineaments into " << num_edges(graph) << " edges." << endl;
    return map;
}


graph_map<> GRAPH::ReadVEC4raster(double transform[8], VECTOR &lines, double map_distance_threshold)
{
    std::vector<line_type> &faults = lines.data;
	GEO g;
	GEOMETRIE geom;
    graph_map<> map(map_distance_threshold, lines.refWKT);
    Graph& graph = map.get_graph();
	geometry::model::multi_linestring<line_type> intersection;
	vector<std::tuple< std::pair<point_type, point_type>, line_type, unsigned int, double >> G;
	std::pair <point_type, point_type> NODES;
	line_type TRACE;
	long double Xmax, Xmin, Ymax, Ymin;
	vertex_type VA, VB;
	unsigned int type;

	const unsigned int FRONT_ENODE = 1 << 0;
	const unsigned int BACK_ENODE  = 1 << 1;
	const double endpoint_threshold = 0.5;
	
	
//we crop the graph 
	polygon_type pl = g.BoundingBox(transform, transform[1]);


//check for edge nodes (crop faults by raster area)---------------------
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
		
		VA = map.add_vertex(NODES.first);
		VB = map.add_vertex(NODES.second);

		//if( (degree(VA, graph) == 0) && (degree(VB, graph) == 0) )
		AddNewEdge(graph, VA, VB, TRACE, fault_length);
				
		if (type & FRONT_ENODE) graph[VA].Enode = true;
		if (type &  BACK_ENODE) graph[VB].Enode = true;
	}

	cout << " Converted " << faults.size() << " faults into " << num_edges(graph) << " edges" << endl;
	
	faults.clear();
	for (auto Eg : make_iterator_range(edges(graph)))
		faults.push_back(graph[Eg].trace);
    return map;
}

Graph GRAPH::ReadVEC4MODEL(VECTOR &lines, box bx, double map_distance_threshold)
{
    std::vector<line_type> &faults = lines.data;
	GEO g;
	Graph graph;
	GEOMETRIE geom;
	graph_map map(graph, map_distance_threshold, lines.refWKT);
	
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

		VA = map.add_vertex(NODES.first);
		VB = map.add_vertex(NODES.second);

		AddNewEdge(graph, VA, VB, TRACE);
				
		if (type & FRONT_ENODE) graph[VA].Enode = true;
		if (type &  BACK_ENODE) graph[VB].Enode = true;
	}

	if (num_edges(graph) > 0)
	{
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
void GRAPH::RemoveSpurs(graph_map<>& map, double minDist)
{
	cout << "Removing spurs" << endl;
    Graph& G = map.get_graph();
	point_type rmP;
	vertex_type U ,u;

	restart:
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
// 					rmP.set<0>((long long).x());
// 					rmP.set<1>((long long)G[U].location.y());
// 					remove_vertex(U, G);
					map.remove_vertex(G[U].location);
				}
				if (degree(u,G) == 0)
				{
// 					cout << "WARNING: Removing edge to vertex " << endl;
// 					rmP.set<0>((long long).x());
// 					rmP.set<1>((long long)G[u].location.y());
// 					remove_vertex(u, G);
					map.remove_vertex(G[u].location);
				}
				goto restart;
			}
		}
	}
	//taking this out entirely. the graph *should* be planar. if it isn't, we should fix what is wrong witht he code.
// 	if (!boyer_myrvold_planarity_test(G))
// 	{
// 		cout << "WARNING: Graph is non-planar!\n"
// 			 << "Attempting to re-split edges with r = " << split_dist/2 << endl;
// 		graph_map<> resplit_map = SplitFaults(map, split_dist/2);
//         map = resplit_map;
// 		if (!boyer_myrvold_planarity_test(G))
// 		{
// 			cout << "ERROR: Could not construct planar graph" << endl;
// 			exit(EXIT_FAILURE);
// 		}
// 	}
	cout << " done \n" << endl;
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
graph_map<> GRAPH::SplitFaults(graph_map<> &map, double minDist)
{
    Graph& graph = map.get_graph();
	cout << "Splitting edges at intersection points" << endl;
	graph_map<> split_map(map.get_dist(), map.get_refWKT()); //the map translates physical coordinates and nodes in the graph.
    Graph &split_graph = split_map.get_graph();
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
		vertex_type prev_vertex = split_map.add_vertex(std::get<1>(cross.front()));
		for (vector<std::tuple<long double, point_type, AttachPoint>>::const_iterator I = cross.begin(); I != cross.end(); I++)
		{
			if (I == cross.begin()) continue;
			point_type intersect = get<1>(*I);
			//bool is_start = geometry::distance(fault.front(), intersect) <= minDist;
			//bool is_end   = geometry::distance(fault.back(),  intersect) <= minDist;
			//if (is_start || is_end) continue;
			NewV = split_map.add_vertex(intersect);
			if (NewV == prev_vertex) continue;
			
			AddNewEdge(split_graph, prev_vertex, NewV, geom.GetSegment(fault, split_graph[prev_vertex].location, split_graph[NewV].location), fault_length);//also remember the length of the original fault
			
			if(boost::edge(prev_vertex, NewV, split_graph).second) //if an edge exists between prev_vertex and NewV
			{
				auto e = boost::edge(prev_vertex, NewV, split_graph).first; //edge descriptor of the above
				graph[e].FaultNb = N; //set the fault number
			}
			prev_vertex = NewV;
		}

		cross.clear();
	}
// 	graph = g;
// 	map = m;
	cout << " done \n" << endl;
    return split_map;
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
void GRAPH::GraphAnalysis(Graph& G, VECTOR lines, int nb, const double angle_param_penalty, string out_filename)
{
	assert (num_vertices(G) != 0 && num_edges(G) != 0);
	cout<< "GRAPH'S ANALYSIS OF NETWORK" << endl;
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
	vector <line_type> Fault_in_component;
	vector < double> Fault_length;
	vector <line_type> lineaments;

	//test whether graph is planer
	if (!boyer_myrvold_planarity_test(G))
		cout << " WARNING: GRAPH IS NOT PLANAR! Proceed with caution" << endl;
	
		ClassifyEdges(G, Inodes, Ynodes, Xnodes, Enodes, II, IC, CC, NbB, totalLength, numK, Area);
		
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
		//string graphFile =  + ;
//         cout << "saving graph stats data with name " << name << ", lines folder: " << lines.folder << " and lines name: " << lines.name << endl;
        string save_name = FGraph::add_prefix_suffix_subdirs(lines.out_path, {graph_subdir}, "graph_statistics", ".tsv", true); //we need to clean this up //lines.folder
//         cout << "the resulting name is " << save_name << endl;
		txtG = FGraph::CreateFileStream(save_name);
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
		
			for (int i = 0; i < numK; i++)
			{
				Inodes = 0, Ynodes = 0, Xnodes = 0, Enodes = 0, NbB = 0, totalLength = 0;
				II = 0, IC =0, CC = 0;
				for (auto Eg : make_iterator_range(edges(G)))
				{
					U = source(Eg, G);
					u = target(Eg, G);
					
					if (G[U].component == i && G[u].component == i)
					{
						Fault_in_component.push_back(G[Eg].trace);
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
				if (NbB > nb)
				{
					cout << "Component no " << i << " containing " << NbB << " branches " << endl;
					string comp = lines.name + "_Graph_component_" + std::to_string(i);// + name
					VECTOR comp_Lineamants;
					comp_Lineamants.folder = lines.folder;
					comp_Lineamants.name = comp;
					comp_Lineamants.refWKT = lines.refWKT;
					comp_Lineamants.data = Fault_in_component ;
                    comp_Lineamants.out_folder = lines.out_folder;
                    comp_Lineamants.out_path = lines.out_path;
                    comp_Lineamants.in_path = lines.in_path;
					STATS stats;
					stats.GetLengthDist(comp_Lineamants);
					stats.KDE_estimation_strikes(comp_Lineamants, angle_param_penalty);

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
						<< " Strike of underlying lineaments" << endl;
				}
				Fault_in_component.clear();
				lineaments.clear();
				Fault_length.clear();
			}
		}
		else 
			cout << "ERROR: FAILED TO WRITE RESULTS!" << endl;

	txtG.close();
}
  
//find shortest path between source and target--------------------------
Graph GRAPH::ShortPath(graph_map<> m, std::string in_filename, std::string out_filename)
{ 
	GEO Gref;
	Graph shortP;
    Graph &G = m.get_graph();
	GEOMETRIE geom;
	line_type SV , TV;
	point_type source, target;
	
	Gref.Get_Source_Target(in_filename.c_str(), source, target);

    bool added_source_vertex, added_target_vertex;
    vertex_type S, T;
    
	std::tie(S, added_source_vertex) = m.add_vertex_isnew(source);
	std::tie(T, added_target_vertex) = m.add_vertex_isnew(target);
// 	BUFFER BS = geom.DefinePointBuffer(source, radius);
// 	BUFFER BT = geom.DefinePointBuffer(target, radius);
    
    double s_dist = std::numeric_limits<double>::infinity(), t_dist = std::numeric_limits<double>::infinity();
    vertex_type s_nearest, t_nearest;
    
	for (auto vd : boost::make_iterator_range(vertices(G))) //these vertices aren't being added properly //they're being attached to theirselves
	{
		if ((S != vd)/* && geometry::within(G[vd].location, BS)*/)
		{
            double dist = geometry::distance(G[S].location, G[vd].location);
            if (dist < s_dist)
            {
                s_dist = dist;
                s_nearest = vd;
            }
		}
		if((T != vd)/* && geometry::within(G[vd].location, BT)*/)
		{
            double dist = geometry::distance(G[T].location, G[vd].location);
            if (dist < t_dist)
            {
                t_dist = dist;
                t_nearest = vd;
            }
		}
// 		geometry::clear(TV);
// 		geometry::clear(SV);
	}
	
    cout << "Adding source and target vertices for shortest path" << endl;
	
    geometry::append(SV, source);
    geometry::append(SV, G[s_nearest].location);
    bool added_source = AddNewEdge(G, S, s_nearest, SV);
    cout << "Added source? " << added_source << endl;
    
    
    geometry::append(TV, target);
    geometry::append(TV, G[t_nearest].location);
    bool added_target = AddNewEdge(G, T, t_nearest, TV);
    cout << "Added target? " << added_target << endl;
	
	vector<double> distances( boost::num_vertices(G));
	vector<edge_type> path;
	vector<vertex_type> predecessors(boost::num_vertices(G));
	
	boost::dijkstra_shortest_paths(G, S,
									boost::weight_map(boost::get(&FEdge::length, G))
									.distance_map(boost::make_iterator_property_map(distances.begin(), boost::get(boost::vertex_index,G)))
									.predecessor_map(boost::make_iterator_property_map(predecessors.begin(), boost::get(boost::vertex_index,G)))
									);

	cout << "Calculating shortest path from " << S << " to " << T << ": "; //predecessors is mapping each node to itself
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
	for(std::vector<edge_type>::size_type i = 0; i != path.size(); i++) //path is only length 1
	{
		vertex_type u_tmp = boost::source(path[i], G);
		vertex_type v_tmp = boost::target(path[i], G);
		edge_type   e_tmp = boost::edge(u_tmp, v_tmp, G).first;

		auto ud = add_vertex(G[u_tmp], shortP);
		auto vd = add_vertex(G[v_tmp], shortP);
        
//         cout << "Added vertices " << ud << " -> " << vd << endl;

		auto added_edge = add_edge(ud, vd, G[e_tmp], shortP); //<- this causes the segmentation fault //u_tmp and v_tmp are the same

		distance += G[e_tmp].length;
	}
	cout <<"Dijkstra shortest paths: " << distance << " ("<<path.size() << ") \n" << endl;

	if (distance > 0)
	{
		vector<line_type> edges_shortPath;
		for (auto Eg : make_iterator_range(edges(shortP)))
			edges_shortPath.push_back(shortP[Eg].trace);
        if (out_filename != "")
        {
            std::string save_filename = FGraph::add_prefix_suffix_subdirs(out_filename, {graph_subdir});
            Gref.WriteSHP_lines(edges_shortPath, m.get_refWKT(), save_filename);
        }
	}
	
	//remove the source and target vertices that we added earlier (if we added them)
	if (added_source_vertex) m.remove_vertex(source);
    if (added_target_vertex) m.remove_vertex(target);
	
	cout << " done " << endl;
	return(shortP);
}

//create minimum spanning tree (backbone of network)--------------------
Graph GRAPH::MinTree (graph_map<> gm, double map_dist_threshold, std::string filename)
{  
	cout << "Generating minimum spanning tree" << endl;
    Graph &G = gm.get_graph();
	Graph min_graph;
	GEO Gref;  
	int i = 0;
	std::string line;
	graph_map map(min_graph, map_dist_threshold);
	line_type fault;
	
	std::vector<edge_type> spanning_tree;
	kruskal_minimum_spanning_tree( G, std::back_inserter(spanning_tree), 
									weight_map( get(&FEdge::length, G)) );

	for (std::vector<edge_type>::iterator ei = spanning_tree.begin();
		ei != spanning_tree.end(); ++ei) 
	{
		vertex_type S = source(*ei, G);
		vertex_type T = target(*ei, G);
		map.add_vertex(G[S].location); 
		map.add_vertex(G[T].location); 
		add_edge(S, T, 
				{geometry::length(G[*ei].trace)/1000, 
					G[*ei].trace},
					min_graph);
		i++;  
	}
	cout <<"minTree Total eges: " << i << endl;
	vector<line_type> edges_minTree;
	for (auto Eg : make_iterator_range(edges(min_graph)))
		edges_minTree.push_back(min_graph[Eg].trace);
    
	if (filename != "")
    {
        std::string save_filename = FGraph::add_prefix_suffix_subdirs(filename, {graph_subdir});
        Gref.WriteSHP_lines(edges_minTree, gm.get_refWKT(), save_filename);
    }
    
	cout << " done \n" << endl;
	return(min_graph);
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

void GRAPH::IntersectionMap(Graph G, VECTOR lines, float cell_size, float search_size)
{
//create intersection density maps
//(one qualitative(intersections per area) and one quantitative(sum of vertex degrees per area))
	GEO georef;
	GEOMETRIE geom;
	vector<pair<point_type, int>> graph_vertices;

	box AOI = geom.ReturnAOI(lines.data);
	polygon_type t_AOI = geom.Return_tigth_AOI(lines.data);
	
	double min_x = geometry::get<geometry::min_corner, 0>(AOI);
	double min_y = geometry::get<geometry::min_corner, 1>(AOI);
	double max_x = geometry::get<geometry::max_corner, 0>(AOI);
	double max_y = geometry::get<geometry::max_corner, 1>(AOI);

	point_type ll(min_x, min_y);
	point_type ul(min_x, max_y);
	point_type lr(max_x, min_y);
	
	double newGeoTransform[6] = {min_x, cell_size, 0 , max_y, 0, (cell_size * (-1))};
	
	int x_size = (long int) ceil((geometry::distance(ll, lr) / cell_size));
	int y_size = (long int) ceil((geometry::distance(ll, ul) / cell_size));
	vector<vector<double> > vec(x_size , vector<double> (y_size, 0));  //count of intersections in area
	vector<vector<double> > vec2(x_size , vector<double> (y_size, 0));  //sum of degree of intersections in area
	
	double cur_y = max_y;
	double cur_x = min_x;
	
	cout << "Calulating intersection density for raster with size \n"
		 << vec.size()<< " x " << vec[0].size() << endl;

	for (auto ve : make_iterator_range(vertices(G)))
	{
		if (degree(ve, G) > 2)
			graph_vertices.push_back(make_pair(G[ve].location, degree(ve, G)));
	}

	typedef pair<box, decltype(graph_vertices)::iterator> box_point; 
	vector<box_point> result;
	geometry::index::rtree<box_point, geometry::index::rstar<16>> IntersecTree;
	
	for (auto it = graph_vertices.begin(); it < graph_vertices.end(); it++)
	{
			BUFFER search_box = geom.DefineSquareBuffer(it->first, 1);
			box bounding_box = boost::geometry::return_envelope<box>(search_box);
			IntersecTree.insert(make_pair(bounding_box, it));
	}

	progress_display * show_progress =  new boost::progress_display(x_size * y_size);
	for (int i = 0; i < x_size; i++)
	{
		cur_y = max_y;
		for (int j = 0; j < y_size; j++)
		{
			point_type cur_pos(cur_x, cur_y); //cur_pos is the top-left corner of the pixel
			point_type minBox(cur_x, (cur_y - cell_size));
			point_type maxBox((cur_x + cell_size), cur_y );
			box pixel(minBox, maxBox);
			if (!geometry::disjoint(pixel,AOI))
			{
				point_type minBox(cur_x, (cur_y - cell_size));
				point_type maxBox((cur_x + cell_size), cur_y );
				box pixel(minBox, maxBox);
				point_type centre;
				boost::geometry::centroid(pixel, centre);
				
				BUFFER cicular_window = geom.DefinePointBuffer(centre, search_size);
				IntersecTree.query(!geometry::index::disjoint(cicular_window), std::back_inserter(result));
				
				int d  = 0;
				int d2 = 0;
				
				for (auto r = result.begin(); r < result.end(); r++)
				{
					if (geometry::within(r->second->first, cicular_window))
					{
						d++;
						d2 += r->second->second;
					}
				}
				
				vec[i][j]  = d;
				vec2[i][j] = d2;
			}
			else
			{
				vec[i][j]  = -256;
				vec2[i][j] = -256;
			}
			cur_y-= cell_size;
			result.clear();
			 ++(*show_progress);
		}
		cur_x += cell_size;
		}

//write the raster file---------------------------------------------
    std::string isec_dens_name = FGraph::add_prefix_suffix_subdirs(lines.out_path, {graph_subdir}, "intersection_density", ".tif");
	georef.WriteRASTER(vec,  lines.refWKT, newGeoTransform, lines, isec_dens_name);
    
    std::string isec_intens_name = FGraph::add_prefix_suffix_subdirs(lines.out_path, {graph_subdir}, "intersection_intensity", ".tif");
	georef.WriteRASTER(vec2, lines.refWKT, newGeoTransform, lines, isec_intens_name);
	cout << " done \n" << endl;
}

void GRAPH::ClassifyLineaments(Graph G, VECTOR &lines, AngleDistribution &angle_dist, float dist, string name)
{
	cout << "Classifying line set based on orientation and intersections" << endl;
	std::clock_t startcputime = std::clock();
	auto wcts = std::chrono::system_clock::now();

	GEOMETRIE geom;
	typedef std::tuple<decltype(lines.data)::iterator, int, int, int, int> deg_v; //storing vertex degrees per line [ line iterator 1 - 3 - 4 - intersections]
	typedef std::pair<point_type, vertex_iter> g_points; 
	
	vector<deg_v> lines_classified;
	geometry::index::rtree<g_points, geometry::index::rstar<16>> point_tree;

	for (auto it = lines.data.begin(); it < lines.data.end(); it++)
		lines_classified.push_back(make_tuple(it, 0, 0, 0, 0));
	
	vertex_iter v, vend;
	for (boost::tie(v, vend) = vertices(G); v != vend; ++v)
		point_tree.insert(std::make_pair(G[*v].location, v));
	
	for (auto l = lines.data.begin(); l < lines.data.end(); l++)
	{
		box bounding_box = boost::geometry::return_envelope<box>(*l);
		std::vector<g_points> candidates;
		line_type line = *l;
		
		auto it = std::find_if(lines_classified.begin(), lines_classified.end(), [&line](const deg_v& e) {return  geometry::equals(*std::get<0>(e), line);});
		if (it != lines_classified.end()) 
		{
			point_tree.query(!geometry::index::disjoint(bounding_box), std::back_inserter(candidates));
				for (auto candidate = candidates.begin(); candidate < candidates.end(); candidate++)
				{
					BUFFER b = geom.DefinePointBuffer(G[*candidate->second].location, dist);
					if (geometry::intersects(line, b))
					{
						if (degree(*candidate->second, G) == 1)
							std::get<1>(*it)++; 

						if (degree(*candidate->second, G) == 3)
							std::get<2>(*it)++; 

						if (degree(*candidate->second, G) == 4)
							std::get<3>(*it)++; 

						if ( degree(*candidate->second, G) == 3 || degree(*candidate->second, G) == 4 )
							std::get<4>(*it)++; 
					}
				}
		}
	}
	
	double cpu_duration = (std::clock() - startcputime) / (double)CLOCKS_PER_SEC;
	std::chrono::duration<double> wctduration = (std::chrono::system_clock::now() - wcts);
	
	std::cout << "Finished in " << cpu_duration << " seconds [CPU Clock] " << std::endl;
	std::cout << "Finished in " << wctduration.count() << " seconds [Wall Clock]" << std::endl;
	
//write a shp file containing the classification------------------------
	STATS stats;
	GDALAllRegister();
	const char* system = lines.refWKT;
	const char *pszDriverName = "ESRI Shapefile";
    GDALDriver *poDriver;
    GDALDataset *poDS;
    OGRLayer *poLayer;
    
    name.append(".shp");
    const char* Name = name.c_str();
    
    cout << Name << endl;
    OGRSpatialReference oSRS;
    oSRS.importFromWkt(&system);

    poDriver = (GDALDriver*) GDALGetDriverByName(pszDriverName );
    if( poDriver == NULL )
    {
        printf( "%s driver not available.\n", pszDriverName );
        exit( 1 );
    }
    
    poDriver->SetDescription("bla");

	poDS = poDriver->Create(Name, 0, 0, 0, GDT_Unknown, NULL );
	if( poDS == NULL )
	{
		printf( "Creation of output file failed.\n" );
		exit( 1 );
	}
	
	poLayer = poDS->CreateLayer( "graph_branches", &oSRS, wkbLineString, NULL );
	if( poLayer == NULL )
	{
		printf( "Layer creation failed.\n" );
		exit( 1 );
	}
	
	OGRFieldDefn oField( "ID", OFTInteger );
	oField.SetWidth(10);
	if( poLayer->CreateField( &oField ) != OGRERR_NONE )
	{
		printf( "Creating 'No' field failed.\n" );
		exit( 1 );
	}
	
	OGRFieldDefn oField0( "Length", OFTReal );
	oField0.SetWidth(50);
	oField0.SetPrecision(5);
	if( poLayer->CreateField( &oField0 ) != OGRERR_NONE )
	{
		printf( "Creating 'Length' field failed.\n" );
		exit( 1 );
	}
	
	OGRFieldDefn oField1( "Angle", OFTReal );
	oField1.SetWidth(10);
	oField1.SetPrecision(3);
	if( poLayer->CreateField( &oField1 ) != OGRERR_NONE )
	{
		printf( "Creating 'Angle' field failed.\n" );
		exit( 1 );
	}
	
	if (angle_dist.gaussians.size() != 0)
	{
		OGRFieldDefn oField2( "Set", OFTInteger );
		oField1.SetWidth(10);
		if( poLayer->CreateField( &oField2 ) != OGRERR_NONE )
		{
			printf( "Creating 'Set' field failed.\n" );
			exit( 1 );
		}
	}
//vertex counts as 
		OGRFieldDefn oField3( "I", OFTInteger );
		oField1.SetWidth(10);
		if( poLayer->CreateField( &oField3 ) != OGRERR_NONE )
		{
			printf( "Creating 'I' field failed.\n" );
			exit( 1 );
		}
	
		OGRFieldDefn oField4( "Y", OFTInteger );
		oField1.SetWidth(10);
		if( poLayer->CreateField( &oField4 ) != OGRERR_NONE )
		{
			printf( "Creating 'Y' field failed.\n" );
			exit( 1 );
		}
		
		OGRFieldDefn oField5( "X", OFTInteger );
		oField1.SetWidth(10);
		if( poLayer->CreateField( &oField5 ) != OGRERR_NONE )
		{
			printf( "Creating 'X' field failed.\n" );
			exit( 1 );
		}
		
		OGRFieldDefn oField6( "totals x y", OFTInteger );
		oField1.SetWidth(10);
		if( poLayer->CreateField( &oField6 ) != OGRERR_NONE )
		{
			printf( "Creating 'totals x y' field failed.\n" );
			exit( 1 );
		}
	
	int id = 0;
	for (auto it = lines_classified.begin(); it != lines_classified.end(); ++it) 
	{
		line_type line = *std::get<0>(*it);
		OGRLineString l;
		OGRFeature *poFeature;
		geometry::unique(line);
		float L = (float) geometry::length(line);
		double strike = (double)(atan2(line.front().x() - line.back().x(), line.front().y() - line.back().y())) 
						* (180 / math::constants::pi<double>());
						
		if (strike  < 0) 
			strike  += 180;

		poFeature = OGRFeature::CreateFeature( poLayer->GetLayerDefn() );
		poFeature->SetField( "ID", id );
		poFeature->SetField( "Length", L);
		poFeature->SetField( "Angle", strike);
		
		if (angle_dist.gaussians.size() != 0)
			poFeature->SetField( "Set", stats.CheckGaussians(angle_dist, strike));
		
		poFeature->SetField( "I",  std::get<1>(*it));
		poFeature->SetField( "Y",  std::get<2>(*it));
		poFeature->SetField( "X",  std::get<3>(*it));
		poFeature->SetField( "totals x y", std::get<4>(*it));
		
		BOOST_FOREACH(point_type P, line) 
			l.addPoint(P.x(), P.y());
		poFeature->SetGeometry( &l );
		
		if( poLayer->CreateFeature( poFeature ) != OGRERR_NONE )
		{
			printf( "Failed to create feature in shapefile.\n" );
			exit( 1 );
		}
		OGRFeature::DestroyFeature( poFeature );
		id++;
	}
    GDALClose( poDS );
    cout << " done \n" << endl;
}

//parameters: graph, pressure1 , pressure2, adn bool whether vertical or horizontal gradient (vertical if true)
void AssignGrad(Graph G, float p1, float p2, bool vert, const char *refWKT)
{
	//creating a gradient raster with a buffer of 5 cell sizes around the AOI
	GEO georef;
	GEOMETRIE geom;

	//float** values[x_dim][y_dim];
	RASTER<float> raster;
	vector<line_type> lines;
	
	for (auto Eg : boost::make_iterator_range(edges(G))) 
		lines.push_back(G[Eg].trace);
	line_type l = geom.ShortestLine(lines);
	double min_l = floor(geometry::length(l)/3);
	cout << "Aiming for cell size of: " << min_l << endl;
	
	box AOI = geom.ReturnAOI(lines);
	
	double min_x = geometry::get<geometry::min_corner, 0>(AOI) + 5*min_l ; 
	double max_y = geometry::get<geometry::max_corner, 1>(AOI) - 5*min_l ; 

	point_type ll ,lr, ur;
	ll.set<0>(geometry::get<geometry::min_corner, 0>(AOI) - 5*min_l );
	ll.set<1>(geometry::get<geometry::min_corner, 1>(AOI) + 5*min_l);
	
	lr.set<0>(geometry::get<geometry::max_corner, 0>(AOI) + 5*min_l );
	lr.set<1>(geometry::get<geometry::min_corner, 1>(AOI) - 5*min_l );
	
	ur.set<0>(geometry::get<geometry::max_corner, 0>(AOI) + 5*min_l );
	ur.set<1>(geometry::get<geometry::max_corner, 1>(AOI) + 5*min_l );
	
	int x_dim = geometry::distance(ll, lr) / min_l ;
	int y_dim = geometry::distance(lr, ur) / min_l ;

	cout << "Building raster with dimesnions: "<< x_dim << " " << y_dim << endl;
	raster.transform[0] = ll.x();
	raster.transform[1] = min_l;
	raster.transform[2] = 0;
	raster.transform[3] = ur.y();
	raster.transform[4] = 0;
	raster.transform[5] = -min_l;
	raster.transform[6] = x_dim;
	raster.transform[7] = y_dim;
	
	//set reference to teh one of the shp file
	raster.refWKT = refWKT;
	
	float change;
	if (p1 > p2)
		change = -(p1 - p2) / x_dim;
	
	if (p1 < p2)
		change = (p2 - p1) / x_dim;
	
	if (p1 == p2)
	{
		cout << "no differntial pressure" << endl;
		change = 0;
	}
	
	raster.values = new float*[x_dim];
	for (int i = 0; i < x_dim ; i++)
		raster.values[i] = new float[y_dim];
		
	float value = p1;
	for (int x = 0; x < x_dim ; x++)
	{
		if (vert)
			value = p1;
		for (int y = 0; y < y_dim ; y++)
		{
			raster.values[x][y] = value;
			if (vert)
				value += change;
		}
		if (!vert)
			value += change;
	}
	for (auto Ve : boost::make_iterator_range(vertices(G)))
	{
		G[Ve].data = georef.getValue(G[Ve].location, raster.transform, raster.values);
		if (std::isnan(G[Ve].data)) 
		{
			cout <<"Error: vertex " << Ve << " read a value of nan" << endl;
				cout << G[Ve].location.x() << " " << G[Ve].location.y() << endl;
			}
	}
	delete raster.values;
}

void GRAPH::MaximumFlow_R(Graph G, string st_filename, string capacity_type, const char *refWKT, std::string out_filename)
{
	GEO georef;
	point_type s, t;
	georef.Get_Source_Target(st_filename.c_str(), s, t);
	cout<< "Maximum flow with raster data." << endl;
	
	DGraph dg = georef.MakeDirectedGraph(G);
	georef.setup_maximum_flow(dg, capacity_type);
	double mf =  georef.maximum_flow(dg, s, t);
	cout << "maximum flow is: " << mf << endl;
	if (out_filename != "")
    {
        std::string name = FGraph::add_prefix_suffix_subdirs(out_filename, {graph_subdir}, "max_flow_R_");//
        georef.WriteSHP_maxFlow(dg, refWKT, name.c_str());
    }
	cout << " done \n" << endl;
}

void GRAPH::MaximumFlow_VG(Graph G, string st_filename, float top, float bottom, string capacity_type, const char *refWKT, std::string out_filename)
{
	GEO georef;
	point_type s, t;
	georef.Get_Source_Target(st_filename.c_str(), s, t);
	cout<< "Maximum flow with vertical gradient: " << top << "-" << bottom << endl;  
	AssignGrad(G, top, bottom, false, refWKT);
	
	DGraph dg = georef.MakeDirectedGraph(G);
	georef.setup_maximum_flow(dg, capacity_type);
	double mf =  georef.maximum_flow(dg, s, t);
	cout << "maximum flow is: " << mf << endl;
    if (out_filename != "")
    {
        std::string name = FGraph::add_prefix_suffix_subdirs(out_filename, {graph_subdir}, "max_flow_VG_");
        georef.WriteSHP_maxFlow(dg, refWKT, name.c_str());
    }
	cout << " done \n" << endl;
}

void GRAPH::MaximumFlow_HG(Graph G, string st_filename, float left, float right, string capacity_type, const char *refWKT, std::string out_filename)
{
	GEO georef;
	point_type s, t;
	georef.Get_Source_Target(st_filename.c_str(), s, t);
	cout<< "Maximum flow with horizontal gradient: " << left << "-" << right << endl;  
	AssignGrad(G, left, right, false, refWKT);
	
	DGraph dg = georef.MakeDirectedGraph(G);
	georef.setup_maximum_flow(dg, capacity_type);
	double mf =  georef.maximum_flow(dg, s, t);
    if (out_filename != "")
    {
        string name = FGraph::add_prefix_suffix_subdirs(out_filename, {graph_subdir}, "max_flow_HG_");
        georef.WriteSHP_maxFlow(dg, refWKT, name.c_str());
    }
	cout << " done \n" << endl;
}
