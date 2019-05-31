#include "graph.h"
#include "geometrie.h"
#include "GeoRef.h"

#include <random>

GRAPH::GRAPH ()

{

}

//add new vertex to graph-----------------------------------------------
vertex_type GRAPH::AddNewVertex(map_vertex_type& map, point_type const& key, Graph& graph)
 {
	const double threshold = 1;
	point_int check;
	check.set<0>((int)key.x());
	check.set<1>((int)key.y());
	typename map_vertex_type::iterator it = map.find(check);
	
	if (it != map.end())
	{
		std::vector<vertex_type> &possible = it->second;
		for (auto pit = possible.begin(); pit != possible.end(); pit++)
		{
			if (geometry::distance(graph[*pit].location, key) <= threshold) return *pit;
		}
		vertex_type new_vertex = add_vertex(FVertex<point_type>(key), graph);
		possible.push_back(new_vertex);
		return new_vertex;
	}
	std::vector<vertex_type> possible(1);
	vertex_type new_vertex = add_vertex(FVertex<point_type>(key), graph);
	possible.push_back(new_vertex);
	map[check] = possible;
	return new_vertex;
	//~ return it->second; 
 }

//add new edge to teh graph---------------------------------------------
 void GRAPH::AddNewEdge(Graph& G, vertex_type S, vertex_type T, line_type FaultSeg)
 {
	 //cout << "Assing edge between " << S << " and " << T << " with number of vertices " << num_vertices(G) << endl;
	 if (!boost::edge(S, T, G).second &&
		!geometry::equals(G[S].location, G[T].location))
		add_edge(S, T, 
			{geometry::length(FaultSeg)/1000, FaultSeg},
				G);
 }
 
 //draw graph as png (only fro samll graphs)----------------------------
 void GRAPH::DrawGraph(Graph G)
 {
 int systemRet;
 std::ofstream dot_file("PrimaryFaults.dot"); 
 
 write_graphviz(dot_file, G, default_writer(), 
	make_label_writer(get(&FEdge::length, G)));

 systemRet = system("dot -Tpng PrimaryFaults.dot -o PrimaryFaults.png"); 
	if(systemRet == -1)
		cout << "Failed to write png!" << endl;
 }
 
 //read vector data into graph------------------------------------------
 void GRAPH::ReadVEC(Graph& graph, map_vertex_type& map, std::vector<line_type> &faults)
 {
 vector <point_type> Intersec;
 vector <std::tuple< std::pair<point_type, point_type>, line_type, int >> G;
 std::pair <point_type, point_type> NODES;
 line_type TRACE;
 long double Xmax, Xmin, Ymax, Ymin;
 polygon_type pl;
 vertex_type VA, VB;
 box bx;
 int type;
 GEO g;
 GEOMETRIE geom;
 
 const double aoi_buffer = 10;
 
 Xmax = GeoTransform[0] + GeoTransform[1] * GeoTransform[6] - aoi_buffer;    // west
 Xmin = GeoTransform[0] + aoi_buffer; 									    // east
 Ymax = GeoTransform[3] - aoi_buffer; 									   // north
 Ymin = GeoTransform[3] + GeoTransform[5] * GeoTransform[7] + aoi_buffer; //south
 
 bx.min_corner().set<0>( Xmin );
 bx.min_corner().set<1>( Ymin );
 bx.max_corner().set<0>( Xmax );
 bx.max_corner().set<1>( Ymax );

//check fore edge nodes (crop faults by raster area)--------------------

	geometry::convert(bx, pl);

	BOOST_FOREACH(line_type const& Fault, faults)
	{
		if (geometry::within(Fault.front(), pl) ||
			geometry::within(Fault.back(), pl))
		{
			geometry::intersection(Fault, pl, Intersec);
			
			if (Intersec.size() > 0)
			{
				if (geometry::within(Fault.front(), pl))
					G.push_back(make_tuple(make_pair(Fault.front(),Intersec.at(0)), 
						geom.GetSegment(Fault, Intersec.at(0), Fault.front()), 1));
				
					if (geometry::within(Fault.back(), pl))
						G.push_back(make_tuple(make_pair(Intersec.at(0),Fault.back()), 
						geom.GetSegment(Fault, Intersec.at(0), Fault.back()), 2));
			}
				
				if (Intersec.size() == 0)
					G.push_back(make_tuple(make_pair(Fault.front(),Fault.back()), Fault, 0));
			Intersec.clear();
		}
	}

//create edges and nods for faults in the graph-------------------------
	for (typename vector <std::tuple< std::pair<point_type, point_type>, line_type, int >>::const_iterator it = G.begin(); it != G.end(); ++it)
	{
	NODES = get<0> (*it);
	TRACE = get<1> (*it);
	type  = get<2> (*it);
	
	VA = AddNewVertex(map, NODES.first, graph);
	VB = AddNewVertex(map, NODES.second, graph);

	//if( (degree(VA, graph) == 0) && (degree(VB, graph) == 0) )
		AddNewEdge(graph, VA, VB, TRACE);
			
		if (type == 1)
			graph[VB].Enode = true;
				
		if (type == 2)
			graph[VA].Enode = true;
	}

	cout << " Converted " << faults.size() << " faults into " << num_edges(graph) << " edges" << endl;
	
	faults.clear();
		for (auto Eg : make_iterator_range(edges(graph)))
			faults.push_back(graph[Eg].trace);
 }
  
  
  void GRAPH::CheckNetwork(Graph& G, map_type& map, double minDist)
 {
 /*GEOMETRIE geom;
 float angle, angle1, angle2;
 vertex_type U, u;
 line_type fault;
 vector<point_type> rmVert;
 vector<edge_type> rmEdge;

//check length of edges-------------------------------------------------
	for (auto eg : boost::make_iterator_range(edges(G))) 
	{
		if (geometry::length(G[eg].trace) < minDist)
		{
		U = source(eg, G);
		u = target(eg, G);

			if (degree(U, G) == 1 || 
				degree(u, G) == 1)
				rmEdge.push_back(eg);
		}
	}
//remove short edges----------------------------------------------------
	sort(rmEdge.begin(), rmEdge.end()); 
	rmEdge.erase(std::unique(rmEdge.begin(), rmEdge.end()), rmEdge.end());  

	BOOST_FOREACH(edge_type E, rmEdge)
	{
	U = source(E, G);
	u = target(E, G);

		if (edge(U, u, G).second)
		{
		cout << "WARNING: Removed edge to vertex! " << endl;
				remove_edge(U, u, G);
		}
	}
//check degrees of vertices---------------------------------------------
	for (auto vd : boost::make_iterator_range(vertices(G))) 
	{
		if (degree(vd, G) == 0)
			rmVert.push_back(G[vd].location);
		
		if (degree(vd, G) == 2)
		{
		auto be = adjacent_vertices(vd, G);
		   for (auto beit = be.first; beit != be.second; ++beit)
			{
				
			angle1 = (float)(atan2(G[vd].location.x() - G[*beit].location.x(), G[vd].location.y() - G[*beit].location.y())) 
							* (180 / math::constants::pi<double>());
							
							
							
				if (angle1 < 0) 
					angle1 += 180;
				
				if (angle != 0)
				
			cout << angle << " " << angle1 << endl;
			angle = angle1;	
			}

		}
		
		if (degree(vd, G) > 4)
		{
		cout << "ERROR: Too may edges at vertex \n"
			<< G[vd].location.x() << " " << G[vd].location.y() << endl;
			exit(EXIT_FAILURE);
		}
	}
	
	for (auto vd : boost::make_iterator_range(vertices(G))) 
	{
		BOOST_FOREACH(point_type P, rmVert)
		{
			if (geometry::equals(P, G[vd].location))
			{
			cout << "WARNING: Removed vertex!" << endl;
				map.erase(G[vd].location);
				remove_vertex(vd, G);
			}
		}
	}*/
 }
 
 //create the graph by insertic vertices at intersections and 
 //split traces into edges----------------------------------------------
void GRAPH::CreateGraph(Graph& graph, map_vertex_type& map, double minDist )
{
 GEOMETRIE geom;
 long double Distance;
 line_type fault, fault2;
 edge_iter Eg, Eg_end, Eg2, Eg2_end;
 vertex_type U, u, NewV, NewV2;
 vector<point_type> Intersec;
 vector<point_type> rmVert;

 vector <std::pair<vertex_type, vertex_type >> removeEdge;
 std::pair<vertex_type, vertex_type > removeS, removeT;
 vector <std::tuple<long double, point_type, edge_iter>> cross;

 typedef vector <std::tuple < edge_iter, vector< std::tuple<long double, point_type, edge_iter>> >> UpGraph;
 UpGraph GraphBuild;
 
 const double threshold = 1.0;
 
 for (tie(Eg, Eg_end) = edges(graph); Eg != Eg_end; ++Eg)
 {
	fault  = graph[*Eg].trace;
	for (tie(Eg2, Eg2_end) = edges(graph); Eg2 != Eg2_end; ++Eg2)
	{
	fault2 = graph[*Eg2].trace;	
	Intersec.clear();
		if(!geometry::equals(fault, fault2)) 
		{
//check for X-node------------------------------------------------------
			if(geometry::intersects(fault, fault2)) 
			{
				geometry::intersection(fault, fault2, Intersec);
				Distance =  geometry::length(geom.GetSegment(fault, fault.front(), Intersec.at(0)));
				cross.push_back(make_tuple(Distance, Intersec.at(0), Eg2));
			}
//check for Y-node------------------------------------------------------
			else if (!geometry::intersects(fault, fault2))
			{ 
				if (geometry::distance( fault, fault2.front()) < minDist )
				{
					Distance =  geometry::length(geom.GetSegment(fault, fault.front(), fault2.front()));
					cross.push_back(make_tuple(Distance, fault2.front(), Eg2));
				}
				if (geometry::distance(fault, fault2.back()) < minDist )
				{
					Distance =  geometry::length(geom.GetSegment(fault, fault.front(), fault2.back()));
					cross.push_back(make_tuple(Distance, fault2.back(), Eg2));
				}
			}
		}
	}
	if (cross.size() > 0)
		GraphBuild.push_back(make_tuple(Eg, cross));
	cross.clear();
 }

//Now update the Graph with intersections-------------------------------	
 for (typename UpGraph::const_iterator i = GraphBuild.begin(); i != GraphBuild.end(); ++i) 
 {
	Eg = get<0>(*i); 
	fault = graph[*Eg].trace;
	U = source(*Eg, graph);
	u = target(*Eg, graph);

	cross = get<1>(*i);

	if (cross.size() > 1)
	{
		if (boost::edge(U, u, graph).second) 
			removeEdge.push_back(make_pair(U, u));
		geom.SortDist(cross);
		for (typename vector <std::tuple<long double, point_type, edge_iter>>::const_iterator I = cross.begin(); I != cross.end(); I++)
		{
			NewV = AddNewVertex(map, get<1>(*I), graph);
			auto nx = std::next(I, 1);
			if (abs(get<0>(*I) - get<0>(cross.at(0))) <= threshold)
				AddNewEdge(graph, U, NewV, geom.GetSegment( fault, graph[U].location, graph[NewV].location));

			if (abs(get<0>(*I) - get<0>(cross.back())) <= threshold)
				AddNewEdge(graph, u, NewV, geom.GetSegment( fault, graph[u].location, graph[NewV].location));

			else
			{
				NewV2 = AddNewVertex(map, get<1>(*nx), graph);
				AddNewEdge(graph, NewV, NewV2, geom.GetSegment( fault, graph[NewV].location, graph[NewV2].location));
			}
		}
	}
	else
	{
		NewV = AddNewVertex(map, get<1>(cross.at(0)), graph);
		AddNewEdge(graph, U, NewV, geom.GetSegment( fault, graph[NewV].location, graph[U].location));
		AddNewEdge(graph, u, NewV, geom.GetSegment( fault, graph[NewV].location, graph[u].location));

		if (!geometry::equals(graph[U].location, graph[NewV].location) &&
			!geometry::equals(graph[u].location, graph[NewV].location) &&
			 boost::edge(U, u, graph).second) 
				removeEdge.push_back(make_pair(U, u));
	}
 }
/*
 //Now remove edges-----------------------------------------------------
 sort(removeEdge.begin(), removeEdge.end()); 
 removeEdge.erase(std::unique(removeEdge.begin(), removeEdge.end()), removeEdge.end());  
 for ( vector < pair<vertex_type, vertex_type> >::const_iterator it = removeEdge.begin();  
	it != removeEdge.end () ;  it++)
		remove_edge(get<0>(*it), get<1>(*it), graph);
		* */
}

void GRAPH::SplitFaults(Graph& graph, map_vertex_type& map, double minDist )
{
 Graph g;
 map_vertex_type m;
 GEOMETRIE geom;
 long double Distance;
 line_type fault, fault2;
 edge_iter Eg, Eg_end, Eg2, Eg2_end;
 vertex_type U, u, NewV, NewV2;
 vector<point_type> Intersec;
 vector<point_type> rmVert;

 vector <std::pair<vertex_type, vertex_type >> removeEdge;
 std::pair<vertex_type, vertex_type > removeS, removeT;
 vector <std::tuple<long double, point_type, edge_iter>> cross;

 typedef vector <std::tuple < edge_iter, vector< std::tuple<long double, point_type, edge_iter>> >> UpGraph;
 UpGraph GraphBuild;
 
 const double threshold = 1.0;
 
 for (tie(Eg, Eg_end) = edges(graph); Eg != Eg_end; ++Eg)
 {
	fault  = graph[*Eg].trace;
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
			cross.push_back(make_tuple(Distance, Intersec.at(0), Eg2));
		}
//check for Y-node------------------------------------------------------
		else
		{ 
			if (geometry::distance( fault, fault2.front()) < minDist )
			{
				Distance =  geometry::length(geom.GetSegment(fault, fault.front(), fault2.front()));
				cross.push_back(make_tuple(Distance, fault2.front(), Eg2));
			}
			if (geometry::distance(fault, fault2.back()) < minDist )
			{
				Distance =  geometry::length(geom.GetSegment(fault, fault.front(), fault2.back()));
				cross.push_back(make_tuple(Distance, fault2.back(), Eg2));
			}
		} 
	}
	
	geom.SortDist(cross);
	U = AddNewVertex(m, fault.front(), g);
	u = AddNewVertex(m, fault.back(), g);
	vertex_type prev_vertex = U;
	
	for (typename vector <std::tuple<long double, point_type, edge_iter>>::const_iterator I = cross.begin(); I != cross.end(); I++)
	{
		point_type intersect = get<1>(*I);
		bool is_start = geometry::distance(fault.front(), intersect) <= minDist;
		bool is_end   = geometry::distance(fault.back(),  intersect) <= minDist;
		if (is_start || is_end) continue;
		NewV = AddNewVertex(m, intersect, g);
		AddNewEdge(g, prev_vertex, NewV, geom.GetSegment(fault, g[prev_vertex].location, g[NewV].location));
		prev_vertex = NewV;
	}
	AddNewEdge(g, prev_vertex, u, geom.GetSegment(fault, g[prev_vertex].location, g[u].location));
	cross.clear();
 }
 
 cout << endl;
 graph = g;
 map = m;
}

//Topolgy analysis of graph---------------------------------------------
void GRAPH::GraphAnalysis(Graph& G, vector<float>& METRIC)
 {
 assert (num_vertices(G) != 0 && num_edges(G) != 0);
 
//1D intesnity P10 = (NE/2pi r) * (pi /2)
//Dimesnionless intesity P21 * Lb

 vector<int> component(num_vertices(G));
 int Inodes = 0,  Xnodes = 0, Ynodes = 0, Enodes = 0, NbB = 0; 
 int Nc, NbN, Nl, numK;
 float Cl, Cb, NCfreq, B20,  P20, P21 , B22;
 double A, totalLength = 0, avLenL, avLenB;
 box envBox;
 point_type Ul;
 vertex_type U, u;
 polygon_type area;
 map_type COUNT;
 typename map_type::const_iterator it;
  
 /***********************************************************************
 * Graph analysis based on:												*
 * Sanderson, D. J., Peacock, D. C., Nixon, C. W., & Rotevatn, A. (2018)* 
 * Graph theory and the analysis of fracture networks.					*
 * Journal of Structural Geology.										*
 **********************************************************************/

//test whether graph is planer
	if (boyer_myrvold_planarity_test(G))
	{
//number of connected components of the graph
	numK = connected_components(G, &component[0]);

	for (size_t i = 0; i < boost::num_vertices (G); ++i)
	{
		for (int ii = 0; ii < numK; ii++)
		{
		if (component[i] == ii)
			G[i].component = ii ;
		}
	}
 
//classify edges
	for (auto Eg : make_iterator_range(edges(G)))
	{
		U = source(Eg, G);
		u = target(Eg, G);
		G[Eg].component = std::to_string(G[U].component) + "-" + std::to_string(G[U].component);
		
		
			it = COUNT.find(G[U].location);
			if (it == COUNT.end())
			{
			COUNT[G[U].location] = U;
				if (G[U].Enode == false)
				{
					if (degree(U, G) == 1 && G[U].Enode == false)
						Inodes++;
					
					if (degree(U, G) == 3 && G[U].Enode == false) 
						Ynodes++;
						
					if (degree(U, G) == 4 && G[U].Enode == false)
						Xnodes++;
				}
				else if (G[U].Enode == true)
					Enodes++;
			}
			
			it = COUNT.find(G[u].location);
			if (it == COUNT.end())
			{
			COUNT[G[u].location] = u;
				if (G[u].Enode == false)
				{
					if (degree(u, G) == 1 && G[u].Enode == false)
						Inodes++;
					if (degree(u, G) == 3 && G[u].Enode == false)
						Ynodes++;
					if (degree(u, G) == 4 && G[u].Enode == false)
						Xnodes++;
				}
				else if (G[u].Enode == true)
					Enodes++;
			}

		if (degree(U, G) == 1 && degree(u, G) == 1) 
			G[Eg].BranchType = "I-I";
	
		if ( degree(U, G) == 3 && degree(u, G) == 3) 
			G[Eg].BranchType = "Y-Y";

		if ( degree(U, G) == 4 && degree(u, G) == 4) 
			G[Eg].BranchType = "X-X";

		if ((degree(U, G) == 1 && degree(u, G) == 3) ||
			(degree(U, G) == 3 && degree(u, G) == 1))
				G[Eg].BranchType = "I-Y";

		if ((degree(U, G) == 1 && degree(u, G) == 4) ||
			(degree(U, G) == 4 && degree(u, G) == 1))
				G[Eg].BranchType = "I-X";

		if ((degree(U, G) == 3 && degree(u, G) == 4) ||
			(degree(U, G) == 4 && degree(u, G) == 3))
				G[Eg].BranchType = "Y-X";
				
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

//total length of fault traces and area to be analysed
	totalLength += G[Eg].length;
	geometry::append(area, G[Eg].trace);
	NbB++;
 }
 envBox = boost::geometry::return_envelope<box>(area); 
 Ul.set<0>(envBox.min_corner().get<0>());
 Ul.set<1>(envBox.max_corner().get<1>());

 A = geometry::area(envBox)/1000 ; 
  
  cout <<"Area: " << A << endl;
//Number of connections
Nc = Ynodes + Xnodes;
 
//Number of branches NB = (I + 3Y + 4X) / 2 
NbN = (Inodes + 3*Ynodes + 4*Xnodes) / 2;

//Number of lines NL = (I + 2Y) / 2 
Nl = (Inodes + 2*Ynodes) / 2;


//Average Length of lines L / NL 
if (Nl != 0)
{
	avLenL = totalLength / Nl;
	
	//Connections per line 2(Y+X) / NL 
	Cl = 2*(Ynodes + Xnodes) / Nl;
}

if (A != 0)
{
	//Connecting node frequency (NC km-2)
	NCfreq = Nc / A;

	//Branch frequency (B20)
	B20 = NbN / A;

	//Line frequency (P20)
	P20 = Nl / A;

	//2D intesity (P21 = L / A 
	P21 = totalLength / A;
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

//write results to txt-------------------------------------------------

	cout<< "\nGRAPH'S ANALYSIS OF NETWORK \n"
		<< "Nodes: " << num_vertices(G) << " EDGES: " << num_edges(G) << "\n"
		<< "Edgenodes: " << Enodes << "\n"
		<< "Number of connected Nodes: " << Nc << " (Xnodes + Ynodes)" << "\n"
		<< "Branches: " <<NbN << " Lines: "<< Nl << " Number of Branches: " << NbB << "\n"
		<< "Number of components (c): " << numK << endl;

	ofstream txtG ("Graph_results.txt");
	if (txtG.is_open())  
	{ 
		txtG <<"*********************************************************************\n"
			 <<"Graph analysis based on:											\n"
			 <<"Sanderson, D. J., Peacock, D. C., Nixon, C. W., & Rotevatn, A. (2018)\n"
			 <<"Graph theory and the analysis of fracture networks.				\n"
			 <<"Journal of Structural Geology.									\n"
			 <<"*********************************************************************\n" << endl;

		txtG << "Nodes: " << num_vertices(G) << " EDGES: " << num_edges(G) << "\n"
			 << "Edgenodes: " << Enodes << "\n"
			 << "Number of connected Nodes: " << Nc << " (Xnodes + Ynodes)" << "\n"
			 << "Branches: " <<NbN << " Lines: "<< Nl <<"\n"
			 << "Number of Branches: " << NbB << "\n"
			 << "Connections per line: " << Cl <<"\n" 
			 << "Connections per Branch: "<< Cb << "\n"
			 << "Average Branch length: " << avLenB << "\n"
			 << "Average Line length: " << avLenL << "\n"
			 << "Connecting node frequency: " << NCfreq << "\n"
			 << "Branch frequency: " << B20 << "\n"
			 << "Line frequency: " << P20 << "\n"
			 << "2D Intesnsity: " << P21 << "\n"
			 << "Dimesnonless intesity: " << B22 << "\n"
			 << "Average degree of network (2e/n): " << (float) 2 * num_edges(G)/ num_vertices(G) << "\n"
			 << "Average number of connections per line: " << (float) 2 * (Xnodes + Ynodes) / num_vertices(G) << "\n"
			 << "Number of components (c): " << numK << "\n"
			 << "number of faces (f): " << num_edges(G) + numK - num_vertices(G) +1 << "\n" << endl;

//Analyse the different component of the network------------------------
		for (int i =0; i < numK; i++)
		{
			Inodes = 0, Ynodes = 0, Xnodes = 0, Enodes = 0, NbB = 0, totalLength = 0;
			for (auto Eg : make_iterator_range(edges(G)))
			{
				U = source(Eg, G);
				u = target(Eg, G);
				
				if (G[U].component == i && G[u].component == i)
				{
					if (G[U].Enode == false || G[u].Enode == false)
					{ 
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
				NbN = (Inodes + 3*Ynodes + 4*Xnodes) / 2;
			}
			if (NbB > 1)
			{
			txtG << "COMPONENT NO. " << i << "\n"
				 << "Branches: " << NbB   << "\n" 
				 << "Branches (calc): "   << NbN << "\n" 
				 << "Inodes: " << Inodes  << "\n" 
				 << "Ynodes: " << Ynodes  << "\n" 
				 << "Xnodes: " << Xnodes  << "\n" 
				 << "Enodes: " << Enodes  << "\n"
				 <<" Connections (Y +X): "<< (Ynodes + Xnodes) << "\n"
				 <<" Cummulative length: "<< totalLength << "\n"
				 <<" Average length: "    << totalLength / NbB<< "\n" << endl;
			 }
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
	METRIC.push_back ((float) 2 *    num_edges(G)   / num_vertices(G)); 
	METRIC.push_back ((float)(2 * Xnodes + Ynodes) / num_vertices(G));
 }
 
 //populate vicinity of traces with potential fracture centres----------
void GRAPH::CreateFractures(Graph& G, map_vertex_type& map, vector<line_type> FAULTS, string const& filename)
 {
	std::clock_t startcputime = std::clock();
	cout << "Distribute fractures: "<< endl;
	GEOMETRIE geom;
	//int Nb;
	//float frac;
	GEO geo;
	//box Box;
	//line_type Fault;

	//BUFFER buffer;
	vector<BUFFER> buffer2;	 

	//long double len;
	long double Xmax, Xmin, Ymax, Ymin;
	polygon_type pl;
	box bx;

	Xmax = GeoTransform[0] + GeoTransform[1] * GeoTransform[6] - 1000;    // west
	Xmin = GeoTransform[0] + 1000; 									  // east
	Ymax = GeoTransform[3] - 1000; 									 // north
	Ymin = GeoTransform[3] + GeoTransform[5] * GeoTransform[7] + 1000;	//south
 
	bx.min_corner().set<0>( Xmin );
	bx.min_corner().set<1>( Ymin );
	bx.max_corner().set<0>( Xmax );
	bx.max_corner().set<1>( Ymax );
	geometry::convert(bx, pl);

	BOOST_FOREACH(line_type Fault, FAULTS)
	{
		buffer2.push_back( geom.DefineLineBuffer(Fault, 5)); 
	}
	int totalF;
	int F = 0;
	#pragma omp parallel for
	for(unsigned int Fau = 0; Fau < FAULTS.size(); Fau++)
	{
	line_type Fault = FAULTS.at(Fau);

		long double len = geometry::length(Fault);
		float frac = len/25;
		
		BUFFER buffer = geom.DefineLineBuffer(Fault, frac);
		box Box = geometry::return_envelope<box>(buffer);
		
		int Nb =  10 * sqrt(2*len) ;

		for (int i = 0; i < Nb; i++)
		{
			point_type ranP;
			int ITER = 0;
			bool found = false, rejec;
			vertex_type NewV;

			while (!found)
			{
				rejec = false;

				double X = (geometry::get<geometry::min_corner, 0>(Box)) + (((long double) rand()) / 
					(long double) RAND_MAX) * (geometry::get<geometry::max_corner, 0>(Box) - (geometry::get<geometry::min_corner, 0>(Box)));      

				double Y = (geometry::get<geometry::min_corner, 1>(Box)) + (((long double) rand()) / 
					(long double) RAND_MAX) * (geometry::get<geometry::max_corner, 1>(Box) - (geometry::get<geometry::min_corner, 1>(Box)));  

				geometry::set<0>(ranP, X);
				geometry::set<1>(ranP, Y); 

				if (geometry::within(ranP, pl))
				{
					BOOST_FOREACH(BUFFER B , buffer2)
					{
						if (geometry::within(ranP, B)){
							rejec = true;
							break;
						}
					}
					if (!rejec)
					{
						if(geometry::within(ranP, buffer)) 
						{
							#pragma omp critical
							{
								NewV = AddNewVertex(map, ranP, G); 
								//assign elevation for fracture
								double hmax = geo.getElevation(G[NewV].location, filename);
								double hmin = -100;

								double h = ((double) rand() / (RAND_MAX)) ;
								double hc = hmin + (h * hmax);
								G[NewV].elevation = hc;
								totalF++;
							}
							//release lock
							if (frac > 100)
							{
								frac -=   2 * exp((1 - (i*10 / floor((0.5*Nb)) )));
								buffer = geom.DefineLineBuffer(Fault, frac);
								Box = geometry::return_envelope<box>(buffer);
							}
							found = true;
						}
					}
				}
				if (ITER >  100)
					found = true;
				ITER ++;
			}
		}
	if(omp_get_thread_num() == 0)
		cout << "fault: "<< F << "/"<< FAULTS.size() << endl;
	F++;
	}
	cout << "Created "<< totalF << " fractures along " <<FAULTS.size() << " faults" << endl;
	cout << " in " << (clock() - startcputime) / (double)CLOCKS_PER_SEC << " seconds [CPU Clock] \n" << endl;
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

  for(vertex_type u = predecessors[T]; u != T; T = u, u = predecessors[T])
  {
    std::pair<edge_type, bool> edge_pair = boost::edge(u,T,G);
    path.push_back( edge_pair.first );
  }

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

if (distance > 0)
Gref.WriteSHP(shortP, "ShortestPath.shp");
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
  
  std::vector < edge_type > spanning_tree;
  kruskal_minimum_spanning_tree( G, std::back_inserter(spanning_tree), 
								 weight_map( get(&FEdge::length, G)) );

  for (std::vector < edge_type  >::iterator ei = spanning_tree.begin();
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
 Gref.WriteSHP(min_graph, "MinTree.shp");
 cout <<"minTree Total eges: " << i << endl;
}
