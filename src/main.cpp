//https://map.sarig.sa.gov.au/ 
#include <iostream>
#include <iomanip>
#include <chrono>
#include <ctime>
#include <thread>

#include "graph.h"
#include "GeoRef.h" 
#include "geometrie.h"
#include "main.h"
#include "stats.h"
#include "model.h"
#include "main.h"

using namespace std; 
using namespace FGraph;

int main(int argc, char *argv[]) 
{ 
	std::clock_t startcputime = std::clock();
	GEO geo;
	GRAPH G; 
	STATS stats;
	GEOMETRIE geom;
 	MODEL m;
	
	Graph graph, r_graph;
	s_t source_target;
	map_vertex_type map; //a map of the vertices in graph G, for quick retrieval of vertices by their location
	point_type source, target; //source and target for shortest path calculation
	
	VECTOR lines = geo.ReadVector(argc, argv[1]);		 // read the first layer of the shape file
	geo.CorrectNetwork(lines.data, 5);					// rejoin faults that are incorrectly split in the data file

	geom.CentreDistanceMap(lines, 2500.0);
	geom.P_Maps(lines, 5000.0);

 	//stats.CreateStats(lines); // statistical analysis of network
 	//stats.GetLengthDist(lines); // test for three distributions of length 
 	//stats.DoBoxCount(lines); // Boxcounting algorithm (NEED to run for several staring points)
 	
 	//stats.KDE_estimation_strikes(lines, "blub");

	G.ReadVEC(graph, map, lines.data); //convert the faults into a graph
	G.SplitFaults(graph, map, 5); //split the faults in the graph into fault segments, according to the intersections of the faults
	G.RemoveSpurs(graph, map, 10); //remove any spurs from the graph network
	
	G.GraphAnalysis(graph, lines, 10); //graph, vector data, minimum number of branches per component to analyse
	geo.WriteGraph(graph, lines, "a", false);
	
	VECTOR exComp = G.ComponentExtract(graph, lines, 0); //graph, vector data, of component to extract
	geo.ReadPoints(lines.folder+"/points.shp", lines, source_target);

	r_graph = geo.RasterGraph(lines, 5, 10, lines.folder+"/dem.tif");
	
	geo.WriteGraph(r_graph, lines, "b", true);
	

	//m.WriteGmsh_2D(true, r_graph, 10, "bla.msh");
	cout << "Finished in " << (clock() - startcputime) / (double)CLOCKS_PER_SEC << " seconds [CPU Clock] \n" << endl;
	return EXIT_SUCCESS; 
} 
