//https://map.sarig.sa.gov.au/ 

#include "graph.h"
#include "GeoRef.h" 
#include "geometrie.h"
#include "main.h"
#include "stats.h"
#include "model.h"

#include <iostream>
#include <iomanip>
#include <chrono>
#include <ctime>
#include <thread>

using namespace std; 
using namespace FGraph;

string ext;

void Check(int argc, char* f1)
{
	if (argc == 1+1 )
	{
		ext =  strrchr(f1,'.') ;

		if (ext == ".txt" ||  ext == ".shp")
			cout << "Fault  data from " << ext << "-file." << endl;
		else
		{
			cout << "ERROR: Wrong vector format provided  (" << ext << ")" << endl;
			exit (EXIT_FAILURE);
		}
	}
	else
	{
		cout << " ERROR: Not enough input argmuments! (Please provide a vector file)" << endl;
		exit (EXIT_FAILURE);
	}
}

int main(int argc, char *argv[]) 
{ 
	cout<< "******************************************** \n"
		<< "* Boost version: " << BOOST_LIB_VERSION  << "                    *\n"
		<< "* Armadillo version: " << arma::arma_version::as_string() << "*\n"
		<< "******************************************** "<< endl;

	GEO geo;
	GRAPH G; 
	STATS stats;
	GEOMETRIE geom;
	MODEL m;
	map_vertex_type map; //a map of the vertices in graph G, for quick retrieval of vertices by their location

	point_type source, target; //source and target for shortest path calculation
	Graph graph1, frac_graph;
	
	double Dist; //distance threshold for considering different points to be the same location
	vector<line_type> faults;  //a vector of faults, which holds the data as it is read in from the shapefile
	
	Check(argc, argv[1]); // check that we have a shapefile for the fault location data, and a raster file for the elevation data
	string vecFile = (string) argv[1];

	ofstream txtF; 
	txtF.open ("Fault_Statistics.csv"); 
	txtF << vecFile <<endl;
	txtF.close();
	
	ofstream txtF2; 
	txtF2.open ("Graph_Statistics.csv"); 
	txtF2 << vecFile <<endl;
	txtF2.close();
	
	cout << "\nEnter minimum distance in m: ";
	std::cin >> Dist;
	
	std::clock_t startcputime = std::clock();
//----------------------------------------------------------------------
	//read in the fault information from the vector file
	if (ext == ".shp")
		geo.read_shp(vecFile, faults);
	if (ext == ".txt")
		geo.read_wkt(vecFile, faults);  
//----------------------------------------------------------------------
	geom.CentreDistanceMap(vecFile, 1000);

	geo.CorrectNetwork(faults, Dist);//rejoin faults that are incorrectly split in the data file
	stats.CreateStats(txtF, faults); // statistical analysis of network
	G.ReadVEC(graph1, map, faults); //convert the faults into a graph
	G.SplitFaults(graph1, map, Dist); //split the faults in the graph into fault segments, according to the intersections of the faults
	G.RemoveSpurs(graph1, map, Dist); //remove any spurs from the graph network
	
	G.GraphAnalysis(graph1, faults, txtF2, 10); 
	geo.WriteSHP(graph1,  "Branches.shp");
	geo.WriteSHP2(graph1, "Vertices.shp");
	G.ComponentExtract(graph1, faults);
	//---------------------------------------------------------------------- 

	source.set<0>(14301336.315685500000);
	source.set<1>(-1800551.207095710000);
	target.set<0>(14260164.968410300000);
	target.set<1>(-1809965.628127270000);
	G.ShortPath(graph1, map, source, target, 500);
	G.MinTree(graph1);
	//G.CreateFractures(frac_graph, map, faults, rasFile ); 
	//geo.WriteTxt(frac_graph, "frac_graph.txt");
	//----------------------------------------------------------------------
	stats.AddData(faults, source, target); //this asks for additional raster files
	

	//m.WriteGmsh_2D(faults, graph1, vecFile);
	cout << "Finished in " << (clock() - startcputime) / (double)CLOCKS_PER_SEC << " seconds [CPU Clock] \n" << endl;
	return EXIT_SUCCESS; 
} 
