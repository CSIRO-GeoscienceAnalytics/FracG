//https://map.sarig.sa.gov.au/ 

#include "graph.h"
#include "GeoRef.h" 
#include "geometrie.h"
#include "main.h"
#include "stats.h"

#include <iostream>
#include <iomanip>
#include <chrono>
#include <ctime>
#include <thread>

using namespace std; 
using namespace FGraph;

string ext, ext2;

void Check(int argc, char* f1, char* f2)
{
	if (argc == 2+1 )
	{
		ext =  strrchr(f1,'.') ;
		ext2 = strrchr(f2,'.') ;
		if (ext == ".txt" ||  ext == ".shp")
			cout << "Fault  data from " << ext << "-file." << endl;
		else
		{
			cout << "ERROR: Wrong vector format provided  (" << ext << ")" << endl;
			exit (EXIT_FAILURE);
		}
		
		if (ext2 == ".tif") 
			cout << "Raster data from " << ext2 << "-file." << endl;
		else
		{
			cout << "ERROR: Wrong raster format (" << ext2 << ")" << endl;
			exit (EXIT_FAILURE);
		}
	}
	else
	{
		cout << " ERROR: Not enough input argmuments! (Please provide vector and raster file)" << endl;
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
	GRAPH G; //graph of fault network, which is to be created by this program
	STATS stats;
	GEOMETRIE geom;
	seg_tree Segment_Tree;
	map_vertex_type map; //a map of the vertices in graph G, for quick retrieval of vertices by their location
	FSTATS  faultStats;
	FSTATS2 faultStats2;
	point_type source, target; //source and target for shortest path calculation
	Graph graph, frac_graph;
	
	double Dist; //distance threshold for considering different points to be the same location
	double** RASTER; //holds the raster data of the file
	vector<line_type> faults;  //a vector of faults, which holds the data as it is read in from the shapefile
	vector<float> METRICS;   //1: degree of network,

	if (argc < 3)
	{
		cout << "ERROR: Only " << argc - 1 << " arguments found, need at least 2" << endl;
		exit(-1);
	}
	
	Check(argc, argv[1], argv[2]); // check that we have a shapefile for the fault location data, and a raster file for the elevation data
	string vecFile = (string) argv[1];
	string rasFile = (string) argv[2];
	
	ofstream txtF; 
	txtF.open ("Fault_Statistics.csv"); 
	txtF << vecFile <<endl;
	txtF.close();
	
	cout << "\nEnter minimum distance in m: ";
	std::cin >> Dist;
	
	std::clock_t startcputime = std::clock();
	//----------------------------------------------------------------------
	//read in the fault information
	if (ext == ".shp")
		geo.read_shp(vecFile, faults);
	if (ext == ".txt")
		geo.read_wkt(vecFile, faults);  
		
	//rejoin faults that are incorrectly split in the data file
	geo.CorrectNetwork(faults, Dist);

	geo.GetRasterProperties(rasFile, RASTER);
	
	stats.CreateStats(txtF, faults);
	G.ReadVEC(graph, map, faults); //convert the faults into a graph
	G.SplitFaults(graph, map, Dist); //split the faults in the graph into fault segments, according to the intersections of the faults
	
	G.RemoveSpurs(graph, map, Dist); //remove any spurs from the graph network
	
	//G.DrawGraph(graph); //draw an image showing the network of fault intersections
		
	geo.AssignValuesAll(graph, rasFile); //assign elevation values to the vertices in the graph
	G.GraphAnalysis(graph, METRICS); 

	stats.DEManalysis(graph, 60, rasFile, RASTER);
		
	geo.WriteSHP(graph, "GRAPH.shp"); //write out the graph data as a file
	geo.WriteTxt(graph, argv[1]); //write out the graph data as a text file

	//---------------------------------------------------------------------- 
	source.set<0>(14301336.315685500000);
	source.set<1>(-1800551.207095710000);
	
	target.set<0>(14260164.968410300000);
	target.set<1>(-1809965.628127270000);

	//G.ShortPath(graph, map, source, target, 500);
	G.MinTree (graph);
		
	//G.CreateFractures(frac_graph, map, faults, rasFile ); 
	//geo.WriteTxt(frac_graph, "frac_graph.txt");
	//----------------------------------------------------------------------
	cout << "Finished in " << (clock() - startcputime) / (double)CLOCKS_PER_SEC << " seconds [CPU Clock] \n" << endl;
	stats.AdditionalData(faults, graph);
	
	seg_tree T = geo.Seg_Tree(faults);
	
	return EXIT_SUCCESS; 
} 
