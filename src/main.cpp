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
 if (argc != 2 )
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
 GRAPH G;
 STATS stats;
 GEOMETRIE geom;
 seg_tree Segment_Tree;
 map_type map;
 FSTATS  faultStats;
 FSTATS2 faultStats2;
 point_type source, target;
 Graph graph, frac_graph;
 
 double Dist;
 double** RASTER;
 vector<line_type> faults;  
 vector<float> METRICS;   //1: degree of network,

 Check(argc, argv[1], argv[2]);
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
	if (ext == ".shp")
		geo.read_shp(vecFile, faults);
	if (ext == ".txt")
		geo.read_wkt(vecFile, faults);  
	geo.CorrectNetwork(faults, Dist);



 geo.GetRasterProperties(rasFile, RASTER);
 
	stats.CreateStats(txtF, faults);

	G.ReadVEC(graph, map,faults);
	//G.CreateGraph(graph, map, Dist);
//	G.CheckNetwork(graph, map, Dist);
 geo.AssignValuesAll(graph, rasFile);
	G.GraphAnalysis(graph, METRICS); 

 stats.DEManalysis(graph, 60, rasFile, RASTER);
	
 geo.WriteSHP(graph, "GRAPH.shp"); 
 geo.WriteTxt(graph, argv[1]);

//---------------------------------------------------------------------- 
   source.set<0>(14301336.315685500000);
   source.set<1>(-1800551.207095710000);
   
   target.set<0>(14260164.968410300000);
   target.set<1>(-1809965.628127270000);

	//G.ShortPath(graph, map, source, target, 500);
	//G.MinTree (graph);
	
	//G.CreateFractures(frac_graph, map, faults, rasFile ); 
 //geo.WriteTxt(frac_graph, "frac_graph.txt");
//----------------------------------------------------------------------
 cout << "Finished in " << (clock() - startcputime) / (double)CLOCKS_PER_SEC << " seconds [CPU Clock] \n" << endl;
 //stats.AdditionalData(faults, graph);
 
 seg_tree T = geo.Seg_Tree(faults);
  
 return EXIT_SUCCESS; 
} 
