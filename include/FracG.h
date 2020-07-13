#include <iostream>
#include <iomanip>
#include <chrono>
#include <ctime>
#include <thread>

#include "../include/graph.h"
#include "../include/GeoRef.h" 
#include "../include/geometrie.h"
#include "../include/fracg.h"
#include "../include/stats.h"
#include "../include/fracg.h"
#include "../include/model.h"
#include "../include/util.h"

using namespace std; 

namespace FGraph
{
	GEO geo;
	GRAPH G;
	STATS stats;
	GEOMETRIE geom;
	MODEL m;
	
	//these are external variables that can be acessed by all functions
	double split_dist;			//save the user defined split distance in case we have to re-split (for building the graph)
	vector<std::tuple<double, double, double>> gauss_params;	//holds the parameters of the gaussian(s) determined by KDE
	const char *refWKT_shp;
}
