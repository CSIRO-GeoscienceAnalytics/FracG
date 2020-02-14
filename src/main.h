#define _DEFAULT_SOURCE 1

#include <omp.h>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstdint>
#include <chrono>
#include <limits>
#include <cmath>
#include <sys/time.h>
#include <ctime>
#include <algorithm>
#include <iterator>
#include <tuple>
#include <list>
#include <assert.h> 

#include <limits>
#include <thread>
 
#include <boost/shared_ptr.hpp>
#include <boost/config.hpp>
#include <boost/foreach.hpp>

#include <boost/random.hpp>
#include <boost/random/uniform_real.hpp>
#include <boost/random/random_device.hpp>

#include <boost/progress.hpp>
#include "boost/multi_array.hpp"

#include <boost/range/adaptors.hpp>
#include <boost/math/distributions.hpp>
#include <boost/math/constants/constants.hpp>
#include <boost/math/interpolators/barycentric_rational.hpp>

#include <boost/graph/copy.hpp>
#include <boost/graph/random.hpp>
#include <boost/graph/graphviz.hpp> 
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/graph_utility.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/adjacency_iterator.hpp>
#include <boost/property_map/property_map.hpp>
#include <boost/graph/connected_components.hpp>
#include <boost/graph/dijkstra_shortest_paths.hpp>
#include <boost/graph/boyer_myrvold_planar_test.hpp>
#include <boost/graph/kruskal_min_spanning_tree.hpp>
#include <boost/graph/prim_minimum_spanning_tree.hpp>

#include <boost/geometry/geometry.hpp>
#include <boost/geometry/algorithms/buffer.hpp> 
#include <boost/geometry/algorithms/equals.hpp>
#include <boost/geometry/algorithms/expand.hpp> 
#include <boost/geometry/algorithms/envelope.hpp> 
#include <boost/geometry/algorithms/intersection.hpp> 
#include <boost/geometry/geometries/geometries.hpp>
#include <boost/geometry/geometries/point_xy.hpp>
#include <boost/geometry/geometries/linestring.hpp>
#include <boost/geometry/multi/geometries/multi_linestring.hpp>
#include <boost/geometry/io/wkt/read.hpp>
#include <boost/geometry/geometries/point.hpp>
#include <boost/geometry/geometries/box.hpp>
#include <boost/geometry/index/rtree.hpp>

#include "/usr/include/gdal/ogr_geometry.h"
#include "/usr/include/gdal/gdal_priv.h"
#include "/usr/include/gdal/ogrsf_frmts.h"
#include "/usr/include/gdal/cpl_conv.h"
#include "/usr/include/gdal/cpl_port.h"
#include "/usr/include/gdal/gdalwarper.h"

#pragma once
using namespace boost;
using namespace std;

namespace FGraph
{
	//======================GEOMETRY======================================== 
	//typedef geometry::model::point<double, 2, geometry::cs::geographic <geometry::degree>> point_type;  
	typedef geometry::model::d2::point_xy<long double>  point_type;
	typedef boost::geometry::model::d2::point_xy<long long> point_int;
	typedef boost::geometry::model::segment<point_type> Segment;
	typedef geometry::model::linestring<point_type> line_type;
	typedef geometry::model::multi_linestring<line_type> multiline_type;
	
	//declare types for buffer
	typedef geometry::model::polygon<point_type> polygon_type;
	typedef geometry::model::multi_polygon<polygon_type> BUFFER;
	typedef geometry::model::box<point_type> box;
	
	struct READ
	{
		string name;
		double** values;
		double transform[8];
		READ(string n, double** v, double t[8])
		{
			name	= n;
			values	= v;
			memcpy(transform, t, 8*sizeof(double));
		}
	}; 
	
	//=======================GRAPH==========================================
	// Fault Vetex Properties
	template <typename Point>
	struct FVertex
	{
		FVertex()
		{
			geometry::assign_zero(location);
		}
		FVertex(Point const& loc)
			: location(loc)
		{
		}
		Point location;
		double elevation;
		double Pressure;
		bool Enode = false;
		int component;
		double data;
	};

	// Fault Edges Properties
	struct FEdge
	{
		long double length; //length of the fault segment that makes up this edge
		line_type trace; //the fault segment that makes up this edge
		std::string BranchType;
		int component;
		double Centre;
		double MeanValue;
		double CentreGrad;
		double CrossGrad;
		double ParalGrad;
		int FaultNb;
		double fault_length; //the length of the entire fault that the fault segment is taken from
		friend std::ostream& operator<<(std::ostream &os, const FEdge &de) {return os << "l " << de.length << ", full length = " << de.fault_length;}
	};
	
		struct FEdge2
	{
		long double length;
		line_type trace;
		int FaultNb;
	};
	
	//now we can define the graph as an adjacency list
	//we also need a vertex descriptor to add vertices to the graph
	//and a map that stores the vertices that have already been added (control)
	typedef adjacency_list <boost::vecS, boost::vecS, boost::undirectedS,
		FVertex<point_type>, FEdge> Graph;  
	
	//now we can define the graph as an adjacency list
	//we also need a vertex descriptor to add vertices to the graph
	//and a map that stores the vertices that have already been added (control)
	typedef adjacency_list <boost::vecS, boost::vecS, boost::undirectedS,
		FVertex<point_type>, FEdge2> MGraph;

	typedef graph_traits<Graph>::vertex_descriptor vertex_type;
	typedef graph_traits<Graph>::vertex_iterator vertex_iter;
	typedef graph_traits<Graph>::adjacency_iterator adj_iter;
	typedef graph_traits<Graph>::edge_descriptor edge_type;
	typedef graph_traits<Graph>::edge_iterator edge_iter;  

	typedef pair<vertex_iter, vertex_iter> VertexPair;
	typedef pair<edge_iter, edge_iter> EdgePair;

	typedef std::vector<edge_type> short_path;
	typedef map<point_type, vertex_type, geometry::less<point_type> > map_type;
	typedef map<point_int, std::vector<vertex_type>, geometry::less<point_int> > map_vertex_type;
	
	
	typedef adjacency_list_traits<boost::vecS, boost::vecS, boost::directedS> DGraphTraits;
	typedef DGraphTraits::edge_descriptor dedge_type;
	typedef DGraphTraits::vertex_descriptor dvertex_type;
	//Structure for a directed vertex, for use in the maximum flow algorithm
	struct DVertex
	{
		point_type location; //the node's physical location
		bool Enode; //whether or not this node is an end node (off the edge of the raster file)
		double elevation; //the elevation of the vertex's location
		long long index;
		dedge_type predecessor; //store the edge to this vertex's predecessor
		boost::default_color_type colour; //colour property used by the maximum flow algorithm
		double distance; //distance from this vertex to the sink
		DVertex(){};
		DVertex(point_type loc, bool end, double elev, long long ind) : location(loc), Enode(end), elevation(elev), index(ind) {};
	};
	
	//Structure for a directed edge, for use in the maximum flow algorithm
	struct DEdge
	{
		double length; //the length of this fault/fracture segment
		double full_length; //the full length of the fault/fracture that this segment belongs to (used to determine the width of the damage zone)
		line_type trace;
		double capacity; //the capacity for this edge
		double residual_capacity; //the unused capacity for this edge
		dedge_type reverse; //holds a link/pointer/whatever to the edge that links the same two vertices, but in the opposite direction
		DEdge(){};
		DEdge(double len, double full_len, line_type tra, double cap) : length(len), full_length(full_len), trace(tra),capacity(cap), residual_capacity(cap) {};
		DEdge(double len, double full_len,line_type tra) : length(len), full_length(full_len), trace(tra), capacity(0), residual_capacity(0) {};
		friend std::ostream& operator<<(std::ostream &os, const DEdge &de) {return os << "l " << de.length << ", full length = " << de.full_length;}
	};
	//<edge list for each vertex, vectex list, un/directed, vertex properties, edge properties, graph properties, (all?) edges list>
	typedef adjacency_list<boost::vecS, boost::vecS, boost::directedS, DVertex, DEdge> DGraph;
	
	typedef std::vector<std::pair<double, double>> FSTATS;
	typedef std::vector<std::tuple<double, double, double>> FSTATS2; 

	typedef boost::multi_array<double, 2> Raster_type;
	
	typedef std::pair<Segment, size_t> seg_index;
	typedef std::pair<point_type, size_t> p_index;
	typedef std::tuple<point_type,size_t, double_t> pl_index;
	typedef geometry::index::rtree<p_index, geometry::index::rstar<16> > p_tree;
	typedef geometry::index::rtree<pl_index, geometry::index::rstar<16> > pl_tree;
	typedef geometry::index::rtree<seg_index, geometry::index::linear<8>> seg_tree;
}
