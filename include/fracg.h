/****************************************************************/
/*				DO NOT MODIFY THIS HEADER						*/
/*					FRACG - FRACture Graph						*/
/*				Network analysis and meshing software			*/
/*																*/
/*						(c) 2021 CSIRO							*/
/*			GNU General Public Licence version 3 (GPLv3)		*/
/*																*/
/*						Prepared by CSIRO						*/
/*																*/
/*					See license for full restrictions 			*/
/****************************************************************/
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

#include <sys/stat.h>
#include <sys/types.h>
#include <dirent.h>
#include <unistd.h>

#include <limits>
#include <thread>

#include <gsl/gsl_vector.h>
#include <gsl/gsl_statistics.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_multifit_nlinear.h>

#include <boost/shared_ptr.hpp>
#include <boost/config.hpp>
#include <boost/foreach.hpp>

#include <boost/random.hpp>
#include <boost/random/uniform_real.hpp>
#include <boost/random/random_device.hpp>

#if (BOOST_VERSION >= 107200)
	#include <boost/timer/progress_display.hpp>
#else
	#include <boost/progress.hpp>
#endif
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
#include <boost/graph/dijkstra_shortest_paths.hpp>
#include <boost/graph/betweenness_centrality.hpp>
#include <boost/graph/connected_components.hpp>
#include <boost/graph/boyer_myrvold_planar_test.hpp>
#include <boost/graph/kruskal_min_spanning_tree.hpp>
#include <boost/graph/prim_minimum_spanning_tree.hpp>
#include <boost/graph/boykov_kolmogorov_max_flow.hpp>

#include <boost/algorithm/clamp.hpp>

#include <boost/filesystem.hpp>

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

#include "ogr_geometry.h"
#include "gdal_priv.h"
#include "ogrsf_frmts.h"
#include "cpl_conv.h"
#include "cpl_port.h"
#include "gdalwarper.h"

#include <boost/variant.hpp>
#include <boost/variant/static_visitor.hpp>
#include <boost/variant/apply_visitor.hpp>
#include <boost/math/tools/roots.hpp>
#include <boost/math/special_functions/sign.hpp>
#include <armadillo>

#pragma once

namespace FracG
{
	//======================GEOMETRY======================================== 
	//typedef geometry::model::point<double, 2, geometry::cs::geographic <geometry::degree>> point_type;  
	typedef boost::geometry::model::d2::point_xy<long double>  point_type;
	typedef boost::geometry::model::d2::point_xy<long long> point_int;
	typedef boost::geometry::model::segment<point_type> Segment;
	typedef boost::geometry::model::linestring<point_type> line_type;
	typedef boost::geometry::model::multi_linestring<line_type> multiline_type;
	
	//declare types for buffer
	typedef boost::geometry::model::polygon<point_type> polygon_type;
	typedef boost::geometry::model::multi_polygon<polygon_type> BUFFER;
	typedef boost::geometry::model::box<point_type> box;

	struct VECTOR
	{
            std::string folder;
            boost::filesystem::path in_path;
            std::string out_folder;
            boost::filesystem::path out_path;
            std::string name;
            std::string refWKT;
            std::vector<line_type> data;
            using LINE_IT = decltype(data)::iterator;
	};
	
	
	template <typename T> 
	struct RASTER
	{
            std::string name;
            std::string refWKT;
            double transform[8];
            T** values;
	}; 
    
    //enum for directions, used for specifiying borders of the area of interest
	enum Direction: char {LEFT='l', RIGHT='r', TOP='t', BOTTOM='b', NONE='\0'};
	//=======================GRAPH==========================================
	// Fault Vetex Properties
	template <typename Point>
	struct FVertex
	{
        
        FVertex()
        {
                boost::geometry::assign_zero(location);
        }
        FVertex(Point const& loc)
                : location(loc)
        {
        }
        Point location; //this vertex's position in 2D space
        bool Enode = false; //whether or not this vertex intersects with the bounding box/area of interest (AOI)
        Direction enode_dir = Direction::NONE; //which edge of bounding box that this vertex intersects
        int component; //the ID of the Connected Component that this vertex is a part of
        double data;
	};

	// Fault Edges Properties
	struct FEdge
	{
		long double length; //length of the fault segment that makes up this edge
		line_type trace; //the fault segment that makes up this edge
		std::string BranchType;
		int FaultNb;
		int component;
		double fault_length; //the length of the entire fault that the fault segment is taken from
		double angle;
		int set;
		int index;
		
		double Centre;
		double MeanValue;
		double CentreGrad;
		double CrossGrad;
		double ParalGrad;
		friend std::ostream& operator<<(std::ostream &os, const FEdge &de) 
			{return os << "l " << de.length << ", full length = " << de.fault_length;}
	};
	
        struct FEdge2
	{
            long double length;
            line_type trace;
            int FaultNb;
	};
	
	typedef std::pair<point_type, point_type> s_t; //source and target for graph
	
	//now we can define the graph as an adjacency list
	//we also need a vertex descriptor to add vertices to the graph
	//and a map that stores the vertices that have already been added (control)
	typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS,
		FVertex<point_type>, FEdge> Graph;  
	
	//now we can define the graph as an adjacency list
	//we also need a vertex descriptor to add vertices to the graph
	//and a map that stores the vertices that have already been added (control)
	typedef boost::adjacency_list <boost::vecS, boost::vecS, boost::undirectedS,
		FVertex<point_type>, FEdge2> MGraph;

	typedef boost::graph_traits<Graph>::vertex_descriptor vertex_type;
	typedef boost::graph_traits<Graph>::vertex_iterator vertex_iter;
	typedef boost::graph_traits<Graph>::adjacency_iterator adj_iter;
	typedef boost::graph_traits<Graph>::edge_descriptor edge_type;
	typedef boost::graph_traits<Graph>::edge_iterator edge_iter;  

	typedef std::pair<vertex_iter, vertex_iter> VertexPair;
	typedef std::pair<edge_iter, edge_iter> EdgePair;

	typedef std::vector<edge_type> short_path;
	typedef std::map<point_type, vertex_type, boost::geometry::less<point_type> > map_type;
	typedef std::map<point_int, std::vector<vertex_type>, boost::geometry::less<point_int> > map_vertex_type;
	
	
	typedef boost::adjacency_list_traits<boost::vecS, boost::vecS, boost::directedS> DGraphTraits;
	typedef DGraphTraits::edge_descriptor dedge_type;
	typedef DGraphTraits::vertex_descriptor dvertex_type;
	//Structure for a directed vertex, for use in the maximum flow algorithm
	struct DVertex
	{
            point_type location; //the node's physical location
            bool Enode; //whether or not this node is an end node (off the edge of the raster file)
            Direction direction = Direction::NONE;
            double data = 0; //the data value at the vertex's location
            long long index = 0;
            dedge_type predecessor; //store the edge to this vertex's predecessor
            boost::default_color_type colour; //colour property used by the maximum flow algorithm
            double distance = 0; //distance from this vertex to the sink
            DVertex(){};
            DVertex(point_type loc, bool end, Direction dir, double data, long long ind) : location(loc), Enode(end), direction(dir), data(data), index(ind) {};
            DVertex(FVertex<point_type> v) : location(v.location), Enode(v.Enode), direction(v.enode_dir), data(v.data) {};
            DVertex(FVertex<point_type> v, long long ind) : location(v.location), Enode(v.Enode), direction(v.enode_dir), data(v.data), index(ind) {};
	};
	
	//Structure for a directed edge, for use in the maximum flow algorithm
	struct DEdge
	{
            double length; //the length of this fault/fracture segment
            double full_length; //the full length of the fault/fracture that this segment belongs to (used to determine the width of the damage zone)
            double angle; //the orientation of this edge
            line_type trace;
            double capacity; //the capacity for this edge
            double residual_capacity; //the unused capacity for this edge
            dedge_type reverse; //holds a link/pointer/whatever to the edge that links the same two vertices, but in the opposite direction
            DEdge() : length(0), full_length(0), capacity(0), residual_capacity(0) {};
            DEdge(double len, double full_len, line_type tra, double cap) : length(len), full_length(full_len), trace(tra), capacity(cap), residual_capacity(cap) {};
            DEdge(double len, double full_len, line_type tra)             : length(len), full_length(full_len), trace(tra), capacity(  0), residual_capacity(  0) {};
            
            DEdge(double cap) : capacity(cap), residual_capacity(cap) {};
            friend std::ostream& operator<<(std::ostream &os, const DEdge &de) {return os << "l " << de.length << ", full length = " << de.full_length;}
	};
	//<edge list for each vertex, vectex list, un/directed, vertex properties, edge properties, graph properties, (all?) edges list>
	typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::directedS, DVertex, DEdge> DGraph; //directed graph
	typedef boost::multi_array<double, 2> Raster_type;
	typedef std::pair<point_type, size_t> p_index;
	typedef std::tuple<point_type,size_t, double_t> pl_index;

	typedef boost::geometry::index::rtree<p_index, boost::geometry::index::rstar<16> > p_tree;
	typedef boost::geometry::index::rtree<pl_index, boost::geometry::index::rstar<16> > pl_tree;
}
