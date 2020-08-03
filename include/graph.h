#ifndef _GRAPH_h
#define _GRAPH_h
#include <type_traits>
#include "../include/fracg.h"



using namespace boost;
using namespace FGraph;



namespace bg  = boost::geometry;
namespace bgi = boost::geometry::index;

//an rtree-backed map which associates a value with a particular point location
template <typename VT = vertex_type, typename PT = point_type>
class point_map//
{
    typedef typename std::pair<PT, VT> element;
    
private:
    bgi::rtree<element, bgi::rstar<16>> rtree; //rtree object which is the core of this data structure
    double default_dist; //the distances used here are to allow for points in the data set and points asked for not being at the exact same location
    //points are considered the same if they are within *dist* distance of each other
    //this value is the default distance, to be used if the user does not choose a different distance threshold
    
    //return the nearest element to a given point
    std::optional<element> get_nearest_element(std::vector<element> &elements, PT &point, double dist = -1)
    {
        if (dist < 0) dist = default_dist;
        double closest_dist = std::numeric_limits<double>::infinity();
        decltype(elements.begin()) closest;
        if (elements.size() <= 0)  return std::nullopt;
        for (auto it = elements.begin(); it < elements.end(); it++)
        {
            double dist = bg::distance(it->first, point);
            if (dist < closest_dist)
            {
                closest_dist = dist;
                closest = it;
            }
        }
        if (closest_dist > dist) return std::nullopt; //nearest object is too far away
        return *closest;
    }
    
public:
    point_map(double dist)
    {
        default_dist = dist;
        rtree = decltype(rtree)();
    }
    
    
    //copy constructor
    point_map(const point_map<VT, PT> &ref):rtree(ref.rtree), default_dist(ref.default_dist) {};
    
    point_map<VT, PT>& operator=(const point_map<VT, PT>& other)
    {
        this->rtree = other.rtree;
        this->default_dist = other.default_dist;
        return *this;
    }
    
    //get the nearest element that is within the specified distance of the given point
    std::optional<element> get_nearest_element(PT &point, double dist=-1)
    {
        std::vector<element> results;
        int values_found = rtree.query(bgi::nearest(point, 1), std::back_inserter(results));
        return get_nearest_element(results, point, dist);
    }
    
    //get the nearest value to the given point, if it is within the specified distance
    std::optional<VT> get_value(PT &point, double dist = -1)
    {
//         if (values_found < 1)
//         {
//             //found no points
//             return std::nullopt;
//         }
        std::optional<element> nearest = get_nearest_element(point, dist);
        if (!nearest) return std::nullopt;
        return nearest->second;
    }
    
    //add a value, and return either the given value or the value that is already in place
    //the second value is true iff the value is new
    std::pair<VT, bool> add_value(PT &point, VT &value, double dist = -1)
    {
        std::optional<VT> result = get_value(point, dist);
        if (result) return std::pair(*result, false);// if a value exists, return it
        element el(point, value);
        rtree.insert(el);
        return std::pair(value, true);
    }
    
    //set this value to the location, replacing an existing value if it exists
    //return true iff the value is new (ie, there is no existing point within that distance)
    bool set_value(PT &point, VT &value, double dist = -1)
    {
        std::optional<element> nearest = get_nearest_element(point, dist);
        element new_element = std::pair(point, value);
        bool is_new = true;
        if (nearest)
        {
            //replace this value
            rtree.remove(nearest);
            is_new = false;
        }
        rtree.insert(new_element);
        return is_new;
    }
    
    //remove a particular value from the map, returning the associated value if it exists
    std::optional<VT> remove_value(PT &point, double dist = -1)
    {
        std::optional<element> candidate = get_nearest_element(point, dist);
        if (!candidate) return std::nullopt; //no candidate found, nothing to remove
        rtree.remove(*candidate);
        return candidate->second;
    }
    
    double get_dist() {return default_dist;}
};

//class that stores vertices in a graph by their point location
template <typename PT = point_type, typename VT = vertex_type, typename GT = Graph>
class graph_map{
private:
    point_map<VT, PT> pm; //map to associate points (locations) with vertices in the graph
    GT &graph; //graph object to which the vertices belong
    GT graph_holder; //if the user doesn't supply a graph reference, use this to hold the data
    
public:
    graph_map(GT &graph_obj, double dist) : pm(dist), graph(graph_obj) { }
    
    graph_map(double dist) : pm(dist), graph(graph_holder) { }
    
    graph_map(const graph_map<PT, VT, GT> &other) : pm(other.pm), graph(other.graph), graph_holder(other.graph_holder) { }
    
    graph_map<PT, VT, GT>& operator=(const graph_map<PT, VT, GT> &other)
    {
        this->pm = other.pm;
        this->graph = other.graph;
        this->graph_holder = other.graph_holder;
        return *this;
    }
    
    std::optional<VT> get_vertex(PT &point, double dist = -1)
    {
        pm.get_value(point, dist);
    }
    
    std::pair<VT, bool> add_vertex_isnew(PT &point, double dist = -1)
    {
        std::optional<VT> existing = pm.get_value(point, dist);
        if (existing) return std::pair(*existing, false);
        VT new_vertex = boost::add_vertex(FVertex<PT>(point), graph);
        std::pair<VT, bool> added_value = pm.add_value(point, new_vertex, dist);
//         if (!added_value->second) std::cerr << "Error: Added a vertex at " << point.x << ", " << point.y << " which was not present at the first check, but did exist at when it was added to the map: " << added_value.first << std::endl;
        return std::pair(added_value.first, true);
    }
    
    VT add_vertex(PT &point, double dist = -1)
    {
        VT vertex;
        bool is_new;
        std::tie(vertex, is_new) = add_vertex_isnew(point, dist);
        return vertex;
    }
    
    std::optional<VT> remove_vertex(PT &point, double dist = -1)
    {
        std::optional<VT> candidate = pm.remove_value(point, dist);
        if (!candidate) return std::nullopt;
        boost::remove_vertex(*candidate, graph);
        return *candidate;
    }
    
    GT &get_graph() { return graph; }
    
    double get_dist() {return pm.get_dist();}
    
};

class GRAPH
{
	public:

	friend class GEOMETRIE;
	friend class STATS;
	friend class GEO;

	GRAPH();
	~GRAPH()
	{}
	;  
	
	void RemoveSpurs(graph_map<> & map, double minDist);

	vertex_type AddNewVertex(map_vertex_type& map, point_type const& key, Graph& graph);
	vertex_type GetVertex(map_vertex_type& map, point_type const& key, Graph& graph);
	bool AddNewEdge(Graph& G, vertex_type S, vertex_type T, line_type FaultSeg);
	bool AddNewEdge(Graph& G, vertex_type S, vertex_type T, line_type FaultSeg, double FaultLength);
	graph_map<point_type, vertex_type, Graph> ConvertLinesToGraph(std::vector<line_type> &faults, double distance_threshold);
	
	graph_map<> ReadVEC4raster(double transform[8], std::vector<line_type> &faults, double distance_threshold);
	
	Graph ReadVEC4MODEL(std::vector<line_type> faults, box bx, double map_distance_threshold);
	void CreateGraph(Graph& graph, map_vertex_type& map, double minDist );
	graph_map<> SplitFaults(graph_map<>& map, double minDist );
	void GraphAnalysis(Graph& G, VECTOR lines, int nb, string name);

	VECTOR ComponentExtract(Graph G, VECTOR lines, int nb);
	void IntersectionMap(Graph G, VECTOR lines, float cell_size, float search_size);
	void ClassifyLineaments(Graph G, VECTOR lines, float dist, string name);
	
	Graph MinTree (Graph G, double map_dist_threshold, std::string out_filename="");
	Graph ShortPath(graph_map<> m, std::string in_filename, std::string out_filename="");
	
	void MaximumFlow_R(Graph G, string st_filename, string type, std::string out_filename);
	void MaximumFlow_VG(Graph G, string st_filename, float top, float bottom, string capacity_type, std::string out_filename);
	void MaximumFlow_HG(Graph G, string st_filename, float left, float rigth, string capacity_type, std::string out_filename);
};
#endif
