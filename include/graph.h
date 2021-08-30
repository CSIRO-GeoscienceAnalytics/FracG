#ifndef _GRAPH_h
#define _GRAPH_h
#include <type_traits>
#include "../include/fracg.h"
#include "../include/stats.h"

using namespace boost;
namespace FracG
{
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
        std::optional<element> GetNearestElement(PT &point, double dist=-1)
        {
            std::vector<element> results;
            int values_found = rtree.query(bgi::nearest(point, 1), std::back_inserter(results));
            return get_nearest_element(results, point, dist);
        }

        //get the nearest value to the given point, if it is within the specified distance
        std::optional<VT> GetValue(PT &point, double dist = -1)
        {
            std::optional<element> nearest = GetNearestElement(point, dist);
            if (!nearest) return std::nullopt;
            return nearest->second;
        }

        //add a value, and return either the given value or the value that is already in place
        //the second value is true iff the value is new
        std::pair<VT, bool> AddValue(PT &point, VT &value, double dist = -1)
        {
            std::optional<VT> result = GetValue(point, dist);
            if (result) return std::pair(*result, false);// if a value exists, return it
            element el(point, value);
            rtree.insert(el);
            return std::pair(value, true);
        }

        //set this value to the location, replacing an existing value if it exists
        //return true iff the value is new (ie, there is no existing point within that distance)
        bool SetValue(PT &point, VT &value, double dist = -1)
        {
            std::optional<element> nearest = GetNearestElement(point, dist);
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
        std::optional<VT> RemoveValue(PT &point, double dist = -1)
        {
            std::optional<element> candidate = GetNearestElement(point, dist);
            if (!candidate) return std::nullopt; //no candidate found, nothing to remove
            rtree.remove(*candidate);
            return candidate->second;
        }

        //return the (default) distance threshold used by this point map
        double GetDist() {return default_dist;}
    };

    //class that stores vertices in a graph by their point location
    template <typename PT = point_type, typename VT = vertex_type, typename GT = Graph>
    class graph_map{
    private:
        point_map<VT, PT> pm; //map to associate points (locations) with vertices in the graph
        GT &graph; //graph object to which the vertices belong
        GT graph_holder; //if the user doesn't supply a graph reference, use this to hold the data
        std::string reference_wkt; //well-known-text which describes the reference

    public:
        graph_map(GT &graph_obj, double dist, std::string ref_str="") : pm(dist), graph(graph_obj), reference_wkt(ref_str)
        {
            typename GT::vertex_iterator vi, vend;
            for (std::tie(vi, vend) = boost::vertices(graph); vi != vend; vi++)
            {
                typename graph_traits<GT>::vertex_descriptor v = *vi;
                pm.AddValue(graph[*vi].location, v, 0); //add vertices with a distance threshold, to ensure that the point map has all the vertices that are actually in the graph, even if they violate the distance threshold given here
            }
        }

        graph_map(double dist, std::string ref_str="") : pm(dist), graph(graph_holder), reference_wkt(ref_str) { }

        graph_map(const graph_map<PT, VT, GT> &other) : pm(other.pm), graph(other.graph), graph_holder(other.graph_holder), reference_wkt(other.reference_wkt) { }

        graph_map<PT, VT, GT>& operator=(const graph_map<PT, VT, GT> &other)
        {
            this->pm = other.pm;
            this->graph = other.graph;
            this->graph_holder = other.graph_holder;
            return *this;
        }

        //return an (optional) vertex for the given point, if a vertex exists within the distance threshold from that point
        std::optional<VT> GetVertex(PT &point, double dist = -1)
        {
            pm.GetValue(point, dist);
        }

        //Add a vertex to the graph at a particular point, or return an existing vertex if one already exists within the distance threshold of the given point
        VT AddVertex(PT &point, double dist = -1)
        {
            VT vertex;
            bool is_new;
            std::tie(vertex, is_new) = AddVertexIsNew(point, dist);
            return vertex;
        }

        //add a new vertex, or return an existing vertex if one is within the distance threshold of the given point, and also return a boolean stating whether or not the returned vertex is newly created
        std::pair<VT, bool> AddVertexIsNew(PT &point, double dist = -1)
        {
            std::optional<VT> existing = pm.GetValue(point, dist);
            if (existing) return std::pair(*existing, false);
            VT new_vertex = boost::add_vertex(FVertex<PT>(point), graph);
            std::pair<VT, bool> added_value = pm.AddValue(point, new_vertex, dist);
    //         if (!added_value->second) std::cerr << "Error: Added a vertex at " << point.x << ", " << point.y << " which was not present at the first check, but did exist at when it was added to the map: " << added_value.first << std::endl;
            return std::pair(added_value.first, true);
        }

        //remove a vertex from the graph and point map, if a vertex exists within the distance threshold if the given point
        std::optional<VT> RemoveVertex(PT &point, double dist = -1)
        {
            std::optional<VT> candidate = pm.RemoveValue(point, dist);
            if (!candidate) return std::nullopt;
            boost::remove_vertex(*candidate, graph);
            return *candidate;
        }

        //return the underlying graph object
        GT &GetGraph() { return graph; }

        //return the (default) distance threshold used by this point map
        double GetDist() {return pm.GetDist();}

        //get the Well Known Text which reference information for this graph
        const char *GetRefWKT() { return reference_wkt.c_str(); }
        
        //calculate and return the bounding box of the edges in the graph
        box GetBoundingBox()
        {
            polygon_type holder;
            typename GT::edge_iterator e, estart, eend;
            boost::tie(estart, eend) = boost::edges(graph);
            for (e = estart; e != eend; e++)
            {
                bg::append(holder, graph[*e].trace);
            }
            return bg::return_envelope<box>(holder);
        }
    };

    void RemoveSpurs(graph_map<> & map, double minDist);

    vertex_type AddNewVertex(map_vertex_type& map, point_type const& key, Graph& graph);
    vertex_type GetVertex(map_vertex_type& map, point_type const& key, Graph& graph);
    bool AddNewEdge(Graph& G, vertex_type S, vertex_type T, line_type FaultSeg);
    bool AddNewEdge(Graph& G, vertex_type S, vertex_type T, line_type FaultSeg, double FaultLength);
    graph_map<point_type, vertex_type, Graph> ConvertLinesToGraph(std::vector<line_type> &faults, std::string refWKT, double distance_threshold);

    graph_map<> ReadVEC4raster(double transform[8], VECTOR &lines, double distance_threshold);

    Graph Read4MODEL(Graph g, box bx, double map_distance_threshold);
    void CreateGraph(Graph& graph, map_vertex_type& map, double minDist );
    graph_map<> SplitFaults(graph_map<>& map, double minDist );
    void GraphAnalysis(Graph& G, VECTOR lines, int nb, AngleDistribution angle_dist, std::string name);

    void ComponentExtract(Graph G, VECTOR lines, int nb);
    void IntersectionMap(Graph G, VECTOR lines, float cell_size, float search_size, bool resample);
    void ClassifyLineaments(Graph G, VECTOR &lines, AngleDistribution &angle_dist, float dist, std::string name);

	void MakeCroppedGraph(Graph &g, box AOI, double crop_amount);
	
	void Betweeness_centrality(Graph G);
    Graph MinTree (graph_map<> gm, std::string out_filename="");
    Graph ShortPath(graph_map<> m, std::string in_filename, std::string out_filename="");

	void SetBoundaryPoints(Graph G, point_type& s, point_type& t, bool vert_grad);
    double MaximumFlowGradient(graph_map<> G, Direction flow_direction, Direction pressure_direction, double start_pressure, double end_pressure, double border_amount, std::string capacity_type, std::string refWKT, std::string out_filename);
	void MaxFlowTensor(graph_map<> G, std::string capacity_type, const char *refWKT, std::string out_filename);
}
#endif
