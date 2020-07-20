#ifndef _GRAPH_h
#define _GRAPH_h
#include "../include/fracg.h"



using namespace boost;
using namespace FGraph;

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
	
	void RemoveSpurs(Graph& G, map_vertex_type& map, double minDist);

	vertex_type AddNewVertex(map_vertex_type& map, point_type const& key, Graph& graph);
	vertex_type GetVertex(map_vertex_type& map, point_type const& key, Graph& graph);
	bool AddNewEdge(Graph& G, vertex_type S, vertex_type T, line_type FaultSeg);
	bool AddNewEdge(Graph& G, vertex_type S, vertex_type T, line_type FaultSeg, double FaultLength);
	void ReadVEC(Graph& graph, map_vertex_type& map, std::vector<line_type> &faults);
	
	void ReadVEC4raster(Graph& graph, double transform[8], map_vertex_type& map, std::vector<line_type> &faults);
	
	Graph ReadVEC4MODEL(std::vector<line_type> faults, box bx);
	void CreateGraph(Graph& graph, map_vertex_type& map, double minDist );
	void SplitFaults(Graph& graph, map_vertex_type& map, double minDist );
	void GraphAnalysis(Graph& G, VECTOR lines, int nb, string name);

	VECTOR ComponentExtract(Graph G, VECTOR lines, int nb);
	void IntersectionMap(Graph G, VECTOR lines, float cell_size, float search_size);
	void ClassifyLineaments(Graph G, VECTOR lines, float dist, string name);
	
	Graph MinTree (Graph G, std::string out_filename="");
	Graph ShortPath(Graph G, map_vertex_type m, std::string in_filename, std::string out_filename="");
	
	void MaximumFlow_R(Graph G, map_vertex_type map, string st_filename, string type, std::string out_filename);
	void MaximumFlow_VG(Graph G, map_vertex_type map, string st_filename, float top, float bottom, string capacity_type, std::string out_filename);
	void MaximumFlow_HG(Graph G, map_vertex_type map, string st_filename, float left, float rigth, string capacity_type, std::string out_filename);
};
#endif
