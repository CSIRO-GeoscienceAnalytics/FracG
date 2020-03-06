#ifndef _GRAPH_h
#define _GRAPH_h
#include "main.h"

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
	void AddNewEdge(Graph& G, vertex_type S, vertex_type T, line_type FaultSeg);
	void AddNewEdge(Graph& G, vertex_type S, vertex_type T, line_type FaultSeg, double FaultLength);
	void ReadVEC(Graph& graph, map_vertex_type& map, std::vector<line_type> &faults);
	void ReadVEC4raster(Graph& graph, RASTER raster, map_vertex_type& map, std::vector<line_type> &faults);
	Graph ReadVEC4MODEL(std::vector<line_type> faults, box bx);
	void CreateGraph(Graph& graph, map_vertex_type& map, double minDist );
	void SplitFaults(Graph& graph, map_vertex_type& map, double minDist );
	void GraphAnalysis(Graph& G, VECTOR lines, int nb);
	void ShortPath(Graph G, map_vertex_type m, point_type source, point_type target, double radius);
	void MinTree (Graph G);
	VECTOR ComponentExtract(Graph G, VECTOR lines, int nb);
};
#endif
