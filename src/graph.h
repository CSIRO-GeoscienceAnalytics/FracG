#ifndef _GRAPH_h
#define _GRAPH_h
#include "main.h"

using namespace boost;
using namespace FGraph;

class GRAPH
{
	public:

	friend class GEOMETRIE;
	friend class GEO;

	GRAPH();
	~GRAPH()
	{}
	;  
	void DrawGraph(Graph G);
	void CheckNetwork(Graph& G, map_type& map, double minDist);

	vertex_type AddNewVertex(map_vertex_type& map, point_type const& key, Graph& graph);
	void AddNewEdge(Graph& G, vertex_type S, vertex_type T, line_type FaultSeg);
	void ReadVEC(Graph& graph, map_vertex_type& map, std::vector<line_type> &faults);
	void CreateGraph(Graph& graph, map_vertex_type& map, double minDist );
	void SplitFaults(Graph& graph, map_vertex_type& map, double minDist );
	void GraphAnalysis(Graph& G, vector<float>& METRIC);
	void CreateFractures(Graph& G, map_vertex_type& map, vector<line_type> FAULTS, string const& filename);
	void ShortPath(Graph G, map_vertex_type m, point_type source, point_type target, double radius);
	void MinTree (Graph G);
 };
#endif
