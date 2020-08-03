#ifndef _MODEL_h
#define _MODEL_h
#include <gmsh.h>
#include "../include/fracg.h"
#include "../include/geometrie.h"
#include "../include/graph.h"

namespace model = gmsh::model;
namespace factory = gmsh::model::occ;

class MODEL
{
	public:
		geometry::index::rtree<p_index, geometry::index::rstar<16>> DistTree;
		MODEL();
		~MODEL()
		{}
		
		
	void BuildPointTree(Graph G);
	
	void addLineament(line_type line, int source, int target, int &p_tag, int &l_tag, float lc);
	void WriteGmsh_2D(bool output, Graph G, int nb_cells, string out_filename);
	void SampleNetwork_2D(bool output, vector<line_type> faults, int nb_cells, int nb_samples, double map_distance_threshold, string out_filename);
	
};
#endif
