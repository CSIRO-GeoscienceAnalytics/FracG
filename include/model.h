#ifndef _MODEL_h
#define _MODEL_h
#include <gmsh.h>
#include "../include/fracg.h"
#include "../include/geometrie.h"
#include "../include/graph.h"

namespace model = gmsh::model;
namespace factory = gmsh::model::occ;

namespace FracG
{
	
	typedef boost::geometry::index::rtree<FracG::p_index, geometry::index::rstar<16>> dist_tree;
		
	dist_tree BuildPointTree(Graph &G);
	
	void AddLineament(line_type line, int source, int target, int &p_tag, int &l_tag, float lc, dist_tree &dtree);
	void WriteGmsh2D(dist_tree &dtree, bool output, Graph G, int nb_cells, string out_filename);
	void SampleNetwork2D(dist_tree &dtree, bool output, VECTOR &lines, int nb_cells, int nb_samples, double map_distance_threshold, string out_filename);
	
}
#endif
