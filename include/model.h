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
	
	p_tree BuildPointTree(Graph &G);
	
	void AddLineament(line_type line, int source, int target, int &p_tag, int &l_tag, float lc, p_tree &dist_tree);
	void WriteGmsh2D(p_tree &dist_tree, bool output, Graph G, int nb_cells, std::string out_filename);
	void SampleNetwork2D(p_tree &dist_tree, bool output, VECTOR &lines, int nb_cells, int nb_samples, double map_distance_threshold, std::string out_filename);
	
}
#endif
