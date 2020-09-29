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

	void WriteGmsh_2D(bool output, Graph G, int nb_cells, int gmsh_min_cl, double gmsh_min_dist, double gmsh_max_dist, std::string out_filename);
	void SampleNetwork_2D(bool output, VECTOR &lines, int nb_cells, int nb_samples, double map_distance_threshold, std::string out_filename);
	
}
#endif
