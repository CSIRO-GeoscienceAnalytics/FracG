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
	void WriteGmsh_2D(bool output, Graph G, int nb_cells, double gmsh_min_cl, double gmsh_min_dist, double gmsh_max_dist, bool gmsh_in_meters, bool name_ss, std::string out_filename);
	void WriteGmsh_3D(FracG::AngleDistribution, bool output, Graph G, int nb_cells, double gmsh_min_cl, double gmsh_min_dist, double gmsh_max_dist, double z,double gmsh_point_tol, bool gmsh_in_meters, bool name_ss, std::string out_filename);
	void SampleNetwork_2D(bool show_output, double gmsh_sw_size, Graph g, int nb_cells, int nb_samples, bool gmsh_in_meters, std::string filename);
}
#endif
