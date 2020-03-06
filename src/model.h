#ifndef _MODEL_h
#define _MODEL_h
#include <gmsh/Gmsh.h>
#include "main.h"
#include "geometrie.h"
#include "graph.h"

namespace model = gmsh::model;
namespace factory = gmsh::model::occ;

class MODEL
{
	public:

		MODEL();
		~MODEL()
		{}
		
	void WriteGmsh_2D(bool output, Graph G, int nb_cells, string filename);
	void SampleNetwork_2D(bool output, vector<line_type> faults, int nb_cells, int nb_samples, string filename);
};
#endif
