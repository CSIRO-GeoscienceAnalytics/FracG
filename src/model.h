/*
 * https://gitlab.onelab.info/gmsh/gmsh/blob/master/demos/api/t1.cpp
 * 
 * https://gitlab.onelab.info/gmsh/gmsh/wikis/Gmsh-compilation
 * 
 * 
* http://git.dev.opencascade.org/gitweb/?p=occt.git;a=snapshot;h=refs/tags/V7_3_0;sf=tgz
* tar zxf occt.tgz
* cd occt-V7_3_0
* mkdir build
* cd build
* cmake -DCMAKE_BUILD_TYPE=Release -DBUILD_MODULE_Draw=0 -DBUILD_MODULE_Visualization=0 -DBUILD_MODULE_ApplicationFramework=0 ..
* make
* sudo make install
* 
* git clone http://gitlab.onelab.info/gmsh/gmsh.git
* cd gmsh
* mkdir build
* cd build
* cmake -DENABLE_BUILD_DYNAMIC=1 -DENABLE_CGNS=0 -DENABLE_MED=0 ..
* make
* http://gmsh.info/doc/texinfo/gmsh.html#Compiling-the-source-code
*/

#ifndef _MODEL_h
#define _MODEL_h
#include "main.h"
#include "geometrie.h"

#include <typeinfo>
#include <armadillo>
//#include <gmsh.h>

using namespace FGraph;
using namespace arma;
using namespace std;

class MODEL
{
	public:

		MODEL();
		~MODEL()
		{}
		

	void WriteGmsh_2D(vector<line_type> faults, Graph G, string filename);
};
#endif
 
 
 
