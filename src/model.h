#ifndef _MODEL_h
#define _MODEL_h
#include "main.h"
#include <typeinfo>
#include <armadillo>

using namespace FGraph;
using namespace arma;
using namespace std;

class MODEL
{
	public:


		MODEL();
		~MODEL()
		{}

	void WriteGeo(vector<line_type> faults, string filename);
};
#endif
 
 
 
