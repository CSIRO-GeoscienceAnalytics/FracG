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

	void GaussianMix(vec data1, int g);

};
#endif
 
 
 
