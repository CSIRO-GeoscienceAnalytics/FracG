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

//Template functions for Armadillo--------------------------------------
template<typename T>
mat PrincipalComponent(T faultStats)
{
 int I = 0;
 mat coeff, score;
 vec latent,tsquared;
 
	mat data0(2, faultStats.size(),fill::zeros );
	 for (typename FSTATS::const_iterator it = faultStats.begin(); it != faultStats.end(); it++)
	 {	 
		data0(0,I) = it->first; 
		data0(1,I) = it->second; 
		I++;
	 }
 princomp(coeff, score, latent, tsquared, data0);
 coeff.print("coeff: ");
 return coeff;
 }
 
void GradientDecent(std::vector<std::pair<double, double>> faultStats);
void Covariance();
void GaussianMix();

 };
 #endif
 
 
 
