#ifndef _STATS_h
#define _STATS_h
#include "main.h"
#include <armadillo>

using namespace FGraph;
using namespace boost;
using namespace std;
using namespace arma;

class STATS
{
	public:
		friend class GEOMETRIE;

		STATS();
		~STATS()
		{}
	;  
	
	vec Bootstrapping(vec Data);
	double PointExtractor(point_type P, double radius, double Transform[8], double** raster);
	double LineExtractor(line_type L, double radius, double Transform[8], double** raster);

	void BoxCountQuadTree(const vector<line_type> &faults, vector<std::tuple<double, long long int>> &counts, const double start_scale, const int max_depth, point_type &offset);
	void DoBoxCount(vector<line_type> &faults);
	void BoxCount(vector<line_type> faults, double &D);
	void AdditionalData(vector<line_type> faults, Graph& G);
	void CreateStats(ofstream& txtF, vector<line_type> faults);
	void DEManalysis(Graph& G, double radius, string filename, double** RASTER);
	FSTATS KMCluster(bool input, int No, FSTATS faultStats);
	double MinVarBuf(line_type L,  double GeoTransform[8], double** raster);
	double ScanLineDensity(vector<line_type> faults);

	void LSFitter(FSTATS length_strike);

	void WriteData(ofstream& txtF);
};
#endif
