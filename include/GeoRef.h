#ifndef _geo_h
#define _geo_h
#include "../include/fracg.h"

#include <boost/filesystem.hpp>


using namespace FGraph;

extern OGRSpatialReference GeoRef;
extern double GeoTransform[8];
extern const char *GeoProj;

class GEO
{
	public:
	
		GEO();
		~GEO()
		{}
		;  
		
	struct different_id_p
	{
		different_id_p(size_t i) : id(i) {}
		bool operator()(std::pair<point_type, size_t> const& v) const { return v.second != id; }
		size_t id;
	};

	struct different_id_pl
	{
		different_id_pl(size_t i) : id(i) {}
		bool operator()(std::tuple<point_type, size_t, double_t> const& v) const { return std::get<1>(v) != id; }
		size_t id;
	};

	struct larger
	{
		larger(double_t l) : len(l) {}
		bool operator()(std::tuple<point_type, size_t, double_t> const& v) const 
		{ 
			return std::get<2>(v) > len; 
		}
		double_t len;
	};
	
	void WriteGraph(Graph g, VECTOR lines, string subF);
	void WriteGraph_R(Graph G, VECTOR lines, string subF);
	
	VECTOR ReadVector(std::string in_file, std::string out_directory="");
	void ReadPoints(std::string const& filename, VECTOR lines, std::pair<point_type, point_type> &source_target);

	void Point_Tree(vector<p_index> pointsF, vector<p_index>& closest);
	void Point_Tree2(vector<pl_index> points, vector<pl_index>& closest, double max_len);
	void Point_Tree3(vector<p_index> points,  vector<p_index>& closest, int nb);

	void WRITE_SHP(VECTOR &lineaments, gauss_params &angle_dist, string name);
	void CorrectNetwork(vector<line_type>&F, double dist);
	
	void WriteSHP_lines(vector<line_type>lineaments, const char* refWKT, string name);
	void WriteSHP_maxFlow(DGraph G, const char* refWKT, string name );
	void WriteSHP(Graph G, string name);
	
//Raster function-------------------------------------------------------
	polygon_type BoundingBox(double transform[8], double buffer);
	template<typename T> void AnalyseRaster(VECTOR lines, double dist, RASTER<T> raster);
	
	Graph BuildRasterGraph(VECTOR lines, double split, double spur, double map_distance_threshold, const double angle_param_penalty, string name);
	template<typename T> T** RasterConvert(int rows, int cols, T **M);
	template<typename T> void AssignValuesGraph(Graph& G, RASTER<T> raster);
	template<typename T> T getValue(point_type p, double transform[8], T** values);
	
//these functions convert the values from rasters into doubles----------
	template<typename T> double LineExtractor(line_type L, RASTER<T> raster);
	template<typename T> double CentreGradient(line_type F, RASTER<T> raster);
	template<typename T> double CrossGradient(line_type F, RASTER<T> raster);
	template<typename T> double ParallelGradient(line_type F, RASTER<T> raster);

	void WriteRASTER(vector<vector<double>> data, char* SpatialRef, double adfGeoTransform[6], VECTOR &input_file, string suffix = "");
	template<typename T> void WriteRASTER_struc(RASTER<T> raster);
	
//Maximum Flow----------------------------------------------------------
	void Get_Source_Target(const char* Name, point_type &Source, point_type &Target);
	DGraph MakeDirectedGraph(Graph &g);
	void setup_maximum_flow(DGraph &dg, string type);
	double maximum_flow(DGraph &dg, point_type source, point_type target);
	double maximum_flow(DGraph &dg, dvertex_type s, dvertex_type t);
};
#endif
