#ifndef _geo_h
#define _geo_h

#include <boost/graph/boykov_kolmogorov_max_flow.hpp>
#include "main.h"
#include <boost/algorithm/clamp.hpp>

using namespace FGraph;

extern double GeoTransform[8];
extern OGRSpatialReference GeoRef;
extern const char *GeoProj;

class GEO
{
	public:
	
		GEO();                           
		~GEO()
		{}
		;  
		
	struct different_id
	{
		different_id(size_t i) : id(i) {}
		bool operator()(std::pair<Segment, size_t> const& v) const { return v.second != id; }
		size_t id;
	};


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
	
	void ReprojectVECTOR2WGS84(std::string const& filename);

	template<typename T> T** RasterConvert(int rows, int cols, T **M);
	
	seg_tree Seg_Tree(vector<line_type>F);
	
	void Point_Tree(vector<p_index> pointsF, vector<p_index>& closest);
	void Point_Tree2(vector<pl_index> points, vector<pl_index>& closest, double max_len);
	void Point_Tree3(vector<p_index> points,  vector<p_index>& closest, int nb);
	
	double getElevation(point_type p, std::string const& filename);
	template<typename T> T getElevationFromArray(point_type p, const T *data);
	void AssignValues(Graph& G, std::string const& filename);
	void AssignValuesAll(Graph& G, std::string const& filename);
	
	template<typename T>void AssignValuesGraph(Graph& G, double transform[8], T** values);
	DGraph MakeDirectedGraph(Graph &g);
	void setup_maximum_flow(DGraph &dg);
	double maximum_flow(DGraph &dg, point_type source, point_type target);
	double maximum_flow(DGraph &dg, dvertex_type s, dvertex_type t);
	polygon_type BoundingBox(double transform[8], double raster_crop_size);
	template<typename T> T getValue(point_type p, double transform[8], T** values);
	double LineExtractor(line_type L, double Transform[8], double** raster);
	template<typename T>double CentreGradient(line_type F, double transform[8], T** values);
	double CrossGradient(line_type F, double transform[8], double** values);
	double ParallelGradient(line_type F, double transform[8], double** values);
	
	template<typename T> int readRaster(std::string const filename, T *&data);
	void GetRasterProperties(std::string const& filename, double**& RASTER);
	void read_wkt(std::string const& filename, std::vector<line_type>& lineString);
	void read_shp(std::string const& filename, std::vector<line_type>& lineString);
	void WriteTxt(Graph& g, string const& filename);
	void WriteSHP(Graph G, const char* Name);
	void WriteSHP2(Graph G, const char* Name);
	void CorrectNetwork(vector<line_type>&F, double dist);
	
	void WriteRASTER(vector<vector<double>> data, char* SpatialRef, double adfGeoTransform[6], const char* Name);
	void WriteSHP_lines(vector<line_type>lineaments, const char* Name);
	void WriteSHP_maxFlow(DGraph G,  const char* Name);

	void Source_Target(const char* Name, point_type &Source, point_type &Target);

};
#endif
