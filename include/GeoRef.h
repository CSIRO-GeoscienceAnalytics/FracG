/****************************************************************/
/*				DO NOT MODIFY THIS HEADER						*/
/*					FRACG - FRACture Graph						*/
/*				Network analysis and meshing software			*/
/*																*/
/*						(c) 2021 CSIRO							*/
/*			GNU General Public Licence version 3 (GPLv3)		*/
/*																*/
/*						Prepared by CSIRO						*/
/*																*/
/*					See license for full restrictions 			*/
/****************************************************************/
#ifndef _geo_h
#define _geo_h
#include "../include/fracg.h"
#include "../include/stats.h"
#include <boost/filesystem.hpp>

namespace FracG
{
	const std::string raster_subdir="raster";

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
	
	void WriteGraph(Graph g, VECTOR lines, std::string subF, bool raster, bool skip_betweenness_centrality);
	
	VECTOR ReadVector(std::string in_file, std::string out_directory="");
	void ReadPoints(std::string const& filename, VECTOR lines, std::pair<point_type, point_type> &source_target);

	void PointTree(std::vector<p_index> pointsF, std::vector<p_index>& closest);
	void PointTree2(std::vector<pl_index> points, std::vector<pl_index>& closest, double max_len);
	void Point_Tree3(std::vector<p_index> points,  std::vector<p_index>& closest, int nb);

	void WriteShapefile(VECTOR &lineaments, AngleDistribution &angle_dist, std::string name);
	void CorrectNetwork(std::vector<line_type>&F, double dist, double angl_threshold, double dfd_thres);
	
	void WriteSHP_lines(std::vector<line_type>lineaments, std::string refWKT, std::string name);
	void WriteSHP_maxFlow(DGraph G, std::string refWKT, std::string name );
	void WriteSHP(Graph G, std::string name);
	
//Raster function-------------------------------------------------------
	polygon_type BoundingBox(double transform[8], double buffer);
	template<typename T> void AnalyseRaster(VECTOR lines, double dist, RASTER<T> raster);
	
	Graph BuildRasterGraph(VECTOR lines, double split, double spur, double map_distance_threshold, double raster_dist, AngleDistribution angle_dist, const double param_penalty, std::string name, bool skip_betweenness_centrality);
	template<typename T> T** RasterConvert(int rows, int cols, T **M);
	template<typename T> void AssignValuesGraph(Graph& G, RASTER<T> raster);
	template<typename T> T GetRasterValue(point_type p, double transform[8], T** values);
	
//these functions convert the values from rasters into doubles----------
	template<typename T> double LineExtractor(line_type L, RASTER<T> raster);
	template<typename T> double CentreGradient(line_type F, RASTER<T> raster);
	template<typename T> double CrossGradient(line_type F, RASTER<T> raster);
	template<typename T> double ParallelGradient(line_type F, RASTER<T> raster);

	template<typename T> void WriteRasterStruct(RASTER<T> raster);
	void WriteRASTER(std::vector<std::vector<double>> data, std::string SpatialRef, double adfGeoTransform[6], std::string suffix = "");
	void gdal_resample(std::string srcfname, std::string dstfname);
//Maximum Flow----------------------------------------------------------
	void GetSourceTarget(const char* Name, point_type &Source, point_type &Target);
	DGraph MakeDirectedGraph(Graph &g);
	void SetupMaximumFlow(DGraph &dg, std::string type, double gradient_angle = 0);
	double MaximumFlow(DGraph &dg, point_type source, point_type target);
	double MaximumFlow(DGraph &dg, dvertex_type s, dvertex_type t);
    double CalculateMaximumFlow(DGraph &dg, dvertex_type s, dvertex_type t);
};
#endif
