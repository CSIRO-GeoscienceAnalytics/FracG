#ifndef _geo_h
#define _geo_h

#include "main.h"

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



 template<typename T> T** RasterConvert(int rows, int cols, T **M);
 
 seg_tree Seg_Tree(vector<line_type>F);
 
 void Point_Tree(vector<p_index> pointsF, vector<p_index>& closest);
 void Point_Tree2(vector<pl_index> points, vector<pl_index>& closest, double max_len);
 
 int getElevation(point_type p, std::string const& filename);
 int getElevationFromArray(point_type p, const int *data);
 void AssignValues(Graph& G, std::string const& filename);
 void AssignValuesAll(Graph& G, std::string const& filename);
 template<typename T> int readRaster(std::string const filename, T *&data);
 void GetRasterProperties(std::string const& filename, double**& RASTER);
 void read_wkt(std::string const& filename, std::vector<line_type>& lineString);
 void read_shp(std::string const& filename, std::vector<line_type>& lineString);
 void WriteTxt(Graph& g, string const& filename);
 void WriteSHP(Graph G, const char* Name);
 void CorrectNetwork(vector<line_type>&F, double dist);

};
#endif
