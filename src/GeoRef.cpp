#include "../include/GeoRef.h"
#include "../include/graph.h"
#include "../include/geometrie.h"
#include "../include/stats.h"
#include "../include/util.h"

namespace fs = boost::filesystem;

const std::string raster_subdir="raster";

GEO::GEO ()
{
}

void CheckReferenceSystem(const char* s1, const char* s2)
{
	OGRSpatialReference oSRS1, oSRS2;
	oSRS1.importFromWkt(&s1);
	oSRS2.importFromWkt(&s2);
	int err = oSRS1.IsSameGeogCS(&oSRS2);
	if (err != 1)
	{
		cout << "differnet reference systems!" << endl;
		exit (EXIT_FAILURE);
	}
}


// create and open filestream in folder "statistics"
string CreateDir(VECTOR &input_file, std::initializer_list<string> folders)
{
	string folder_name(input_file.folder);
	while (folder_name.back() == '/' ) folder_name.pop_back(); //gdal/ESRI shapefiles don't like double slashes at the start
	while (folder_name.back() == '\\') folder_name.pop_back();
	if (folder_name.size() <= 0) folder_name = ".";
	for (auto it = folders.begin(); it != folders.end(); it++) folder_name += "/"+*it;
	const char* folder_cstr = folder_name.c_str();
	boost::filesystem::create_directories(folder_name);
	return folder_name;
}

//read vector data from shp file---------------------------------------
void read_shp(std::string const& filename, VECTOR& data)
{
	point_type P;
	line_type Line;
	const char * name = filename.c_str();
	
	GDALAllRegister();
	GDALDataset *poDS = static_cast<GDALDataset*>
	(
		GDALOpenEx( name, GDAL_OF_VECTOR, NULL, NULL, NULL )
	);
	if( poDS == NULL )
	{
		printf( " Opening shapefile \"%s\" failed.\n", name );
		exit( 1 );
	}  
	
	OGRSpatialReference * pOrigSrs = poDS->GetLayer( 0 )-> GetSpatialRef();
	if ( pOrigSrs )
	{
		if (!pOrigSrs->IsProjected() )
		{
			cout << "ERROR: vector data without spatial reference /not a projected reference system" << endl;
			exit(EXIT_FAILURE);
		}
		if ( pOrigSrs->IsProjected() )
			 pOrigSrs->exportToWkt( &data.refWKT );
	}

	OGRLayer  *poLayer = poDS->GetLayer( 0 );
	
	poLayer->ResetReading();
	OGRFeature *poFeature;
	while( (poFeature = poLayer->GetNextFeature()) != NULL )
	{
		OGRGeometry *poGeometry = poFeature->GetGeometryRef();

//check for multilinestring and split into individula lines-------------
		if ( poGeometry != NULL && wkbFlatten(poGeometry->getGeometryType()) == wkbMultiLineString)	
		{
			cout <<"MultiLineString detected, splitting into individual linestrings" << endl;
			OGRMultiLineString *poMuLine = (OGRMultiLineString *) poGeometry;
			
			for (auto line : poMuLine)
			{
				OGRLineString *poLine = (OGRLineString *) line;
				for (int i = 0; i < poLine->getNumPoints(); i++)
				{
					P.set<0>(poLine->getX(i));
					P.set<1>(poLine->getY(i));
					geometry::append(Line, P);
				}
//correction for double points and digitized circles--------------------
				boost::geometry::unique(Line);
				if (!geometry::equals(Line.front(), Line.back()))
					data.data.push_back(Line);
				else
					cout <<"!!! DETECTED CIRCULAR LINE AT " << setprecision(10) << Line.front().x() << " " << Line.front().y() << endl;
				Line.clear();
			}
		}
		
//reading in the lines in case of a simple linestring-------------------
		else if( poGeometry != NULL
				&& wkbFlatten(poGeometry->getGeometryType()) == wkbLineString)
		{
			OGRLineString *poLine = (OGRLineString *) poGeometry;
			for (int i = 0; i < poLine->getNumPoints(); i++)
			{
				P.set<0>(poLine->getX(i));
				P.set<1>(poLine->getY(i));
				geometry::append(Line, P);
			}
//correction for double points and digitized circles--------------------
			boost::geometry::unique(Line);
			if (!geometry::equals(Line.front(), Line.back()))
				data.data.push_back(Line);
			else
				cout <<"!!! DETECTED CIRCULAR LINE AT " << setprecision(10) << Line.front().x() << " " << Line.front().y() << endl;
			Line.clear();
		}
		OGRFeature::DestroyFeature( poFeature );
	}
	GDALClose( poDS );
	cout << "read " << data.data.size() << " faults from shp" << endl;
}
 
//Convert c++ template type to GDAL Datatype
template <typename T> struct GetGDALDataTypeTraits
{
	static const GDALDataType datatype;
};
//note that gdal has no 64bit integer datatypes
template<> struct GetGDALDataTypeTraits<unsigned char> {static const GDALDataType datatype = GDT_Byte;};
template<> struct GetGDALDataTypeTraits<char> {static const GDALDataType datatype = GDT_Byte;};
template<> struct GetGDALDataTypeTraits<unsigned short> {static const GDALDataType datatype = GDT_UInt16;};
template<> struct GetGDALDataTypeTraits<short> {static const GDALDataType datatype = GDT_Int16;};
template<> struct GetGDALDataTypeTraits<unsigned int> {static const GDALDataType datatype = GDT_UInt32;};
template<> struct GetGDALDataTypeTraits<int> {static const GDALDataType datatype = GDT_Int32;};
template<> struct GetGDALDataTypeTraits<float> {static const GDALDataType datatype = GDT_Float32;};
template<> struct GetGDALDataTypeTraits<double> {static const GDALDataType datatype = GDT_Float64;};

//call this function with the given template type to get the corresponding gdal datatype
template <typename T> static GDALDataType GetGDALDataType()
{
	return GetGDALDataTypeTraits<T>::datatype;
}

//create array with size of raster
template<typename T> 
T** RasterConvert(int rows, int cols, T **M)
{
    M = new T*[rows];
    for (int i = 0; i < rows; i++){
        M[i] = new T[cols];
    }
    return M;
}

//get mean raster value around line
template <typename T>double GEO::LineExtractor(line_type L, RASTER<T> raster)
{
	//note that raster values are converted into double 
	GEOMETRIE geom;

	polygon_type pl = BoundingBox(raster.transform, 1);
	box AOI;
	BUFFER envelop;
	point_type point;
	int maxX, minX, maxY, minY;
	vector<double> D;
    double radius = 1.5 * (( abs(raster.transform[1]) + abs(raster.transform[5]))/2) ;

	envelop = geom.DefineLineBuffer(L, radius);
	geometry::envelope(envelop, AOI);

	maxX = (int)(geometry::get<geometry::max_corner, 0>(AOI) - raster.transform[0]) / raster.transform[1];
	minX = (int)(geometry::get<geometry::min_corner, 0>(AOI) - raster.transform[0]) / raster.transform[1]; 
	maxY = (int)abs((geometry::get<geometry::max_corner, 1>(AOI) - raster.transform[3]) / raster.transform[5]);
	minY = (int)abs((geometry::get<geometry::min_corner, 1>(AOI) - raster.transform[3]) / raster.transform[5]);
	
	for (int x = minX; x < maxX; x++)
	{
		point.set<0>(raster.transform[0] + x * raster.transform[1]);
		for (int y = maxY; y < minY; y++)
		{
			point.set<1>(raster.transform[3] + y * raster.transform[5]);

			if (!geometry::disjoint(point, envelop) && geometry::within(point, pl))
				D.push_back((double)getValue(point, raster.transform, raster.values));
		}
	}

	if (D.size() > 0)
		return(accumulate(D.begin(), D.end(), 0.0)/D.size());
	
	else
		return(0);
 }
 
 template <typename T> double GEO::CentreGradient(line_type F, RASTER<T> raster)
{
// gradinet across teh centre of lineament
	GEO georef;
	line_type cross;
	point_type point, p1, p2;
	polygon_type pl = georef.BoundingBox(raster.transform, 1);
	double len =  2 * (abs(raster.transform[1]) + abs(raster.transform[5]) / 2);
	geometry::centroid(F, point);
	
	double rx = F.back().x() - F.front().x();
	double ry = F.back().y() - F.front().y();
	
	double l = sqrt(rx*rx + ry*ry);
	
	p1.set<0>((point.x() + ( ry/l) * len ));
	p1.set<1>((point.y() + (-rx/l) * len ));
	
	p2.set<0>((point.x() + (-ry/l) * len ));
	p2.set<1>((point.y() + ( rx/l) * len ));
	
	geometry::append(cross, p1);
	geometry::append(cross, point);
	geometry::append(cross, p2);

	if (!geometry::disjoint(p1, pl) && !geometry::disjoint(p2, pl))
	{	
		double p_v1 = (double)getValue(p1, raster.transform, raster.values);
		double p_v2 = (double)getValue(p2, raster.transform, raster.values);
		return(abs(p_v1 - p_v2)/len);
	}
	else
		return(0);
}

template <typename T>
double GEO::CrossGradient(line_type F, RASTER<T> raster)
{
//mean gradient across all segmetns of lineament
//create a functor that boost can use in for_each_segment
	GEOMETRIE::Perpencicular <geometry::model::referring_segment<point_type>> functor;
	functor.d =  2 * (abs(raster.transform[1]) + abs(raster.transform[5]) / 2);	 //set the length of the crossing lines to twice the mean pixel resolution
	functor  = geometry::for_each_segment(F, functor );
	
	polygon_type pl = BoundingBox(raster.transform, 1);
	vector<double> D;
	
	BOOST_FOREACH(line_type cross, functor.all)
	{
		if (geometry::within(cross.front(), pl) && geometry::within(cross.back(),pl))
		{
			double v2 = (double) getValue(cross.front(), raster.transform, raster.values);
			double v1 = (double) getValue(cross.back(), raster.transform, raster.values);
			D.push_back(abs(v1-v2)/geometry::length(cross));
		}
	}
	
	if (D.size() != 0)
		return(accumulate(D.begin(), D.end(), 0.0)/D.size() );
	else
		return(0);
}

template <typename T>
double GEO::ParallelGradient(line_type F, RASTER<T> raster)
{
//slope between lineametn tips
	polygon_type pl = BoundingBox(raster.transform, 1);
	
	if (geometry::within(F.front(), pl) && geometry::within(F.back(), pl))
	{
		double v2 = (double) getValue(F.front(), raster.transform, raster.values);
		double v1 = (double) getValue(F.back(), raster.transform, raster.values);
		return( abs(v1-v2) / geometry::length(F));
	}
	else
		return(0);
}


//make a directed graph from an undirected graph
//assumes that AssignValuesGraph() has already been called on the graph, and an appropriate raster file, to initialise the height values (which are stored in FVertex.data, not FVertex.elevation)
DGraph GEO::MakeDirectedGraph(Graph &g)
{
	DGraph dg; //this is the directed graph object that will be made, then returned
	
	std::map<vertex_type, dvertex_type> vertex_map; //use this to keep track of which nodes in the new graph correspond to particular nodes in the existing graph
	
	//start by making the indices
	Graph::vertex_iterator v, vstart, vend;
	boost::tie(vstart, vend) = boost::vertices(g);
    int index = 0;
	for (v = vstart; v < vend; v++)
	{
		FVertex<point_type> fv = g[*v];
		DVertex dv(fv.location, fv.Enode, fv.data, v - vstart);
		dvertex_type dvd = boost::add_vertex(dv, dg);
		vertex_map[*v] = dvd;
        dg[dvd].index = index++; //initialise vertex index list
	}
	
	Graph::edge_iterator e, estart, eend;
	boost::tie(estart, eend) = boost::edges(g);
	for (e = estart; e != eend; e++) //can't do e < eend here, I suppose these iterators are unordered
	{
		vertex_type s = source(*e, g); //source and target from the unidirectinal graph
		vertex_type t = target(*e, g);
		dvertex_type ds = vertex_map[s]; //get the corresponding source and target for the new graph
		dvertex_type dt = vertex_map[t];
		FEdge fe = g[*e];
		
		DEdge de(fe.length, fe.fault_length, fe.trace); //currently one set of edge properties. will fill out the other values inside the maximum flow function.
		dedge_type frwd, back;
		bool fadded, badded;
		boost::tie(frwd, fadded) = boost::add_edge(ds, dt, de, dg);
		boost::tie(back, badded) = boost::add_edge(dt, ds, de, dg);
		if (!fadded || !badded) cout << "Warning: Edge from " << dg[ds].index << " to " << dg[dt].index << " already exists (forward " << !fadded << ", back " << !badded << ")" << endl;
		dg[frwd].reverse = back;
		dg[back].reverse = frwd;
	}
	return dg;
}

//setting the capacity
double length_capacity(DEdge de)
{
	/*we set matrix permeability(p_m) to 0.5 and fracture permeability (p_f) to 1
	* fracture flow capcaity is calucalted as:
	* F = (p_f * w_f) / (p_m * l/2)
	* for the width of the fracuture (w_f) we assume a square-root scaling
	*/
	double p_f = 10;
	double p_m = 5;
	double w_f = sqrt(de.full_length); 
	double flow = (p_f * w_f) / (p_m * (de.full_length / 2));
	return(flow);
}

double orientation_capacity(DEdge de, float sig1)
{
	//we assume a guassian function for capacity depnding on roentation
	double angle = (float)(atan2(de.trace.front().x() - de.trace.back().x(), de.trace.front().y() - de.trace.back().y())) 
							* (180 / math::constants::pi<double>());
	if (angle  < 0) 
		angle  += 180;

	double stddev = 25;
	double mean = sig1 + 90;
	double flow = 1/ (stddev * sqrt(2*math::constants::pi<double>())) * exp(-0.5 * (std::pow( ((angle - mean) / stddev),2)));
	flow *= 1000;
	cout << "angle " << angle << " " << flow << endl;
	return(flow);
}

double length_orientation_capacity(DEdge de, float sig1)
{
	/*we set matrix permeability(p_m) to 0.5 and fracture permeability (p_f) to 1
	* fracture flow capcaity is calucalted as:
	* F = (p_f * w_f) / (p_m * l/2)
	* for the width of the fracuture (w_f) we assume a square-root scaling
	* we assume a guassian function for capacity depnding on roentation
	*/

	double angle = (float)(atan2(de.trace.front().x() - de.trace.back().x(), de.trace.front().y() - de.trace.back().y())) 
							* (180 / math::constants::pi<double>());
	if (angle  < 0) 
		angle  += 180;

	double stddev = 25;
	double mean = sig1 + 90;
	double p_f = 10;
	double p_m = 5;
	double w_f = sqrt(de.full_length); 
	double cap1 = (p_f * w_f) / (p_m * (de.full_length / 2));
	double flow = cap1* (1/ (stddev * sqrt(2*math::constants::pi<double>())) * exp(-0.5 * (std::pow( ((angle - mean) / stddev),2))));
	flow *= 1000;
	return(flow);
}

//the capacity values for the maximum flow algorithm
void GEO::setup_maximum_flow(DGraph &dg, string cap_type)
{
	float sig1;
	enum {ERROR, l, lo, o};
	static std::map<string, int> c_type;
	if (c_type.empty() )
	{
		c_type["l"] = l;			//capacity depending on length
		c_type["lo"] = lo;			//capacity depending on length and orientation
		c_type["o"] = o;			//capacity depending on orientation
	}
	int type = c_type[cap_type];
	
	switch (type) 
	{
			case o:
			case lo:
			{
				cout << "Enter maximum stress orientation: " ;
				cin >> sig1;
				cout << endl;
			} break;
			default: 
			{
			} break;
	}

	DGraph::edge_iterator e, estart, eend;
	boost::tie(estart, eend) = boost::edges(dg);
	for (e = estart; e != eend; e++)
	{
		double flow;
		switch (type) 
		{
			case l:
			{
				flow = length_capacity(dg[*e]);
			} break;
			
			case lo:
			{
				flow = length_orientation_capacity(dg[*e], sig1);
			} break;
			
			case o:
			{
				flow =  orientation_capacity(dg[*e], sig1);
			} break;
			
			default: 
			{
				cout << "ERROR: Wrong capacity type defined" << endl;
				exit(EXIT_FAILURE);
			} break;
		}
		dg[*e].capacity = flow;
		dg[*e].residual_capacity = flow;
	}
}


//calculate the maximum flow from source s to sink/target t
double c_maximum_flow(DGraph &dg, dvertex_type s, dvertex_type t)
{
	auto capacity_map = boost::get(&DEdge::capacity, dg);
	auto res_cap_map  = boost::get(&DEdge::residual_capacity, dg);
	auto rev_map      = boost::get(&DEdge::reverse, dg);
	auto pred_map     = boost::get(&DVertex::predecessor, dg);
	auto colour_map   = boost::get(&DVertex::colour, dg);
	auto distance_map = boost::get(&DVertex::distance, dg);//not using a distance map currently. this should be distance to the target vertex, and improves the speed of the algorithm
	auto index_map    = boost::get(&DVertex::index, dg); //no idea what the type is, but it will do for now
                                                                                    
	const double flow = boost::boykov_kolmogorov_max_flow(dg, capacity_map, res_cap_map, rev_map, pred_map, colour_map, distance_map, index_map, s, t);
	return flow;//0
}

//perform the maximum flow from the source location to the target location
double GEO::maximum_flow(DGraph &dg, point_type source, point_type target)
{
// 	DGraph::edge_iterator e, eend, eprev;
// 	boost::tie(e, eend) = boost::edges(dg);
	DGraph::vertex_iterator v, vstart, vend;
	boost::tie(vstart, vend) = boost::vertices(dg);
	dvertex_type s = *vstart, t = *vstart;
	double sdist, tdist;
	sdist = geometry::distance(source, dg[*vstart].location);
	tdist = geometry::distance(target, dg[*vstart].location);
	for (v = vstart; v < vend; v++)// eprev = e;
	{
		const double ds = geometry::distance(source, dg[*v].location);
		if (ds < sdist)
		{
			sdist = ds;
			s = *v;
		}
		const double dt = geometry::distance(target, dg[*v].location);
		if (dt < tdist)
		{
			tdist = dt;
			t = *v;
		}
	}
	if (dg[t].data > dg[s].data){
		std::swap(s, t);
		cout << "Swapping source and target to find maximum flow from higher location to lower location" << endl;
	}
	
	for (v = vstart; v < vend; v++)
    {
        //I'm not sure if this should be distance to the target, or if I should calculate the shortest paths from the target, or perhaps I set it to zero
        dg[*v].distance = geometry::distance(dg[*v].location, dg[t].location); //set distance map
    }
	
	//sanity check, make sure we're not calculating the max flow from a vertex to itself
	if (s == t){
		cout << "Warning: The source and target points for the maximum flow (" << s << ") are the same" << endl;
		return std::numeric_limits<double>::infinity();
	}
	
	cout << "Calculating maximum flow from vertex " << s << " (" << dg[s].data << "m) to vertex " << t << " (" << dg[t].data << ")" << endl;
	
	//clear capacities for verticies that are above the source
	const double start_data = dg[s].data;
	vector<std::pair<dedge_type, double>> saved_capacities;
	DGraph::edge_iterator e, estart, eend;
	boost::tie(estart, eend) = boost::edges(dg);
	for (e = estart; e != eend; e++)
	{
		dvertex_type cs = boost::source(*e, dg); //check source
		dvertex_type ct = boost::target(*e, dg); //check target
		if (dg[cs].data > start_data || dg[ct].data > start_data)
		{
			saved_capacities.push_back(std::make_pair(*e, dg[*e].capacity));
			dg[*e].capacity = 0;
			dg[*e].residual_capacity = 0;
		}
	}
// 	s = boost::source(*e, dg); //still testing, so pick two arbitrary locations
// 	e++;
// 	eprev = e;
// 	t = boost::target(*eprev, dg);
	double max_flow = c_maximum_flow(dg, s, t);
	//now restore capacities
	//or don't restore capacities, and instead use them to mark the faults as unusuable
// 	for (auto p = saved_capacities.begin(); p != saved_capacities.end(); p++)
// 	{
// 		dedge_type e = p->first;
// 		double cap = p->second;
// 		dg[e].capacity = cap;
// 	}
	return max_flow;
}



polygon_type GEO::BoundingBox(double transform[8], double buffer)
{
	long double Xmax, Xmin, Ymax, Ymin;
	polygon_type pl;
	box bx;

	//create new box of raster
	Xmax = transform[0] + transform[1] * transform[6] - buffer;   // west
	Xmin = transform[0] + buffer; 								  // east
	Ymax = transform[3] + buffer; 								 // north
	Ymin = transform[3] + transform[5] * transform[7] - buffer ;	// south
	
	bx.min_corner().set<0>( Xmin );
	bx.min_corner().set<1>( Ymin );
	bx.max_corner().set<0>( Xmax );
	bx.max_corner().set<1>( Ymax );
	geometry::convert(bx, pl);
	
	return(pl);
}

VECTOR GEO::ReadVector(std::string in_filename, std::string out_directory)
{
    //this needs to be re-written
//     if (argc < 1)
//     {
//         cerr << "ERROR: need at least one argument (shapefile to read)" << endl;
//         exit (EXIT_FAILURE);
//     }
	VECTOR lineaments;
    fs::path in_file(in_filename);
// 	string file = f; //this crashes if there are no arrguments
	string ext, name;
	
    
    
// 	if (argc == 1+1 )
// 	{
    ext =  in_file.extension().string();
//     name = file.substr(0, file.find_last_of("."));
    
    if (ext == ".txt" ||  ext == ".shp") //we don't read from text files
        cout << "Fault  data from " << ext << "-file." << endl;
    else
    {
        cout << "ERROR: Wrong vector format provided  (" << ext << ")" << endl;
        exit (EXIT_FAILURE);
    }
// 	}
// 	else
// 	{
// 		cout << " ERROR: Not enough input argmuments! (Please provide a vector file)" << endl;
// 		exit (EXIT_FAILURE);
// 	}
	
	
// 	size_t first_folder = name.find_first_of("/\\");
// 	while (first_folder != std::string::npos && first_folder + 1 < name.size() && (name[first_folder+1] == '/' || name[first_folder+1] == '\\')) name.erase(first_folder+1, 1);
// 	
// 	
// 	size_t folder_index = name.find_last_of("/\\"); //this should work on both linux (/) and windows (\\) 
// 	string folder = (folder_index != std::string::npos) ? name.substr(0               , folder_index + 1 ) : "./";
// 	name          = (folder_index != std::string::npos) ? name.substr(folder_index + 1, std::string::npos) : name;
	
	lineaments.folder = in_file.parent_path().string();
	lineaments.name = in_file.stem().string();
    lineaments.out_folder = (out_directory == "") ? lineaments.folder : out_directory;
    lineaments.in_path = in_file.parent_path();
    lineaments.out_path = fs::path(lineaments.out_folder);
    
//     cout << " folder = " << lineaments.folder << ", name = " << lineaments.name << ", out_folder = " << lineaments.out_folder << ", in_path = " << lineaments.in_path << ", out_path = " << lineaments.out_path << std::endl;

	//read information from the vector file
	if (ext == ".shp")
		read_shp(in_filename, lineaments);
	else
		cerr << "ERROR: no shape file definend" << endl;
	
	//set variable in FracG namespace
	refWKT_shp = lineaments.refWKT; //TODO: store this in the local variables, not a global. or rather, read it from local varaibles
	return(lineaments);
}

//assign raster values (primarily elevation, from a digital elevation map (DEM)) to the elements of a graph
//same effect as AssignValues() above, but more efficient
//read the entire file at once, rather than one at a time. gives faster IO. This assumes that we can fit the entire file in memory at once.
template <typename T>
void GEO::AssignValuesGraph(Graph& G, RASTER<T> raster)
{
	STATS stats;
	line_type curEdge;
	point_type centre;
	cout << " Assigning raster data to graph." << endl;
	
	for (auto Ve : boost::make_iterator_range(vertices(G)))
	{
		G[Ve].data = getValue(G[Ve].location, raster.transform, raster.values);
		if (std::isnan(G[Ve].data)) cout <<"Error: vertex " << Ve << " read a value of nan" << endl;
	}

	for (auto Eg : boost::make_iterator_range(edges(G))) 
	{
		curEdge = G[Eg].trace;
		geometry::centroid(G[Eg].trace, centre);

		G[Eg].Centre     = getValue(centre, raster.transform, raster.values);
		G[Eg].MeanValue  = LineExtractor<T>(G[Eg].trace, raster) ;
		G[Eg].CentreGrad = CentreGradient(G[Eg].trace, raster);
		G[Eg].CrossGrad  = CrossGradient(G[Eg].trace, raster);
		G[Eg].ParalGrad  = ParallelGradient(G[Eg].trace, raster) ;
	}
	cout << " done" << endl;
}
//make these instatiations so that other cpp files can find them
template void GEO::AssignValuesGraph<int>(Graph& G, RASTER<int> raster);
template void GEO::AssignValuesGraph<float>(Graph& G, RASTER<float> raster);
template void GEO::AssignValuesGraph<double>(Graph& G, RASTER<double> raster);

template <typename T>
RASTER <T>ReadRaster(string filename)
{
	
	/*
	[0]  top left x 
	[1]  w-e pixel resolution 
	[2]  0 
	[3]  top left y 
	[4]  0 
	[5]  n-s pixel resolution (negative value) 
	[6]  x-size
	[7]  y-size
	*/
	cout << "Reading raster " << filename;
	T** raster_data;
	const char * name = filename.c_str();
	RASTER <T>raster;
	double adfGeoTransform[6];
	
	
	GDALAllRegister();
	GDALDataset  *poDataset;
	poDataset = (GDALDataset *) GDALOpen( name, GA_ReadOnly );
    if( poDataset == NULL )
    {
        cout << "\n ERROR: cannot open raster file " << endl;
		exit(EXIT_FAILURE);
    }
    else
    {
			std::size_t pos = filename.find("."); 
			raster.name = filename.substr(0, pos);
			
			if( poDataset->GetProjectionRef()  != NULL )
				raster.refWKT = poDataset->GetProjectionRef() ;
				
				//check consitency of refernce system with shp file
				CheckReferenceSystem(raster.refWKT, refWKT_shp);
				
			if( poDataset->GetGeoTransform( adfGeoTransform ) == CE_None )
				memcpy ( &raster.transform, &adfGeoTransform, sizeof(adfGeoTransform) );
				
					
			GDALRasterBand *band = poDataset -> GetRasterBand(1);  
			
			poDataset = (GDALDataset *) GDALOpen( name, GA_ReadOnly );
			
			if( poDataset == NULL )
				std::cout<<"no file"<< std::endl;

			int xSize = poDataset->GetRasterXSize();
			int ySize = poDataset->GetRasterYSize();
			raster.transform[6] = xSize;
			raster.transform[7] = ySize;
			
			cout << " of size: "<<  xSize << " x " << ySize << endl;
			
			raster.values = (T**)RasterConvert(poDataset->GetRasterXSize(), poDataset->GetRasterYSize(), raster.values);
			
			for(int row = 0; row < poDataset->GetRasterXSize(); row++) 
			{
				int err = band->RasterIO(GF_Read, row, 0, 1, ySize, &(raster.values[row][0]), 1, ySize, GetGDALDataType<T>(), 0, 0, nullptr);  //read coloun of raster (switch datatype depending on raster's data type)
				if (err != 0)
				{
					cerr << " ERROR reading from raster" << endl;
					exit (EXIT_FAILURE);
				}
			}

			GDALClose( poDataset );
	}
	return(raster);
}

template<typename T> 
Graph RasterGraph(VECTOR lines, int split, int spurs, RASTER<T> raster)
{
	GRAPH g;
	Graph raster_graph;
	map_vertex_type map;
	
	cout << "Building graph for raster " << raster.name << endl;
	g.ReadVEC4raster(raster_graph, raster.transform, map, lines.data);
	g.SplitFaults(raster_graph, map, split); //split the faults in the graph into fault segments, according to the intersections of the faults
	g.RemoveSpurs(raster_graph, map, spurs); //remove any spurs from the graph network
	cout << "done" << endl;
	return(raster_graph);
}

template<typename T>
void WriteSHP_r(VECTOR lines, int dist, RASTER<T> raster)
{
	GEO georef;
	std::string save_name = FGraph::add_prefix_suffix_subdirs(lines.out_path, {"raster_shp"}, lines.name, ".shp", true);
	FGraph::CreateDir(save_name);
    
	GDALAllRegister();
	const char* ref = raster.refWKT;
	const char *pszDriverName = "ESRI Shapefile";
	GDALDriver *poDriver;
	GDALDataset *poDS;
	OGRLayer *poLayer;
	OGRSpatialReference oSRS;
	oSRS.importFromWkt(&ref);

	poDriver = (GDALDriver*) GDALGetDriverByName(pszDriverName );
	if( poDriver == NULL )
	{
		printf( "%s driver not available.\n", pszDriverName );
		exit( 1 );
	}
	
	poDriver->SetDescription("RasterData");
	
	poDS = poDriver->Create(save_name.c_str(), 0, 0, 0, GDT_Unknown, NULL );
	if( poDS == NULL )
	{
		printf( "Creation of output file failed.\n" );
		exit( 1 );
	}

	poLayer = poDS->CreateLayer( (raster.name).c_str(), &oSRS, wkbLineString, NULL );
	if( poLayer == NULL )
	{
		printf( "Layer creation failed.\n" );
		exit( 1 );
	}
	
	OGRFieldDefn oField( "ID", OFTInteger );
	oField.SetWidth(10);
	if( poLayer->CreateField( &oField ) != OGRERR_NONE )
	{
		printf( "Creating 'No' field failed.\n" );
		exit( 1 );
	}
	
	OGRFieldDefn oField0( "Length", OFTReal );
	oField0.SetWidth(20);
	oField0.SetPrecision(10);
	if( poLayer->CreateField( &oField0 ) != OGRERR_NONE )
	{
		printf( "Creating 'Length' field failed.\n" );
		exit( 1 );
	}

	OGRFieldDefn oField1( "Angle", OFTReal );
	oField1.SetWidth(10);
	oField1.SetPrecision(3);
	if( poLayer->CreateField( &oField1 ) != OGRERR_NONE )
	{
		printf( "Creating 'Angle' field failed.\n" );
		exit( 1 );
	}

// raster---------------------------------------------------------------
		OGRFieldDefn oField2( "Centre", OFTReal );
		oField2.SetWidth(20);
		oField2.SetPrecision(5);
		if( poLayer->CreateField( &oField2 ) != OGRERR_NONE )
		{
			printf( "Creating 'Centre' field failed.\n" );
			exit( 1 );
		}
		
		OGRFieldDefn oField3( "MeanValue", OFTReal );
		oField3.SetWidth(20);
		oField3.SetPrecision(5);
		if( poLayer->CreateField( &oField3 ) != OGRERR_NONE )
		{
			printf( "Creating 'MeanValue' field failed.\n" );
			exit( 1 );
		}
		
		OGRFieldDefn oField4( "CentreGrad", OFTReal );
		oField4.SetWidth(20);
		oField4.SetPrecision(5);
		if( poLayer->CreateField( &oField4 ) != OGRERR_NONE )
		{
			printf( "Creating 'CentreGrad' field failed.\n" );
			exit( 1 );
		}

		OGRFieldDefn oField5( "CrossGrad", OFTReal );
		oField5.SetWidth(10);
		oField5.SetPrecision(5);
		if( poLayer->CreateField( &oField5 ) != OGRERR_NONE )
		{
			printf( "Creating 'CrossGrad' field failed.\n" );
			exit( 1 );
		}

		OGRFieldDefn oField6( "PGrad", OFTReal );
		oField6.SetWidth(10);
		oField6.SetPrecision(5);
		if( poLayer->CreateField( &oField6 ) != OGRERR_NONE )
		{
			printf( "Creating 'PGrad' field failed.\n" );
			exit( 1 );
		}

	int id = 0;
	BOOST_FOREACH(line_type line, lines.data)
	{
		OGRLineString l;
		OGRFeature *poFeature;
		geometry::unique(line);
		
		double strike = (double)(atan2(line.front().x() - line.back().x(), line.front().y() - line.back().y())) 
						* (180 / math::constants::pi<double>());
		if (strike  < 0) 
			strike  += 180;
			
		point_type cent; 
		geometry::centroid(line, cent);
		double c_val   = (double) georef.getValue<T>(cent, raster.transform, raster.values);
		double mean_v  = (double) georef.LineExtractor<T>(line, raster);
		double grad_c  = (double) georef.CentreGradient<T>(line, raster);
		double grad_cc = (double) georef.CrossGradient<T>(line, raster);
		double grad_p  = (double) georef.ParallelGradient<T>(line, raster);
						
		poFeature = OGRFeature::CreateFeature( poLayer->GetLayerDefn() );
		poFeature->SetField( "ID", id );
		poFeature->SetField( "Length", (double)geometry::length(line));
		poFeature->SetField( "Angle", strike);
//add raster data to graph shp file-------------------------------------
		poFeature->SetField( "Centre", c_val);
		poFeature->SetField( "MeanValue", mean_v);
		poFeature->SetField( "CentreGrad", grad_c);
		poFeature->SetField( "CrossGrad", grad_cc);
		poFeature->SetField( "PGrad", grad_p);

		BOOST_FOREACH(point_type P, line) 
			l.addPoint(P.x(), P.y());
		poFeature->SetGeometry( &l );
		
		if( poLayer->CreateFeature( poFeature ) != OGRERR_NONE )
		{
			printf( "Failed to create feature in shapefile.\n" );
			exit( 1 );
		}
		OGRFeature::DestroyFeature( poFeature );
		id++;
	}
	GDALClose( poDS );
}


//sort the cells by the distance from a point (closest to farthest)
template<typename T>
struct CellSorter {
CellSorter(point_type origin) { this->origin = origin; }
	bool operator () (std::pair<box, T>& lhs, std::pair<box, T>& rhs) 
	{
		point_type p1;
		point_type p2;
		geometry::centroid(lhs.first, p1);
		geometry::centroid(rhs.first, p2);

		if (geometry::distance(p1, origin) > geometry::distance(p2, origin))
			return true;
		else
			return false;
	}
	point_type origin;
};

template<typename T>
void WriteTXT(VECTOR lines, int dist, RASTER<T> raster, string filename)
{
    std::string raster_name = fs::path(raster.name).stem().string();
	string tsv_filename = FGraph::add_prefix_suffix(lines.out_path, lines.name + "_" + raster_name, ".tsv", true);
	ofstream txtF = FGraph::CreateFileStream(tsv_filename);
	
	GEO georef;
	GEOMETRIE geom;
	GEOMETRIE::GetMidSegment <geometry::model::referring_segment<point_type>> functor;
	
	int x_size   = raster.transform[6];
	int y_size   = raster.transform[7];
	
	double ul_x  = raster.transform[0];
	double ul_y  = raster.transform[3];
	double x_res = raster.transform[1];
	double y_res = raster.transform[5];
	double pdfMin = 0, pdfMax = 0, pdfMean = 0, pdfStdDev = 0;
	
	//first we build a searchabel tree from the raster data
	typedef std::pair<box, T> grid_cell;
	geometry::index::rtree< grid_cell, geometry::index::quadratic<16> > rtree;
	
	double cur_x = ul_x;
	for (int x = 0; x < x_size; x++)
	{
		double cur_y = ul_y;
		for (int y =0; y< y_size; y++)
		{
			box cur_box;
			point_type value;
			geometry::set<geometry::min_corner, 0>(cur_box, cur_x); 
			geometry::set<geometry::min_corner, 1>(cur_box, cur_y + y_res);
			geometry::set<geometry::max_corner, 0>(cur_box, cur_x + x_res); 
			geometry::set<geometry::max_corner, 1>(cur_box, cur_y);
			geometry::centroid(cur_box, value);
			rtree.insert(std::make_pair(cur_box, georef.getValue(value, raster.transform, raster.values)));
			cur_y += y_res;
		}
		cur_x += x_res;
	}

//get raster statistics-------------------------------------------------
	GDALAllRegister();
	GDALDataset  *poDataset;
	poDataset = (GDALDataset *) GDALOpen(filename.c_str(), GA_ReadOnly );
	if( poDataset == NULL )
	{
		cout << "\n ERROR: cannot open raster file " << endl;
		exit(EXIT_FAILURE);
	}
	else
	{
			GDALRasterBand *band = poDataset -> GetRasterBand(1);  
			 if (band->ComputeStatistics(false, &pdfMin, &pdfMax, &pdfMean,  &pdfStdDev, NULL, NULL))
					cout << "WARNING: cannot compute raster statistics" << endl;
	}
	GDALClose( poDataset );
	txtF	<< "Statiscis " << raster.name << "\n" 
			<< "Min \t" 	<< pdfMin << "\n" 
			<< "Max \t" 	<< pdfMax << "\n" 
			<< "Mean \t" 	<< pdfMean << "\n" 
			<< "StdDev \t" 	<< pdfStdDev << endl;
			
	txtF << "Pixel values \t" << dist << endl;
	txtF << "No \t Parallel Profile \t Cross Profile" << endl;
	
//now loop over lines to obtain profile and crossing--------------------
	int nb = 0;
	BOOST_FOREACH(line_type l, lines.data)
	{
		txtF << nb << "\t";
		std::vector<grid_cell> result_p; //results of profile
		std::vector<grid_cell> result_c; //re3sults of cross
		
//get profile values----------------------------------------------------
		BUFFER profile = geom.DefineLineBuffer(l, dist+1);
		rtree.query(geometry::index::intersects(profile), std::back_inserter(result_p));
		std::sort(result_p.begin(), result_p.end(), CellSorter<T>(l.front()));
		
		for(auto it = result_p.begin(); it != result_p.end(); ++it) 
		{
			if (it < result_p.end()-1)
				txtF << it->second << ", ";
			else
				txtF << it->second << "\t";
		}

//get crossing values---------------------------------------------------
//first get a line though the centre of the lineament
		line_type cross, m_seg;
		point_type centre, p1, p2;
		double len =  geometry::length(l);
		double rx = l.back().x() - l.front().x();
		double ry = l.back().y() - l.front().y();
		double le = sqrt(rx*rx + ry*ry);
		
		geometry::centroid(l, centre);
		
		p1.set<0>((centre.x() + ( ry/le) * len ));
		p1.set<1>((centre.y() + (-rx/le) * len ));
		p2.set<0>((centre.x() + (-ry/le) * len ));
		p2.set<1>((centre.y() + ( rx/le) * len ));
	
		geometry::append(cross, p1);
		geometry::append(cross, centre);
		geometry::append(cross, p2);
		
//here we get the segment cut by the line 
		functor.cross = cross;
		functor  = geometry::for_each_segment(l, functor );
		
//now we calculate the actual line through teh centre with orientation obtained from teh segmetn
		rx = functor.midSeg.back().x() - functor.midSeg.front().x();
		ry = functor.midSeg.back().y() - functor.midSeg.front().y();
		len =  dist * x_res ;
		le = sqrt(rx*rx + ry*ry);
		
		p1.set<0>((centre.x() + ( ry/le) * len ));
		p1.set<1>((centre.y() + (-rx/le) * len ));
		p2.set<0>((centre.x() + (-ry/le) * len ));
		p2.set<1>((centre.y() + ( rx/le) * len ));
	
		geometry::append(cross, p1);
		geometry::append(cross, centre);
		geometry::append(cross, p2);
		
		BUFFER crossing = geom.DefineLineBuffer(cross, dist+1);
		rtree.query(geometry::index::intersects(crossing), std::back_inserter(result_c));
		std::sort(result_c.begin(), result_c.end(), CellSorter<T>(cross.front()));

		for(auto it = result_c.begin(); it != result_c.end(); ++it) 
		{
			if (it < result_c.end()-1)
				txtF << it->second << ", ";
			else
				txtF << it->second << endl;;
		}
		nb++;
	}
}

template<typename T>
void GEO::AnalyseRaster(VECTOR lines, int dist, RASTER<T> raster)
{
	
	//two things are happening here:
	//1.) building a shape file 
	//2.) creating a txt file containing profiles and cross profiles
	GEO georef;
	GEOMETRIE geom;
	string filename = raster.name;
	raster = ReadRaster<T>(raster.name);
	WriteSHP_r<T>(lines, dist, raster);
	WriteTXT(lines, dist, raster, filename);
 }
template void GEO::AnalyseRaster<int>(VECTOR lines, int dist, RASTER<int> raster);
template void GEO::AnalyseRaster<float>(VECTOR lines, int dist, RASTER<float> raster);
template void GEO::AnalyseRaster<double>(VECTOR lines, int dist, RASTER<double> raster);

Graph GEO::BuildRasterGraph(VECTOR lines, int split, int spur, string raster_filename)
{
	cout << "Generating graph linked to raster " << raster_filename << endl;
	GRAPH G;
	Graph raster_graph;
	//switch case based on raster type
	const char * name = raster_filename.c_str();

	//build a map to use switch case for differnt data types
	enum {ERROR, Byte, UInt16, Int16, UInt32, Int32, Float32, Float64 };
	static std::map<string, int> d_type;
	if ( d_type.empty() )
	{
		d_type["Byte"] = Byte;
		d_type["UInt16"] = UInt16;
		d_type["Int16"] = Int16;
		d_type["UInt32"] = UInt32;
		d_type["Int32"] = Int32;
		d_type["Float32"] = Float32;
		d_type["Float64"] = Float64;
	}
	
	GDALDataset  *poDataset;
	GDALAllRegister();
	poDataset = (GDALDataset *) GDALOpen( name, GA_ReadOnly );
	GDALRasterBand  *poBand;
	poBand = poDataset->GetRasterBand( 1 );
	std::string datatype (GDALGetDataTypeName(poBand->GetRasterDataType()));
	GDALClose( poDataset );
	
	cout << "Raster datatype is " ;
	int type = d_type[datatype];
	 switch (type) 
	 {
		case Byte:
		{
			cout	<< "byte (will not perform graph analysis)" << endl;
			RASTER<char> R;
			R = ReadRaster<char>(raster_filename);
			raster_graph = RasterGraph(lines, split, spur, R);
			WriteGraph_R(raster_graph, lines, "raster");
		} break;
		 
		case UInt16:
		{
			cout << "uint16" << endl;
			RASTER<int> R;
			R = ReadRaster<int>(raster_filename);
			raster_graph = RasterGraph(lines, split, spur, R);
			AssignValuesGraph<int>(raster_graph, R);
			G.GraphAnalysis(raster_graph, lines, 10, raster_filename);
			WriteGraph_R(raster_graph, lines, "raster");
		} break;
			
		case Int16:
		{
			cout << "int16" << endl;
			RASTER<int> R;
			R = ReadRaster<int>(raster_filename);
			raster_graph = RasterGraph(lines, split, spur, R);
			AssignValuesGraph<int>(raster_graph, R);
			G.GraphAnalysis(raster_graph, lines, 10, raster_filename);
			WriteGraph_R(raster_graph, lines, "raster");
		} break;
			
		case UInt32:
		{
			cout << "uint32" << endl;
			RASTER<int> R;
			R = ReadRaster<int>(raster_filename);
			raster_graph = RasterGraph(lines, split, spur, R);
			AssignValuesGraph<int>(raster_graph, R);
			G.GraphAnalysis(raster_graph, lines, 10, raster_filename);
			WriteGraph_R(raster_graph, lines, "raster");
		} break;
			
		case Int32:
		{
			cout << "int32" << endl;
			RASTER<int> R;
			R = ReadRaster<int>(raster_filename);
			raster_graph = RasterGraph(lines, split, spur, R);
			AssignValuesGraph<int>(raster_graph, R);
			G.GraphAnalysis(raster_graph, lines, 10, raster_filename);
			WriteGraph_R(raster_graph, lines, "raster");
		} break;
		
		case Float32:
		{
			cout << "Float32" << endl;
			RASTER<float> R;
			R = ReadRaster<float>(raster_filename);
			raster_graph = RasterGraph(lines, split, spur, R);
			AssignValuesGraph<float>(raster_graph, R);
			G.GraphAnalysis(raster_graph, lines, 10, raster_filename);
			WriteGraph_R(raster_graph, lines, "raster");
		} break;
			
		case Float64:
		{
			cout << "Float64" << endl;
			RASTER<double> R;
			R = ReadRaster<double>(raster_filename);
			raster_graph = RasterGraph(lines, split, spur, R);
			AssignValuesGraph<double>(raster_graph, R);
			G.GraphAnalysis(raster_graph, lines, 10, raster_filename);
			WriteGraph_R(raster_graph, lines, "raster");

		} break;
			
		default: 
		{
			cout << " unknown" << endl;
			cout << "ERROR: Could not determine dataype of raster" << endl;
		}
		break;
	 }
	cout << " done \n" << endl;
	return(raster_graph);
}

void GEO::ReadPoints(std::string const& filename, VECTOR lines, std::pair<point_type, point_type> &source_target)
{
	point_type P;
	line_type Line;
	const char* refWKT = lines.refWKT;
	const char * name = filename.c_str();
	vector<point_type> points;
	
	OGRSpatialReference init_ref;
	init_ref.importFromWkt(&refWKT);
	const char* datum = init_ref.GetAttrValue("Datum") ;
	
	GDALAllRegister();
    GDALDataset *poDS = static_cast<GDALDataset*>(
        GDALOpenEx( name, GDAL_OF_VECTOR, NULL, NULL, NULL ));
    if( poDS == NULL )
    {
        printf( "Open failed.\n" );
        exit( 1 );
    }
    
    OGRSpatialReference * pOrigSrs = poDS->GetLayer( 0 )-> GetSpatialRef();
	const char *datum2 = pOrigSrs->GetAttrValue("Datum") ;
	if ( pOrigSrs )
	{
		if (!pOrigSrs->IsProjected() )
		{
			cout << "ERROR: vector data without spatial reference" << endl;
			exit(EXIT_FAILURE);
		}
		if ( pOrigSrs->IsProjected() )
		{
			string d1 = to_string(datum);
			if(to_string(datum)!= to_string(datum2))
			{
				cout << " ERROR: datum of point data not consitent with initial datum" << endl;
				return;
			}
		}
	}

	OGRLayer  *poLayer = poDS->GetLayer ( 0 );
	OGRFeatureDefn *poFDefn = poLayer->GetLayerDefn();
	poLayer->ResetReading();
	OGRFeature *poFeature;
	while( (poFeature = poLayer->GetNextFeature()) != NULL )
    {
        for( int iField = 0; iField < poFDefn->GetFieldCount(); iField++ )
        {
            OGRFieldDefn *poFieldDefn = poFDefn->GetFieldDefn( iField );

            switch( poFieldDefn->GetType() )
            {
                case OFTInteger:
                    poFeature->GetFieldAsInteger( iField );
                    break;
                case OFTInteger64:
                    poFeature->GetFieldAsInteger64( iField );
                    break;
                case OFTReal:
                    poFeature->GetFieldAsDouble(iField);
                    break;
                case OFTString:
                    poFeature->GetFieldAsString(iField);
                    break;
                default:
                     poFeature->GetFieldAsString(iField);
                    break;
            }
        }

        OGRGeometry *poGeometry = poFeature->GetGeometryRef();
        if( poGeometry != NULL
                && wkbFlatten(poGeometry->getGeometryType()) == wkbPoint )
        {
            OGRPoint *poPoint = (OGRPoint *) poGeometry;
            point_type p;
			geometry::set<0>(p, poPoint->getX());
			geometry::set<1>(p, poPoint->getY());
			points.push_back(p);
        }
        else
        {
            printf( "no point geometry\n" );
        }
        OGRFeature::DestroyFeature( poFeature );
    }
    GDALClose( poDS );

    if (points.size() == 2)
    {
		source_target.first  = points[0];
		source_target.second = points[1];
	}
	else
		cout << " ERROR: inconsitent point number in source - target file" << endl;
}

//read a value from teh raster at a point-------------------------------
template <typename T>
T GEO::getValue(point_type p, double transform[8], T** values)
{
	polygon_type pl = BoundingBox(transform, 0);
	if(geometry::within(p,pl))
	{
		int x = (int) (p.x() - transform[0]) / transform[1];
		int y = (int) abs((p.y() - transform[3]) / transform[5]);
		x = boost::algorithm::clamp(x, 0, (int)transform[6]);
		y = boost::algorithm::clamp(y, 0, (int)transform[7]);
		return(values[x][y]);
	}
	else 
		return std::numeric_limits<T>::quiet_NaN();
}

//this is a function to write a raster file from teh datastructure RASTER
//mainly for debug purposes but migth be useful
template<typename T>
void GEO::WriteRASTER_struc(RASTER<T> raster)
{
	GDALDataset *poDstDS;
	GDALDriver *poDriver;
	char *pszSRS_WKT = NULL;
	char **papszMetadata;
	int x = raster.transform[6];
	int y = raster.transform[7];
	const char *pszFormat = "GTiff";

	poDriver = GetGDALDriverManager()->GetDriverByName(pszFormat);
	if( poDriver == NULL )
	{
		cout << "ERROR: Could not load GDAL driver ( GTiff)" << endl;
		EXIT_FAILURE;
	}
	
	papszMetadata = poDriver->GetMetadata();

	poDstDS = poDriver->Create("test.tif", x,y, 1, GDT_Float32, NULL);
	
	poDstDS->SetGeoTransform( raster.transform );

	poDstDS->SetProjection( raster.refWKT);

	GDALRasterBand *poBand = poDstDS->GetRasterBand(1);

	float *rowBuff = (float*) CPLMalloc(sizeof(float)*x);
	cout << "Writing raster (struc) " << endl;
	progress_display * show_progress =  new boost::progress_display(y * x);
	
	for(int row = 0; row < y; row++) 
	{
		for(int col = 0; col < x; col++) 
		{
			rowBuff[col] = (float) raster.values[col][row];
				++(*show_progress);
		}
		int err = poBand->RasterIO(GF_Write, 0, row, x, 1, rowBuff, x, 1, GDT_Float32, 0, 0);
	}
	poBand->SetNoDataValue(-256);
	GDALClose( (GDALDatasetH) poDstDS );
	cout << " done \n" << endl;
}
template void GEO::WriteRASTER_struc<int>(RASTER<int> raster);
template void GEO::WriteRASTER_struc<float>(RASTER<float> raster);
template void GEO::WriteRASTER_struc<double>(RASTER<double> raster);

 void GEO::WriteRASTER(vector<vector<double>> data, char* SpatialRef, double adfGeoTransform[6], VECTOR &input_file, string suffix)
{
  //get the path of the current directory
	char cur_path[256];
	getcwd(cur_path, 255);

	GDALDataset *poDstDS;
	GDALDriver *poDriver;
	char *pszSRS_WKT = NULL;
	char **papszMetadata;
	int x = data.size();
	int y = data[0].size();
	const char *pszFormat = "GTiff";

    std::string out_name = FGraph::add_prefix_suffix_subdirs(input_file.out_path, {raster_subdir}, input_file.name, suffix, true);
	FGraph::CreateDir(out_name);

	poDriver = GetGDALDriverManager()->GetDriverByName(pszFormat);
	if( poDriver == NULL )
	{
		cout << "ERROR: Could not load GDAL driver ( GTiff)" << endl;
		EXIT_FAILURE;
	}
	
	papszMetadata = poDriver->GetMetadata();
	poDstDS = poDriver->Create(out_name.c_str(), x,y, 1, GDT_Float32, NULL);
	poDstDS->SetGeoTransform( adfGeoTransform );
	poDstDS->SetProjection( SpatialRef );
	
	GDALRasterBand *poBand = poDstDS->GetRasterBand(1);

	float *rowBuff = (float*) CPLMalloc(sizeof(float)*x);
	cout << "Writing raster " << endl;

	for(int row = 0; row < y; row++) 
	{
		for(int col = 0; col < x; col++) 
		{
			rowBuff[col] = (float) data[col][row];

		}
		int err = poBand->RasterIO(GF_Write, 0, row, x, 1, rowBuff, x, 1, GDT_Float32, 0, 0);
	}
	poBand->SetNoDataValue(-256);
	GDALClose( (GDALDatasetH) poDstDS );
	chdir(cur_path);
}
  void GEO::WriteSHP_lines(vector<line_type>lineaments, string name)
 {
 	const char *pszDriverName = "ESRI Shapefile";
    GDALDriver *poDriver;
    GDALAllRegister();
    
    GDALDataset *poDS;
    OGRLayer *poLayer;
    
    name = FGraph::add_prefix_suffix(name, "", ".shp");
    FGraph::CreateDir(name);
    const char* Name = name.c_str();
    const char* refWKT = refWKT_shp;
    
	OGRSpatialReference oSRS;
    oSRS.importFromWkt(&refWKT);
    poDriver = (GDALDriver*) GDALGetDriverByName(pszDriverName );
    if( poDriver == NULL )
    {
        printf( "%s driver not available.\n", pszDriverName );
        exit( 1 );
    }
    
	poDS = poDriver->Create(Name, 0, 0, 0, GDT_Unknown, NULL );
	if( poDS == NULL )
	{
		printf( "Creation of output file failed.\n" );
		exit( 1 );
	}
	
	poLayer = poDS->CreateLayer(Name, NULL, wkbLineString, NULL );
	if( poLayer == NULL )
	{
		printf( "Layer creation failed.\n" );
		exit( 1 );
	}
	
	OGRFieldDefn oField( "ID", OFTString );

	oField.SetWidth(32);

	if( poLayer->CreateField( &oField ) != OGRERR_NONE )
	{
		printf( "Creating ID field failed.\n" );
		exit( 1 );
	}
	
	int id = 0;
	BOOST_FOREACH(line_type line, lineaments)
	{
		geometry::unique(line);
		OGRFeature *poFeature;
        poFeature = OGRFeature::CreateFeature( poLayer->GetLayerDefn() );
        poFeature->SetField( "ID", id );

		OGRLineString l;
		BOOST_FOREACH(point_type P, line) 
			l.addPoint(P.x(), P.y());
		
		poFeature->SetGeometry( &l );

		if( poLayer->CreateFeature( poFeature ) != OGRERR_NONE )
		{
			printf( "Failed to create feature in shapefile.\n" );
			exit( 1 );
		}

		OGRFeature::DestroyFeature( poFeature );
		id++;
	}
    GDALClose( poDS );
 }
 
 void GEO::WRITE_SHP(VECTOR lineaments, string name)
 {
	STATS stats;
	GDALAllRegister();
	const char* refWKT = lineaments.refWKT;
	const char *pszDriverName = "ESRI Shapefile";
    GDALDriver *poDriver;
    GDALDataset *poDS;
    OGRLayer *poLayer;
    
    name = FGraph::add_prefix_suffix(name, "", ".shp");
    const char* Name = name.c_str();

    OGRSpatialReference oSRS;
    oSRS.importFromWkt(&refWKT);

    poDriver = (GDALDriver*) GDALGetDriverByName(pszDriverName );
    if( poDriver == NULL )
    {
        printf( "%s driver not available.\n", pszDriverName );
        exit( 1 );
    }
    
    poDriver->SetDescription(name.c_str());

	poDS = poDriver->Create(Name, 0, 0, 0, GDT_Unknown, NULL );
	if( poDS == NULL )
	{
		printf( "Creation of output file failed.\n" );
		exit( 1 );
	}
	
	poLayer = poDS->CreateLayer( (lineaments.name).c_str(), &oSRS, wkbLineString, NULL );
	
	if( poLayer == NULL )
	{
		printf( "Layer creation failed.\n" );
		exit( 1 );
	}
	
	OGRFieldDefn oField( "ID", OFTInteger );
	oField.SetWidth(10);
	if( poLayer->CreateField( &oField ) != OGRERR_NONE )
	{
		printf( "Creating 'No' field failed.\n" );
		exit( 1 );
	}
	
	OGRFieldDefn oField0( "Length", OFTReal );
	oField0.SetWidth(50);
	oField0.SetPrecision(5);
	if( poLayer->CreateField( &oField0 ) != OGRERR_NONE )
	{
		printf( "Creating 'Length' field failed.\n" );
		exit( 1 );
	}
	
	OGRFieldDefn oField1( "Angle", OFTReal );
	oField1.SetWidth(10);
	oField1.SetPrecision(3);
	if( poLayer->CreateField( &oField1 ) != OGRERR_NONE )
	{
		printf( "Creating 'Angle' field failed.\n" );
		exit( 1 );
	}
	
	if (gauss_params.size() != 0)
	{
		OGRFieldDefn oField2( "Set", OFTInteger );
		oField1.SetWidth(10);
		if( poLayer->CreateField( &oField2 ) != OGRERR_NONE )
		{
			printf( "Creating 'Set' field failed.\n" );
			exit( 1 );
		}
	}

	int id = 0;
	BOOST_FOREACH(line_type line, lineaments.data) 
	{
		OGRLineString l;
		OGRFeature *poFeature;
		geometry::unique(line);
		float L = (float) geometry::length(line);
		
		double strike = (double)(atan2(line.front().x() - line.back().x(), line.front().y() - line.back().y())) 
							* (180 / math::constants::pi<double>());
		if (strike  < 0) 
			strike  += 180;

		poFeature = OGRFeature::CreateFeature( poLayer->GetLayerDefn() );
		poFeature->SetField( "ID", id );
		poFeature->SetField( "Length", L);
		poFeature->SetField( "Angle", strike);
		if (gauss_params.size() != 0)
			poFeature->SetField( "Set", stats.CheckGaussians(strike));

		BOOST_FOREACH(point_type P, line) 
			l.addPoint(P.x(), P.y());
		poFeature->SetGeometry( &l );
		
		if( poLayer->CreateFeature( poFeature ) != OGRERR_NONE )
		{
			printf( "Failed to create feature in shapefile.\n" );
			exit( 1 );
		}
		OGRFeature::DestroyFeature( poFeature );
		id++;
	}
    GDALClose( poDS );
    cout << "created shp-file " << name << endl << endl;
 }
 
//convert graph edges to shp-file---------------------------------------
void WriteSHP_g_lines(Graph G, char *refWKT, const char* Name)
{
	GDALAllRegister();
	const char* ref = refWKT;
	const char *pszDriverName = "ESRI Shapefile";
    GDALDriver *poDriver;
    GDALDataset *poDS;
    OGRLayer *poLayer;

    OGRSpatialReference oSRS;
    oSRS.importFromWkt(&ref);

    poDriver = (GDALDriver*) GDALGetDriverByName(pszDriverName );
    if( poDriver == NULL )
    {
        printf( "%s driver not available.\n", pszDriverName );
        exit( 1 );
    }
    
    poDriver->SetDescription("graph_edges");

	poDS = poDriver->Create(Name, 0, 0, 0, GDT_Unknown, NULL );
	if( poDS == NULL )
	{
		printf( "Creation of output file failed.\n" );
		exit( 1 );
	}
	
	poLayer = poDS->CreateLayer( "graph_branches", &oSRS, wkbLineString, NULL );
	if( poLayer == NULL )
	{
		printf( "Layer creation failed.\n" );
		exit( 1 );
	}
	
	OGRFieldDefn oField( "ID", OFTInteger );
	oField.SetWidth(10);
	if( poLayer->CreateField( &oField ) != OGRERR_NONE )
	{
		printf( "Creating 'No' field failed.\n" );
		exit( 1 );
	}
	
	OGRFieldDefn oField0( "Length", OFTReal );
	oField0.SetWidth(10);
	oField0.SetPrecision(5);
	if( poLayer->CreateField( &oField0 ) != OGRERR_NONE )
	{
		printf( "Creating 'Length' field failed.\n" );
		exit( 1 );
	}

	OGRFieldDefn oField1( "BranchType", OFTString );
	oField1.SetWidth(10);
	if( poLayer->CreateField( &oField1 ) != OGRERR_NONE )
	{
		printf( "Creating 'BranchType' field failed.\n" );
		exit( 1 );
	}
	
	OGRFieldDefn oField2( "Component", OFTInteger );
	oField1.SetWidth(20);
	if( poLayer->CreateField( &oField2 ) != OGRERR_NONE )
	{
		printf( "Creating 'Component' field failed.\n" );
		exit( 1 );
	}
	
	OGRFieldDefn oField3( "Angle", OFTReal );
	oField2.SetWidth(10);
	oField2.SetPrecision(3);
	if( poLayer->CreateField( &oField3 ) != OGRERR_NONE )
	{
		printf( "Creating 'Angle' field failed.\n" );
		exit( 1 );
	}

	int id = 0;
	for (auto Eg : boost::make_iterator_range(edges(G))) 
	{
		OGRLineString l;
		OGRFeature *poFeature;
		line_type line = G[Eg].trace;
		geometry::unique(line);
		float L = (float) G[Eg].length;
		double strike = (double)(atan2(line.front().x() - line.back().x(), line.front().y() - line.back().y())) 
						* (180 / math::constants::pi<double>());
		if (strike  < 0) 
			strike  += 180;
		int C = G[Eg].component;
		
        poFeature = OGRFeature::CreateFeature( poLayer->GetLayerDefn() );
        poFeature->SetField( "ID", id );
        poFeature->SetField( "Length", L);
		poFeature->SetField( "BranchType", G[Eg].BranchType.c_str());
		poFeature->SetField( "Component", C);
		poFeature->SetField( "Angle", strike);
        
		BOOST_FOREACH(point_type P, line) 
			l.addPoint(P.x(), P.y());
		poFeature->SetGeometry( &l );
		
		if( poLayer->CreateFeature( poFeature ) != OGRERR_NONE )
		{
			printf( "Failed to create feature in shapefile.\n" );
			exit( 1 );
		}

		OGRFeature::DestroyFeature( poFeature );
		id++;
	}
	GDALClose( poDS );
}

//convert graph edges to shp-file---------------------------------------

void WriteSHP_g_lines_R(Graph G, char *refWKT, const char* Name)
{
		GDALAllRegister();
		const char* ref = refWKT;
		const char *pszDriverName = "ESRI Shapefile";
		GDALDriver *poDriver;
		GDALDataset *poDS;
		OGRLayer *poLayer;

		OGRSpatialReference oSRS;
		oSRS.importFromWkt(&ref);

		poDriver = (GDALDriver*) GDALGetDriverByName(pszDriverName );
		if( poDriver == NULL )
		{
			printf( "%s driver not available.\n", pszDriverName );
			exit( 1 );
		}
		
		poDriver->SetDescription("graph_edges_raster");

		poDS = poDriver->Create(Name, 0, 0, 0, GDT_Unknown, NULL );
		if( poDS == NULL )
		{
			printf( "Creation of output file failed.\n" );
			exit( 1 );
		}
		
		poLayer = poDS->CreateLayer( "graph_branches", &oSRS, wkbLineString, NULL );
		if( poLayer == NULL )
		{
			printf( "Layer creation failed.\n" );
			exit( 1 );
		}
		
		OGRFieldDefn oField( "ID", OFTInteger );
		oField.SetWidth(10);
		if( poLayer->CreateField( &oField ) != OGRERR_NONE )
		{
			printf( "Creating 'No' field failed.\n" );
			exit( 1 );
		}
		
		OGRFieldDefn oField0( "Length", OFTReal );
		oField0.SetWidth(10);
		oField0.SetPrecision(5);
		if( poLayer->CreateField( &oField0 ) != OGRERR_NONE )
		{
			printf( "Creating 'Length' field failed.\n" );
			exit( 1 );
		}

		OGRFieldDefn oField1( "BranchType", OFTString );
		oField1.SetWidth(10);
		if( poLayer->CreateField( &oField1 ) != OGRERR_NONE )
		{
			printf( "Creating 'BranchType' field failed.\n" );
			exit( 1 );
		}
		
		OGRFieldDefn oField2( "Component", OFTInteger );
		oField1.SetWidth(20);
		if( poLayer->CreateField( &oField2 ) != OGRERR_NONE )
		{
			printf( "Creating 'Component' field failed.\n" );
			exit( 1 );
		}
		
		OGRFieldDefn oField3( "Angle", OFTReal );
		oField3.SetWidth(10);
		oField3.SetPrecision(3);
		if( poLayer->CreateField( &oField3 ) != OGRERR_NONE )
		{
			printf( "Creating 'Angle' field failed.\n" );
			exit( 1 );
		}

	// raster---------------------------------------------------------------
			OGRFieldDefn oField4( "Centre", OFTReal );
			oField4.SetWidth(25);
			oField4.SetPrecision(5);
			if( poLayer->CreateField( &oField4 ) != OGRERR_NONE )
			{
				printf( "Creating 'Centre' field failed.\n" );
				exit( 1 );
			}
			
			OGRFieldDefn oField5( "MeanValue", OFTReal );
			oField5.SetWidth(25);
			oField5.SetPrecision(5);
			if( poLayer->CreateField( &oField5 ) != OGRERR_NONE )
			{
				printf( "Creating 'MeanValue' field failed.\n" );
				exit( 1 );
			}
			
			OGRFieldDefn oField6( "CentreGrad", OFTReal );
			oField6.SetWidth(25);
			oField6.SetPrecision(5);
			if( poLayer->CreateField( &oField6 ) != OGRERR_NONE )
			{
				printf( "Creating 'CentreGrad' field failed.\n" );
				exit( 1 );
			}

			OGRFieldDefn oField7( "CrossGrad", OFTReal );
			oField7.SetWidth(25);
			oField7.SetPrecision(5);
			if( poLayer->CreateField( &oField7 ) != OGRERR_NONE )
			{
				printf( "Creating 'CrossGrad' field failed.\n" );
				exit( 1 );
			}

			OGRFieldDefn oField8( "PGrad", OFTReal );
			oField8.SetWidth(10);
			if( poLayer->CreateField( &oField8 ) != OGRERR_NONE )
			{
				printf( "Creating 'PGrad' field failed.\n" );
				exit( 1 );
			}

		int id = 0;
		for (auto Eg : boost::make_iterator_range(edges(G))) 
		{
			OGRLineString l;
			OGRFeature *poFeature;
			line_type line = G[Eg].trace;
			geometry::unique(line);
			float L = (float) G[Eg].length;
			double strike = (double)(atan2(line.front().x() - line.back().x(), line.front().y() - line.back().y())) 
							* (180 / math::constants::pi<double>());
			if (strike  < 0) 
				strike  += 180;
							
			int C = G[Eg].component;
			poFeature = OGRFeature::CreateFeature( poLayer->GetLayerDefn() );
			poFeature->SetField( "ID", id );
			poFeature->SetField( "Length", L);
			poFeature->SetField( "BranchType", G[Eg].BranchType.c_str());
			poFeature->SetField( "Component", C);
			poFeature->SetField( "Angle", strike);
	//add raster data to graph shp file-------------------------------------
			poFeature->SetField( "Centre", 		(float) G[Eg].Centre);
			poFeature->SetField( "MeanValue", 	(float) G[Eg].MeanValue);
			poFeature->SetField( "CentreGrad", 	(float) G[Eg].CentreGrad);
			poFeature->SetField( "CrossGrad", 	(float) G[Eg].CrossGrad);
			poFeature->SetField( "PGrad", 		(float) G[Eg].ParalGrad );

			BOOST_FOREACH(point_type P, line) 
				l.addPoint(P.x(), P.y());
			poFeature->SetGeometry( &l );
			
			if( poLayer->CreateFeature( poFeature ) != OGRERR_NONE )
			{
				printf( "Failed to create feature in shapefile.\n" );
				exit( 1 );
			}

			OGRFeature::DestroyFeature( poFeature );
			id++;
		}
		GDALClose( poDS );
}

void GEO::WriteSHP_maxFlow(DGraph G, string name)
{
	const char *pszDriverName = "ESRI Shapefile";
    GDALDriver *poDriver;
    GDALAllRegister();
    
    GDALDataset *poDS;
    OGRLayer *poLayer;
    
    name = FGraph::add_prefix_suffix(name, "", ".shp");
    FGraph::CreateDir(name);
    const char* Name = name.c_str();
    const char* refWKT = refWKT_shp;
    
	OGRSpatialReference oSRS;
    oSRS.importFromWkt(&refWKT);
    poDriver = (GDALDriver*) GDALGetDriverByName(pszDriverName );
    if( poDriver == NULL )
    {
        printf( "%s driver not available.\n", pszDriverName );
        exit( 1 );
    }
    
	poDS = poDriver->Create(Name, 0, 0, 0, GDT_Unknown, NULL );
	if( poDS == NULL )
	{
		printf( "Creation of output file failed.\n" );
		exit( 1 );
	}
	
	poLayer = poDS->CreateLayer(Name, NULL, wkbLineString, NULL );
	if( poLayer == NULL )
	{
		printf( "Layer creation failed.\n" );
		exit( 1 );
	}

	OGRFieldDefn oField( "No", OFTInteger );
	oField.SetWidth(10);
	if( poLayer->CreateField( &oField ) != OGRERR_NONE )
	{
		printf( "Creating 'No' field failed.\n" );
		exit( 1 );
	}

	OGRFieldDefn oField1( "Capacity", OFTReal );
	oField1.SetWidth(10);
	oField1.SetPrecision(5);
	if( poLayer->CreateField( &oField1 ) != OGRERR_NONE )
	{
		printf( "Creating 'Capacity' field failed.\n" );
		exit( 1 );
	}

	OGRFieldDefn oField2( "Residual", OFTReal );
	oField2.SetWidth(10);
	oField2.SetPrecision(5);
	if( poLayer->CreateField( &oField2 ) != OGRERR_NONE )
	{
		printf( "Creating 'Data' field failed.\n" );
		exit( 1 );
	}
	
	OGRFieldDefn oField3( "Fill", OFTReal );
	oField3.SetWidth(10);
	oField3.SetPrecision(5);
	if( poLayer->CreateField( &oField3 ) != OGRERR_NONE )
	{
		printf( "Creating 'Data' field failed.\n" );
		exit( 1 );
	}
	
	OGRFieldDefn oField4( "CapUsed", OFTReal );
	oField4.SetWidth(10);
	oField4.SetPrecision(5);
	if( poLayer->CreateField( &oField4 ) != OGRERR_NONE )
	{
		printf( "Creating 'Capacity Used' field failed.\n" );
		exit( 1 );
	}
	
	OGRFieldDefn oField5( "RelCapUsed", OFTReal );
	oField5.SetWidth(10);
	oField5.SetPrecision(5);
	if( poLayer->CreateField( &oField5 ) != OGRERR_NONE )
	{
		printf( "Creating 'Relative Capacity Used' field failed.\n" );
		exit( 1 );
	}
	
	OGRFieldDefn oField6( "AbsCapUsed", OFTReal );
	oField6.SetWidth(10);
	oField6.SetPrecision(5);
	if( poLayer->CreateField( &oField6 ) != OGRERR_NONE )
	{
		printf( "Creating 'Absolute Capacity Used' field failed.\n" );
		exit( 1 );
	}

	int NO = 0;
	OGRFeature *poFeature;
	poFeature = OGRFeature::CreateFeature( poLayer->GetLayerDefn() );
	//write a WKT file and shp file-----------------------------------------
	for (auto Eg : boost::make_iterator_range(edges(G))) 
	{
		line_type fault = G[Eg].trace;
		poFeature->SetField( "No", NO  );
		poFeature->SetField( "Capacity", G[Eg].capacity);
		poFeature->SetField( "Residual", G[Eg].residual_capacity);

		//convenience names
		const double  cap = G[Eg].capacity;
		const double rcap = G[Eg].residual_capacity;
		
		poFeature->SetField( "Fill", (rcap/cap));
		poFeature->SetField( "CapUsed", (cap-rcap));
		poFeature->SetField("RelCapUsed", (0 < cap) && (rcap <= cap) ? 1 - rcap/cap : 0); //relative amount of capacity used
		poFeature->SetField("AbsCapUsed", (0 < cap) && (rcap <= cap) ? cap - rcap : 0); //absolute amount of capacity used

		OGRLineString l;
		BOOST_FOREACH(point_type P,fault) 
			l.addPoint(P.x(), P.y());
		
		poFeature->SetGeometry( &l );
		if( poLayer->CreateFeature( poFeature ) != OGRERR_NONE )
		{
			printf( "Failed to create feature in shapefile.\n" );
			exit( 1 );
		}
		NO++;
	}
	OGRFeature::DestroyFeature( poFeature );
	GDALClose( poDS );
}

//convert graph edges to shp-file---------------------------------------
void WriteSHP_g_points(Graph G, char* refWKT, const char* Name)
{
	point_type point;
	std::string Point;
	GDALAllRegister();
	GDALDataset *poDS;
	GDALDriver *poDriver;
	OGRLayer *poLayer;
	OGRFeature *poFeature;
	const char* ref = refWKT;
	const char* pszDriverName = "ESRI Shapefile";
	poDriver = GetGDALDriverManager()->GetDriverByName(pszDriverName );
	poDS = poDriver->Create(Name, 0, 0, 0, GDT_Unknown, NULL );
	
	OGRSpatialReference oSRS;
	oSRS.importFromWkt(&ref);
	
	if( poDS == NULL )
	{
		printf( "Creation of shp-file failed.\n" );
		exit( 1 );
	}

	poLayer = poDS->CreateLayer(Name, &oSRS, wkbPoint, NULL );
	if( poLayer == NULL )
	{
		printf( "Layer creation failed.\n" );
		exit( 1 );
	}
	
	OGRFieldDefn oField( "No", OFTInteger );
	oField.SetWidth(10);
	if( poLayer->CreateField( &oField ) != OGRERR_NONE )
	{
		printf( "Creating 'No' field failed.\n" );
		exit( 1 );
	}
	
	OGRFieldDefn oField1( "Degree", OFTInteger );
	oField1.SetWidth(10);
	if( poLayer->CreateField( &oField1 ) != OGRERR_NONE )
	{
		printf( "Creating 'Degree' field failed.\n" );
		exit( 1 );
	}
	
	OGRFieldDefn oField2( "Component", OFTInteger );
	oField2.SetWidth(10);
	if( poLayer->CreateField( &oField2 ) != OGRERR_NONE )
	{
		printf( "Creating 'Component' field failed.\n" );
		exit( 1 );
	}

	int NO = 0;
	poFeature = OGRFeature::CreateFeature( poLayer->GetLayerDefn() );
	//write shp file----------------------------------------------------
	for (auto Ve : boost::make_iterator_range(vertices(G))) 
	{
		point = G[Ve].location;
		int de = degree(Ve, G);
		int co = G[Ve].component;
		float da = G[Ve].data;

		poFeature->SetField( "No", NO);
		poFeature->SetField( "Degree", de);
		poFeature->SetField( "Component", co);

		OGRPoint PO;
		PO.setX(point.x());
		PO.setY(point.y());
		poFeature->SetGeometry( &PO );
		if( poLayer->CreateFeature( poFeature ) != OGRERR_NONE )
		{
			printf( "Failed to create feature in shapefile.\n" );
			exit( 1 );
		}
		Point.clear();
		NO++;
	}
	OGRFeature::DestroyFeature( poFeature );
	GDALClose( poDS );
}

void GEO::WriteGraph(Graph G, VECTOR lines, string subF)
{
	cout << "starting writegraph" << endl;
	assert (num_vertices(G) != 0 && num_edges(G) != 0);
	char cur_path[256];
	char* reference = lines.refWKT;
	getcwd(cur_path, 255);
	
    string subdir_name = FGraph::add_prefix_suffix_subdirs(lines.out_path, {"graph_shp", subF}, "", "/");
    FGraph::CreateDir(subdir_name);
	
	string n_b =  FGraph::add_prefix_suffix(subdir_name, "", lines.name + "_branches.shp");
	string n_v =  FGraph::add_prefix_suffix(subdir_name, "", lines.name + "_vertices.shp");
	const char* Name_b = n_b.c_str();
	const char* Name_v = n_v.c_str();

	WriteSHP_g_lines(G, reference, Name_b);
	WriteSHP_g_points(G, reference, Name_v);
	cout << " done" << endl;
}

//TODO: Why is this a separate function from the above?
void GEO::WriteGraph_R(Graph G, VECTOR lines, string subF)
{
	cout << "starting writegraph (raster)" << endl;
	assert (num_vertices(G) != 0 && num_edges(G) != 0);
	char cur_path[256];
	char* reference = lines.refWKT;
	getcwd(cur_path, 255);
	string subdir_name = FGraph::add_prefix_suffix_subdirs(lines.out_path, {"graph_shp", subF}, "", "/");
	FGraph::CreateDir(subdir_name);
	
	string n_b =  FGraph::add_prefix_suffix(subdir_name, "", lines.name + "_branches.shp");
	string n_v =  FGraph::add_prefix_suffix(subdir_name, "", lines.name + "_vertices.shp");
	const char* Name_b = n_b.c_str();
	const char* Name_v = n_v.c_str();

	WriteSHP_g_lines_R(G, reference, Name_b);
	WriteSHP_g_points(G, reference, Name_v);
	cout << " done" << endl;
}

void GEO::Point_Tree(vector<p_index> points,  vector<p_index>& closest)
{
	p_tree rt(points.begin(), points.end());

	for ( size_t i = 0 ; i < points.size() ; ++i )
		rt.query(geometry::index::nearest(points[i].first, 1) && geometry::index::satisfies(different_id_p(i)), std::back_inserter(closest) );        
}

void GEO::Point_Tree2(vector<pl_index> points, vector<pl_index>& closest, double max_len)
{
 pl_tree rt(points.begin(), points.end());
 
 for ( size_t i = 0 ; i < points.size() ; ++i )
 {
	if (get<2>(points[i]) >= max_len)
	{
		closest.push_back(make_tuple(get<0>(points[i]),i,0));
	}
	else
	{
		rt.query(geometry::index::nearest(get<0>(points[i]), 1) && geometry::index::satisfies(different_id_pl(get<1>(points[i])))
			&& geometry::index::satisfies(larger(get<2>(points[i]))),
				std::back_inserter(closest) );  
	} 
 }         
}



//calculate the angle (in radians) between the connecting end points of two line_type's that sharea  common point
//the x_front boolean values are true if the user is checking the front of that line segment, and false if the user wants to know the angle from the end of that line segment
//returns close to 0 if the lines are close to having the same angle, or larger angles if the line_types are not in line with each other
double CalculateAngleDifference(line_type &a, bool a_front, line_type &b, bool b_front)
{
	//get the two relevant end points of the line segments, as it is the first/last two points that define the ending segment of the line_type
	point_type a1, a2, b1, b2;
	//normalise the points so that they are in the order (b2 -> b1) <-> (a1 -> a2)
	//(matching to the "front" of a and the "back" of b, relative to this diagram)
	if (a_front)
	{
		a1 = a[0];
		a2 = a[1];
	} else {
		const int s = a.size();
		a1 = a[s-1];
		a2 = a[s-2];
	}
	if (b_front)
	{
		b1 = b[0];
		b2 = b[1];
	} else {
		const int s = b.size();
		b1 = b[s-1];
		b2 = b[s-2];
	}

	//calculate the segments
	//remember, these are relative to the diagram above
	const double dxa = a2.x() - a1.x();
	const double dya = a2.y() - a1.y();
	const double dxb = b1.x() - b2.x();
	const double dyb = b1.y() - b2.y();
	//dot product: a . b == a.x * b.x + a.y * b.y == |a| |b| cos(theta) -> theta = acos(a . b / (|a| |b|))
	const double la = std::sqrt(std::pow(dxa, 2) + std::pow(dya, 2));
	const double lb = std::sqrt(std::pow(dxb, 2) + std::pow(dyb, 2));
	const double dotprod = dxa*dxb + dya*dyb;
	const double theta = acos(dotprod / (la*lb)); //calculate angle, and subract from 1pi radians/180 degrees to get the angle difference from being a straight line
	return theta;
}

//   typedef the rtree here so both functions can see the same type. not putting it in the header because no other function uses it
typedef list<line_type> unmerged_type; //might want to change this to a set or map
typedef std::tuple<point_type, unmerged_type::iterator, bool> endpoint_value_type;
typedef geometry::index::rtree<endpoint_value_type, geometry::index::rstar<16>> endpoint_rtree_type;

//convenience function to remove enfpoints from the endpoint rtree, given an iterator
void remove_endpoints(endpoint_rtree_type &rtree, list<line_type>::iterator it)
{
	rtree.remove(std::make_tuple(it->front(), it, true ));
	rtree.remove(std::make_tuple(it->back (), it, false));
}

//sometimes the data has faults that are split into different entries in the shapefile
//this function checks the vector of input line_types for faults that were split, and merges them back together
//this checks for any and all fault segments that should be merged back with base
bool MergeConnections(unmerged_type &faults, endpoint_rtree_type &endpoints, line_type &base, double distance_threshold)
{
	//our thresholds for merging fault segments together are:
	//1) Their endpoints are within threshold distance of each other
	//2) Only two fault segments meet at the intersection point of the two candidate fault segments
	//3) The angle between the two candidate fault segments is within angle_threshold (in radians)
	//remember that a single fault could be split into more than two fault segments
	//change: the damage zone distance depends on the length of the fault, so we need to ensure that they are all joined properly
	//therefore, we still want to merge fault segments when there are multiple other intersecting faults
	//and we choose to merge the segment with the smallest angle difference
	
	bool changed;
	const double angle_threshold = 45 * math::constants::pi<double>() / 180; //25 degrees, in radians
	
	const double dt = distance_threshold; //convenience name
	
	do {
		//a front match is one that matches to the front of the base fault/line segment, a back match is one that matches to the back of the base fault/line segment
		//use faults.end() as a sentinel value to mark that a match hasn't been found
		unmerged_type::iterator front_match = faults.end(), back_match = faults.end(), candidate;
		bool front_match_loc = false, back_match_loc = false; //true iff we're matching to the front of the other linestring
		std::vector<std::tuple<unmerged_type::iterator, bool>> front_matches, back_matches;
		bool match_loc;
		//use boxes, because boost doesn't have a circle geometry object (just a polygon/linestring that approximates a circle)
		box front_box(point_type(base.front().x() - dt, base.front().y() - dt), point_type(base.front().x() + dt, base.front().y() + dt));
		box  back_box(point_type(base.back ().x() - dt, base.back ().y() - dt), point_type(base.back ().x() + dt, base. back().y() + dt));
		vector<endpoint_value_type> candidates;
		endpoints.query(geometry::index::intersects(front_box), back_inserter(candidates));
		for (auto cand_it = candidates.begin(); cand_it != candidates.end(); cand_it++)
		{
			unmerged_type::iterator comp_it = std::get<1>(*cand_it);
			//the line segments could be in any order, so we need to check all four pairs of endpoints
			if (geometry::distance(base.front(), comp_it->front()) <= distance_threshold) front_matches.push_back(make_tuple(comp_it,  true));
			if (geometry::distance(base.front(), comp_it-> back()) <= distance_threshold) front_matches.push_back(make_tuple(comp_it, false));
		}
		
		//now check for possible matches against the back of this (base) fault
		candidates.clear();
		endpoints.query(geometry::index::intersects(back_box), back_inserter(candidates));
		for (auto cand_it = candidates.begin(); cand_it != candidates.end(); cand_it++)
		{
			unmerged_type::iterator comp_it = std::get<1>(*cand_it);
			if (geometry::distance(base.back() , comp_it->front()) <= distance_threshold)  back_matches.push_back(make_tuple(comp_it,  true));
			if (geometry::distance(base.back() , comp_it-> back()) <= distance_threshold)  back_matches.push_back(make_tuple(comp_it, false));
		}
			
		changed = false;
		//check for matches to the front of base
		double best_angle = std::numeric_limits<double>::infinity();
		for (auto match_it = front_matches.begin(); match_it != front_matches.end(); match_it++)
		{
			std::tie(candidate, match_loc) = *match_it;
			const double theta = CalculateAngleDifference(base, true, *candidate, match_loc);
			const double angle_diff = abs(theta);
			if (angle_diff <= angle_threshold && angle_diff < best_angle)
			{
				front_match = candidate;
				best_angle = angle_diff;
				front_match_loc = match_loc;
			}
		}
		//if there is more than one candidate, check to see if they should be merged into each other instead
		if (front_match != faults.end() && front_matches.size() >= 2)
		{
			for (auto it1 = front_matches.begin(); it1 != front_matches.end(); it1++)
			{
				unmerged_type::iterator f1 = std::get<0>(*it1);
				bool front1 = std::get<1>(*it1);
				for (auto it2 = it1; it2 != front_matches.end(); it2++)
				{
					if (it1 == it2) continue;
					unmerged_type::iterator f2 = std::get<0>(*it2);
					bool front2 = std::get<1>(*it2);
					double angle_diff = abs(CalculateAngleDifference(*f1, front1, *f2, front2));
					if (angle_diff < best_angle)
					{
						front_match = faults.end();
						break;
					}
				}
				if (front_match == faults.end()) break;
			}
		}
		//check for matches to the back of base
		best_angle = std::numeric_limits<double>::infinity(); //reset the best angle varaible
		for (auto match_it = back_matches.begin(); match_it != back_matches.end(); match_it++)
		{
			std::tie(candidate, match_loc) = *match_it;
			double theta = CalculateAngleDifference(base, false, *candidate, match_loc);
			const double angle_diff = abs(theta);
			if (angle_diff <= angle_threshold && angle_diff < best_angle)
			{
				back_match = candidate;
				best_angle = angle_diff;
				back_match_loc = match_loc;
			}
		}
		
		//if there is more than one candidate, check to see if they should be merged into each other instead
		if (back_match != faults.end() && back_matches.size() >= 2)
		{
			for (auto it1 = back_matches.begin(); it1 != back_matches.end(); it1++)
			{
				unmerged_type::iterator f1 = std::get<0>(*it1);
				bool front1 = std::get<1>(*it1);
				for (auto it2 = it1; it2 != back_matches.end(); it2++)
				{
					if (it1 == it2) continue;
					unmerged_type::iterator f2 = std::get<0>(*it2);
					bool front2 = std::get<1>(*it2);
					double angle_diff = abs(CalculateAngleDifference(*f1, front1, *f2, front2));
					if (angle_diff < best_angle)
					{
						back_match = faults.end();
						break;
					}
				}
				if (back_match == faults.end()) break;
			}
		}
		
		
		
// 		this is clunky, but we need to wait to erase the front match, because we need the iterator to be valid long enough to compare it to the back match
		if ( (front_match != faults.end()) && (front_match == back_match) ) 
		{
			//in this case, the same other fault is considered a match for either end of the fault segment currently being considered
			//so figure out which end it other fault should be matched to

			double front_angle = abs(CalculateAngleDifference(base, true, *front_match, front_match_loc));
			double back_angle  = abs(CalculateAngleDifference(base, false, *back_match,  back_match_loc));
			//mark the one with a worse angle as the one that doesn't match
			if (front_angle <= back_angle) back_match = faults.end();
			else front_match = faults.end();
		}
		
		if (front_match != faults.end())
		{
			line_type prepend = *front_match;
			if (front_match_loc) geometry::reverse(prepend);
			geometry::append(prepend, base);
			base = prepend;
			remove_endpoints(endpoints, front_match);
			faults.erase(front_match);
			changed = true;
		}
		
		if (back_match != faults.end()){
			line_type append = *back_match;
			if (!back_match_loc) geometry::reverse(append);
			geometry::append(base, append);
			remove_endpoints(endpoints, back_match);
			faults.erase(back_match);
			changed = true;
		}
		
	} while (changed);
	return true;
}

//sort by distance from point
struct LineSorter_orig {
LineSorter_orig(point_type origin) { this->origin = origin; }
    bool operator () (line_type& lhs, line_type& rhs) 
    {
		point_type center_lhs, center_rhs;
		geometry::centroid(lhs, center_lhs);
		geometry::centroid(rhs, center_rhs);
        if (geometry::distance(center_lhs, origin) < geometry::distance(center_rhs, origin))
          return true;
        else
          return false;
	}
	point_type origin;
};

//compare two line types after sorting
struct LineCompare {
	bool operator()(const line_type& lhs, const line_type& rhs){
		 if (geometry::equals(lhs, rhs))
          return true;
        else
          return false;
	}
};

//sometimes the data has faults that are split into different entries in the shapefile
//this function checks the vector of input line_types for faults that were split, and merges them back together
void GEO::CorrectNetwork(vector<line_type>&F, double dist)
{
	//first we need to remove dublicates. We sort the lines by distance to a reference point (from line centroid)
	//and then remove every dublicates
	GEOMETRIE geom;
	box AOI = geom.ReturnAOI(F);
	point_type origin(geometry::get<geometry::min_corner, 0>(AOI), geometry::get<geometry::min_corner, 1>(AOI));
	bool found_dublicate = true;
	
		do{
			int size = F.size();
			std::sort(F.begin(), F.end(), LineSorter_orig(origin));
			std::vector<line_type>::iterator new_lines_3;
			new_lines_3 = std::unique(F.begin(), F.end(), LineCompare()); 
			F.erase(new_lines_3, F.end());
			
			if (size == F.size())
				found_dublicate = false;
			
		}	while (found_dublicate);
	cout << F.size() << " lines remaining after removing dublicates \n" << endl;
//----------------------------------------------------------------------
	
	vector<line_type> merged_faults; //the new list of merged faults, to be built by this function
	unmerged_type unmerged_faults; //a list of faults that have not been merged yet. the get removed from this list once they've been added to merged_faults, either by being merged with another fault segment, or on its own
	unmerged_faults.insert(std::end(unmerged_faults), std::begin(F), std::end(F));
	
	endpoint_rtree_type endpoint_tree;
	for (auto it = unmerged_faults.begin(); it != unmerged_faults.end(); it++)
	{
		endpoint_tree.insert(std::make_tuple(it->front(), it, true ));
		endpoint_tree.insert(std::make_tuple(it->back (), it, false));
	}
	
	cout <<"merging split faults - starting" <<endl;
	while (!unmerged_faults.empty())
	{
		line_type base = unmerged_faults.front();
		remove_endpoints(endpoint_tree, unmerged_faults.begin());
		unmerged_faults.pop_front(); //I don't know why this doesn't also return the value being pop'd
		MergeConnections(unmerged_faults, endpoint_tree, base, dist);
		merged_faults.push_back(base);
	}
	cout <<"merging split faults - finished" <<endl;
	F = merged_faults;
	cout << F.size() << " lines remaining after merging" << endl;
}



void GEO::Get_Source_Target(const char* Name, point_type &Source, point_type &Target)
{
	char* refWKT;
	GDALAllRegister();
	GDALDataset *poDS = static_cast<GDALDataset*>
	(
		GDALOpenEx( Name, GDAL_OF_VECTOR, NULL, NULL, NULL )
	);
	if( poDS == NULL )
	{
		printf( " Opening source/target shapefile \"%s\" failed.\n", Name );
		exit( 1 );
	}  
	
	OGRSpatialReference * pOrigSrs = poDS->GetLayer( 0 )-> GetSpatialRef();
	if ( pOrigSrs )
	{
		if (!pOrigSrs->IsProjected() )
		{
			cout << "ERROR: vector data without spatial reference /not a projected reference system" << endl;
			exit(EXIT_FAILURE);
		}
		if ( pOrigSrs->IsProjected() )
			 pOrigSrs->exportToWkt( &refWKT );
	}
	OGRLayer  *poLayer = poDS->GetLayer( 0 );
	
	poLayer->ResetReading();
	OGRFeature *poFeature;
	int p = 0;
	while( (poFeature = poLayer->GetNextFeature()) != NULL )
	{
		if (p > 1)
		{
			cout << "more than two points in shape file" << endl;
			break;
		}
        OGRGeometry *poGeometry = poFeature->GetGeometryRef();
        if( poGeometry != NULL
                && wkbFlatten(poGeometry->getGeometryType()) == wkbPoint )
        {
            OGRPoint *poPoint = (OGRPoint *) poGeometry;
            if ( p == 0)
            {
				geometry::set<0>(Source, poPoint->getX());
				geometry::set<1>(Source, poPoint->getY());
			}
			else
			{
				geometry::set<0>(Target, poPoint->getX());
				geometry::set<1>(Target, poPoint->getY());
			}
        }
        else
        {
            printf( "no point geometry\n" );
        }
        OGRFeature::DestroyFeature( poFeature );
        p++;
    }
	GDALClose( poDS );
	cout << setprecision(10) << "Source: " << Source.x() << ", " << Source.y() << "\n"
		 << "Target: " << Target.x() << ", " << Target.y() << endl;
}
 
