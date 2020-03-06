#include "GeoRef.h"
#include "graph.h"
#include "geometrie.h"
#include "stats.h"

GEO::GEO ()
{
}

double GeoTransform[8];


// create and open filestream in folder "statistics"
string CreateDir(VECTOR &input_file, std::initializer_list<string> folders)
{
	
	string folder_name(input_file.folder);
	while (folder_name.back() == '/' ) folder_name.pop_back(); //gdal/ESRI shapefiles don't like double slashes at the start
	while (folder_name.back() == '\\') folder_name.pop_back();
	if (folder_name.size() <= 0) folder_name = ".";
	for (auto it = folders.begin(); it != folders.end(); it++) folder_name += "/"+*it;
	const char* folder_cstr = folder_name.c_str();
	//if(!opendir(folder_cstr))
		//mkdir(folder_cstr, 0777);
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
			cout << "ERROR: vector data without spatial reference" << endl;
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
		if( poGeometry != NULL
				&& wkbFlatten(poGeometry->getGeometryType()) == wkbLineString)
		{
			OGRLineString *poLine = (OGRLineString *) poGeometry;
			for (int i = 0; i < poLine->getNumPoints(); i++)
			{
				P.set<0>(poLine->getX(i));
				P.set<1>(poLine->getY(i));
				geometry::append(Line, P);
			}
//correction for double points and digitized circles---------------------
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

//assign raster values (primarily elevation, from a digital elevation map (DEM)) to the elements of a graph
//same effect as AssignValues() above, but more efficient
//read the entire file at once, rather than one at a time. gives faster IO. This assumes that we can fit the entire file in memory at once.

template <typename T>
void GEO::AssignValuesGraph(Graph& G, RASTER raster)
{
	STATS stats;
	line_type curEdge;
	point_type centre;
	
	for (auto Ve : boost::make_iterator_range(vertices(G)))
	{
		G[Ve].data = getValue(G[Ve].location, raster.transform, raster.values);
		if (std::isnan(G[Ve].data)) cout <<"Error: vertex " << Ve << " read a value of nan" << endl;
	}

	cout << "edge" << endl;

	for (auto Eg : boost::make_iterator_range(edges(G))) 
	{
		curEdge = G[Eg].trace;
		geometry::centroid(G[Eg].trace, centre);

	    G[Eg].Centre     = getValue(centre, raster.transform, raster.values);
	//	G[Eg].MeanValue  = LineExtractor(G[Eg].trace, transform, values) ;
		//G[Eg].CentreGrad = CentreGradient(G[Eg].trace, raster.transform, raster.values);
// 		cout << G[Eg].CentreGrad << endl;
	//	G[Eg].CrossGrad  = CrossGradient(G[Eg].trace, transform, values);
	//	G[Eg].ParalGrad  = ParallelGradient(G[Eg].trace, transform, values) ;
	}
	cout <<"done" << endl;
}
//make these instatiations so that other cpp files can find them
template void GEO::AssignValuesGraph<double>(Graph& G, RASTER raster);
template void GEO::AssignValuesGraph<float>(Graph& G, RASTER raster);
template void GEO::AssignValuesGraph<int>(Graph& G, RASTER raster);

//make a directed graph from an undirected graph
//assumes that AssignValuesGraph() has already been called on the graph, and an appropriate raster file, to initialise the height values (which are stored in FVertex.data, not FVertex.elevation)
DGraph GEO::MakeDirectedGraph(Graph &g)
{
	DGraph dg; //this is the directed graph object that will be made, then returned
	
	std::map<vertex_type, dvertex_type> vertex_map; //use this to keep track of which nodes in the new graph correspond to particular nodes in the existing graph
	
	//start by making the indices
	Graph::vertex_iterator v, vstart, vend;
	boost::tie(vstart, vend) = boost::vertices(g);
	for (v = vstart; v < vend; v++)
	{
		FVertex<point_type> fv = g[*v];
		DVertex dv(fv.location, fv.Enode, fv.data, v - vstart);
		dvertex_type dvd = boost::add_vertex(dv, dg);
		vertex_map[*v] = dvd;
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
		
// 		cout << "added edge (" << ds f<< ", " << dt << ") from " << *e << endl;
	}
	return dg;
}

//set up the capacity values for the maximum flow algorithm
void GEO::setup_maximum_flow(DGraph &dg)
{
	DGraph::edge_iterator e, estart, eend;
	boost::tie(estart, eend) = boost::edges(dg);
	for (e = estart; e != eend; e++)
	{
// 		dvertex_type s = source(*e, dg), t = target(*e, dg);
// 		double hs = dg[s].elevation, ht = dg[t].elevation;

		DEdge de = dg[*e];
		const double area = (1 - std::exp(-de.full_length/100)); //this will become a function later
		//const double area = (1 - std::exp(-de.full_length));
		const double flow = area;// * (std::max<double>(hs - ht, 0)) / de.length;
		dg[*e].capacity = flow;
		dg[*e].residual_capacity = flow;
		
// 		cout << "capacity of edge " << *e << " from " << dg[s].index << " to " << dg[t].index << " is " << dg[*e].capacity << endl;
	}
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
	if (dg[t].elevation > dg[s].elevation){
		std::swap(s, t);
		cout << "Swapping source and target to find maximum flow from higher location to lower location" << endl;
	}
	
	//sanity check, make sure we're not calculating the max flow from a vertex to itself
	if (s == t){
		cout << "Warning: The source and target points for the maximum flow (" << s << ") are the same" << endl;
		return std::numeric_limits<double>::infinity();
	}
	
	cout << "Calculating maximum flow from vertex " << s << " (" << dg[s].elevation << "m) to vertex " << t << " (" << dg[t].elevation << ")" << endl;
	
	//clear capacities for verticies that are above the source
	const double start_elevation = dg[s].elevation;
	vector<std::pair<dedge_type, double>> saved_capacities;
	DGraph::edge_iterator e, estart, eend;
	boost::tie(estart, eend) = boost::edges(dg);
	for (e = estart; e != eend; e++)
	{
		dvertex_type cs = boost::source(*e, dg); //check source
		dvertex_type ct = boost::target(*e, dg); //check target
		if (dg[cs].elevation > start_elevation || dg[ct].elevation > start_elevation)
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
	double max_flow = maximum_flow(dg, s, t);
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

//calculate the maximum flow from source s to sink/target t
double GEO::maximum_flow(DGraph &dg, dvertex_type s, dvertex_type t)
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

polygon_type GEO::BoundingBox(double transform[8])
{
	long double Xmax, Xmin, Ymax, Ymin;
	polygon_type pl;
	box bx;

	//create new box of raster
	Xmax = transform[0] + transform[1] * transform[6] ;    // west
	Xmin = transform[0] ; 								  // east
	Ymax = transform[3] ; 								 // north
	Ymin = transform[3] + transform[5] * transform[7] ;	// south
	
	bx.min_corner().set<0>( Xmin );
	bx.min_corner().set<1>( Ymin );
	bx.max_corner().set<0>( Xmax );
	bx.max_corner().set<1>( Ymax );
	geometry::convert(bx, pl);
	
	return(pl);
}

double GEO::CrossGradient(line_type F, double transform[8], double** values)
{
	GEOMETRIE::Perpencicular <geometry::model::referring_segment<point_type>> functor;
	polygon_type pl = BoundingBox(transform);
	vector<double> D;

	functor  = geometry::for_each_segment(F, functor );

	BOOST_FOREACH(line_type cross, functor.all)
	{
		if (geometry::within(cross.front(), pl) && geometry::within(cross.back(),pl))
		{
			if(getValue(cross.front(), transform, values) > getValue(cross.back(), transform, values))
				D.push_back(
				(getValue(cross.front(), transform, values) - getValue(cross.back(), transform, values))/geometry::length(cross)
				);
			else
				D.push_back( (
				getValue(cross.back(), transform, values) - getValue(cross.front(), transform, values))/geometry::length(cross)
				);
		}
	}

	vec DATA(D);
	if (DATA.size() != 0)
		return( arma::mean(DATA) );
	else
		return(0);
}

double GEO::ParallelGradient(line_type F, double transform[8], double** values)
{
	GEOMETRIE::SegPoints	 <geometry::model::referring_segment<point_type>> functor2;
	polygon_type pl = BoundingBox(transform);

	functor2 = geometry::for_each_segment(F, functor2);
	vector<double>D;

	for (vector<std::tuple<point_type, point_type, point_type>>::
	const_iterator I = functor2.Points.begin(); I != functor2.Points.end(); ++I) 
	{
// 		int no = 0;
// 		double grad1 = 0, grad2 = 0; //unused
		if (geometry::within(get<0>(*I), pl) && geometry::within(get<1>(*I), pl) && geometry::within(get<2>(*I), pl))
		{
			D.push_back((getValue(get<0>(*I), transform, values) - getValue(get<1>(*I), transform, values))  +
						(getValue(get<1>(*I), transform, values) - getValue(get<2>(*I), transform, values)));
		}
	}
	vec DATA(D);
	if (DATA.size() > 0)
		return(arma::mean(DATA));
	else
		return(0);
}

template <typename T>
double GEO::CentreGradient(line_type F, double transform[8], T** values)
{
	//Convert from map to pixel coordinates and read values at centre of line
	//Only works for geotransforms with no rotation.
	line_type cross;
	point_type point, p1, p2;
	polygon_type pl = BoundingBox(transform);
	double len =  ceil((abs(transform[1]) + abs(transform[5])) / 2);
	
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

	if (geometry::within(p1, pl) && geometry::within(p2, pl))
	{	
		/*
		if (getValue(p1, transform, values) > getValue(p2, transform, values))
			return( (getValue(p1, transform, values) - getValue(p2,transform, values)) / geometry::length(cross) );
		else
			return( (getValue(p2, transform, values) - getValue(p1,transform, values)) / geometry::length(cross));
			*/
		return (abs(getValue(p1, transform, values) - getValue(p2,transform, values))) / geometry::length(cross);
	}
	else
		return(0);
}

VECTOR GEO::ReadVector(int argc, char* f)
{
	VECTOR lineaments;
	string file = f;
	string ext, name;
	
	if (argc == 1+1 )
	{
		ext =  strrchr(f,'.') ;
		name = file.substr(0, file.find_last_of("."));
		
		if (ext == ".txt" ||  ext == ".shp")
			cout << "Fault  data from " << ext << "-file." << endl;
		else
		{
			cout << "ERROR: Wrong vector format provided  (" << ext << ")" << endl;
			exit (EXIT_FAILURE);
		}
	}
	else
	{
		cout << " ERROR: Not enough input argmuments! (Please provide a vector file)" << endl;
		exit (EXIT_FAILURE);
	}
	
	
	size_t first_folder = name.find_first_of("/\\");
	while (first_folder != std::string::npos && first_folder + 1 < name.size() && (name[first_folder+1] == '/' || name[first_folder+1] == '\\')) name.erase(first_folder+1, 1);
	
	
	size_t folder_index = name.find_last_of("/\\"); //this should work on both linux (/) and windows (\\) 
	string folder = (folder_index != std::string::npos) ? name.substr(0               , folder_index + 1 ) : "./";
	name          = (folder_index != std::string::npos) ? name.substr(folder_index + 1, std::string::npos) : name;
	
	lineaments.folder = folder;
	lineaments.name = name;

	//read information from the vector file
	if (ext == ".shp")
		read_shp(file, lineaments);
	else
		cerr << "ERROR: no shape file definend" << endl;
		
	return(lineaments);
}

template<typename T>
RASTER ReadRaster(const string filename)
{
	/*
	[0]  top left x 
	[1]  w-e pixel resolution 
	[2]  0 
	[3]  top left y 
	[4]  0 
	[5]  n-s pixel resolution (negative value) 
	*/
	cout << "Reading raster " << filename;
	float** raster_data;
	const char * name = filename.c_str();
	RASTER raster;
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
			
			raster.values = RasterConvert(poDataset->GetRasterXSize(), poDataset->GetRasterYSize(), raster.values);
			
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


Graph GEO::RasterGraph(VECTOR lines, int split, int spurs, string filename)
{
	const char * name = filename.c_str();
	Graph raster_graph;
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
	
	GRAPH g;
	RASTER R;
	map_vertex_type map;
	int type = d_type[datatype];
	 switch (type) 
	 {
		case Byte:
		{
			cout << "byte" << endl;
			R = ReadRaster<char>(filename);
			cout << " " << endl;
		} break;
		 
		case UInt16:
		{
			cout << "uint16" << endl;
			R = ReadRaster<unsigned short>(filename);
		} break;
			
		case Int16:
		{
			cout << "int16" << endl;
			R = ReadRaster<short>(filename);
		} break;
			
		case UInt32:
		{
			cout << "uint32" << endl;
			R = ReadRaster<unsigned int>(filename);
		} break;
			
		case Int32:
		{
			cout << "int32" << endl;
			R = ReadRaster<int>(filename);
		} break;
		
		case Float32:
		{
			cout << "Float32" << endl;
			R = ReadRaster<float>(filename);
		} break;
			
		case Float64:
		{
			cout << "Float64" << endl;
			R = ReadRaster<double>(filename);
		} break;
			
		default: 
			cout << "ERROR: Could not determine dataype of raster" << endl;
			break;
	 }
	 
	 cout << "Building graph for raster " << R.name << endl;
	 g.ReadVEC4raster(raster_graph, R, map, lines.data);
	 g.SplitFaults(raster_graph, map, split); //split the faults in the graph into fault segments, according to the intersections of the faults
	 g.RemoveSpurs(raster_graph, map, spurs); //remove any spurs from the graph network
	 cout << "done" << endl;
	 
	 AssignValuesGraph<float>(raster_graph, R);
	 
	 
	 
	 cout << "finished" << endl;

	return(raster_graph);
}

void GEO::ReadPoints(std::string const& filename, VECTOR lines, std::pair<point_type, point_type> &source_target)
{
	point_type P;
	line_type Line;
	const char * name = filename.c_str();
	vector<point_type> points;
	
	OGRSpatialReference init_ref;
	init_ref.importFromWkt(&lines.refWKT);
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




template <typename T>
T GEO::getValue(point_type p, double transform[8], T** values)
{
	polygon_type pl = BoundingBox(transform);
	if(geometry::within(p,pl))
	{
		int x = (int) (p.x() - transform[0]) / transform[1];
		int y = (int) abs((p.y() - transform[3]) / transform[5]);
		x = boost::algorithm::clamp(x, 0, (int)transform[6]);
		y = boost::algorithm::clamp(y, 0, (int)transform[7]);
		
		cout << values[x][y] << endl;
		
		return(values[x][y]);
	}
	else 
		return std::numeric_limits<double>::quiet_NaN();
}

double GEO::LineExtractor(line_type L, double Transform[8], double** raster)
{
	polygon_type pl = BoundingBox(Transform);
	polygon_type pl2;
	GEOMETRIE geom;
	box AOI;
	BUFFER envelop;
	
	point_type point;
	int maxX, minX, maxY, minY;
	double radius = 1; 
	vector<double> D;
	
    radius = Transform[1] * Transform[5] /2;
	envelop = geom.DefineLineBuffer(L, radius);
	geometry::envelope(envelop, AOI);
	geometry::convert(AOI, pl2);
	maxX = (int)(geometry::get<geometry::max_corner, 0>(AOI) - Transform[0]) / Transform[1];
	minX = (int)(geometry::get<geometry::min_corner, 0>(AOI) - Transform[0]) / Transform[1]; 
	maxY = (int)abs((geometry::get<geometry::max_corner, 1>(AOI) - Transform[3]) / Transform[5]);
	minY = (int)abs((geometry::get<geometry::min_corner, 1>(AOI) - Transform[3]) / Transform[5]);
	
	if (geometry::within(pl2, pl))
	{
		for (int x = minX; x < maxX; x++)
		{
			point.set<0>(Transform[0] + x * Transform[1]);
			for (int y = maxY; y < minY; y++)
			{
				point.set<1>(Transform[3] + y * Transform[5]);
				
				if (geometry::within(point, envelop))
					D.push_back(raster[x][y]);
			}
		}
	}
	else
		return(0);
		
	vec DATA(D);
	if (DATA.size() > 0)
		return(arma::mean(DATA));
	else
	return(0);
 }
  
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

	string subdir = CreateDir(input_file, {"raster"});

	poDriver = GetGDALDriverManager()->GetDriverByName(pszFormat);
	if( poDriver == NULL )
	{
		cout << "ERROR: Could not load GDAL driver ( GTiff)" << endl;
		EXIT_FAILURE;
	}
	
	papszMetadata = poDriver->GetMetadata();
	poDstDS = poDriver->Create((subdir + "/" + input_file.name + suffix).c_str(), x,y, 1, GDT_Float32, NULL);
	poDstDS->SetGeoTransform( adfGeoTransform );
	poDstDS->SetProjection( SpatialRef );
	
	GDALRasterBand *poBand = poDstDS->GetRasterBand(1);

	float *rowBuff = (float*) CPLMalloc(sizeof(float)*x);
	cout << "Writing raster: " << endl;
	progress_display * show_progress =  new boost::progress_display(y * x);
	for(int row = 0; row < y; row++) 
	{
		for(int col = 0; col < x; col++) 
		{
			rowBuff[col] = (float) data[col][row];
				++(*show_progress);
		}
		poBand->RasterIO(GF_Write, 0, row, x, 1, rowBuff, x, 1, GDT_Float32, 0, 0);
	}
	GDALClose( (GDALDatasetH) poDstDS );
	chdir(cur_path);
}
 
 void GEO::WriteSHP_lines(vector<line_type>lineaments, const char* Name)
 {
 	const char *pszDriverName = "ESRI Shapefile";
    GDALDriver *poDriver;
    GDALAllRegister();
    
    GDALDataset *poDS;
    OGRLayer *poLayer;

    poDriver = (GDALDriver*) GDALGetDriverByName(pszDriverName );
    if( poDriver == NULL )
    {
        printf( "%s driver not available.\n", pszDriverName );
        exit( 1 );
    }
    
	poDS = poDriver->Create( "component.shp", 0, 0, 0, GDT_Unknown, NULL );
	if( poDS == NULL )
	{
		printf( "Creation of output file failed.\n" );
		exit( 1 );
	}
	
	poLayer = poDS->CreateLayer( "point_out", NULL, wkbLineString, NULL );
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
 
 
//convert graph edges to shp-file---------------------------------------
void WriteSHP_g_lines(Graph G, char *refWKT, const char* Name, bool raster)
{
	GDALAllRegister();
	const char *pszDriverName = "ESRI Shapefile";
    GDALDriver *poDriver;
    GDALDataset *poDS;
    OGRLayer *poLayer;

    OGRSpatialReference oSRS;
    oSRS.importFromWkt(&refWKT);

    poDriver = (GDALDriver*) GDALGetDriverByName(pszDriverName );
    if( poDriver == NULL )
    {
        printf( "%s driver not available.\n", pszDriverName );
        exit( 1 );
    }
    
    poDriver->SetDescription("bla");

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
	oField1.SetWidth(10);
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

//only if raster file was loaded----------------------------------------
	if (raster)
	{
		OGRFieldDefn oField4( "Component", OFTReal );
		oField2.SetWidth(10);
		if( poLayer->CreateField( &oField4 ) != OGRERR_NONE )
		{
			printf( "Creating 'Component' field failed.\n" );
			exit( 1 );
		}

		OGRFieldDefn oField5( "Centre", OFTReal );
		oField3.SetWidth(10);
		oField3.SetPrecision(5);
		if( poLayer->CreateField( &oField5 ) != OGRERR_NONE )
		{
			printf( "Creating 'Centre' field failed.\n" );
			exit( 1 );
		}
		
		OGRFieldDefn oField6( "MeanValue", OFTReal );
		oField4.SetWidth(10);
		oField4.SetPrecision(5);
		if( poLayer->CreateField( &oField6 ) != OGRERR_NONE )
		{
			printf( "Creating 'MeanValue' field failed.\n" );
			exit( 1 );
		}
		
		OGRFieldDefn oField7( "CentreGrad", OFTReal );
		oField5.SetWidth(10);
		oField5.SetPrecision(5);
		if( poLayer->CreateField( &oField7 ) != OGRERR_NONE )
		{
			printf( "Creating 'CentreGrad' field failed.\n" );
			exit( 1 );
		}

		OGRFieldDefn oField8( "CrossGrad", OFTReal );
		oField6.SetWidth(10);
		oField6.SetPrecision(5);
		if( poLayer->CreateField( &oField8 ) != OGRERR_NONE )
		{
			printf( "Creating 'CrossGrad' field failed.\n" );
			exit( 1 );
		}

		OGRFieldDefn oField9( "ParalGrad", OFTReal );
		oField7.SetWidth(10);
		if( poLayer->CreateField( &oField9 ) != OGRERR_NONE )
		{
			printf( "Creating 'Data' field failed.\n" );
			exit( 1 );
		}
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

void GEO::WriteSHP_maxFlow(DGraph G, const char* Name)
{
	
	line_type fault;
 	std::string line;

	GDALAllRegister();
	GDALDataset *poDS;
	GDALDriver *poDriver;
	OGRLayer *poLayer;
	OGRFeature *poFeature;
	const char *pszDriverName = "ESRI Shapefile";

	poDriver = GetGDALDriverManager()->GetDriverByName(pszDriverName );
	poDS = poDriver->Create(Name, 0, 0, 0, GDT_Unknown, NULL );
	if( poDS == NULL )
	{
		printf( "Creation of shp-file failed.\n" );
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

	int NO = 1;
	poFeature = OGRFeature::CreateFeature( poLayer->GetLayerDefn() );
	//write a WKT file and shp file-----------------------------------------
	for (auto Eg : boost::make_iterator_range(edges(G))) 
	{
		fault = G[Eg].trace;
		poFeature->SetField( "No", NO  );
		poFeature->SetField( "Capacity", G[Eg].capacity);
		poFeature->SetField( "Residual", G[Eg].residual_capacity);
		
		//cout << G[Eg].residual_capacity << endl;
		
		//convenience names
		const double  cap = G[Eg].capacity;
		const double rcap = G[Eg].residual_capacity;
		
		poFeature->SetField( "Fill", (rcap/cap));
		
		poFeature->SetField( "CapUsed", (cap-rcap));
		
		poFeature->SetField("RelCapUsed", (0 < cap) && (rcap <= cap) ? 1 - rcap/cap : 0); //relative amount of capacity used
		
		poFeature->SetField("AbsCapUsed", (0 < cap) && (rcap <= cap) ? cap - rcap : 0); //absolute amount of capacity used

		OGRLineString l;
//  		line.append("LINESTRING(");
		BOOST_FOREACH(point_type P,fault) 
		{
// 			line.append(std::to_string(P.x()) + " " + std::to_string(P.y()));
//  			if (!geometry::equals(P, fault.back()))
//  				line.append(", ");
			l.addPoint(P.x(), P.y());
		}
// 		line.append( ")");
//  		const char* branch = (const char*) line.c_str();
//  		l.importFromWkt(&branch);
		poFeature->SetGeometry( &l );
		if( poLayer->CreateFeature( poFeature ) != OGRERR_NONE )
		{
			printf( "Failed to create feature in shapefile.\n" );
			exit( 1 );
		}
 		line.clear();
		NO++;
	}
	OGRFeature::DestroyFeature( poFeature );
	GDALClose( poDS );
}



//convert graph edges to shp-file---------------------------------------
void WriteSHP_g_points(Graph G, char* refWKT, const char* Name, bool raster)
{
	point_type point;
	std::string Point;
	GDALAllRegister();
	GDALDataset *poDS;
	GDALDriver *poDriver;
	OGRLayer *poLayer;
	OGRFeature *poFeature;
	
	const char *pszDriverName = "ESRI Shapefile";
	poDriver = GetGDALDriverManager()->GetDriverByName(pszDriverName );
	poDS = poDriver->Create(Name, 0, 0, 0, GDT_Unknown, NULL );
	
	OGRSpatialReference oSRS;
    oSRS.importFromWkt(&refWKT);
	
	
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
	
	OGRFieldDefn oField3( "Data", OFTReal );
	oField3.SetWidth(10);
	oField3.SetPrecision(5);
	if( poLayer->CreateField( &oField3 ) != OGRERR_NONE )
	{
		printf( "Creating 'Data' field failed.\n" );
		exit( 1 );
	}
	
	cout << "here" << endl;
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
		poFeature->SetField( "Data", da);
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

void GEO::WriteGraph(Graph G, VECTOR lines, string subF, bool raster)
{
	cout << "starting writegraph" << endl;
	assert (num_vertices(G) != 0 && num_edges(G) != 0);
	char cur_path[256];
	char* reference = lines.refWKT;
	getcwd(cur_path, 255);
	string subdir_name = CreateDir(lines, {"Vector", subF});
	
	string n_b =  subdir_name + "/" + lines.name + "_branches.shp";
	string n_v =  subdir_name + "/" + lines.name + "_vertices.shp";
	const char* Name_b = n_b.c_str();
	const char* Name_v = n_v.c_str();

	WriteSHP_g_lines(G, reference, Name_b, false);
	WriteSHP_g_points(G, reference, Name_v, false);
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
	if (get<2>(points[i]) > max_len)
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

void GEO::Point_Tree3(vector<p_index> points,  vector<p_index>& closest, int nb)
{
	p_tree rt(points.begin(), points.end());
	
	for (size_t i = 0 ; i < points.size() ; ++i )
	if (i < nb)
		continue;
	else
	{
		rt.query(geometry::index::nearest(points[i].first, 1) && geometry::index::satisfies(different_id_p(i)),
		std::back_inserter(closest) );   
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
	const double angle_threshold = 25 * math::constants::pi<double>() / 180; //25 degrees, in radians
	
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

//sometimes the data has faults that are split into different entries in the shapefile
//this function checks the vector of input line_types for faults that were split, and merges them back together
void GEO::CorrectNetwork(vector<line_type>&F, double dist)
{
	vector<line_type> merged_faults; //the new list of merged faults, to be built by this function
	unmerged_type unmerged_faults; //a list of faults that have not been merged yet. thye get removed from this list once they've been added to merged_faults, either by being merged with another fault segment, or on its own
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
}

void GEO::Source_Target(const char* Name, point_type &Source, point_type &Target)
{
	
}
 
