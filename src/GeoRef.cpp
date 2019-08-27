#include "GeoRef.h"
#include "geometrie.h"
#include "stats.h"

GEO::GEO ()

{
	
}

double GeoTransform[8];
OGRSpatialReference GeoRef;
const char *GeoProj;

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
T** GEO::RasterConvert(int rows, int cols, T **M)
{
    M = new T*[rows];
    for (int i = 0; i < rows; i++){
        M[i] = new T[cols];
    }
    return M;
}

//assign raster values (primarily elevation, from a digital elevation map (DEM)) to the elements of a graph
void GEO::AssignValues(Graph& G, std::string const& filename)
{
	for (auto vd : boost::make_iterator_range(vertices(G))) 
		G[vd].elevation = getElevation(G[vd].location, filename);
}

//assign raster values (primarily elevation, from a digital elevation map (DEM)) to the elements of a graph
//same effect as AssignValues() above, but more efficient
//read the entire file at once, rather than one at a time. gives faster IO. This assumes that we can fit the entire file in memory at once.
void GEO::AssignValuesAll(Graph& G, std::string const& filename)
{
	int *data = nullptr;
	readRaster<int>(filename, data);
	for (auto vd : boost::make_iterator_range(vertices(G))) 
	{       
		point_type p = G[vd].location;
		G[vd].elevation = getElevationFromArray(p, data);
		G[vd].Pressure = G[vd].elevation * 9.81 * 1;
	}
	CPLFree(data);
}

void GEO::AssignValuesGraph(Graph& G, double transform[8], double** values)
{
	STATS stats;
	line_type curEdge;
	point_type centre;

	for (auto Ve : boost::make_iterator_range(vertices(G)))
		G[Ve].data = getValue(G[Ve].location, transform, values);
	
	cout << "edge" << endl;
	
	for (auto Eg : boost::make_iterator_range(edges(G))) 
	{
		curEdge = G[Eg].trace;
		geometry::centroid(G[Eg].trace, centre);
		vertex_type U = source(Eg, G);
		vertex_type u = target(Eg, G);
		
		G[Eg].Centre = getValue(centre, transform, values);
		G[Eg].MeanValue = LineExtractor(G[Eg].trace, transform, values) ;
		G[Eg].CentreGrad = CentreGradient(G[Eg].trace, transform, values);
		G[Eg].CrossGrad = CrossGradient(G[Eg].trace, transform, values);
		G[Eg].ParalGrad =  ParallelGradient(G[Eg].trace, transform, values) ;
	}
	cout <<"done" << endl;
}

polygon_type GEO::BoundingBox(double transform[8], double raster_crop_size)
{
	long double Xmax, Xmin, Ymax, Ymin;
	polygon_type pl;
	box bx;
	
	//create new box of raster
	Xmax = transform[0] + transform[1] * transform[6] - raster_crop_size;      // west
	Xmin = transform[0] + raster_crop_size; 								  // east
	Ymax = transform[3] - raster_crop_size; 								 // north
	Ymin = transform[3] + transform[5] * transform[7] + raster_crop_size;	// south
		
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
	polygon_type pl = BoundingBox(transform, 10);
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
	polygon_type pl = BoundingBox(transform, 10);

	functor2 = geometry::for_each_segment(F, functor2);
	vector<double>D;

	for (vector<std::tuple<point_type, point_type, point_type>>::
	const_iterator I = functor2.Points.begin(); I != functor2.Points.end(); ++I) 
	{
		int no = 0;
		double grad1 = 0, grad2 = 0;
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
double GEO::CentreGradient(line_type F, double transform[8], double** values)
{
	const double raster_crop_size = 10;
	//Convert from map to pixel coordinates and read values at centre of line
	//Only works for geotransforms with no rotation.
	box bx;
	line_type cross;
	point_type point, p1, p2;
	long double Xmax, Xmin, Ymax, Ymin;
	polygon_type pl = BoundingBox(transform, raster_crop_size);
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
		if (getValue(p1, transform, values) > getValue(p2, transform, values))
			return( (getValue(p1, transform, values) - getValue(p2,transform, values)) / geometry::length(cross) );
		else
			return( (getValue(p2, transform, values) - getValue(p1,transform, values)) / geometry::length(cross));
	}
	else
		return(0);
}

//convenience function to read a value from a 2D array, using a point as an index
int GEO::getElevationFromArray(point_type p, const int *data){
	int xind = (p.x() - GeoTransform[0]) / GeoTransform[1];
	int yind = abs((p.y() - GeoTransform[3]) / GeoTransform[5]);
	return data[yind*((int)GeoTransform[6])+xind];
}

//read in the contents of a raster file to a user-specified one-dimensional array, for a given datatype
//Free the data with CPLFree(data)
template <typename T>
int GEO::readRaster(std::string const filename, T *&data){
	const char * name = filename.c_str();
	GDALDataset  *poDataset;
	GDALAllRegister();

	poDataset = (GDALDataset *) GDALOpen( name, GA_ReadOnly );
	if( poDataset == NULL )
		std::cout<<"no file"<< std::endl;

	GDALRasterBand *band = poDataset -> GetRasterBand(1);
	cout << "xSize = " << GeoTransform[6] << ", XSize = " << band->GetXSize() << ", ySize = " << GeoTransform[7] << ", YSize = " << band->GetYSize() << endl;
	int xSize = GeoTransform[6];
	int ySize = GeoTransform[7];
	data = (T *) CPLMalloc(sizeof(*data) * xSize * ySize); //the buffer that gets filled up
	
	if (data == NULL)
	{
		cerr << "ERROR: unable to allocate " << sizeof(*data) << " x " << xSize << " x " << ySize << " = " << sizeof(*data) * xSize * ySize << " bytes to read in raster file" << endl;
		return -1;
	}
	
	int err = band->RasterIO(GF_Read,0,0,xSize,ySize,data,xSize,ySize,GetGDALDataType<T>(),0,0,nullptr); //GDT_Int32

	if (err != 0)
		cout << "ERROR: cannot read elevation data!" << endl;
	
	GDALClose( poDataset );
	
	return err;
}

 //read pixel value of DEM at coodinate
int GEO::getElevation(point_type p, std::string const& filename)
{
	const char * name = filename.c_str();
	GDALDataset  *poDataset;
	GDALAllRegister();
	int value = 0, X, Y, Data ;

	poDataset = (GDALDataset *) GDALOpen( name, GA_ReadOnly );
	if( poDataset == NULL ){
		std::cout<<"no file \"" << filename << "\"" << " to get elevation data from" << std::endl;
		exit(-1);
	}

	GDALRasterBand *band = poDataset -> GetRasterBand(1);       

	X = (p.x() - GeoTransform[0]) / GeoTransform[1]; 
	Y = abs((p.y() - GeoTransform[3]) / GeoTransform[5]);

	Data = band->RasterIO(GF_Read,X,Y,1,1,&value,1,1,band->GetRasterDataType(),1,1,nullptr);

	if (Data != 0)
		cout << "ERROR: cannot read elevation data!" << endl;

	GDALClose( poDataset );
	return (value);
}

double GEO::getValue(point_type p, double transform[8], double** values)
{
	polygon_type pl = BoundingBox(transform, 10);
	if(geometry::within(p,pl))
	{
		int x = (int) (p.x() - transform[0]) / transform[1];
		int y = (int) abs((p.y() - transform[3]) / transform[5]);
		return(values[x][y]);
	}
	else 
		return std::numeric_limits<double>::quiet_NaN();
}

double GEO::LineExtractor(line_type L, double Transform[8], double** raster)
{
	polygon_type pl = BoundingBox(Transform, 10);
	polygon_type pl2;
	GEOMETRIE geom;
	box AOI;
	BUFFER envelop;
	
	point_type point;
	int maxX, minX, maxY, minY;
	double M = 0, radius;
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



//obtain information about raster's projection--------------------------
void GEO::GetRasterProperties(std::string const& filename, double**& RASTER)
{
	/*
	[0]  top left x 
	[1]  w-e pixel resolution 
	[2]  0 
	[3]  top left y 
	[4]  0 
	[5]  n-s pixel resolution (negative value) 
	[6]  X-size
	[7]  Y-size
	*/
	
	GDALDataset  *poDataset;
	GDALAllRegister();
	double adfGeoTransform[6];
	int data;
	int value = 0;

	poDataset = (GDALDataset *) GDALOpen( filename.c_str(), GA_ReadOnly );
	if( poDataset == NULL ){
		std::cout<<"no file \"" << filename << "\" to read raster properties from" << std::endl;
		exit(-1);
	}

	printf( "Driver: %s/%s\n",
		poDataset->GetDriver()->GetDescription(),
		poDataset->GetDriver()->GetMetadataItem( GDAL_DMD_LONGNAME ) );

	printf( "Size is %dx%dx%d\n",
		poDataset->GetRasterXSize(), poDataset->GetRasterYSize(),
		poDataset->GetRasterCount() );

	if( poDataset->GetProjectionRef()  != NULL )
	{
		GeoProj = poDataset->GetProjectionRef() ;
		GeoRef.SetProjCS(poDataset->GetProjectionRef());
	}
	
	if( poDataset->GetGeoTransform( adfGeoTransform ) == CE_None )
		for (int i = 0; i < 6; i++)
			GeoTransform[i] = adfGeoTransform[i];
			
	printf( "Origin = (%.6f,%.6f)\n",
		adfGeoTransform[0], adfGeoTransform[3] );
		
	printf( "Pixel Size = (%.6f,%.6f)\n",
		adfGeoTransform[1], adfGeoTransform[5] ); 
		
	for (int i = 0; i < 6; i++)
		GeoTransform[i] = adfGeoTransform[i]; 
				
	GeoTransform[6] = poDataset->GetRasterXSize();
	GeoTransform[7] = poDataset->GetRasterYSize();

	//create multi-array from data------------------------------------------

	GDALRasterBand *band = poDataset -> GetRasterBand(1);        

	RASTER = RasterConvert(poDataset->GetRasterXSize(), poDataset->GetRasterYSize(), RASTER);
	for (int i = 0; i < poDataset->GetRasterXSize(); i++)
	{
		for (int j = 0; j < poDataset->GetRasterYSize(); j++)
		{
			data = band->RasterIO(GF_Read,i,j,1,1,&value,1,1,band->GetRasterDataType(),1,1,nullptr);
			if (data == 0)
				RASTER[i][j] = value;
			if (data != 0)
				RASTER[i][j] = 0;
		}
	}
	cout << "Converted input raster to array. \n" << endl;
	GDALClose( poDataset );
}
 
//read vector data from txt file (WKT format)--------------------------
void GEO::read_wkt(std::string const& filename, std::vector<line_type>& lineString)
{
	string name;	 
	line_type geometry;
	ifstream cpp_file(filename.c_str());
	if (cpp_file.is_open())
	{
		while ( !cpp_file.eof() )
		{
			string line;
			getline(cpp_file, line);
			boost::trim(line);
			if ( !line.empty() && !boost::starts_with(line, "#") )
			{
				string::size_type pos = line.find(";");
				if (pos != std::string::npos)
				{
					name = line.substr(pos + 1);
					line.erase(pos);
					trim(line);
					trim(name);
				}
				geometry::read_wkt(line, geometry);
				lineString.push_back(geometry);
				geometry.clear();
			}
		}
	} else {
		cout << "ERROR: Failure reading WKT vector file: " << filename << endl;
		exit(EXIT_FAILURE);
	}
	cout << "read " << lineString.size() << " faults from txt" << endl;
}
 
//read vector data from shp file---------------------------------------
void GEO::read_shp(std::string const& filename, std::vector<line_type>& lineString)
{
	const char * name = filename.c_str();
// 	const char * layer;
	line_type Line;
	point_type P;
	string f = filename.substr(filename.find_last_of('/')+1); 
// 	cout << "f.find(\".\") = " << f.find(".") << " of " << f.size() << ", f.substr = " << f.substr(0, f.find(".")) << endl;
// 	layer = ;
	std::string layer_name = f.substr(0, f.find("."));
// 	cout << "layer_name = \"" << layer_name << "\"" << endl;;

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

// 	cout << "f = " << f << ", layer_name = " << layer_name << ", layer_name.c_str = " << layer_name.c_str() <<endl;

	OGRLayer  *poLayer = poDS->GetLayerByName( layer_name.c_str() );
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
			lineString.push_back(Line);
			Line.clear();
		}
		OGRFeature::DestroyFeature( poFeature );
	}
	GDALClose( poDS );
	cout << "read " << lineString.size() << " faults from shp" << endl;
}
 
 //this creates an output txt-file---------------------------------------
void GEO::WriteTxt(Graph& g, string const& filename)
{
	string output = (string) filename.c_str() + "_results.txt";
	ofstream txtF (output);	

	if (txtF.is_open())  
	{ 
		for (auto vd : boost::make_iterator_range(vertices(g))) 
			txtF << fixed << setprecision(12) << vd << " " << g[vd].location.x() << " " << g[vd].location.y() <<" " <<g[vd].elevation << " " << degree(vd, g) << " " << g[vd].data << "\n";
		txtF.close();
	}
	else cout << "Ne results written" << endl;
}
 

//convert graph edges to shp-file---------------------------------------
void GEO::WriteSHP(Graph G, const char* Name)
{
	line_type fault;
	std::string line;

	GDALAllRegister();
	GDALDataset *poDS;
	GDALDriver *poDriver;
	OGRLayer *poLayer;
	OGRFeature *poFeature;
	OGRLineString l;
	const char *pszDriverName = "ESRI Shapefile";

	char * pszWKT;
	GeoRef.exportToWkt( &pszWKT );

	poDriver = GetGDALDriverManager()->GetDriverByName(pszDriverName );
	poDS = poDriver->Create(Name, 0, 0, 0, GDT_Unknown, NULL );
	if( poDS == NULL )
	{
		printf( "Creation of shp-file failed.\n" );
		exit( 1 );
	}

	poLayer = poDS->CreateLayer(Name, &GeoRef, wkbLineString, NULL );
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

	OGRFieldDefn oField0( "Length", OFTReal );
	oField0.SetWidth(10);
	oField0.SetPrecision(5);
	if( poLayer->CreateField( &oField0 ) != OGRERR_NONE )
	{
		printf( "Creating 'Length' field failed.\n" );
		exit( 1 );
	}

	OGRFieldDefn oField1( "BranchType", OFTString );
	oField1.SetWidth(5);
	if( poLayer->CreateField( &oField1 ) != OGRERR_NONE )
	{
		printf( "Creating 'BranchType' field failed.\n" );
		exit( 1 );
	}

	OGRFieldDefn oField2( "Component", OFTString );
	oField2.SetWidth(10);
	if( poLayer->CreateField( &oField2 ) != OGRERR_NONE )
	{
		printf( "Creating 'Component' field failed.\n" );
		exit( 1 );
	}

	OGRFieldDefn oField3( "Centre", OFTReal );
	oField3.SetWidth(10);
	oField3.SetPrecision(5);
	if( poLayer->CreateField( &oField3 ) != OGRERR_NONE )
	{
		printf( "Creating 'Centre' field failed.\n" );
		exit( 1 );
	}
	
	OGRFieldDefn oField4( "MeanValue", OFTReal );
	oField4.SetWidth(10);
	oField4.SetPrecision(5);
	if( poLayer->CreateField( &oField4 ) != OGRERR_NONE )
	{
		printf( "Creating 'MeanValue' field failed.\n" );
		exit( 1 );
	}
	
	OGRFieldDefn oField5( "CentreGrad", OFTReal );
	oField5.SetWidth(10);
	oField5.SetPrecision(5);
	if( poLayer->CreateField( &oField5 ) != OGRERR_NONE )
	{
		printf( "Creating 'CentreGrad' field failed.\n" );
		exit( 1 );
	}

	OGRFieldDefn oField6( "CrossGrad", OFTReal );
	oField6.SetWidth(10);
	oField6.SetPrecision(5);
	if( poLayer->CreateField( &oField6 ) != OGRERR_NONE )
	{
		printf( "Creating 'CrossGrad' field failed.\n" );
		exit( 1 );
	}

	OGRFieldDefn oField7( "ParalGrad", OFTReal );
	oField7.SetWidth(10);
	if( poLayer->CreateField( &oField7 ) != OGRERR_NONE )
	{
		printf( "Creating 'Data' field failed.\n" );
		exit( 1 );
	}

	int NO = 1;
	poFeature = OGRFeature::CreateFeature( poLayer->GetLayerDefn() );
	//write a WKT file and shp file-----------------------------------------
	for (auto Eg : boost::make_iterator_range(edges(G))) 
	{
		fault = G[Eg].trace;
		float L = (float) G[Eg].length;
		const char *T = G[Eg].BranchType.c_str();
		const char *C = G[Eg].component.c_str();
		
		poFeature->SetField( "No", NO  );
		poFeature->SetField( "Length", L);
		poFeature->SetField( "BranchType", T);
		poFeature->SetField( "Component", C);
		poFeature->SetField( "Centre", G[Eg].Centre);
		poFeature->SetField( "MeanValue", G[Eg].MeanValue);
		poFeature->SetField( "CentreGrad",G[Eg].CentreGrad);
		poFeature->SetField( "CrossGrad", G[Eg].CrossGrad);
		poFeature->SetField( "Paralgrad", G[Eg].ParalGrad);

		line.append("LINESTRING(");
		BOOST_FOREACH(point_type P,fault) 
		{
			line.append(std::to_string(P.x()) + " " + std::to_string(P.y()));
			if (!geometry::equals(P, fault.back()))
				line.append(", ");
		}
		line.append( ")");
		const char* branch = (const char*) line.c_str();
		l.importFromWkt(&branch);
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
void GEO::WriteSHP2(Graph G, const char* Name)
{
	assert (num_vertices(G) != 0 && num_edges(G) != 0);
	vector<int> component(num_vertices(G));
	
	point_type point;
	std::string Point;
	
	GDALAllRegister();
	GDALDataset *poDS;
	GDALDriver *poDriver;
	OGRLayer *poLayer;
	OGRFeature *poFeature;
	OGRPoint PO;
	
	const char *pszDriverName = "ESRI Shapefile";
	char * pszWKT;
	
	GeoRef.exportToWkt( &pszWKT );

	poDriver = GetGDALDriverManager()->GetDriverByName(pszDriverName );
	poDS = poDriver->Create(Name, 0, 0, 0, GDT_Unknown, NULL );
	if( poDS == NULL )
	{
		printf( "Creation of shp-file failed.\n" );
		exit( 1 );
	}

	poLayer = poDS->CreateLayer(Name, &GeoRef, wkbPoint, NULL );
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
	

	int NO = 1;
	poFeature = OGRFeature::CreateFeature( poLayer->GetLayerDefn() );
	
	//write shp file----------------------------------------------------
	for (auto Ve : boost::make_iterator_range(vertices(G))) 
	{
		point = G[Ve].location;
		int de = degree(Ve, G);
		int co = connected_components(G, &component[0]);

		poFeature->SetField( "No", NO);
		poFeature->SetField( "Degree", de);
		poFeature->SetField( "Component", co);
		Point.append("POINT(");
		Point.append(std::to_string(point.x()) + " " + std::to_string(point.y()));
		Point.append( ")");
		const char* p = (const char*) Point.c_str();
		
		PO.importFromWkt(&p);
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




 
//
seg_tree GEO::Seg_Tree(vector<line_type>F)
{
	size_t index = 0;
	std::vector<seg_index> segments;
	
	BOOST_FOREACH(line_type l, F)
	{
		Segment curSeg(l.front(), l.back());
		segments.push_back(std::make_pair(curSeg, index));
		index++;
	}
	
    seg_tree rt(segments.begin(), segments.end());
    vector< vector<seg_index> > closest(segments.size());
    return(rt);
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

//calculate the angle (in radians) between the connecting end points of two line_type's that sharea  common point
//the x_front boolean values are true if the user is checking the front of that line segment, and false if the user wants to know the angle from the end of that line segment
//returns close to 0 if the lines are close to having the same angle, or larger angles if the line_types are not in line with each other
double CalculateAngleDifference(line_type &a, bool a_front, line_type &b, bool b_front)
{
	//get the two relevant end points of the oine segments, as it is the first/last two points that define the ending segment of the line_type
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

//sometimes the data has faults that are split into different entries in the shapefile
//this function checks the vector of input line_types for faults that were split, and merges them back together
//this checks for any and all fault segments that should be merged back with base
bool MergeConnections(std::list<line_type> &faults, line_type &base, double distance_threshold)
{
	//our thresholds for merging fault segments together are:
	//1) Their endpoints are within threshold distance of each other
	//2) Only two fault segments meet at the intersection point of the two candidate fault segments
	//3) The angle between the two candidate fault segments is within angle_threshold (in radians)
	//remember that a single fault could be split into more than two fault segments
	//change: the damage zone distance depends on the length of the fault, so we need to ensure that they are all joined properly
	//therefore, we still want to merge fault segments when there are multiple other intersecting faults
	//and we choose to merge the segment with the smallest angle difference
	const double angle_threshold = 25 * math::constants::pi<double>() / 180; //25 degrees, in radians
	bool changed;
	
	do {
		std::list<line_type>::iterator front_match = faults.end(), back_match = faults.end(), candidate;
		bool front_match_loc = false, back_match_loc = false; //true iff we're matching to the front of the other linestring
		std::vector<std::tuple<std::list<line_type>::iterator, bool>> front_matches, back_matches;
		bool match_loc;
		for (auto comp_it = faults.begin(); comp_it != faults.end(); comp_it++)
		{
			//the line segments could be in any order, so we need to check all four pairs of endpoints
			if (geometry::distance(base.front(), comp_it->front()) <= distance_threshold) front_matches.push_back(make_tuple(comp_it,  true));
			if (geometry::distance(base.front(), comp_it-> back()) <= distance_threshold) front_matches.push_back(make_tuple(comp_it, false));
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
				std::list<line_type>::iterator f1 = std::get<0>(*it1);
				bool front1 = std::get<1>(*it1);
				for (auto it2 = it1; it2 != front_matches.end(); it2++)
				{
					if (it1 == it2) continue;
					std::list<line_type>::iterator f2 = std::get<0>(*it2);
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
		if (front_match != faults.end())
		{
			line_type prepend = *front_match;
			if (front_match_loc) geometry::reverse(prepend);
			geometry::append(prepend, base);
			base = prepend;
			faults.erase(front_match);
			changed = true;
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
				std::list<line_type>::iterator f1 = std::get<0>(*it1);
				bool front1 = std::get<1>(*it1);
				for (auto it2 = it1; it2 != back_matches.end(); it2++)
				{
					if (it1 == it2) continue;
					std::list<line_type>::iterator f2 = std::get<0>(*it2);
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
		if (back_match != faults.end()){
			line_type append = *back_match;
			if (!back_match_loc) geometry::reverse(append);
			geometry::append(base, append);
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
	list<line_type> unmerged_faults; //a list of faults that have not been merged yet. thye get removed from this list once they've been added to merged_faults, either by being merged with another fault segment, or on its own
	unmerged_faults.insert(std::end(unmerged_faults), std::begin(F), std::end(F));
	
	while (!unmerged_faults.empty())
	{
		line_type base = unmerged_faults.front();
		unmerged_faults.pop_front();
		MergeConnections(unmerged_faults, base, dist);
		merged_faults.push_back(base);
	}
	F = merged_faults;
}
 
