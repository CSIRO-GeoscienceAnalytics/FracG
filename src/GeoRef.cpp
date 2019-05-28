#include "GeoRef.h"
#include "geometrie.h"

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
 template<> struct GetGDALDataTypeTraits<unsigned char> {static const GDALDataType datatype = GDT_Byte;};
 template<> struct GetGDALDataTypeTraits<char> {static const GDALDataType datatype = GDT_Byte;};
 template<> struct GetGDALDataTypeTraits<unsigned short> {static const GDALDataType datatype = GDT_UInt16;};
 template<> struct GetGDALDataTypeTraits<short> {static const GDALDataType datatype = GDT_Int16;};
 template<> struct GetGDALDataTypeTraits<unsigned int> {static const GDALDataType datatype = GDT_UInt32;};
 template<> struct GetGDALDataTypeTraits<int> {static const GDALDataType datatype = GDT_Int32;};
 template<> struct GetGDALDataTypeTraits<float> {static const GDALDataType datatype = GDT_Float32;};
 template<> struct GetGDALDataTypeTraits<double> {static const GDALDataType datatype = GDT_Float64;};
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

 void GEO::AssignValues(Graph& G, std::string const& filename)
{
 for (auto vd : boost::make_iterator_range(vertices(G))) 
   G[vd].elevation = getElevation(G[vd].location, filename);
}

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

//convenience function to read a value using a point as an index
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

 //read pixel vaule of DEM at coodinate
 int GEO::getElevation(point_type p, std::string const& filename)
 {
 const char * name = filename.c_str();
 GDALDataset  *poDataset;
 GDALAllRegister();
 int value = 0, X, Y, Data ;

 poDataset = (GDALDataset *) GDALOpen( name, GA_ReadOnly );
	if( poDataset == NULL )
		std::cout<<"no file"<< std::endl;

 GDALRasterBand *band = poDataset -> GetRasterBand(1);       

 X = (p.x() - GeoTransform[0]) / GeoTransform[1]; 
 Y = abs((p.y() - GeoTransform[3]) / GeoTransform[5]);

 Data = band->RasterIO(GF_Read,X,Y,1,1,&value,1,1,band->GetRasterDataType(),1,1,nullptr);

 if (Data != 0)
	cout << "ERROR: cannot read elevation data!" << endl;
 
 GDALClose( poDataset );
 return (value);
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
		if( poDataset == NULL )
			std::cout<<"no file"<< std::endl;
			

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
		while (! cpp_file.eof() )
		{
			string line;
			getline(cpp_file, line);
			boost::trim(line);
			if (! line.empty() && ! boost::starts_with(line, "#"))
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
	}
 else 
 {
 cout<<"ERROR: Failure reading file: " << filename << endl;
 exit(EXIT_FAILURE);
 }
 cout << "read " << lineString.size() << " faults from txt" << endl;
 }
 
 //read vector data from shp file---------------------------------------
 void GEO::read_shp(std::string const& filename, std::vector<line_type>& lineString)
 {
	const char * name = filename.c_str();
	const char * layer;
	line_type Line;
	point_type P;
	string f = filename.substr(filename.find("/")+1); 
	layer =  f.substr(0, f.find(".")).c_str();

	GDALAllRegister();
    GDALDataset *poDS = static_cast<GDALDataset*>
    (
        GDALOpenEx( name, GDAL_OF_VECTOR, NULL, NULL, NULL )
    );
    if( poDS == NULL )
    {
        printf( " Open failed.\n" );
        exit( 1 );
    }   
 
	cout << "f = " << f << ", layer = " << layer <<endl;
 
    OGRLayer  *poLayer = poDS->GetLayerByName( layer );
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

 poLayer = poDS->CreateLayer( "Branches", &GeoRef, wkbLineString, NULL );
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
	
		OGRFieldDefn oField3( "Offset", OFTReal );
 oField2.SetWidth(10);
 oField0.SetPrecision(5);
	if( poLayer->CreateField( &oField3 ) != OGRERR_NONE )
	{
		printf( "Creating 'Offset' field failed.\n" );
		exit( 1 );
	}
	
			OGRFieldDefn oField4( "Data", OFTReal );
 oField2.SetWidth(10);
	if( poLayer->CreateField( &oField4 ) != OGRERR_NONE )
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
		poFeature->SetField( "Offset", G[Eg].offset);
		poFeature->SetField( "Data", G[Eg].data);

		line.append("LINESTRING(");
			BOOST_FOREACH(point_type P,fault) 
			{
				line.append(std::to_string(P.x()) + " " + std::to_string(P.y()));
				if (!geometry::equals(P, fault.back()))
					line.append(", ");
			}
			line.append( ")");
			char* branch = (char*) line.c_str();
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



void GEO::CorrectNetwork(vector<line_type>&F, double dist)
{
 GEOMETRIE geom;
 seg_tree SG = Seg_Tree(F);
 vector<seg_index> front, back;

int count = 0;
 BOOST_FOREACH(line_type f, F)
 {
	front.clear();
	back.clear();
	
	box g ;
	BUFFER y = geom.DefinePointBuffer(f.front(), dist);
	geometry::convert(g, y);

	//SG.query(geometry::index::nearest(f.front(),1), std::back_inserter(front));   
		//if (geometry::disjoint(front.at(0), geom.DefinePointBuffer(f.back(), dist));
               

//     cout << F.at(gg.second).front().x() << endl;
            
    if (front.size() > 0)
		count++;
		
    if (back.size() > 0)
		count++;
      
            

 }
 cout << "detected " << count << " bad vertices" << endl;
}
 
