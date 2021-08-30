#include "../include/GeoRef.h"
#include "../include/graph.h"
#include "../include/geometrie.h"
#include "../include/stats.h"
#include "../include/util.h"

namespace FracG
{
	namespace fs = boost::filesystem;

	const std::string raster_subdir="raster";

	//compare two well known text strings, to check if they are equivalent
	void CheckReferenceSystem(std::string wkt1, std::string wkt2)
	{
		OGRSpatialReference oSRS1, oSRS2;
		oSRS1.importFromWkt(wkt1.c_str());
		oSRS2.importFromWkt(wkt2.c_str());
		int err = oSRS1.IsSame(&oSRS2);
		if (err != 1)
		{
			std::cout << "differnet reference systems!" << std::endl;
			exit (EXIT_FAILURE);
		}
	}

	// create and open filestream in folder "statistics"
	std::string CreateDir(VECTOR &input_file, std::initializer_list<std::string> folders)
	{
		std::string folder_name(input_file.folder);
		while (folder_name.back() == '/' ) folder_name.pop_back(); //gdal/ESRI shapefiles don't like double slashes at the start
		while (folder_name.back() == '\\') folder_name.pop_back();
		if (folder_name.size() <= 0) folder_name = ".";
		for (auto it = folders.begin(); it != folders.end(); it++) folder_name += "/"+*it;
		const char* folder_cstr = folder_name.c_str();
		boost::filesystem::create_directories(folder_name);
		return folder_name;
	}

	//read vector data from shp file---------------------------------------
	void ReadShapefile(std::string const& filename, VECTOR& data)
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
				std::cout << "ERROR: vector data without spatial reference /not a projected reference system" << std::endl;
				exit(EXIT_FAILURE);
			}
			if ( pOrigSrs->IsProjected() )
			{
				 char *ref_well_known_text;
				 pOrigSrs->exportToWkt( & ref_well_known_text);
				 data.refWKT = ref_well_known_text;
			}
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
				std::cout <<"MultiLineString detected, splitting into individual linestrings" << std::endl;
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
						std::cout <<"!!! DETECTED CIRCULAR LINE AT " << std::setprecision(10) << Line.front().x() << " " << Line.front().y() << std::endl;
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
					std::cout <<"!!! DETECTED CIRCULAR LINE AT " << std::setprecision(10) << Line.front().x() << " " << Line.front().y() << std::endl;
				Line.clear();
			}
			OGRFeature::DestroyFeature( poFeature );
		}
		GDALClose( poDS );
		std::cout << "read " << data.data.size() << " faults from shp" << std::endl;
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
	template <typename T>double LineExtractor(line_type L, RASTER<T> raster, double dist)
	{
		//note that raster values are converted into double 
		//ode sno use dist anymore
		polygon_type pl = BoundingBox(raster.transform, 0);
		box AOI;
		BUFFER envelop;
		point_type point;
		int maxX, minX, maxY, minY;
		std::vector<double> D;
		double radius = 1.5*(( abs(raster.transform[1]) + abs(raster.transform[5]))/2) ;

		envelop = DefineLineBuffer(L, radius);
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
					D.push_back((double)GetRasterValue(point, raster.transform, raster.values));
			}
			geometry::clear(point);
		}

		if (D.size() > 0)
			return(accumulate(D.begin(), D.end(), 0.0)/D.size());
		else
			return(0);
	 }

	 template <typename T> double CentreGradient(line_type F, RASTER<T> raster, double dist)
	{
	// gradient across the centre of lineament
		line_type cross;
		point_type point, p1, p2;
		polygon_type pl = BoundingBox(raster.transform, 0);
		double len =  dist * (abs(raster.transform[1]) + abs(raster.transform[5]) / 2);
		
		point = GetCentre(F);

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
			double p_v1 = (double)GetRasterValue(p1, raster.transform, raster.values);
			double p_v2 = (double)GetRasterValue(p2, raster.transform, raster.values);
			return(abs(p_v1 - p_v2)/len);
		}
		else
			return(0);
	}

	template <typename T>
	double CrossGradient(line_type F, RASTER<T> raster, double dist)
	{
	//mean gradient across all segmetns of lineament
	//create a functor that boost can use in for_each_segment
		Perpencicular <geometry::model::referring_segment<point_type>> functor;
		functor.d =  dist * (abs(raster.transform[1]) + abs(raster.transform[5]) / 2);	 //set the length of the crossing lines to twice the mean pixel resolution
		functor  = geometry::for_each_segment(F, functor );

		polygon_type pl = BoundingBox(raster.transform, 0);
		std::vector<double> D;

		BOOST_FOREACH(line_type cross, functor.all)
		{
			if (geometry::within(cross.front(), pl) && geometry::within(cross.back(),pl))
			{
				double v2 = (double) GetRasterValue(cross.front(), raster.transform, raster.values);
				double v1 = (double) GetRasterValue(cross.back(), raster.transform, raster.values);
				D.push_back(abs(v1-v2)/geometry::length(cross));
			}
		}

		if (D.size() != 0)
			return(accumulate(D.begin(), D.end(), 0.0)/D.size() );
		else
			return(0);
	}

	template <typename T>
	double ParallelGradient(line_type F, RASTER<T> raster)
	{
	//slope between lineametn tips
		polygon_type pl = BoundingBox(raster.transform, 0);

		if (geometry::within(F.front(), pl) && geometry::within(F.back(), pl))
		{
			double v2 = (double) GetRasterValue(F.front(), raster.transform, raster.values);
			double v1 = (double) GetRasterValue(F.back(), raster.transform, raster.values);
			return( abs(v1-v2) / geometry::length(F));
		}
		else
			return(0);
	}

	//make a directed graph from an undirected graph
	//assumes that AssignValuesGraph() has already been called on the graph, and an appropriate raster file, to initialise the height values (which are stored in FVertex.data, not FVertex.elevation)
	DGraph MakeDirectedGraph(Graph &g)
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
			DVertex dv(fv.location, fv.Enode, fv.enode_dir, fv.data, v - vstart);
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
			if (!fadded || !badded) std::cout << "Warning: Edge from " << dg[ds].index << " to " << dg[dt].index << " already exists (forward " << !fadded << ", back " << !badded << ")" << std::endl;
			dg[frwd].reverse = back;
			dg[back].reverse = frwd;
		}
		return dg;
	}

	//setting the capacity
	double LengthCapacity(DEdge de)
	{
		//double flow = sqrt(de.full_length / 1e3); 
		double flow = 1;
		return(flow);
	}

	double OrientationCapacity(DEdge de, float sig1)
	{
		//we assume a guassian function for capacity depending on orientation
		double angle = (double)(atan2(de.trace.front().y() - de.trace.back().y(), de.trace.front().x() - de.trace.back().x())) 
								* (180 / math::constants::pi<double>());
		if (angle  < 0) 
			angle  += 180;

		double stddev = 25;
		double mean = sig1 ;
		double flow = 1 / (stddev * sqrt(2*math::constants::pi<double>())) * exp(-0.5 * (std::pow( ((angle - mean) / stddev),2)));
		flow *= 1e3;
		return(flow);
	}

	double LengthOrientationCapacity(DEdge de, float sig1)
	{
		double angle = (double)(atan2(de.trace.front().y() - de.trace.back().y(), de.trace.front().x() - de.trace.back().x())) 
								* (180 / math::constants::pi<double>());
		if (angle  < 0) 
			angle  += 180;

		double stddev = 25;
		double mean = sig1 ;
		double cap1 = sqrt(de.full_length / 1e3); 
		double flow = cap1 * (1/ (stddev * sqrt(2*math::constants::pi<double>())) * exp(-0.5 * (std::pow( ((angle - mean) / stddev),2))));
		flow *= 1e5;
		return(flow);
	}

	//the capacity values for the maximum flow algorithm
	void SetupMaximumFlow(DGraph &dg, std::string cap_type, double gradient_angle)
	{
		float sig1;
		enum {ERROR, l, lo, o};
		static std::map<std::string, int> c_type;
		if (c_type.empty() )
		{
			c_type["l"] = l;			//capacity depending on length
			c_type["lo"] = lo;		//capacity depending on length and orientation
			c_type["o"] = o;		//capacity depending on orientation
		}
		int type = c_type[cap_type];

		std::vector<double> all_length;
		switch (type) 
		{
			case l:
			{
				for (auto Eg : make_iterator_range(edges(dg)))	
					all_length.push_back(geometry::length(dg[Eg].trace));
							
			} break;

			case o:
			{
				
			} break;
				
			case lo:
			{
					
			} break;


			default: 
			{
				std::cout << "ERROR: Wrong capacity type defined" << std::endl;
				exit(EXIT_FAILURE);
			} break;
				
		}
		

		std::vector<double> capacities; //need to check for maximum value
		DGraph::edge_iterator e, estart, eend;
		boost::tie(estart, eend) = boost::edges(dg);
		for (e = estart; e != eend; e++)
		{
			double flow;
			switch (type) 
			{
				case l:
				{
					flow = LengthCapacity(dg[*e]);
				} break;

				case o:
				{
					flow =  OrientationCapacity(dg[*e], 122);
				} break;
				
				case lo:
				{
					flow = LengthOrientationCapacity(dg[*e], 111);
				} break;


				default: 
				{
					std::cout << "ERROR: Wrong capacity type defined" << std::endl;
					exit(EXIT_FAILURE);
				} break;
			}
			vertex_type source_vertex = boost::source(*e, dg);
            vertex_type target_vertex = boost::target(*e, dg);
            
			//std::cout << flow << " " << ( (dg[source_vertex].data - dg[target_vertex].data)) / (geometry::length(dg[*e].trace)/1e3) << std::endl;
			
            flow *=  ( (dg[source_vertex].data - dg[target_vertex].data)) / (geometry::length(dg[*e].trace)/1e3);
				
			//std::cout << " " << flow << std::endl;
			
			           		
			dg[*e].capacity = flow;
			dg[*e].residual_capacity = flow;
			capacities.push_back(flow);
		}
		
		double max_cap = *std::max_element(capacities.begin(), capacities.end());
	
		for (e = estart; e != eend; e++)
		{
			dg[*e].capacity /= max_cap;
			dg[*e].residual_capacity /= max_cap;
		}
	}

	//calculate the maximum flow from source s to sink/target t
	//sets up the parameters needed for boost's max flow algorithm, and executes it
	double CalculateMaximumFlow(DGraph &dg, dvertex_type s, dvertex_type t)
	{
		auto capacity_map = boost::get(&DEdge::capacity, dg);
		auto res_cap_map  = boost::get(&DEdge::residual_capacity, dg);
		auto rev_map      = boost::get(&DEdge::reverse, dg);
		auto pred_map     = boost::get(&DVertex::predecessor, dg);
		auto colour_map   = boost::get(&DVertex::colour, dg);
		auto distance_map = boost::get(&DVertex::distance, dg);//not using a distance map currently. this should be distance to the target vertex, and improves the speed of the algorithm
		auto index_map    = boost::get(&DVertex::index, dg); 

		const double flow = boost::boykov_kolmogorov_max_flow(dg, capacity_map, res_cap_map, rev_map, pred_map, colour_map, distance_map, index_map, s, t);
		return flow;
	}

	//perform the maximum flow from the source location to the target location
	double MaximumFlow(DGraph &dg, point_type source, point_type target)
	{
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
			std::cout << "Swapping source and target to find maximum flow from higher location to lower location" << std::endl;
		}

		for (v = vstart; v < vend; v++)
		{
			//I'm not sure if this should be distance to the target, or if I should calculate the shortest paths from the target, or perhaps I set it to zero
			dg[*v].distance = geometry::distance(dg[*v].location, dg[t].location); //set distance map
		}

		//sanity check, make sure we're not calculating the max flow from a vertex to itself
		if (s == t){
			std::cout << "Warning: The source and target points for the maximum flow (" << s << ") are the same" << std::endl;
			return std::numeric_limits<double>::infinity();
		}

		std::cout << "Calculating maximum flow from vertex " << s <<  " to vertex " << t << std::endl;
		//clear capacities for verticies that are above the source
		const double start_data = dg[s].data;
		std::vector<std::pair<dedge_type, double>> saved_capacities;
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

		double max_flow = CalculateMaximumFlow(dg, s, t);
		return max_flow;
	}

	polygon_type BoundingBox(double transform[8], double buffer)
	{
		long double Xmax, Xmin, Ymax, Ymin;
		polygon_type pl;
		box bx;

		//create new box of raster
		Xmax = (transform[0] + transform[1] * transform[6]) - buffer;   // west
		Xmin = transform[0] + buffer; 								  // east
		Ymax = transform[3] + buffer; 								 // north
		Ymin = (transform[3] + transform[5] * transform[7]) - buffer ;	// south

		bx.min_corner().set<0>( Xmin );
		bx.min_corner().set<1>( Ymin );
		bx.max_corner().set<0>( Xmax );
		bx.max_corner().set<1>( Ymax );
		geometry::convert(bx, pl);
		return(pl);
	}

	VECTOR ReadVector(std::string in_filename, std::string out_directory)
	{
		VECTOR lineaments;
		fs::path in_file(in_filename);
		std::string ext, name;
		ext = in_file.extension().string();

		if (ext == ".txt" ||  ext == ".shp") //we don't read from text files
			std::cout << "Reading fault data from a " << ext << "-file." << std::endl;
		else
		{
			std::cout << "ERROR: Wrong vector format provided  (" << ext << ")" << std::endl;
			exit (EXIT_FAILURE);
		}
		lineaments.folder = in_file.parent_path().string();
		lineaments.name = in_file.stem().string();
		lineaments.out_folder = (out_directory == "") ? lineaments.folder : out_directory;
		lineaments.in_path = in_file.parent_path();
		lineaments.out_path = fs::path(lineaments.out_folder);

		//read information from the vector file
		if (ext == ".shp")
			ReadShapefile(in_filename, lineaments);
		else
			std::cerr << "ERROR: no shape file definend" << std::endl;
			
		if (!lineaments.data.size() > 2)
			std::cerr << "ERROR: input data contains less than two lineaments! " << std::endl;
		
		return(lineaments);
	}

	//assign raster values (primarily elevation, from a digital elevation map (DEM)) to the elements of a graph
	//same effect as AssignValues() above, but more efficient
	//read the entire file at once, rather than one at a time. gives faster IO. This assumes that we can fit the entire file in memory at once.
	template <typename T>
	void AssignValuesGraph(Graph& G, RASTER<T> raster, double dist)
	{
		line_type curEdge;
		point_type centre;
		std::cout << " Assigning raster data to graph." << std::endl;

		for (auto Ve : boost::make_iterator_range(vertices(G)))
		{
			G[Ve].data = GetRasterValue(G[Ve].location, raster.transform, raster.values);
			if (std::isnan(G[Ve].data)) std::cout <<"Error: vertex " << Ve << " read a value of nan" << std::endl;
		}

		for (auto Eg : boost::make_iterator_range(edges(G))) 
		{
			curEdge = G[Eg].trace;
			geometry::centroid(G[Eg].trace, centre);

			G[Eg].Centre     = GetRasterValue(centre, raster.transform, raster.values);
			G[Eg].MeanValue  = LineExtractor<T>(G[Eg].trace, raster, dist) ;
			G[Eg].CentreGrad = CentreGradient(G[Eg].trace, raster, dist);
			G[Eg].CrossGrad  = CrossGradient(G[Eg].trace, raster, dist);
			G[Eg].ParalGrad  = ParallelGradient(G[Eg].trace, raster) ;
		}
		std::cout << " done" << std::endl;
	}
	//make these instatiations so that other cpp files can find them
	template void AssignValuesGraph<int>(Graph& G, RASTER<int> raster, double dist);
	template void AssignValuesGraph<float>(Graph& G, RASTER<float> raster, double dist);
	template void AssignValuesGraph<double>(Graph& G, RASTER<double> raster, double dist);



	template <typename T>
	RASTER <T>ReadRaster(std::string in_filename, std::string &comparison_refWKT)
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
		
		std::cout << "Reading raster " << in_filename;
		T** raster_data;
		const char * name = in_filename.c_str();
		RASTER <T>raster;
		double adfGeoTransform[6];

		GDALAllRegister();
		GDALDataset  *poDataset;
		poDataset = (GDALDataset *) GDALOpen( name, GA_ReadOnly );
		if( poDataset == NULL )
		{
			std::cout << "\n ERROR: cannot open raster file " << std::endl;
			exit(EXIT_FAILURE);
		}
		else
		{
				raster.name = in_filename;

				if( poDataset->GetProjectionRef()  != NULL )
					raster.refWKT = poDataset->GetProjectionRef() ;

					//check consitency of refernce system with shp file
					CheckReferenceSystem(raster.refWKT, comparison_refWKT);

				if( poDataset->GetGeoTransform( adfGeoTransform ) == CE_None )
					memcpy ( &raster.transform, &adfGeoTransform, sizeof(adfGeoTransform) );

				GDALRasterBand *band = poDataset -> GetRasterBand(1); 
				double no_data =  band->GetNoDataValue();

				if( poDataset == NULL )
					std::cout<<"no file"<< std::endl;

				raster.transform[6] = poDataset->GetRasterXSize();
				raster.transform[7] = poDataset->GetRasterYSize();
				std::cout << " of size: "<<  poDataset->GetRasterXSize() << " x " << poDataset->GetRasterYSize()<< std::endl;

				raster.values = (T**)RasterConvert(poDataset->GetRasterXSize(), poDataset->GetRasterYSize(), raster.values);
				for(int row = 0; row < poDataset->GetRasterXSize(); row++) 
				{
					int err = band->RasterIO(GF_Read, row, 0, 1, poDataset->GetRasterYSize(), &(raster.values[row][0]), 1, poDataset->GetRasterYSize(), 
								GetGDALDataType<T>(), 0, 0, nullptr);  //read column of raster (switch datatype depending on raster's data type)
					if (err != 0)
					{
						std::cerr << " ERROR reading from raster" << std::endl;
						exit (EXIT_FAILURE);
					}
				}
//set no data values to 0-----------------------------------------------
				int noData_count = 0;
				for (int i = 0; i < poDataset->GetRasterXSize(); i++)
				{
					for (int j = 0; j < poDataset->GetRasterYSize(); j++)
					{
						if (raster.values[i][j] == no_data)
						{
							raster.values[i][j] = 0;
							noData_count++;
						}
					}
				}
				std:: cout << "Found " << noData_count << " NoData cells in raster (" << no_data << ")" <<std::endl;
			GDALClose( poDataset );
		}
		return(raster);
	}
	
	
	

	//TODO: these are ints, and not descriptive. split seems to be the distance threshold for splitting the graph into fault segments, and spurs seems to be the distance threshold for removing spurs from the graph
	template<typename T> 
	graph_map<> RasterGraph(VECTOR lines, int split, int spurs, RASTER<T> raster, double distance_threshold)
	{
		std::cout << "Building graph for raster " << raster.name << std::endl;
		graph_map<> map = ReadVEC4raster(raster.transform, lines, distance_threshold);
		graph_map<> split_map = SplitFaults(map, split); //split the faults in the graph into fault segments, according to the intersections of the faults
		RemoveSpurs(split_map, spurs); //remove any spurs from the graph network
		std::cout << "done" << std::endl;
		return split_map;
	}

	//write raster-augmented shapefile
	//TODO: We should separate the calculates and writing into different functions
	//This doesn't use "dist"
	template<typename T>
	void WriteSHP_r(VECTOR lines, double dist, RASTER<T> raster, std::string save_name)
	{
		GDALAllRegister();
		FracG::CreateDir(save_name);
		const char *pszDriverName = "ESRI Shapefile";
		GDALDriver *poDriver;
		GDALDataset *poDS;
		OGRLayer *poLayer;
		OGRSpatialReference oSRS;
		oSRS.importFromWkt(lines.refWKT.c_str());

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
			
			std::cout <<"created shp file" <<std::endl;

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

			point_type cent = GetCentre(line);
			double c_val   = (double) GetRasterValue<T>(cent, raster.transform, raster.values);
			double mean_v  = (double) LineExtractor<T>(line, raster, dist);
			double grad_c  = (double) CentreGradient<T>(line, raster, dist);
			double grad_cc = (double) CrossGradient<T>(line, raster, dist);
			double grad_p  = (double) ParallelGradient<T>(line, raster);

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
	void WriteTXT(VECTOR lines, double dist, RASTER<T> raster, std::string tsv_filename)
	{
		std::ofstream txtF = FracG::CreateFileStream(tsv_filename);
		
		int x_size   = raster.transform[6];
		int y_size   = raster.transform[7];

		double ul_x  = raster.transform[0];
		double ul_y  = raster.transform[3];
		double x_res = raster.transform[1];
		double y_res = raster.transform[5];
		double pdfMin = 0, pdfMax = 0, pdfMean = 0, pdfStdDev = 0;

		//get raster statistics-------------------------------------------------
		GDALAllRegister();
		GDALDataset  *poDataset;
		poDataset = (GDALDataset *) GDALOpen(raster.name.c_str(), GA_ReadOnly );
		if( poDataset == NULL )
		{
			std::cout << "\n ERROR: cannot open raster file " << std::endl;
			exit(EXIT_FAILURE);
		}
		else
		{
				GDALRasterBand *band = poDataset -> GetRasterBand(1);  
				 if (band->ComputeStatistics(false, &pdfMin, &pdfMax, &pdfMean,  &pdfStdDev, NULL, NULL))
						std::cout << "WARNING: cannot compute raster statistics" << std::endl;
		}
		GDALClose( poDataset );
		txtF	<< "Statistics " << raster.name << "\n" 
				<< "Min \t" 	<< pdfMin << "\n" 
				<< "Max \t" 	<< pdfMax << "\n" 
				<< "Mean \t" 	<< pdfMean << "\n" 
				<< "StdDev \t" 	<< pdfStdDev << std::endl;

		txtF << "Pixel values \t" << dist << std::endl;
		txtF << "No \t Parallel Profile \t Cross Profile" << std::endl;

		int nb = 0;
		BOOST_FOREACH(line_type l, lines.data)
		{
			txtF << nb << "\t";
			
			//exract values along the line
			BOOST_FOREACH(point_type p, l)
			{
			double V = (double)GetRasterValue(p, raster.transform, raster.values);
				if (!geometry::equals(p, l.back()) )
				{
					if(!std::isnan(V) )
						txtF << V << "," ;
				}
				else
					if(!std::isnan(V) )
						txtF << V << "\t" ;
			}

			//exract values along the centre profile
			line_type cross, m_seg;
			point_type centre, p1, p2;
			double len = dist * ((abs(raster.transform[1]) + abs(raster.transform[5]))/2);
			
			double rx = l.back().x() - l.front().x();
			double ry = l.back().y() - l.front().y();
			double le = sqrt(rx*rx + ry*ry);
			
			centre = FracG::GetCentre(l);
			
			p1.set<0>((centre.x() + ( ry/le) * len ));
			p1.set<1>((centre.y() + (-rx/le) * len ));
			p2.set<0>((centre.x() + (-ry/le) * len ));
			p2.set<1>((centre.y() + ( rx/le) * len ));
			geometry::append(cross, p1);
			geometry::append(cross, centre);
			geometry::append(cross, p2);
			
			//create point along line
			double length = geometry::length(cross);
			double raster_res = (abs(raster.transform[1]) + abs(raster.transform[5]))/2;
				
			BOOST_FOREACH(point_type p, cross)
			{
				double V = (double)GetRasterValue(p, raster.transform, raster.values);
				//std::cout << std::setprecision(15) << p.x() << " " << p.y() << std::endl;
				if (!geometry::equals(p, cross.back()) )
				{
					if(!std::isnan(V) )
						txtF << V << "," ;
				}
				else
					if(!std::isnan(V) )
						txtF << V ;
			}
			txtF << "\n";
			nb++;
		}
	}
	template void WriteTXT<int>(VECTOR lines, double dist, RASTER<int> raster, std::string tsv_filename);
	template void WriteTXT<float>(VECTOR lines, double dist, RASTER<float> raster, std::string tsv_filename);
	template void WriteTXT<double>(VECTOR lines, double dist, RASTER<double> raster, std::string tsv_filename);
	
	template<typename T>
	void AnalyseRaster(VECTOR lines, double dist, RASTER<T> raster)
	{
		//two things are happening here:
		//1.) building a shape file 
		//2.) creating a txt file containing profiles and cross profiles
		std::string filename = raster.name;
		raster = ReadRaster<T>(raster.name, lines.refWKT);
		std::string raster_shp_filename = FracG::AddPrefixSuffixSubdirs(lines.out_path, {"raster_shp"}, "raster_augmented_shapefile", ".shp", true);
		WriteSHP_r<T>(lines, dist, raster, raster_shp_filename);
		std::string raster_tsv_filename = FracG::AddPrefixSuffixSubdirs(lines.out_path, {raster_subdir}, "raster_parallel_cross_profiles", ".csv", true);
		WriteTXT(lines, dist, raster, raster_tsv_filename);
	 }
	template void AnalyseRaster<int>(VECTOR lines, double dist, RASTER<int> raster);
	template void AnalyseRaster<float>(VECTOR lines, double dist, RASTER<float> raster);
	template void AnalyseRaster<double>(VECTOR lines, double dist, RASTER<double> raster);

	Graph BuildRasterGraph(VECTOR lines, double split, double spur, double map_distance_threshold, double raster_stats_dist, AngleDistribution angle_dist, std::string raster_filename, bool skip_betweenness_centrality)
	{
		Graph raster_graph;
		if (!raster_filename.empty())
		{
			std::cout << "Generating graph linked to raster " << raster_filename << std::endl;

			//switch case based on raster type
			const char * name = raster_filename.c_str();

			//build a map to use switch case for differnt data types
			enum {ERROR, Byte, UInt16, Int16, UInt32, Int32, Float32, Float64 };
			static std::map<std::string, int> d_type;
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

			std::cout << "Raster datatype is " ;
			int type = d_type[datatype];
			 switch (type) 
			 {
				case Byte:
				{
					std::cout	<< "byte (will not perform analysis)" << std::endl;
					RASTER<char> R;
					R = ReadRaster<char>(raster_filename, lines.refWKT);
					graph_map<> map = RasterGraph(lines, split, spur, R, map_distance_threshold);
					raster_graph = map.GetGraph();
					WriteGraph(raster_graph, lines, "raster", true, skip_betweenness_centrality);
				} break;

				case UInt16:
				{
					std::cout << "uint16" << std::endl;
					RASTER<int> R;
					R = ReadRaster<int>(raster_filename, lines.refWKT);
					graph_map<> map = RasterGraph(lines, split, spur, R, map_distance_threshold);
					raster_graph = map.GetGraph();
					AssignValuesGraph<int>(raster_graph, R, raster_stats_dist);
					GraphAnalysis(raster_graph, lines, 10, angle_dist, raster_filename);
					WriteGraph(raster_graph, lines, "raster", true, skip_betweenness_centrality);
				} break;

				case Int16:
				{
					std::cout << "int16" << std::endl;
					RASTER<int> R;
					R = ReadRaster<int>(raster_filename, lines.refWKT);
					graph_map<> map = RasterGraph(lines, split, spur, R, map_distance_threshold);
					raster_graph = map.GetGraph();
					AssignValuesGraph<int>(raster_graph, R, raster_stats_dist);
					GraphAnalysis(raster_graph, lines, 10, angle_dist, raster_filename);
					WriteGraph(raster_graph, lines, "raster", true, skip_betweenness_centrality);
				} break;

				case UInt32:
				{
					std::cout << "uint32" << std::endl;
					RASTER<int> R;
					R = ReadRaster<int>(raster_filename, lines.refWKT);
					graph_map<> map = RasterGraph(lines, split, spur, R, map_distance_threshold);
					raster_graph = map.GetGraph();
					AssignValuesGraph<int>(raster_graph, R, raster_stats_dist);
					GraphAnalysis(raster_graph, lines, 10, angle_dist, raster_filename);
					WriteGraph(raster_graph, lines, "raster", true, skip_betweenness_centrality);
				} break;

				case Int32:
				{
					std::cout << "int32" << std::endl;
					RASTER<int> R;
					R = ReadRaster<int>(raster_filename, lines.refWKT);
					graph_map<> map = RasterGraph(lines, split, spur, R, map_distance_threshold);
					raster_graph = map.GetGraph();
					AssignValuesGraph<int>(raster_graph, R, raster_stats_dist);
					GraphAnalysis(raster_graph, lines, 10, angle_dist, raster_filename);
					WriteGraph(raster_graph, lines, "raster", true, skip_betweenness_centrality);
				} break;

				case Float32:
				{
					std::cout << "Float32" << std::endl;
					RASTER<float> R;
					R = ReadRaster<float>(raster_filename, lines.refWKT);
					graph_map<> map = RasterGraph(lines, split, spur, R, map_distance_threshold);
					raster_graph = map.GetGraph();
					AssignValuesGraph<float>(raster_graph, R, raster_stats_dist);
					GraphAnalysis(raster_graph, lines, 10, angle_dist, raster_filename);
					WriteGraph(raster_graph, lines, "raster", true, skip_betweenness_centrality);
				} break;

				case Float64:
				{
					std::cout << "Float64" << std::endl;
					RASTER<double> R;
					R = ReadRaster<double>(raster_filename, lines.refWKT);
					graph_map<> map = RasterGraph(lines, split, spur, R, map_distance_threshold);
					raster_graph = map.GetGraph();
					AssignValuesGraph<double>(raster_graph, R, raster_stats_dist);
					GraphAnalysis(raster_graph, lines, 10, angle_dist, raster_filename);
					WriteGraph(raster_graph, lines, "raster", true, skip_betweenness_centrality);

				} break;

				default: 
				{
					std::cout << " unknown" << std::endl;
					std::cout << "ERROR: Could not determine dataype of raster" << std::endl;
				}
				break;
			 }
			std::cout << " done \n" << std::endl;
		}
		return(raster_graph);
	}

	void ReadPoints(std::string const& filename, VECTOR lines, std::pair<point_type, point_type> &source_target)
	{
		point_type P;
		line_type Line;
		const char* name = filename.c_str();
		std::vector<point_type> points;

		OGRSpatialReference init_ref;
		init_ref.importFromWkt(lines.refWKT.c_str());
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
				std::cout << "ERROR: vector data without spatial reference" << std::endl;
				exit(EXIT_FAILURE);
			}
			if ( pOrigSrs->IsProjected() )
			{
				std::string d1 = to_string(datum);
				if(to_string(datum)!= to_string(datum2))
				{
					std::cout << " ERROR: datum of point data not consitent with initial datum" << std::endl;
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
			std::cout << " ERROR: inconsitent point number in source - target file" << std::endl;
	}

	//read a value from the raster at a point-------------------------------
	template <typename T>
	T GetRasterValue(point_type p, double transform[8], T** values)
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

	//this is a function to write a raster file from the datastructure RASTER
	//mainly for debug purposes but might be useful
	template<typename T>
	void WriteRasterStruct(RASTER<T> raster)
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
			std::cout << "ERROR: Could not load GDAL driver ( GTiff)" << std::endl;
			EXIT_FAILURE;
		}

		papszMetadata = poDriver->GetMetadata();

		poDstDS = poDriver->Create("test.tif", x,y, 1, GDT_Float32, NULL);

		poDstDS->SetGeoTransform( raster.transform );

		poDstDS->SetProjection( raster.refWKT.c_str());

		GDALRasterBand *poBand = poDstDS->GetRasterBand(1);

		float *rowBuff = (float*) CPLMalloc(sizeof(float)*x);
		std::cout << "Writing raster (struc) " << std::endl;
		progress_display* show_progress =  new boost::progress_display(y * x);

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
		std::cout << " done \n" << std::endl;
	}
	template void WriteRasterStruct<int>(RASTER<int> raster);
	template void WriteRasterStruct<float>(RASTER<float> raster);
	template void WriteRasterStruct<double>(RASTER<double> raster);

	 void WriteRASTER(std::vector<std::vector<double>> data, std::string SpatialRef, double adfGeoTransform[6], std::string out_filename)
	{
		GDALDataset *poDstDS;
		GDALDriver *poDriver;
		char *pszSRS_WKT = NULL;
		char **papszMetadata;
		int x = data.size();
		int y = data[0].size();
		const char *pszFormat = "GTiff";
		FracG::CreateDir(out_filename);

		poDriver = GetGDALDriverManager()->GetDriverByName(pszFormat);
		if( poDriver == NULL )
		{
			std::cout << "ERROR: Could not load GDAL driver ( GTiff)" << std::endl;
			EXIT_FAILURE;
		}

		papszMetadata = poDriver->GetMetadata();
		poDstDS = poDriver->Create(out_filename.c_str(), x,y, 1, GDT_Float32, NULL);
		poDstDS->SetGeoTransform( adfGeoTransform );
		poDstDS->SetProjection( SpatialRef.c_str() );

		GDALRasterBand *poBand = poDstDS->GetRasterBand(1);

		float *rowBuff = (float*) CPLMalloc(sizeof(float)*x);
		std::cout << "Writing raster " << std::endl;

		for(int row = 0; row < y; row++) 
		{
			for(int col = 0; col < x; col++) 
			{
				rowBuff[col] = (float) data[col][row];

			}
			int err = poBand->RasterIO(GF_Write, 0, row, x, 1, rowBuff, x, 1, GDT_Float32, 0, 0);
		}
		CPLFree(rowBuff);
		poBand->SetNoDataValue(-256);
		GDALClose( (GDALDatasetH) poDstDS );
	}
	
	  void WriteSHP_lines(std::vector<line_type> lineaments, std::string refWKT, std::string name)
	 {
		GDALAllRegister();
		const char *pszDriverName = "ESRI Shapefile";
		GDALDriver *poDriver;
		GDALDataset *poDS;
		OGRLayer *poLayer;
		OGRSpatialReference oSRS;
		oSRS.importFromWkt(refWKT.c_str());

		name = FracG::AddPrefixSuffix(name, "", ".shp");
		FracG::CreateDir(name);
		const char* Name = name.c_str();

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

		poLayer = poDS->CreateLayer(Name, &oSRS, wkbLineString, NULL );
		if( poLayer == NULL )
		{
			printf( "Layer creation failed.\n" );
			exit( 1 );
		}

		OGRFieldDefn oField( "ID", OFTString );
		oField.SetWidth(10);
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

	 //write vector data to disk
	 //TODO: this can write gauss parameter-associated values, but is only called before those values are calculated
	 void WriteShapefile(VECTOR &lineaments, AngleDistribution &angle_dist, std::string name)
	 {
		 
		GDALAllRegister();
		const char *pszDriverName = "ESRI Shapefile";
		GDALDriver *poDriver;
		GDALDataset *poDS;
		OGRLayer *poLayer;

		name = FracG::AddPrefixSuffix(name, "", ".shp");
		const char* Name = name.c_str();

		OGRSpatialReference oSRS;
		oSRS.importFromWkt(lineaments.refWKT.c_str());

		poDriver = (GDALDriver*) GDALGetDriverByName(pszDriverName);
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

		if (angle_dist.gaussians.size() != 0)
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
			if (angle_dist.gaussians.size() != 0)
				poFeature->SetField( "Set", CheckGaussians(angle_dist, strike));

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
		std::cout << "created shp-file " << name << std::endl << std::endl;
	 }

	//convert graph edges to shp-file---------------------------------------
	void WriteSHP_g_lines(Graph G, std::string refWKT, std::string Name, bool raster)
	{
		GDALAllRegister();
		const char *pszDriverName = "ESRI Shapefile";
		GDALDriver *poDriver;
		GDALDataset *poDS;
		OGRLayer *poLayer;

		OGRSpatialReference oSRS;
		oSRS.importFromWkt(refWKT.c_str());

		poDriver = (GDALDriver*) GDALGetDriverByName(pszDriverName );
		if( poDriver == NULL )
		{
			printf( "%s driver not available.\n", pszDriverName );
			exit( 1 );
		}

		poDriver->SetDescription("graph_edges");

		poDS = poDriver->Create(Name.c_str(), 0, 0, 0, GDT_Unknown, NULL );
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
		oField2.SetWidth(20);
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
		
		if (raster)
		{
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
			}

		int id = 0;
		for (auto Eg : boost::make_iterator_range(edges(G))) 
		{
			OGRLineString l;
			OGRFeature *poFeature;
			line_type line = G[Eg].trace;
			geometry::unique(line);
			float L = (float) G[Eg].length;
			double strike = (double)(atan2(line.front().y() - line.back().y(), line.front().x() - line.back().x())) 
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
		
		if (raster)
		{
		//add raster data to graph shp file-------------------------------------
			poFeature->SetField( "Centre", 		(float) G[Eg].Centre);
			poFeature->SetField( "MeanValue", 	(float) G[Eg].MeanValue);
			poFeature->SetField( "CentreGrad", 	(float) G[Eg].CentreGrad);
			poFeature->SetField( "CrossGrad", 	(float) G[Eg].CrossGrad);
			poFeature->SetField( "PGrad", 		(float) G[Eg].ParalGrad );
		}
			
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

	void WriteSHP_maxFlow(DGraph G, std::string refWKT, std::string name)
	{
		GDALAllRegister();
		const char *pszDriverName = "ESRI Shapefile";
		GDALDriver *poDriver;
		GDALDataset *poDS;
		OGRLayer *poLayer;

		OGRSpatialReference oSRS;
		oSRS.importFromWkt(refWKT.c_str());

		name = FracG::AddPrefixSuffix(name, "", ".shp");
		FracG::CreateDir(name);
		const char* Name = name.c_str();

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

		poLayer = poDS->CreateLayer(Name, &oSRS, wkbLineString, NULL );
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
	void WriteSHP_g_points(Graph G, const char* refWKT, const char* Name, bool raster, bool skip_betweenness_centrality)
	{
		GDALAllRegister();
		point_type point;
		std::string Point;
		GDALDataset *poDS;
		GDALDriver *poDriver;
		OGRLayer *poLayer;
		OGRFeature *poFeature;
		const char* pszDriverName = "ESRI Shapefile";
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
		
		OGRFieldDefn oField3( "Rel_bc", OFTReal );
		oField3.SetWidth(10);
		if( poLayer->CreateField( &oField3 ) != OGRERR_NONE )
		{
			printf( "Creating betweeness centrality field failed.\n" );
			exit( 1 );
		}
		
		OGRFieldDefn oField4( "Abs_bc", OFTInteger );
		oField4.SetWidth(10);
		if( poLayer->CreateField( &oField4 ) != OGRERR_NONE )
		{
			printf( "Creating betweeness centrality field failed.\n" );
			exit( 1 );
		}
		
		if (raster)
		{
			OGRFieldDefn oField5( "data", OFTReal );
			oField5.SetWidth(10);
			if( poLayer->CreateField( &oField5 ) != OGRERR_NONE )
			{
				printf( "Creating data field failed.\n" );
				exit( 1 );
			}
		}

		std::cout << "centrality started" << std::endl;

		//Betweeness centrality----------------------------------------
		boost::property_map<Graph, boost::vertex_index_t>::type  v_index = get(boost::vertex_index, G);
		std::vector< double > vertex_property_vec(boost::num_vertices(G), 0.0);
		iterator_property_map< std::vector< double >::iterator, boost::property_map<Graph, boost::vertex_index_t>::type> vertex_property_map(vertex_property_vec.begin(), v_index);
		float factor = 0;
		if (!skip_betweenness_centrality)
		{
			brandes_betweenness_centrality(G, vertex_property_map);
			factor = 2 / (num_vertices(G)*num_vertices(G) - 3*num_vertices(G) +2);
		}
		
		std::cout << "centrality done" << std::endl;
		
		int NO = 0;
		poFeature = OGRFeature::CreateFeature( poLayer->GetLayerDefn() );
		//write shp file----------------------------------------------------
		for (auto Ve : boost::make_iterator_range(vertices(G))) 
		{
			point = G[Ve].location;
			float da = G[Ve].data;

			poFeature->SetField( "No", (int) NO );
			poFeature->SetField( "Degree", (int) degree(Ve, G));
			poFeature->SetField( "Component", (int) G[Ve].component);
			
			poFeature->SetField( "Abs_bc", (double) vertex_property_map[Ve]);
			if (factor != 0)
				poFeature->SetField( "Rel_bc",  (float) vertex_property_map[Ve]/factor);
			else
				poFeature->SetField( "Rel_bc",  0);
			
			if (raster)
				poFeature->SetField( "data", (double) G[Ve].data);
				
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

	void WriteGraph(Graph G, VECTOR lines, std::string subF, bool raster, bool skip_betweenness_centrality)
	{
		std::cout << "starting writegraph" << std::endl;
		assert (num_vertices(G) != 0 && num_edges(G) != 0);

		const char* reference = lines.refWKT.c_str();

		std::string subdir_name = FracG::AddPrefixSuffixSubdirs(lines.out_path, {"graph_shp", subF}, "", "/");
		FracG::CreateDir(subdir_name);

		std::string n_b =  FracG::AddPrefixSuffix(subdir_name, "graph_branches", ".shp");
		std::string n_v =  FracG::AddPrefixSuffix(subdir_name, "graph_vertices", ".shp");
		const char* Name_b = n_b.c_str();
		const char* Name_v = n_v.c_str();

		WriteSHP_g_lines(G, reference, Name_b, raster);//false
		WriteSHP_g_points(G, reference, Name_v, raster, skip_betweenness_centrality);//false
		std::cout << " done" << std::endl;
	}

// 	//TODO: Why is this a separate function from the above?
// 	void WriteGraph_R(Graph G, VECTOR lines, std::string subF)
// 	{
// 		std::cout << "starting writegraph (raster)" << std::endl;
// 		assert (num_vertices(G) != 0 && num_edges(G) != 0);
// 
// 		const char* reference = lines.refWKT.c_str();
// 
// 		std::string subdir_name = FracG::AddPrefixSuffixSubdirs(lines.out_path, {"graph_shp", subF}, "", "/");
// 		FracG::CreateDir(subdir_name);
// 
// 		std::string n_b =  FracG::AddPrefixSuffix(subdir_name, "graph_branches", ".shp");
// 		std::string n_v =  FracG::AddPrefixSuffix(subdir_name, "graph_vertices", ".shp");
// 		const char* Name_b = n_b.c_str();
// 		const char* Name_v = n_v.c_str();
// 
// 		WriteSHP_g_lines(G, reference, Name_b, true);
// 		WriteSHP_g_points(G, reference, Name_v, true);
// 		std::cout << " done" << std::endl;
// 	}

	//TODO: Give these more indicative names
	//distance from centre to centre of features
	void PointTree(std::vector<p_index> points,  std::vector<p_index>& closest)
	{
		p_tree rt(points.begin(), points.end());

		for ( size_t i = 0 ; i < points.size() ; ++i )
			rt.query(geometry::index::nearest(points[i].first, 1) && geometry::index::satisfies(different_id_p(i)), std::back_inserter(closest) );        
	}

	//distance from centre to centre of nearest larger feature
	void PointTree2(std::vector<pl_index> points, std::vector<pl_index>& closest, double max_len)
	{
	 pl_tree rt(points.begin(), points.end());

	 for ( size_t i = 0 ; i < points.size() ; ++i )
	 {
		if (get<2>(points[i]) >= max_len)
		{
			closest.push_back(std::make_tuple(get<0>(points[i]),i,0));
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
        const double capped_cosine = std::min(dotprod / (la*lb), 1.0); //this value can go over one, I think it is a floating point rounding error
		const double theta = acos(capped_cosine); //calculate angle, and subract from 1pi radians/180 degrees to get the angle difference from being a straight line
		return theta;
	}

	//   typedef the rtree here so both functions can see the same type. not putting it in the header because no other function uses it
	typedef std::list<line_type> unmerged_type; //might want to change this to a set or map
	typedef std::tuple<point_type, unmerged_type::iterator, bool> endpoint_value_type;
	typedef geometry::index::rtree<endpoint_value_type, geometry::index::rstar<16>> endpoint_rtree_type;

	//convenience function to remove endpoints from the endpoint rtree, given an iterator
	void RemoveEndpoints(endpoint_rtree_type &rtree, std::list<line_type>::iterator it)
	{
		long r1 = rtree.remove(std::make_tuple(it->front(), it, true ));
		long r2 = rtree.remove(std::make_tuple(it->back (), it, false));
	}

	std::tuple<unmerged_type::iterator, bool> ChooseBestMatch(unmerged_type::iterator base_it, std::list<std::tuple<unmerged_type::iterator, bool>> matches, bool is_front, double angle_threshold_radians, unmerged_type::iterator unmatched_marker)
    {
        bool this_match_loc, other_match_loc;
        decltype(matches)::iterator base_match_it = matches.insert(matches.end(), std::make_tuple(base_it, is_front));
        decltype(matches)::iterator candidate1, candidate2;
        
        //need at least two elements to match together
        while (matches.size() >= 2)
        {
            double best_angle = std::numeric_limits<double>::infinity();
            candidate1 = matches.end();
            candidate2 = matches.end();
            for (auto it1 = matches.begin(); it1 != matches.end(); it1++)
            {
                unmerged_type::iterator f1 = std::get<0>(*it1);
                bool front1 = std::get<1>(*it1);
                for (auto it2 = it1; it2 != matches.end(); it2++)
                {
                    if (it1 == it2) continue;
                    unmerged_type::iterator f2 = std::get<0>(*it2);
                    bool front2 = std::get<1>(*it2);
                    double angle_diff = abs(CalculateAngleDifference(*f1, front1, *f2, front2));
                    if (angle_diff < best_angle)
                    {
                        candidate1 = it1;
                        candidate2 = it2;
                        best_angle = angle_diff;
                    }
                }
            }
            //no suitable match found
            if (best_angle > angle_threshold_radians || candidate1 == matches.end())
            {
                return std::make_tuple(unmatched_marker, false);
            }
            if (candidate1 == base_match_it || candidate2 == base_match_it)
            {
                //we have a match for this line segment, so return it
                decltype(matches)::iterator return_it = (candidate1 == base_match_it) ? candidate2 : candidate1;
                return *return_it;//std::make_tuple(this_match, this_match_loc)
            }
            else if (candidate1 != matches.end() && candidate2 != matches.end())
            {
                //two other segments are best matches for each other, so remove them from consideration here
                matches.erase(candidate1);
                matches.erase(candidate2);
            }
        }
        return std::make_tuple(unmatched_marker, false);
    }

	//sometimes the data has faults that are split into different entries in the shapefile
	//this function checks the vector of input line_types for faults that were split, and merges them back together
	//this checks for any and all fault segments that should be merged back with base
	bool MergeConnections(unmerged_type &faults, endpoint_rtree_type &endpoints, unmerged_type::iterator base_it, double distance_threshold, double angl)
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
		const double angle_threshold_radians = angl * math::constants::pi<double>() / 180; 

		const double dt = distance_threshold; //convenience name
		line_type &base = *base_it;

		do {
			//a front match is one that matches to the front of the base fault/line segment, a back match is one that matches to the back of the base fault/line segment
			//use faults.end() as a sentinel value to mark that a match hasn't been found
			unmerged_type::iterator front_match = faults.end(), back_match = faults.end();
			bool front_match_loc = false, back_match_loc = false; //true iff we're matching to the front of the other linestring
			std::list<std::tuple<unmerged_type::iterator, bool>> front_matches, back_matches;
			bool match_loc;
			//use boxes, because boost doesn't have a circle geometry object (just a polygon/linestring that approximates a circle)
			box front_box(point_type(base.front().x() - dt, base.front().y() - dt), point_type(base.front().x() + dt, base.front().y() + dt));
			box  back_box(point_type(base.back ().x() - dt, base.back ().y() - dt), point_type(base.back ().x() + dt, base. back().y() + dt));
			std::vector<endpoint_value_type> candidates;
			endpoints.query(boost::geometry::index::intersects(front_box), back_inserter(candidates));
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
            std::tie(front_match, front_match_loc) = ChooseBestMatch(base_it, front_matches, true, angle_threshold_radians, faults.end());
			
			//check for matches to the back of base
            std::tie(back_match, back_match_loc) = ChooseBestMatch(base_it, back_matches, false, angle_threshold_radians, faults.end());

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
				RemoveEndpoints(endpoints, front_match);
				faults.erase(front_match);
				if (front_match_loc) geometry::reverse(prepend);
				geometry::append(prepend, base);
				base = prepend;
				changed = true;
			}

			if (back_match != faults.end()){
				line_type append = *back_match;
				RemoveEndpoints(endpoints, back_match);
				faults.erase(back_match);
				if (!back_match_loc) geometry::reverse(append);
				geometry::append(base, append);
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
	void CorrectNetwork(std::vector<line_type>&F, double dist, double angl_threshold, double dfd_thres)
	{
		typedef std::pair<box, std::vector<line_type>::iterator> box_line; //a bounding box around the linestring and its iterator
		//first we need to remove duplicates. We sort the lines by distance to a reference point (from line centroid)
		//and then remove duplicates
		box AOI = ReturnAOI(F);
		point_type origin(geometry::get<geometry::min_corner, 0>(AOI), geometry::get<geometry::min_corner, 1>(AOI));
		bool found_duplicate = true;

			do{
				int size = F.size();
				std::sort(F.begin(), F.end(), LineSorter_orig(origin));
				std::vector<line_type>::iterator new_lines_3;
				new_lines_3 = std::unique(F.begin(), F.end(), LineCompare()); 
				F.erase(new_lines_3, F.end());

				if (size == F.size())
					found_duplicate = false;

			}	while (found_duplicate);
		std::cout << F.size() << " lines remaining after removing duplicates \n" << std::endl;
	//----------------------------------------------------------------------

		std::vector<line_type> merged_faults; //the new list of merged faults, to be built by this function
		unmerged_type unmerged_faults; //a list of faults that have not been merged yet. the get removed from this list once they've been added to merged_faults, either by being merged with another fault segment, or on its own
		unmerged_faults.insert(std::end(unmerged_faults), std::begin(F), std::end(F));

		endpoint_rtree_type endpoint_tree;
		for (auto it = unmerged_faults.begin(); it != unmerged_faults.end(); it++)
		{
			endpoint_tree.insert(std::make_tuple(it->front(), it, true ));
			endpoint_tree.insert(std::make_tuple(it->back (), it, false));
		}

		std::cout << "merging split faults - starting" << std::endl;
		while (!unmerged_faults.empty())
		{
			RemoveEndpoints(endpoint_tree, unmerged_faults.begin());
			MergeConnections(unmerged_faults, endpoint_tree, unmerged_faults.begin(), dist, angl_threshold);
			merged_faults.push_back(unmerged_faults.front());//base
			unmerged_faults.pop_front(); //I don't know why this doesn't also return the value being pop'd
		}
		std::cout <<"merging split faults - finished" << std::endl;
		F = merged_faults;
		std::cout << F.size() << " lines remaining after merging" << std::endl;
	}

	void GetSourceTarget(const char* Name, point_type &Source, point_type &Target)
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
					std::cout << "ERROR: vector data without spatial reference /not a projected reference system" << std::endl;
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
					std::cout << "more than two points in shape file" << std::endl;
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
			std::cout << std::setprecision(10) << "Source: " << Source.x() << ", " << Source.y() << "\n"
				 << "Target: " << Target.x() << ", " << Target.y() << std::endl;
		}


	void MutiplyRasterbyCoefficient(const char* input, double coeff)
	{
		CPLErr eErr;
		GDALAllRegister();
		GDALDataset  *poDataset;
		poDataset = (GDALDataset *) GDALOpen( input, GA_Update );
		CPLAssert( eErr == CE_None) ;
		
		int X  = GDALGetRasterXSize( poDataset );
		int Y = GDALGetRasterYSize( poDataset );
		GDALRasterBand* poBand = poDataset ->GetRasterBand(1);
		GDALDataType gdalType = poBand->GetRasterDataType();
		
		for (int x =0; x < X; x++)
		{
			for (int y = 0; y < Y; y++)
			{
			float* value = (float *) CPLMalloc(sizeof(float));
			eErr = poBand->RasterIO( GF_Read, x, y, 
								1, 1,
								value, 
								1, 
								1, 
								gdalType ,
								0, 0 );

			
			value[0] = (float) coeff * value[0];

			eErr = poBand->RasterIO( GF_Write, x, y, 
								1, 1,
								value, 
								1, 
								1, 
								gdalType ,
								0, 0 );
								
			}
		}
		GDALClose( poDataset ); 
	}


	void gdal_resample( std::string srcfname, std::string dstfname)
	{
		GDALAllRegister();
		GDALDatasetH srcDataset;
		GDALDatasetH dstDataset;
		CPLErr eErr;
		
		int dstnrows, dstncols;
		int outnrows = 0;
		int outncols = 0;
		
		const char* srcProjection;
		const char* dstProjection;
		double srcGeoTransform[6];
		double dstGeotransform[6];

		const char *input = srcfname.c_str();
		const char *output = dstfname.c_str();

		srcDataset = GDALOpen( input, GA_ReadOnly );
		CPLAssert( srcDataset != NULL );
		 
		srcProjection = GDALGetProjectionRef( srcDataset ) ;
		CPLAssert( srcProjection != NULL && strlen( srcProjection ) > 0) ;
		GDALDataType sourceDatatype;
		sourceDatatype = GDALGetRasterDataType( GDALGetRasterBand(srcDataset,1) );

		dstProjection = GDALGetProjectionRef( srcDataset ); // set destination projection to source projection
		GDALGetGeoTransform(srcDataset, srcGeoTransform
);		//get the intitial geotransfrom
		//define new geotransform with a 10th of the intitial pixel size
		dstGeotransform[0] = srcGeoTransform[0];
		dstGeotransform[1] = srcGeoTransform[1] / 10;
		dstGeotransform[2] = srcGeoTransform[2];
		dstGeotransform[3] = srcGeoTransform[3];
		dstGeotransform[4] = srcGeoTransform[4];
		dstGeotransform[5] = srcGeoTransform[5] / 10;

		dstncols  = (int)( GDALGetRasterXSize( srcDataset ) * srcGeoTransform[1] ) / dstGeotransform[1];
		dstnrows = (int)( GDALGetRasterYSize( srcDataset ) * srcGeoTransform[5] ) / dstGeotransform[5];

		std::cout << "resampling from " << GDALGetRasterXSize( srcDataset ) << " x " << GDALGetRasterYSize( srcDataset ) << " to " <<
		dstncols << " x " << dstnrows << std::endl;

	  void *handleTransformArg;
	  handleTransformArg = GDALCreateGenImgProjTransformer( srcDataset, srcProjection,
	  NULL, dstProjection, FALSE, 0, 1);
	  CPLAssert( handleTransformArg != NULL);

	  eErr = GDALSuggestedWarpOutput( srcDataset,
		GDALGenImgProjTransform,
		handleTransformArg,
		srcGeoTransform,
		&outncols ,
		&outnrows );
	  CPLAssert( eErr == CE_None);
	  GDALDestroyGenImgProjTransformer( handleTransformArg );
	 
	  GDALDriverH outHandleDriver;
	  GDALDatasetH outDataset;
	  outHandleDriver = GDALGetDriverByName("GTiff");
	  outDataset = GDALCreate( outHandleDriver ,
	  output,
	  dstncols, dstnrows , 1 ,
	  sourceDatatype, NULL);
	  
	  GDALSetProjection( outDataset , dstProjection) ;
	  GDALSetGeoTransform( outDataset, dstGeotransform) ;
	  GDALWarpOptions *options;

	  eErr = GDALReprojectImage( srcDataset ,
		srcProjection,
		outDataset ,
		dstProjection ,
	  GRA_CubicSpline, 0.0, 0.0, NULL,NULL,NULL);
	  CPLAssert( eErr == CE_None) ;
	  GDALClose(outDataset); 
	  
	  //now rescale to initial max values--------------------------------
	  	GDALAllRegister();
		GDALDataset  *poDataset, *poDataset2;
		poDataset = (GDALDataset *) GDALOpen(input, GA_ReadOnly );
		poDataset2 = (GDALDataset *) GDALOpen(output, GA_ReadOnly );
		
		double pdfMin, pdfMax, pdfMean, pdfStdDev, pdfMin2, pdfMax2, pdfMean2, pdfStdDev2;
		if( poDataset == NULL || poDataset2 == NULL)
		{
			std::cout << "\n ERROR: cannot open raster file " << std::endl;
			exit(EXIT_FAILURE);
		}
		else
		{
				GDALRasterBand *band = poDataset -> GetRasterBand(1);  
				GDALRasterBand *band2 = poDataset2 -> GetRasterBand(1);  
				 if (band->ComputeStatistics(false, &pdfMin, &pdfMax, &pdfMean,  &pdfStdDev, NULL, NULL))
						std::cout << "WARNING: cannot compute raster statistics" << std::endl;
				if (band2->ComputeStatistics(false, &pdfMin2, &pdfMax2, &pdfMean2,  &pdfStdDev2, NULL, NULL))
						std::cout << "WARNING: cannot compute raster statistics" << std::endl;
		}
		GDALClose( poDataset );
		
		double scale_coef = 1;
		scale_coef = pdfMax / pdfMax2;
		MutiplyRasterbyCoefficient(output, scale_coef);
	}
}
