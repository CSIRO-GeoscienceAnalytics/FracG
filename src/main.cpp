#include <boost/program_options.hpp>

#include "../include/FracG.h"

using namespace FGraph;

namespace po = boost::program_options;
namespace fs = boost::filesystem;

//string constants for parsing input arguments
const char *SHAPEFILE="shapefile";
const char *OUT_DIR="out_dir";
//const char *OUT_SUBDIR="out_subdir";
const char *DIST_THRESH="dist_thresh";
const char *SCANLINE_COUNT="scanline_count";
const char *RASTER_SPACING="raster_spacing";

const char *PRINT_KMEANS="print_kmeans_progress";
const char *GRAPH_MIN_BRANCHES="graph_min_branches";

const char *CLASSIFY_LINEAMENTS_DIST="classify_lineaments_dist";
const char *RASTER_STATS_DIST="raster_stats_dist";

const char *MAX_FLOW_CAP_TYPE="max_flow_cap_type";

const char *GRAPH_RESULTS_FILE="graph_results_file";
const char *GRAPH_RESULTS_FOLDER="graph_results_folder";

const char *GMSH_CELL_COUNT="gmsh_cell_count";
const char *GMSH_SHOW_OUTPUT="gmsh_show_output";

const char *GMSH_SAMPLE_CELL_COUNT="gmsh_sample_cell_count";
const char *GMSH_SAMPLE_COUNT="gmsh_sample_count";
const char *GMSH_SAMPLE_SHOW_OUTPUT="gmsh_sample_show_output";

int main(int argc, char *argv[])
{ 
	/*BUGS
	* WriteSHp_line creates nosence for the file name
	* WriteSHP_R cannot overwite the file once it is created
	* Shortespath Segmentation fault for test case 1 <- solved
	* 
	* These functions need some attention:
	* 1.) Build raster graph 	(slow)
	* 2.) Intersection map		(slow)
	* 3.) Maximum Flow			(general check)
	* 4.) Please check the graph algorithms and the meshing as well
	* 5.) The folder structure migth need to be cleaned
    * 6.) Had a crash when testing. Probably using undefinded data/bad memory accesses somewhere. Test with valgrind. the gmsh code also had a crash/failure when I tried to view the output. This may be related to a file/other data already exising ("layer .shp already exists"), or a variable name problem.
	* */
	 
    //parse input arguments
    po::options_description desc("FracG Options");
    desc.add_options()
        ("help", "write help message")
        (SHAPEFILE, po::value<string>(), "fault/fracture vector data, as a shapefile")
        (OUT_DIR, po::value<string>()->default_value(""), "Write output to this directory")
        //(OUT_SUBDIR, po::value<string>()->default_value("fracg_output"), "Write output to this subdirectory")
        (DIST_THRESH, po::value<double>()->default_value(1), "Distances under this distance threshold will be considered the same location")
        (SCANLINE_COUNT, po::value<int>()->default_value(50), "Number of scalines to check for orientations (?)")
        (RASTER_SPACING, po::value<double>()->default_value(3000.0), "Pixel size of output raster maps")
        (GRAPH_MIN_BRANCHES, po::value<int>()->default_value(100), "Minimum number of branches(?) required in a component, in order to apply graph(?) analysis")
        (CLASSIFY_LINEAMENTS_DIST, po::value<int>()->default_value(1), "Distance used in ClassifyLineaments")
        (RASTER_STATS_DIST, po::value<int>()->default_value(2), "Distance used in RasterStatistics")
        
        (MAX_FLOW_CAP_TYPE, po::value<string>()->default_value("l"), "Type of capacity to use in maximum flow calculations, l for length, o for orientation, lo for both")
        
        (PRINT_KMEANS, po::bool_switch(), "Print the progress results for kmeans clustering")
        (GRAPH_RESULTS_FILE, po::value<string>()->default_value("first"), "Filename to save graph analysis results to")
        (GRAPH_RESULTS_FOLDER, po::value<string>()->default_value("g"), "Save the resulting graph to this folder")
        
        (GMSH_CELL_COUNT, po::value<int>()->default_value(15), "GMSH Cell Count")
        (GMSH_SHOW_OUTPUT, po::bool_switch(), "Show GMSH output in the GMSH viewer")
        
        (GMSH_SAMPLE_CELL_COUNT, po::value<int>()->default_value(15), "GMSH Sample Network cell count")
        (GMSH_SAMPLE_COUNT, po::value<int>()->default_value(2), "GMSH Sample Network sample count")
        (GMSH_SAMPLE_SHOW_OUTPUT, po::bool_switch(), "GMSH show Sample Network output")
    ;
    po::positional_options_description pdesc;
    pdesc.add(SHAPEFILE, 1); //first positional argument is the shapefile
    po::variables_map vm;
    po::store(po::command_line_parser(argc, argv).options(desc).positional(pdesc).run(),  vm);
    po:notify(vm);
    
    if (vm.count("help"))
    {
        cout << desc;
    }
    
    if (vm.count(SHAPEFILE))
    {
        cout << "The shapefile is given as " << vm[SHAPEFILE].as<string>() << endl;
    } else {
        cerr << "ERROR: Need to give a shapefile for the vector data" << endl;
        return EXIT_FAILURE;
    }
    
	Graph graph;			//graph data structure 
	map_vertex_type map; //a map of the vertices in graph G, for quick retrieval of vertices by their location
	
	//TODO: Need to put in some variables here to control whether or not to do each of these functions
	
	const string shapefile_name = vm[SHAPEFILE].as<string>();
    const string out_dir = vm[OUT_DIR].as<string>();
//     const string out_subdir = vm[OUT_SUBDIR].as<string>();
	
	const double dist_threshold = vm[DIST_THRESH].as<double>(); //NOTE: different things use different distance thresholds, we need to sort out whether or not they should be made consistent, or use separate thresholds
    const double raster_spacing = vm[RASTER_SPACING].as<double>();
    
    const int scanline_count = vm[SCANLINE_COUNT].as<int>();
    const int graph_min_branches = vm[GRAPH_MIN_BRANCHES].as<int>();
    
    //some of these distances are int's, maybe they should be doubles
    const int classify_lineaments_dist = vm[CLASSIFY_LINEAMENTS_DIST].as<int>();
    const int raster_stats_dist = vm[RASTER_STATS_DIST].as<int>();
    
    const string max_flow_cap_type = vm[MAX_FLOW_CAP_TYPE].as<string>();
    
    const bool print_kmeans = vm[PRINT_KMEANS].as<bool>();
    const bool save_kde_params = true; //the code that uses this needs to be improved - it should use local variables, not a global
	
	const string graph_results_filename = vm[GRAPH_RESULTS_FILE].as<string>();
    const string graph_results_folder = vm[GRAPH_RESULTS_FOLDER].as<string>(); //we should allow for this to be empty/null, and use that to signify not saving these results
    
    const int gmsh_cell_count = vm[GMSH_CELL_COUNT].as<int>();
    const int gmsh_show_output = vm[GMSH_SHOW_OUTPUT].as<bool>();
    
    const int gmsh_sample_cell_count = vm[GMSH_SAMPLE_CELL_COUNT].as<int>();
    const int gmsh_sample_count = vm[GMSH_SAMPLE_COUNT].as<int>();
    const bool gmsh_sample_show_output = vm[GMSH_SAMPLE_SHOW_OUTPUT].as<bool>();
	
    //input and output files and directories
    fs::path vector_file(shapefile_name);
    fs::path in_stem = vector_file.stem();
    fs::path source_dir = vector_file.parent_path();
    std::string output_dir_string = out_dir;
    if (out_dir == "") output_dir_string = (vector_file.parent_path() / ("fracg_output_" + in_stem.string())).string();
    fs::path out_path( output_dir_string + "/" );
//     cout << "in file is " << vector_file << ", in stem=" << in_stem << ", source_dir="<<source_dir<<",our_parent="<<out_parent<<", out_path="<<out_path<<endl;
    
	// this is the correction of the network
	VECTOR lines = geo.ReadVector(shapefile_name, out_path.string());		  // read the first layer of the shape file
	geo.CorrectNetwork(lines.data, dist_threshold);					 // rejoin faults that are incorrectly split in the data file (number is the critical search radius)

	geo.WRITE_SHP(lines, FGraph::add_prefix_suffix(shapefile_name, "corrected_"));					// this writes the shp file after correction

	// the following functions analyse staatistical properties of the network
 	stats.GetLengthDist(lines); 							     // test for three distributions of length 
 	stats.DoBoxCount(lines); 								    // Boxcounting algorithm 
 	stats.CreateStats(lines); 								   // statistical analysis of network
 	stats.KDE_estimation_strikes(lines, save_kde_params); 				  //kernel density estimation of orientations (fits gaussians to the data; boolean whether to write in namespace variable)
 	stats.KMCluster(print_kmeans, lines);							 // KM clustering
 	stats.ScanLine(lines, scanline_count);								// sanline analysis of density and spacing (number is number of scalines to generate)
 
 	// Here we create some raster files that characterize the spatial arangement
	geom.CentreDistanceMap(lines, raster_spacing);   //fault centre to fault centre distance (second argument is the pixel resolution)
	geom.P_Maps(lines, raster_spacing); 			//create P20 and P21 map (second argument is the pixel resolution)

	//this creates a georeferences graph, analyses it, and writes two shp files containing edges and vertices of the graph
	G.ReadVEC(graph, map, lines.data); 					  //convert the faults into a graph
	G.SplitFaults(graph, map, dist_threshold);//50 						 //split the faults in the graph into fault segments, according to the intersections of the  (number is merging radsius around line tips)
	G.RemoveSpurs(graph, map, dist_threshold);//100 					 //remove any spurs from the graph network (number is the minimum length of lineamants; everything below will be removed)
	G.GraphAnalysis(graph, lines, graph_min_branches, (out_path / graph_results_filename).string());		//graph, vector data, minimum number of branches per component to analyse
	geo.WriteGraph(graph, lines, graph_results_folder);		//write a point-and line-shapefile containing the elements of the graph (string is subfolder name)
	
	//simple graph algorithms
	Graph m_tree = G.MinTree(graph, (out_path / "min_tree").string());							 //just a minimum spanning tree
	Graph s_path = G.ShortPath(graph, map, (source_dir/"S_T.shp").string(), (out_path / "short_path").string());	//shortest path between points provited by shp-file. number is the merging radius to teh existing graph
	
	//create a shp file with lineaments classified based on orientation (from KDE) and intersecions (from graph)
	G.ClassifyLineaments(graph, lines, classify_lineaments_dist, (out_path/"classified").string());  // number is the vritical distance between lineamnt and intersection point and the string is the filename

	stats.RasterStatistics(lines, raster_stats_dist, (source_dir/"DEM.tif").string());		//parameters are the lineament set , the pixel size for the cross gradinet and the name of the raster file

	//building a graph with raster values assigned to elemnets. Numbers are splitting distance and minimum length
	Graph r_graph = geo.BuildRasterGraph(lines, dist_threshold, dist_threshold, (source_dir/"DEM.tif").string());//5 5, another distance threshold to check
	
	G.MaximumFlow_R(r_graph, map, (source_dir/"S_T.shp").string(), max_flow_cap_type, (out_path/in_stem).string());				  //maximum flow with raster data, capacity derived from length
//	G.MaximumFlow_HG(graph, map, "S_T.shp", 1, 0, "o");			 //maximum flow with horizontal gradient, capacity derived from orientation
//	G.MaximumFlow_VG(graph, map, "S_T.shp", 1, 0, "l");			//maximum flow with vertical gradient, capacity derived from length and orientation 
	
	//create a intersection density map with circular sampling window.
	//First number is pixel size and second number is the search radius.(this is quite slow at the moment; ?smth wrong with the tree?)
	G.IntersectionMap(graph, lines, raster_spacing, raster_spacing);//2000 2500 need to check what values to use here, also need to check the function itself
	
    fs::path mesh_dir = out_path / "mesh";
    
	m.WriteGmsh_2D(gmsh_show_output, graph, gmsh_cell_count, ( mesh_dir / "a_mesh").string());						 //create a 2D mesh. Number is the target elemnt number in x and y and string is the filename
	m.SampleNetwork_2D(gmsh_sample_show_output, lines.data, gmsh_sample_cell_count, gmsh_sample_count, (mesh_dir / "a_messample").string());	//sample the network and create random subnetworks. First number is target elemnt number in x and y and second number is the number of samples.
	return EXIT_SUCCESS;
} 
