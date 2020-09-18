#include <boost/program_options.hpp>


#include "../include/graph.h"
#include "../include/GeoRef.h" 
#include "../include/geometrie.h"
#include "../include/fracg.h"
#include "../include/stats.h"
#include "../include/fracg.h"
#include "../include/model.h"
#include "../include/util.h"

//using namespace FracG;

namespace po = boost::program_options;
namespace fs = boost::filesystem;

//string constants for parsing input arguments
const char *SHAPEFILE="shapefile";
const char *RASTER_FILE="raster_file";
const char *OUT_DIR="out_dir";
//const char *OUT_SUBDIR="out_subdir";
const char *DIST_THRESH="dist_thresh";
const char *MAP_DIST_THRESH="map_dist_thesh";
const char *SPLIT_DIST_THRESH="split_dist_thresh";
const char *SPUR_DIST_THRESH="spur_dist_thresh";
const char *CLASSIFY_LINEAMENTS_DIST="classify_lineaments_dist";
const char *RASTER_STATS_DIST="raster_stats_dist";

const char *RASTER_SPACING="raster_spacing";
const char *ISECT_SEARCH_SIZE="isect_search_size";

const char *SCANLINE_COUNT="scanline_count";

const char *ANGLE_PARAM_PENALTY="angle_param_penalty";

const char *PRINT_KMEANS="print_kmeans_progress";
const char *GRAPH_MIN_BRANCHES="graph_min_branches";

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
        (RASTER_FILE, po::value<string>()->default_value(""), "Raster file to use in analysis")
        (OUT_DIR, po::value<string>()->default_value(""), "Write output to this directory")
        //(OUT_SUBDIR, po::value<string>()->default_value("fracg_output"), "Write output to this subdirectory")
        (DIST_THRESH, po::value<double>()->default_value(1), "Distances under this distance threshold will be considered the same location")
        (MAP_DIST_THRESH, po::value<double>()->default_value(-1), "Distance threshold to consider points to be the same in a point-indexed map")
        (SPLIT_DIST_THRESH, po::value<double>()->default_value(-1), "Distance threshold to use in splitting faults into segments, for considering nearby but separate faults to actually overlap")
        (SPUR_DIST_THRESH, po::value<double>()->default_value(-1), "Distance threshold to use in removing spurs, remove spurs which are shorter than this distance")
        (CLASSIFY_LINEAMENTS_DIST, po::value<double>()->default_value(-1), "Distance used in ClassifyLineaments") //1
        (RASTER_STATS_DIST, po::value<double>()->default_value(-1), "Distance used in RasterStatistics") //2
        
        (RASTER_SPACING, po::value<double>()->default_value(3000.0), "Pixel size of output raster maps")
        (ISECT_SEARCH_SIZE, po::value<double>()->default_value(-1), "Search for intersections within this distance")
        
        (SCANLINE_COUNT, po::value<int>()->default_value(50), "Number of scalines to check for orientations (?)")
        
        (ANGLE_PARAM_PENALTY, po::value<double>()->default_value(2), "Penalty per parameter, when fitting Gaussians to the angle distribution")
        
        (GRAPH_MIN_BRANCHES, po::value<int>()->default_value(100), "Minimum number of branches(?) required in a component, in order to apply graph(?) analysis")
        
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
    pdesc.add(RASTER_FILE, 2); //second positional argument is the raster file
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
    
	FracG::Graph graph;			//graph data structure 
	
	//TODO: Need to put in some variables here to control whether or not to do each of these functions
	
	const string shapefile_name = vm[SHAPEFILE].as<string>();
    string raster_name = vm[RASTER_FILE].as<string>();
    const string out_dir = vm[OUT_DIR].as<string>();
//     const string out_subdir = vm[OUT_SUBDIR].as<string>();
	
	const double dist_threshold = vm[DIST_THRESH].as<double>(); //NOTE: different things use different distance thresholds, we need to sort out whether or not they should be made consistent, or use separate thresholds
    double map_dist_thresh = vm[MAP_DIST_THRESH].as<double>();
    if (map_dist_thresh < 0) map_dist_thresh = dist_threshold; //if not separately specified, use the default distance threshold
    double split_dist_thresh = vm[SPLIT_DIST_THRESH].as<double>();
    if (split_dist_thresh < 0) split_dist_thresh = dist_threshold;
    double spur_dist_thresh = vm[SPUR_DIST_THRESH].as<double>();
    if (spur_dist_thresh < 0) spur_dist_thresh = dist_threshold;
    
    const double raster_spacing = vm[RASTER_SPACING].as<double>();
    double isect_search_size = vm[ISECT_SEARCH_SIZE].as<double>();
    if (isect_search_size < 0) isect_search_size = raster_spacing;
    
    const int scanline_count = vm[SCANLINE_COUNT].as<int>();
    
    const double angle_param_penalty = vm[ANGLE_PARAM_PENALTY].as<double>();
    
    const int graph_min_branches = vm[GRAPH_MIN_BRANCHES].as<int>();
    
    //some of these distances are int's, maybe they should be doubles
    double classify_lineaments_dist = vm[CLASSIFY_LINEAMENTS_DIST].as<double>();
    if (classify_lineaments_dist < 0) classify_lineaments_dist = dist_threshold;
    double raster_stats_dist = vm[RASTER_STATS_DIST].as<double>();
    if (raster_stats_dist < 0) raster_stats_dist = dist_threshold;
    
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
    
    if (raster_name == "") raster_name = (source_dir/"DEM.tif").string();
    
	// this is the correction of the network
	FracG::VECTOR lines = FracG::ReadVector(shapefile_name, out_path.string());		  // read the first layer of the shape file
	FracG::CorrectNetwork(lines.data, dist_threshold);					 // rejoin faults that are incorrectly split in the data file (number is the critical search radius)

	// the following functions analyse staatistical properties of the network
 	FracG::GetLengthDist(lines); 							     // test for three distributions of length 
 	FracG::DoBoxCount(lines); 								    // Boxcounting algorithm 
 	FracG::AngleDistribution angle_distribution = FracG::KdeEstimationStrikes(lines, angle_param_penalty); 				  //kernel density estimation

	FracG::WriteShapefile(lines, angle_distribution, FracG::AddPrefixSuffix(shapefile_name, "corrected_")); // this writes the shp file after correction of orientations (fits gaussians to the data) 
	
 	FracG::CreateStats(lines, angle_distribution); 								   // statistical analysis of network
 	
 	FracG::KMCluster(print_kmeans, lines, angle_distribution);							 // KM clustering
 	FracG::ScanLine(lines, scanline_count, angle_distribution);								// sanline analysis of density and spacing (number is number of scalines to generate)
 
 	// Here we create some raster files that characterize the spatial arrangement
	FracG::CentreDistanceMap(lines, raster_spacing);   //fault centre to fault centre distance (second argument is the pixel resolution)
	FracG::PMaps(lines, raster_spacing); 			//create P20 and P21 map (second argument is the pixel resolution)

	//this creates a geo-referenced graph, analyses it, and writes two shp files containing edges and vertices of the graph
// 	map_vertex_type map;
    FracG::graph_map<FracG::point_type, FracG::vertex_type, FracG::Graph> gm = FracG::ConvertLinesToGraph(lines.data, lines.refWKT, map_dist_thresh); 	  //convert the faults into a graph
    graph = gm.GetGraph();
	FracG::graph_map<> split_map = FracG::SplitFaults(gm, split_dist_thresh);//50 						 //split the faults in the graph into fault segments, according to the intersections of the  (number is merging radsius around line tips)
	FracG::RemoveSpurs(split_map, spur_dist_thresh);//100 					 //remove any spurs from the graph network (number is the minimum length of lineamants; everything below will be removed)
	graph = split_map.GetGraph();
	
	FracG::GraphAnalysis(graph, lines, graph_min_branches, angle_param_penalty, (out_path / graph_results_filename).string());		//graph, vector data, minimum number of branches per component to analyse
	FracG::WriteGraph(graph, lines, graph_results_folder);		//write a point-and line-shapefile containing the elements of the graph (string is subfolder name)
	
	//simple graph algorithms
	FracG::Graph m_tree = FracG::MinTree(split_map, map_dist_thresh, (out_path / "minimum_spanning_tree").string());							 //just a minimum spanning tree
	FracG::Graph s_path = FracG::ShortPath(split_map, (source_dir / "S_T.shp").string(), (out_path / "shortest_path").string());	//shortest path between points provited by shp-file. number is the merging radius to teh existing graph
	
	//create a shp file with lineaments classified based on orientation (from KDE) and intersecions (from graph)
	FracG::ClassifyLineaments(graph, lines, angle_distribution, classify_lineaments_dist, (out_path / "classified").string());  // number is the vritical distance between lineamnt and intersection point and the string is the filename

	FracG::RasterStatistics(lines, raster_stats_dist, raster_name);		//parameters are the lineament set , the pixel size for the cross gradinet and the name of the raster file

	//building a graph with raster values assigned to elemnets. Numbers are splitting distance and minimum length
	FracG::Graph r_graph = FracG::BuildRasterGraph(lines, split_dist_thresh, spur_dist_thresh, map_dist_thresh, angle_param_penalty, raster_name);//5 5, another distance threshold to check
	
	FracG::MaximumFlow_R(r_graph, (source_dir / "S_T.shp").string(), max_flow_cap_type, lines.refWKT, (out_path/in_stem).string());				  //maximum flow with raster data, capacity derived from length
//	MaximumFlow_HG(graph, "S_T.shp", 1, 0, "o");			 //maximum flow with horizontal gradient, capacity derived from orientation
//	MaximumFlow_VG(graph, "S_T.shp", 1, 0, "l");			//maximum flow with vertical gradient, capacity derived from length and orientation 
	
	//create a intersection density map with circular sampling window.
	//First number is pixel size and second number is the search radius.(this is quite slow at the moment; ?smth wrong with the tree?)
	FracG::IntersectionMap(graph, lines, raster_spacing, isect_search_size);//2000 2500 need to check what values to use here, also need to check the function itself
	
    fs::path mesh_dir = out_path / "mesh/";
    
	FracG::dist_tree dtree = FracG::BuildPointTree(graph);
	FracG::WriteGmsh2D(dtree, gmsh_show_output, graph, gmsh_cell_count, ( mesh_dir / "a_mesh").string());						 //create a 2D mesh. Number is the target elemnt number in x and y and string is the filename
	FracG::SampleNetwork2D(dtree, gmsh_sample_show_output, lines, gmsh_sample_cell_count, gmsh_sample_count, map_dist_thresh, (mesh_dir / "a_messample").string());	//sample the network and create random subnetworks. First number is target elemnt number in x and y and second number is the number of samples.
	return EXIT_SUCCESS;
} 
