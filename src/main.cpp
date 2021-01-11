#include <boost/program_options.hpp>
#include "../include/graph.h"
#include "../include/GeoRef.h" 
#include "../include/geometrie.h"
#include "../include/fracg.h"
#include "../include/stats.h"
#include "../include/fracg.h"
#include "../include/model.h"
#include "../include/util.h"

namespace po = boost::program_options;
namespace fs = boost::filesystem;

//string constants for parsing input arguments
const char *SHAPEFILE="shapefile";
const char *RASTER_FILE="raster_file";
const char *SOURCE_FILE="source_file";
const char *OUT_DIR="out_dir";

const char *DIST_THRESH="dist_thresh";
const char *ANGL_THRESH="angl_threshold";
const char *DFD_THRESH="dfd_threshold";

const char *MAP_DIST_THRESH="map_dist_thesh";
const char *SPLIT_DIST_THRESH="split_dist_thresh";
const char *SPUR_DIST_THRESH="spur_dist_thresh";
const char *CLASSIFY_LINEAMENTS_DIST="classify_lineaments_dist";
const char *RASTER_STATS_DIST="raster_stats_dist";

const char *RASTER_SPACING="raster_spacing";
const char *RASTER_SPACING2="raster_spacing2";
const char *ISECT_SEARCH_SIZE="isect_search_size";
const char *RESAMPLE="resample";

const char *ANGLE_PARAM_PENALTY="angle_param_penalty";

const char *SCANLINE_COUNT="scanline_count";
const char *SCANLINE_SPACE="scanline_spaceing";

const char *GRAPH_MIN_BRANCHES="graph_min_branches";

const char *MAX_FLOW_CAP_TYPE="max_flow_cap_type";

const char *GRAPH_RESULTS_FILE="graph_results_file";
const char *GRAPH_RESULTS_FOLDER="graph_results_folder";

const char *GMSH_CELL_COUNT="gmsh_cell_count";
const char *GMSH_SHOW_OUTPUT="gmsh_show_output";

const char *GMSH_MIN_CL="gmsh_min_cl";
const char *GMSH_MAX_DIST="gmsh_max_dist";
const char *GMSH_MIN_DIST="gmsh_min_dist";

const char *GMSH_SAMPLE_CELL_COUNT="gmsh_sample_cell_count";
const char *GMSH_SAMPLE_COUNT="gmsh_sample_count";
const char *GMSH_SAMPLE_SHOW_OUTPUT="gmsh_sample_show_output";

/* TODO:
 * 1.) kde with length
 * 2.) attached intersection boundaries with boundaryies to directed graph (MAXIMUM FLOW)
 * 3.) check capacity calculation(maximum flow)
 * 4.) check shp file generation for maximum flow(unreferenced)
 * */
 
int main(int argc, char *argv[])
{ 	 
	bool raster = true;

    //parse input arguments
    po::options_description desc("FracG Options");
    desc.add_options()
		("help", "write help message")
		(SHAPEFILE, po::value<std::string>(), "fault/fracture vector data, as a shapefile")
		(RASTER_FILE, po::value<std::string>()->default_value(""), "Raster file to use in analysis")
		(SOURCE_FILE, po::value<std::string>()->default_value(""), "source-target file to use in max flow")
		(OUT_DIR, po::value<std::string>()->default_value(""), "Write output to this directory")

		(DIST_THRESH, po::value<double>()->default_value(1), "Distances under this distance threshold will be considered the same location")
		(ANGL_THRESH, po::value<double>()->default_value(25), "Maximum difference in orientation fro merging")
		(DFD_THRESH, po::value<double>()->default_value(1), "Threshold for dicrete frechet distance (line similarity)")
		
		(MAP_DIST_THRESH, po::value<double>()->default_value(-1), "Distance threshold to consider points to be the same in a point-indexed map")
		(SPLIT_DIST_THRESH, po::value<double>()->default_value(-1), "Distance threshold to use in splitting faults into segments, for considering nearby but separate faults to actually overlap")
		(SPUR_DIST_THRESH, po::value<double>()->default_value(-1), "Distance threshold to use in removing spurs, remove spurs which are shorter than this distance")
		(CLASSIFY_LINEAMENTS_DIST, po::value<double>()->default_value(-1), "Distance used in ClassifyLineaments") //1
		(RASTER_STATS_DIST, po::value<double>()->default_value(-1), "Distance used in RasterStatistics") //2

		(RASTER_SPACING, po::value<double>()->default_value(1000.0), "Pixel size of output density/intesity maps")
		(RASTER_SPACING2, po::value<double>()->default_value(500.0), "Pixel size of output distance maps")
		(ISECT_SEARCH_SIZE, po::value<double>()->default_value(-1), "Search for intersections within this distance")
		(RESAMPLE, po::value<bool>()->default_value(false), "Resample raster filed to 1/10 of pixel size (bicubic spline)")

		(SCANLINE_COUNT, po::value<int>()->default_value(50), "Number of scalines for determining intesity and spacing")
		(SCANLINE_SPACE, po::value<double>()->default_value(10), "Minimum spacing of scanlines")

		(ANGLE_PARAM_PENALTY, po::value<double>()->default_value(2), "Penalty per parameter, when fitting Gaussians to the angle distribution")

		(GRAPH_MIN_BRANCHES, po::value<int>()->default_value(100), "Minimum number of branches(?) required in a component, in order to apply graph(?) analysis")

		(MAX_FLOW_CAP_TYPE, po::value<std::string>()->default_value("l"), "Type of capacity to use in maximum flow calculations, l for length, o for orientation, lo for both")

		(GRAPH_RESULTS_FILE, po::value<std::string>()->default_value("first"), "Filename to save graph analysis results to")
		(GRAPH_RESULTS_FOLDER, po::value<std::string>()->default_value("g"), "Save the resulting graph to this folder")

		(GMSH_CELL_COUNT, po::value<int>()->default_value(10), "GMSH Cell Count")
		(GMSH_SHOW_OUTPUT, po::bool_switch(), "Show GMSH output in the GMSH viewer")

		(GMSH_MIN_CL, po::value<double>()->default_value(-1), "Minimum characteristic length for mesh")
		(GMSH_MIN_DIST, po::value<double>()->default_value(-1), "Minimum distance for mesh refinemnt around lineament")
		(GMSH_MAX_DIST, po::value<double>()->default_value(-1), "Maximum distance for mesh refinemnt around lineament")

		(GMSH_SAMPLE_CELL_COUNT, po::value<int>()->default_value(15), "GMSH Sample Network cell count")
		(GMSH_SAMPLE_COUNT, po::value<int>()->default_value(2), "GMSH Sample Network sample count")
		(GMSH_SAMPLE_SHOW_OUTPUT, po::bool_switch(), "GMSH show Sample Network output")
    ;
    po::positional_options_description pdesc;
    pdesc.add(SHAPEFILE, 1); //first positional argument is the shapefile
    pdesc.add(RASTER_FILE, 2); //second positional argument is the raster file
    pdesc.add(SOURCE_FILE, 3); //third positional argument is the source-targhet fro maximum flow
    po::variables_map vm;
    po::store(po::command_line_parser(argc, argv).options(desc).positional(pdesc).run(),  vm);
    po:notify(vm);
    
	if (vm.count("help"))
	{
		std::cout << desc;
		return EXIT_SUCCESS;
	}
    
    if (vm.count(SHAPEFILE))
    {
        std::cout << "The shapefile is given as " << vm[SHAPEFILE].as<std::string>() << std::endl;
    } else {
        std::cerr << "ERROR: Need to give a shapefile for the vector data" << std::endl;
        return EXIT_FAILURE;
    }
    
	FracG::Graph graph;			//graph data structure 
	//TODO: Need to put in some variables here to control whether or not to do each of these functions
	const std::string shapefile_name = vm[SHAPEFILE].as<std::string>();
    std::string raster_name = vm[RASTER_FILE].as<std::string>();
    std::string source_name = vm[SOURCE_FILE].as<std::string>();
    const std::string out_dir = vm[OUT_DIR].as<std::string>();

	const double dist_threshold = vm[DIST_THRESH].as<double>(); //NOTE: different things use different distance thresholds, we need to sort out whether or not they should be made consistent, or use separate thresholds
	const double angl_threshold = vm[ANGL_THRESH].as<double>(); 
	const double dfd_threshold = vm[DFD_THRESH].as<double>();
	
    double map_dist_thresh = vm[MAP_DIST_THRESH].as<double>();
    if (map_dist_thresh < 0) map_dist_thresh = dist_threshold; //if not separately specified, use the default distance threshold
    double split_dist_thresh = vm[SPLIT_DIST_THRESH].as<double>();
    if (split_dist_thresh < 0) split_dist_thresh = dist_threshold;
    double spur_dist_thresh = vm[SPUR_DIST_THRESH].as<double>();
    if (spur_dist_thresh < 0) spur_dist_thresh = dist_threshold;
    
    const double raster_spacing = vm[RASTER_SPACING].as<double>();
    const double raster_spacing2 = vm[RASTER_SPACING2].as<double>();
    double isect_search_size = vm[ISECT_SEARCH_SIZE].as<double>();
    if (isect_search_size < 0) isect_search_size = raster_spacing;
	const bool resample = vm[RESAMPLE].as<bool>();
    
    const int scanline_count = vm[SCANLINE_COUNT].as<int>();
    const double scanline_spaceing = vm[SCANLINE_SPACE].as<double>();
    
    const double angle_param_penalty = vm[ANGLE_PARAM_PENALTY].as<double>();
    
    const int graph_min_branches = vm[GRAPH_MIN_BRANCHES].as<int>();
    
    //some of these distances are int's, maybe they should be doubles
    double classify_lineaments_dist = vm[CLASSIFY_LINEAMENTS_DIST].as<double>();
    if (classify_lineaments_dist < 0) classify_lineaments_dist = dist_threshold;
    double raster_stats_dist = vm[RASTER_STATS_DIST].as<double>();
    if (raster_stats_dist < 0) raster_stats_dist = dist_threshold;
    
    const std::string max_flow_cap_type = vm[MAX_FLOW_CAP_TYPE].as<std::string>();

    const bool save_kde_params = true; //the code that uses this needs to be improved - it should use local variables, not a global
	
	const std::string graph_results_filename = vm[GRAPH_RESULTS_FILE].as<std::string>();
    const std::string graph_results_folder = vm[GRAPH_RESULTS_FOLDER].as<std::string>(); //we should allow for this to be empty/null, and use that to signify not saving these results
    
    const int gmsh_cell_count = vm[GMSH_CELL_COUNT].as<int>();
    const bool gmsh_show_output = vm[GMSH_SHOW_OUTPUT].as<bool>();
    
	double gmsh_min_cl = vm[GMSH_MIN_CL].as<double>();
    double gmsh_min_dist = vm[GMSH_MIN_DIST].as<double>();
    double gmsh_max_dist = vm[GMSH_MAX_DIST].as<double>();
    
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

	//if (raster_name == "") raster_name = (source_dir/"DEM.tif").string();
	if (raster_name == "") raster = false;
	if (source_name == "") source_name  = (source_dir).string();

	// this is the correction of the network
	FracG::VECTOR lines = FracG::ReadVector(shapefile_name, out_path.string());		  // read the first layer of the shape file
	FracG::CorrectNetwork(lines.data, dist_threshold, angl_threshold, dfd_threshold);					 // rejoin faults that are incorrectly split in the data file (number is the critical search radius)
	
	// the following functions analyse statistical properties of the network
 	FracG::GetLengthDist(lines); 							     // test for three distributions of length 
 	FracG::DoBoxCount(lines); 								    // Boxcounting algorithm 
 	FracG::AngleDistribution angle_distribution = FracG::KdeEstimationStrikes(lines, angle_param_penalty); 				  //kernel density estimation
	
	FracG::WriteShapefile(lines, angle_distribution, FracG::AddPrefixSuffix(shapefile_name, "corrected_")); // this writes the shp file after correction of orientations (fits gaussians to the data) 
	FracG::ScanLine(lines, scanline_count, angle_distribution, scanline_spaceing);								// sanline analysis of density and spacing (number is number of scalines to generate)
	FracG::CreateStats(lines, angle_distribution); 								   // statistical analysis of network	

	FracG::P_Maps(lines, raster_spacing, resample); 			//create P20 and P21 map (second argument is the pixel resolution)
	FracG::D_Maps(lines, raster_spacing2, resample);		//create distance maps (second argument is the pixel resolution)

    FracG::graph_map<FracG::point_type, FracG::vertex_type, FracG::Graph> gm = FracG::ConvertLinesToGraph(lines.data, lines.refWKT, map_dist_thresh); 	  //convert the faults into a graph
    graph = gm.GetGraph();
	FracG::graph_map<> split_map = FracG::SplitFaults(gm, split_dist_thresh);//50 						 //split the faults in the graph into fault segments, according to the intersections of the  (number is merging radsius around line tips)
	FracG::RemoveSpurs(split_map, spur_dist_thresh); 					 //remove any spurs from the graph network (number is the minimum length of lineamants; everything below will be removed)
	graph = split_map.GetGraph();
	
	FracG::GraphAnalysis(graph, lines, graph_min_branches, angle_param_penalty, (out_path / graph_results_filename).string());		//graph, vector data, minimum number of branches per component to analyse
	FracG::WriteGraph(graph, lines, graph_results_folder);		//write a point-and line-shapefile containing the elements of the graph (string is subfolder name)
	
	//FracG::ClassifyLineaments(graph, lines, angle_distribution, classify_lineaments_dist, (out_path / "classified").string());  // number is the vritical distance between lineamnt and intersection point and the string is the filename
	FracG::IntersectionMap(graph, lines, raster_spacing, isect_search_size, resample);//2000 2500 need to check what values to use here, also need to check the function itself
	
	//simple graph algorithms
	FracG::Graph m_tree = FracG::MinTree(split_map, map_dist_thresh, (out_path / "minimum_spanning_tree").string());							 //just a minimum spanning tree
	FracG::Graph s_path = FracG::ShortPath(split_map, source_name, (out_path / "shortest_path").string());	//shortest path between points provited by shp-file. number is the merging radius to teh existing graph
	FracG::MaximumFlow_HG(graph,  source_name, 10, 0, "l", lines.refWKT, (out_path/in_stem).string());			 //maximum flow with horizontal gradient, capacity derived from orientation
	FracG::MaximumFlow_VG(graph,  source_name, 10, 0, "l", lines.refWKT, (out_path/in_stem).string());			//maximum flow with vertical gradient, capacity derived from length and orientation 

/*
	FracG::RasterStatistics(lines, raster_stats_dist, raster_name);		//parameters are the lineament set , the pixel size for the cross gradinet and the name of the raster file
	//building a graph with raster values assigned to elemnets. Numbers are splitting distance and minimum length
	FracG::Graph r_graph = FracG::BuildRasterGraph(lines, split_dist_thresh, spur_dist_thresh, map_dist_thresh, angle_param_penalty, raster_name);//5 5, another distance threshold to check
	FracG::MaximumFlow_R(r_graph,  source_name, max_flow_cap_type, lines.refWKT, (out_path/in_stem).string());				  //maximum flow with raster data, capacity derived from length

	fs::path mesh_dir = out_path / "mesh/";
	FracG::WriteGmsh_2D(gmsh_show_output, graph, gmsh_cell_count, gmsh_min_cl, gmsh_min_dist, gmsh_max_dist, ( mesh_dir / "2D_mesh").string());						 //create a 2D mesh. Number is the target elemnt number in x and y and string is the filename
	FracG::SampleNetwork_2D(gmsh_sample_show_output, lines, gmsh_sample_cell_count, gmsh_sample_count, map_dist_thresh, (mesh_dir / "a_messample").string());	//sample the network and create random subnetworks. First number is target elemnt number in x and y and second number is the number of samples.
	FracG::WriteGmsh_3D(gmsh_show_output, graph, gmsh_cell_count, gmsh_min_cl, gmsh_min_dist, gmsh_max_dist, 5000, ( mesh_dir / "3D_mesh").string());						 //create a 3D mesh. Number is the target elemnt number in x and y and string is the filename
	*/
	return EXIT_SUCCESS;
} 
