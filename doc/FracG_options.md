# FracG options

Command line options and default values are listed below.

## Output locations

**out_dir** | _\< std :: string \>_ <br>
_default:_ fracg\_output\_ + \<name-of-shapefile\> <br>
The main directory in which all results will be written.

**graph_results_file** | _\< std :: string \>_ <br>
_default:_ graph\_vertices & graph\_branches<br>
Filenames to save graph analysis.

**graph_results_folder** | _\< std :: string \>_ <br>
_default:_ graph<br>
Folder to save graph analysis.

## Correction parametres and distances

**dist_thresh** | _\< double \>_ <br>
_default:_ 1<br>
Distances under this threshold will be considered the same location. Used for merging fault segments during initial correction procedure, and as distance in the point index map of the graph. The units are metres.

**angl_threshold** | _\< double \>_ <br>
_default:_ 25<br>
Maximum orientation difference in degrees for merging two segments whose tips are within **dist_thresh**.

**split_dist_thresh** | _\< double \>_ <br>
_default:_ dist_thresh<br>
Maximum distance at which nearby but separate faults will be considered to intersect in the graph representation (metres).

**spur_dist_thresh** | _\< double \>_ <br>
_default:_ dist_thresh<br>
Distance threshold to use in removing spurs; spurs which are shorter than this distance will be removed from the graph representation (metres).

**classify_lineaments_dist** | _\< double \>_ <br>
_default:_ dist_thresh<br>
Distance used to classify faults in terms of intersection number along their trace (metres). This distance represents the buffer width around the fault; intersections within this distance are counted for the classification.

**raster_stats_dist** | _\< double \>_ <br>
_default:_ 1.25<br>
Distance used for analysing raster data. The distance is the buffer width around the faults for computing mean values, and the length of the profile lines for computing cross gradients, parallel gradients and cumulative cross-gradients along the faults. The unit is the number of pixels of the input raster. Note that the corresponding distance in metres is the mean of the x- and y-cellsize multiplied by this factor.

**di_raster_spacing** | _\< double \>_ <br>
_default:_ 1000<br>
Pixel size of output density/intensity maps (metres).

**dist_raster_spacing** | _\< double \>_ <br>
_default:_ 500<br>
Pixel size of output distance maps (metres).

**isect_search_size** | _\< double \>_ <br>
_default:_ di\_raster\_spacing <br>
Search for intersections within this distance (metres). Using a circular sampling window this is the radius of the window.

**resample** | _\< bool \>_ <br>
_default:_ false<br>
Resampling all created raster files to a 10th of the initial cell size using cubic spline interpolation.

## Statistical parametres

**angle_param_penalty** | _\< double \>_ <br>
_default:_ 2<br>
Penalty per parameter, when fitting Gaussian distributions to the angle distribution.

**scanline_count** | _\< int \>_ <br>
_default:_ 100<br>
Number of scanlines to construct.

**scanline_spaceing** | _\< double \>_ <br>
_default:_ 10<br>
Minimum spacing of scanlines (metres).

**component** | _\< int \>_ <br>
_default:_ -1<br>
If greater than zero extract this connected component from the graph and build a  shapefile from it.

## Maximum flow options

**max_flow_cap_type** | _\< std :: string \>_ <br>
_default:_ l<br>
Type of capacity to use in maximum flow calculations, l for length, o for orientation, lo for both.

**max_flow_gradient_flow_direction** | _\< std :: string \>_ <br>
_default:_ right<br>
Direction of maximum flow (towards left, right, top, or bottom).

**max_flow_gradient_pressure_direction** | _\< std :: string \>_ <br>
_default:_ right<br>
Direction of pressure gradient for calculating maximum flow (towards left, right, top, or bottom).

**max_flow_gradient_border_amount** | _\< double \>_ <br>
_default:_ 0.05<br>
For maximum flow calculation, the border features are those that intersect with the bounding box that is reduced by this amount (0 to 1, with 0 meaning no reduction and 1 meaning the box is reduced to a point at the centre of the graph).

## Meshing options

**gmsh_cell_count** | _\< int \>_ <br>
_default:_ 10<br>
Target element count in x and y direction. This will determine the usual characteristic length (cl) of the mesh.

**gmsh_crop** | _\< double \>_ <br>
_default:_ 0<br>
Amount by which to crop the borders for mesh generation (percent).

**gmsh_show_output** | _\< bool \>_ <br>
_default:_ false<br>
Show output of gmsh while meshing and the final mesh in the gmsh GUI.

**gmsh_point_tol** | _\< double \>_ <br>
_default:_ 0.1<br>
Point tolerance used by gmsh.

**gmsh_min_cl** | _\< double \>_ <br>
_default:_ cl/10<br>
Minimum characteristic length.

**gmsh_max_dist** | _\< double \>_ <br>
_default:_ cl/2<br>
Maximum distance for refinement around side-sets in 2D.

**gmsh_min_dist** | _\< double \>_ <br>
_default:_ cl/4<br>
Minimum distance for refinement around side-sets in 2D.

**gmsh_name_ss** | _\< bool \>_ <br>
_default:_ false<br>
Name sidesets individually.

**gmsh_sample_cell_count** | _\< int \>_ <br>
_default:_ 2<br>
Number of sampling windows from which 2D-meshes should be generated.

**gmsh_sample_count** | _\< int \>_ <br>
_default:_ 10<br>
Target element count in x and y direction for the sampling windows. 
For rectangular domains this will be the target mean element size along x and y. 
This will determine the usual characteristic length (cl) of the mesh in the sampling window.

**gmsh_sw_size** | _\< double \>_ <br>
_default:_ -1<br>
If positive, this is the fixed size for all sampling windows.

**gmsh_in_metres** | _\< bool \>_ <br>
If true, convert map units into metres.

**gmsh_sample_show_output** | _\< bool \>_ <br>
_default:_ false<br>
Show output of gmsh while meshing and the final mesh in the gmsh GUI for every sampling window.

**gmsh_in_show_metres** | _\< bool \>_ <br>
_default:_ false<br>
Convert coordinates into metres. This may be necessary for small-scale models.

## Skip certain analyses

**skip_length_distribution** | _\< bool \>_ <br>
Skip calculating the length distribution analysis.

**skip_betweenness_centrality** | _\< bool \>_ <br>
Skip calculating the betweenness centrality.

**skip_meshing** | _\< bool \>_ <br>
Skip meshing.


