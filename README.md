# FracG - Fault and fracture analysis and meshing software

FracG is a command line Debian GNU/Linux based software that performs analysis of discontinuity data (e.g. faults and fractures). Currently, the framework is only tested for Ubunutu 20.04 LTS.
The input is provided as line-vector data in shape-file format. Prior to every analysis FracG will try to clean up the given vector data by removing duplicate entries and 
flaws in the topology. The analysis can be
performed with or without obtaining the best-fit models for length distribution, meshing or including raster
data. 

The analyses are as follows:
* Models of fault length distribution (exponential, log-normal, or power-law) and principal orientations
(Gaussians fitted to the Kernel density estimation) obtained by automated fitting
* Statistical parameters including scan line analysis (exported as csv)
* Box-counting to obtain fractal dimension
* Density maps are generated (frequency, intersections, distances)
* Georeferenced graph is created
* Graph algorithms:
  * Betweenness Centrality
  * Minimum Spanning tree
  * Shortest path
  * Maximum Flow
* Classify discontinuities (orientation, intersection numbers)
* Analysis can be extended by including raster data (e.g. DEM)
  * Extract raster values for line geometries
  * Build georeferenced graph including raster data as edge and vertex weights
* Generation of 2D Finite Element(FE) conforming meshes (with discontinuities defined as side sets)
  * Generate series of randomly sampled FE-meshes (rectangular windows)

## Install from Source
FracG is designed for Debian GNU/Linux and requires third-party libraries as outlined below.

First, get the latest GDAL/OGR version, add the PPA to your sources:
```bash
sudo add-apt-repository ppa:ubuntugis/ppa && sudo apt-get update
```

Now the necessary libraries can be installed from the terminal:
```bash
sudo apt-get install \
build-essential \
libgdal-dev \
libboost-all-dev \
liblapack-dev \
libblas-dev \
libgsl-dev \
libgmsh-dev \
libarmadillo-dev
```
Export the environmental variables for gdal
```bash
export CPLUS_INCLUDE_PATH=/usr/include/gdal
export C_INCLUDE_PATH=/usr/include/gdal
```
To obtain gmsh visit http://gmsh.info/.
Note that gmsh's open cascade engine is used for meshing.

## Build
Obtain cmake first.
```bash
sudo apt-get install cmake
```
In the FracG directory, type;
```bash
mkdir build \
cd build \
cmake .. \
sudo make install
```
FracG can now be executed from the command line.

## Licence

FracG is copyright (c) Commonwealth Scientifc and Industrial Research Organisation (CSIRO) ABN 41 687 119 230. Except
where otherwise indicated, including in the Supplementary Licence, CSIRO grants you a licence to the Software on the terms of
the GNU General Public Licence version 3 (GPLv3), distributed at: http://www.gnu.org/licenses/gpl.html.

[The following additional terms apply under clause 7 of GPLv3 to the licence of the Software that is granted by CSIRO.](licence.txt)

## Executing FracG

After installation with cmake FracG will be set as a global executable in your environment. In a terminal choose your current directory in which the files you wish to analyse are located and enter the command:

`FracG <.shp> <.tif> <.shp>`

The arguments are:
1. (required) path/name of the shapefile containing the faults, including extension
1. (optional) path/name of the raster-file in GeoTIFF format including extension
1. (optional) path/name of the point-shape-file including extension. This file needs to contain two points that will be used for computing the shortest path and maximum flow between them.

All input files must be in the same coordinate reference system.

## Options

FracG has several parameters that can defined by the user. To see what options are available type:

`FracG --help`

The optional parameters can be set after the input file(s). They can be given in any order and you can add as many as you want:

`FracG <.shp> --option1 --option2 ...`

### Output file/folder names

**out_dir** | _\< std :: string \>_ <br>
_default:_ fracg\_output\_ + \<name-of-shapefile\> <br>
The main directory in which all results will be written.

**graph_results_file** | _\< std :: string \>_ <br>
_default:_ graph\_vertices & graph\_branches<br>
Filenames to save graph analysis results to.

**graph_results_folder** | _\< std :: string \>_ <br>
_default:_ graph<br>
Folder to save graph analysis results to.

### Correction parameters and distances

**dist_thresh** | _\< double \>_ <br>
_default:_ 1<br>
Distances under this distance threshold will be considered the same location. Used for merging line segments and as distance in the point index map of the graph. The units are meters.

**angl_threshold** | _\< double \>_ <br>
_default:_ 25<br>
Maximum orientation difference in degrees for merging two segments whose tips are with the critical distance.

**split_dist_thresh** | _\< double \>_ <br>
_default:_ dist_thresh<br>
Distance threshold to use in splitting faults into segments, for considering nearby but separate faults to actually overlap. Used for fixing flaws in digitisation leading to false intersection classification. The units are in meters.

**spur_dist_thresh** | _\< double \>_ <br>
_default:_ dist_thresh<br>
Distance threshold to use in removing spurs, remove spurs which are shorter than this distance. Used to correct for false intersection classification. The units are in meters.

**classify_lineaments_dist** | _\< double \>_ <br>
_default:_ dist_thresh<br>
Distance used in to classify lineaments in terms of intersection number along their trace. This distance represents the buffer width around the line and intersections within this distance are counted for the classification. The units are in meters.

**raster_stats_dist** | _\< double \>_ <br>
_default:_ 1.25<br>
Distance used for analysing raster data for the line-strings. The distance is the buffer width around the traces for computing mean values, the length of the profile lines for computing cross gradient, parallel gradients and cumulative cross-gradients along the line-strings. The unit is the number of pixels of the raster of the input raster. Note that the distance in meters is derived as the mean of the x- and y-cellsize multiplied by this factor.

**di_raster_spacing** | _\< double \>_ <br>
_default:_ 1000<br>
Pixel size of output density/intensity maps. The units are in meters.

**dist_raster_spacing** | _\< double \>_ <br>
_default:_ 500<br>
Pixel size of output distance maps. The units are in meters.

**isect_search_size** | _\< double \>_ <br>
_default:_ di\_raster\_spacing <br>
Search for intersections within this distance. Using a circular sampling window this is the radius of the window. The units are in meters.

**resample** | _\< bool \>_ <br>
_default:_ false<br>
Resampling all created raster files to a 10th of the initial cell size using cubic spline interpolation.

### Statistical parameters

**angle_param_penalty** | _\< double \>_ <br>
_default:_ 2<br>
Penalty per parameter, when fitting Gaussian distributions to the angle distribution.

**scanline_count** | _\< int \>_ <br>
_default:_ 100<br>
Number of scanlines to construct.

**scanline_spaceing** | _\< double \>_ <br>
_default:_ 10<br>
Minimum spacing of scanlines in meters.

**component** | _\< int \>_ <br>
_default:_ -1<br>
If greater than zero extract this connected component from the graph and build a line shape-file from it.

### Maximum flow options

**max_flow_cap_type** | _\< std :: string \>_ <br>
_default:_ l<br>
Type of capacity to use in maximum flow calculations, l for length, o for orientation, lo for both.

**max_flow_gradient_flow_direction** | _\< std :: string \>_ <br>
_default:_ right<br>
Target direction of the gradient-based maximum flow (towards left, right, top, or bottom).

**max_flow_gradient_pressure_direction** | _\< std :: string \>_ <br>
_default:_ right<br>
Target direction of the gradient-based maximum flow pressure (towards left, right, top, or bottom).

**max_flow_gradient_border_amount** | _\< double \>_ <br>
_default:_ 0.05<br>
For gradient-based maximum flow, the border features are those that intersect with the bounding box that is reduced by this amount (0 to 1).

### Gmsh options

**gmsh_cell_count** | _\< int \>_ <br>
_default:_ 10<br>
Target element count in x and y direction. For rectangular domains this will be the target mean element size along x and y. This will yield the usual characteristic length (cl) of the model

**gmsh_crop** | _\< double \>_ <br>
_default:_ 0<br>
Amount of borders to cut in percent.

**gmsh_show_output** | _\< bool \>_ <br>
_default:_ false<br>
Show output of gmsh while meshing and the final mesh in the gmsh GUI.

**gmsh_point_tol** | _\< double \>_ <br>
_default:_ 0.1<br>
Point tolerance used by gmsh.

**gmsh_min_cl** | _\< double \>_ <br>
_default:_ cl/10<br>
Minimum characteristic length. By default this is will be a 10th of the usual characteristic length (cl).

**gmsh_max_dist** | _\< double \>_ <br>
_default:_ cl/2<br>
Maximum distance for refinement around side-sets in 2D.

**gmsh_min_dist** | _\< double \>_ <br>
_default:_ cl/4<br>
Minimum distance for refinement around side-sets in 2D.

**gmsh_name_ss** | _\< bool \>_ <br>
_default:_ false<br>
Name sideset individually.

**gmsh_sample_cell_count** | _\< int \>_ <br>
_default:_ 2<br>
Number of sampling windows from which 2D-meshes should be generated.

**gmsh_sample_count** | _\< int \>_ <br>
_default:_ 10<br>
Target element count in x and y direction for the sampling windows. For rectangular domains this will be the target mean element size along x and y. This will yield the usual characteristic length (cl) of the model.

**gmsh_sw_size** | _\< double \>_ <br>
_default:_ -1<br>
If set positive, use a fixed size for all sampling windows.

**gmsh_in_meters** | _\< bool \>_ <br>
If set true, convert map units into meters.

**gmsh_sample_show_output** | _\< bool \>_ <br>
_default:_ false<br>
Show output of gmsh while meshing and the final mesh in the gmsh GUI for every sampling window.

**gmsh_in_show_meters** | _\< bool \>_ <br>
_default:_ false<br>
Convert coordinates into meters. This can be necessary for small scale models.

### Skip certain analyses

**skip_length_distribution** | _\< bool \>_ <br>
Skip calculating the length distribution analysis.

**skip_betweenness_centrality** | _\< bool \>_ <br>
Skip calculating the betweenness centrality.

**skip_meshing** | _\< bool \>_ <br>
Skip meshing.
