# FracG analyses

FracG provides tools to characterize the geometry, spatial arrangement, and topology of a 2D network of geological lineaments (hereafter termed "faults", although they could be fractures, shear zones, dykes or any other linear feature). 
Geometric parameters are investigated in a statistical manner with maximum likelihood fitting to identify the length distribution, and kernel density estimation of the principal orientations. 
The spatial arrangement of the network is investigated by applying window sampling to obtain densities, intensities and spacing. The fractal dimension is obtained via box counting.
A graph representation of the network is used to quantify its connectivity. After converting the network into a graph, classical graph algorithms such as shortest path, betweenness centrality and maximum flow can be applied to further assess the fluid flow properties of the network. 
In addition to these analyses, FracG can generate a conforming mesh of the fault network and host rock suitable for finite element fluid flow simulations.

Further details of the analyses are provided below.

## Correction/cleaning of fault network

 1. Split multi-lines into individual features
 1. Merge lines that are closer than a threshold distance and within a threshold angle
 1. Remove duplicates
 1. Write corrected shapefile to input directory (corrected\_\<original shapefile name\>)

## Geometric properties

### Fault length distribution (optional)

* Determine whether the fault length distribution follows an exponential, log-normal, 
or power-law distribution.
* Output: statistics/length\_distributions.csv

### Principal orientations

* Kernel density estimation of the principal orientations
* Outputs:
  * statistics/angle\_distribution\_kde.csv
  * Weighted by length: statistics/angle\_distribution\_kde\_length\_weights.csv

## Spatial arrangement

### Density, intensity and distance maps

* Scanline sampling is used to obtain average spacing (S) and fault intensity (P10) along a series of parallel lines
  * Scanlines constructed perpendicular to the principal orientations
  * The user can define the total number of scanlines and the desired spacing between lines of the same orientation
  * Output: statistics/scanline\_analysis.csv
* Rectangular window sampling is used to obtain density, intensity and distance maps. Each sampling window becomes a pixel in the resulting raster file. Window/pixel size can be defined by the user. The outputs are:
  * Fault density (total number of faults in each pixel; geometry/P20\_map.tif)
  * Fault intensity (cumulative length of faults in each pixel; geometry/P21\_map.tif)
  * Minimum distance of each pixel centre to the nearest fault (geometry/Distance\_map.tif)
  * Distance of each pixel centre to nearest fault centre (geometry/Center\_distance\_map.tif)

### Fractal dimension

* Obtained by box counting
* Start with boxes with side length equal to the length of the longest fault
* Count number of boxes that contain a fault
* For each box that contains a fault, divide into quadrants and repeat the count at the new size
* Continue reducing box size to user-defined number of recursions
* The fractal dimension is the slope of a best-fit line through log(N) vs. log(S), where N is the number of boxes containing a fault and S is the size of the box. 
* Output: statistics/box\_count.csv

## Topology

### Graph representation

* Topology metrics are obtained from a georeferenced graph, in which the fault network is
represented as a series of vertices (points) and edges (lines). 
* Vertices are classified depending on the number of edges that intersect with them, known as the degree of the vertex. There are 3 classifications:
  * Isolated tips (I type; degree 1)
  * Junctions (Y type; degree 3)
  * Crossings (X type; degree 4)
* The graph is constructed as follows:
   1. Edges and vertices are created from line strings and points in the corrected fault network.
   1. Intersections are identified and vertices added at these locations. Edges within a user-specified distance of each other are considered to intersect.
   1. Short "spurs" are removed, to ensure that Y vertices are not incorrectly classified as X vertices (user can specify the maximum length of spurs to be removed).
* Connected components are identified as collections of edges/vertices that are connected with each other, whereby different connected components cannot be reached from one another
* Outputs: 
  * Full graph: graph\_shp/graph/graph\_branches.shp and graph\_vertices.shp (note: directory can be changed using the option "graph\_results\_folder")
  * Connected components: graph/test.shp
* The following parameters derived from the graph are reported in statistics/graph\_statistics.csv:
   * (add list here)

### Graph algorithms

In the following graph algorithms the edges are assigned weights (called edge weights) which can be either the length of the edge or values derived from [raster analysis](#raster_analysis) (the centre value, the mean value of the edge, or the profile along the edge).

The following algorithms are applied to the graph representation of the fault network:

 * Minimum spanning tree
   * A spanning tree is a tree structure that maintains the overall connectivity of the vertices. A graph can have several of such structures whereas the minimum spanning tree is the one with the lowest sum of edge weights.
   * Output: graph/minimum\_spanning\_tree.shp
 * Shortest path (optional - only if source/target shapefile specified on command line)
   * Calculates shortest path between given source and target vertices (identified in optional second shapefile specified on the command line) using edge weights
   * Output: graph/shortest\_path.shp
 * Betweenness centrality (optional)
   * Betweenness centrality (BC) is a measure of centrality in a (weighted) graph based on shortest paths between pairs of vertices
   * BC of a vertex = number of shortest paths between pairs of vertices that pass through that vertex
   * Output:

### Maximum flow

* The maximum flow algorithm determines which edges would contribute to flow in a given direction with a specified pressure gradient direction.
* The calculation is performed on a directed version of the undirected graph, where each edge in the undirected graph is replaced by 2 directed edges in the directed graph (one in each direction for each undirected edge).
* A fluid pressure is assigned to each vertex based on its position in a pressure gradient of a specified direction (i.e. fluid pressure decreasing towards top, bottom, left or right of a bounding box that encloses the entire fault network and is oriented parallel to the coordinate system). For example, if the fluid pressure is decreasing towards the top, vertices are assigned a pressure between 1 and 0 decreasing linearly with distance from bottom to top of the bounding box.
* Each edge in the graph has a flow direction and capacity:
  * The flow direction for a given edge is determined by the pressure at either end, whereby flow is always from higher to lower pressure.
  * The flow capacity for a given edge is based on its length, orientation, or a combination of length and orientation.
* Outputs:
  * A shapefile showing the amount of flow along each edge for a specified pressure gradient direction and flow direction: graph/max\_flow\_Gradient\_\<original input filname\>.shp (only written if maximum flow \> 0)
  * A matrix showing the total flow across the border of the bounding box for each combination of pressure gradient direction and flow direction: output to console. See [here](doc/FracG_outputs.md) for ordering of the matrix.

## Classification

* Based on orientation and number of intersections
* Output: classified.shp

## Raster analysis (optional)

* If a raster file (.tif) is specified on the command line the topological analyses (described above) will be performed only for the area contained within the raster.
* Faults that cross the borders of the raster are cut to the raster area and their ends marked as edge-nodes (E-nodes). Faults containing an edge node are excluded from the statistical analysis of fault length.
* General statistics of the raster data (e.g. minima, maxima, mean, standard deviation) are reported in a csv file together with the pixel values along strike and perpendicular to strike at the centre of the lines.
* Values in the raster file are used to determine edge and vertex weights. Raster data is collected at every vertex, at the centre of each edge, and along each edge.
  * The weights derived from the raster file can be used in some of the analyses described above. For instance, if the raster data is a digital elevation model (DEM), the path with the maximum topographic gradient can be determined based on the maximum flow algorithm.
* In addition, the gradient across the centre of the edges, the mean gradient along the segments (profile) and the gradient along the branches is calculated.
* Outputs:
  * Raster statistics: raster\_shp/raster\_augmented\_shapefile.shp and raster\_parallel\_cross\_profiles.csv
  * Raster-augmented graph: graph\_shape/raster/graph\_branches.shp and graph\_vertices.shp |

## Mesh generation (optional)

 * The fault network is used to generate a 2D conforming mesh using the open-source mesh generator Gmsh
 * The faults are defined as side sets which are meshed with lower-dimensional (1D) elements, and the intervening space is meshed with triangular elements that are refined around the faults.
   * Mesh refinement is a function of distance to the sidesets
 * Generation of a series of randomly sampled meshes representing rectangular windows within the modelled area
