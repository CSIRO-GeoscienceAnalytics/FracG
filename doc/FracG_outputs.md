# FracG outputs

All outputs are written to **out_dir** (_default:_ fracg\_output\_ + \<name-of-shapefile\>) except the corrected shapefile, which is written to the input file directory.

| **Output** | **Description** | **Location** |
| --- | --- | --- |
| Corrected Shapefile | Shapefile in which faults closer than the threshold distance and within a threshold angle are merged | corrected\_\<original shapefile name\>.shp |
| Length distribution | Fitting parameters of the best-fit length distribution | statistics/length\_distributions.csv |
| Angle Distribution | Output of kernel density estimation of principal orientations | statistics/angle\_distribution\_kde.csv |
| Angle Distribution, weighted by length | Output of kernel density estimation of principal orientations, weighted by length | statistics/angle\_distribution\_kde\_length\_weights.csv |
| Fractal Dimension | Results of box counting. Fractal dimension is the slope of the line log N vs. log S. | statistics/box\_count.csv |
| Vector properties | general stats and properties for each vertex | statistics/vector\_properties.csv |
| Scanline statistics | Fault spacing and intensity along scanlines | statistics/scanline\_analysis.csv |
| Fault Density | Raster map of fault density | geometry/P20\_map.tif |
| Fault Intensity | Raster map of fault intensity | geometry/P21\_map.tif |
| Distance maps | Raster maps of shortest distance to nearest fault, and distance to centre of nearest fault | geometry/Distance\_map.tif and Center\_distance\_map.tif |
| Graph shapefile | Graph representation of fault network | graph\_shp/graph (can be changed using the option "graph\_results\_folder")/graph\_branches.shp and graph\_vertices.shp |
| Graph analysis | Parameters derived from the graph e.g. number and classification of edges, area of graph, average degree of vertices | statistics/graph\_statistics.csv |
| Fault classifications | Based on orientation and intersections | classified.shp |
| Intersection Density | Raster map of intersection density | geometry/interserction\_density.tif |
| Intersection Intensity | Raster map of intersection density, weighted by the degree of the intersection | geometry/intersection\_intensity.tif |
| Connected Components | Collections of edges/vertices that are connected with each other, different connected components cannot be reached from one another | graph/test.shp |
| Minimum Spanning Tree | Spanning tree in the graph with lowest sum of edge weights | graph/minimum\_spanning\_tree.shp |
| Shortest Path | Shortest path (using edge weights) between specified Source and Target vertices, if given | graph/shortest\_path.shp |
| Maximum Flow | Maximum flow in a specified direction (towards top/bottom/left/right) using a specified pressure gradient direction (pressure decreasing towards top/bottom/left/right) | graph/max\_flow\_Gradient\_\<original input filname\>.shp (only written if maximum flow \> 0) |
| Maximum Flow matrix | Maximum flow across the network for each combination of flow direction and pressure direction. See [below](#maximum-flow-matrix) for ordering of the flow matrix. | Printed out to the console |
| Raster statistics | Statistics of the raster file (if raster file is specified) | raster\_shp/raster\_augmented\_shapefile.shp and raster\_parallel\_cross\_profiles.csv |
| Raster-augmented graph | Graph augmented with edge and vertex weights derived from raster file (if raster file is specified) | graph\_shape/raster/graph\_branches.shp and graph\_vertices.shp |
| Mesh | Mesh(es) generated from fault network | mesh/2D\_mesh\_ |

## Maximum flow matrix

The following table indicates the order of components in the maximum flow matrix. Each entry in the matrix represents the "total flow" for a given combination of pressure and flow direction, which is the sum of flow along the border of the network.
Values a to f in the matrix may be non-zero. The directions left, right, top and bottom indicate flow towards or pressure decreasing towards that boundary. Hence, flow right with pressure left results in zero flow. Flow right with pressure top may result in non-zero flow depending on the orientation of the lineaments and their intersections with the edges. Swapping the flow and pressure directions gives the same maximum flow (e.g. flow right with pressure right is the same as flow left with pressure left).

|   | Pressure Right  | Pressure Top  | Pressure Left | Pressure Bottom |
| --- | --- | --- | --- | --- |
| Flow Right  | a | b  | 0 | c  |
| Flow Top  | d | e | f | 0 |
| Flow Left | 0 | c | a | b |
| Flow Bottom  | f | 0 | d | e |