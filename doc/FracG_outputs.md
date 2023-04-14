# FracG outputs

| **Output** | **Description** | **Location** |
| --- | --- | --- |
| Length distribution | | out\_dir / statistics / length\_distributions.csv |
| Box counting/Fractal Dimension | | out\_dir / statistics / bo\_count.csv |
| Angle Distribution | | out\_dir / statistics / angle\_distribution\_kde.csv |
| Angle Distribution, weighted by length | | out\_dir / statistics / angle\_distribution\_kde\_length\_weights.csv |
| vector properties | general stats and properties for each vertex | out\_dir / statistics / vector\_properties.csv |
| Scanline statistics | | out\_dir / statistics / scanline\_analysis.csv |
| Corrected Shapefile | merged shapefile, based on distance threshold - anything closer than the threshold distance is merged together | input path / corrected\_\<original shapefile name\> |
| P20 point density and P21 length density maps | | out\_dir / geometry / P20\_map.tif and P21\_map.tif |
| Distance maps | distance to nearest point, and to nearest centre of a feature | out\_dir / geometry / Distance\_map.tif and Center\_distance\_map.tif |
| Graph properties analysis | | out\_dir / statistics / graph\_statistics.csv |
| Graph shapefile | writing the graph to a file | out\_dir / graph\_shp / graph (can be changed using the option "graph\_results\_folder") / graph\_branches.shp and graph\_vertices.shp |
| Lineament Classifications | based on orientation and intersections | out\_dir / classified.shp |
| Intersection Density | number of intersections by area | out\_dir / geometry / interserction\_density.tif |
| Intersection Intensity | number of intersections per area, weighted by the degree of the vertex | out\_dir / geometry / intersection\_intensity.tif |
| Connected Components | collections of edges/vertices that are connected with each other, different connected components cannot be reached from one another | out\_dir / graph / test.shp |
| Minimum Spanning Tree | | out\_dir / graph / minimum\_spanning\_tree.shp |
| Shortest Path | between given Source and Target vertices, if given | out\_dir / graph / shortest\_path.shp |
| Maximum Flow | Maximum flow in a specified direction (max\_flow\_gradient\_flow\_direction) using a specified pressure gradient direction (max\_flow\_gradient\_pressure\_direction) to limit which edges can carry flow in which direction. These directions indicate flow towards or pressure decreasing towards "left", "right", "top", or "bottom". It identifies the edges of the graph as anything intersecting a box that is a shrunk-down version of the network's overall bounding box. The extent to which is it shrunk is given by max\_flow\_gradient\_border\_amount, with a value from 0 (no shrinking) to 1 (effectively a point at the centre of the graph). | out\_dir / graph / max\_flow\_Gradient\_\<original input filname\>.shp (only written if maximum flow \> 0) |
| Maximum Flow matrix | Maximum flow across the network for each combination of flow direction and pressure direction. See [below](#maximum-flow-matrix) for ordering of the flow matrix. | Printed out to the console |
| Raster statistics if raster file is specified | | out\_dir / raster\_shp / raster\_augmented\_shapefile.shp and raster\_parallel\_cross\_profiles.csv |
| Raster-augmented graph if raster file is specified | | out\_dir / graph\_shape / raster / graph\_branches.shp and graph\_vertices.shp |
| Mesh | | out\_dir / mesh / 2D\_mesh\_ |

## Maximum flow matrix

The following table indicates the order of components in the maximum flow matrix. Entries a to f in the matrix may be non-zero. The directions left, right, top and bottom indicate flow towards or pressure decreasing towards that boundary of the graph. Hence, flow right with pressure left results in zero flow. Flow right with pressure top may result in non-zero flow depending on the orientation of the lineaments and their intersections with the edges. Swapping the flow and pressure directions gives the same maximum flow (e.g. flow right with pressure right is the same as flow left with pressure left).

|   | Pressure Right  | Pressure Top  | Pressure Left | Pressure Bottom |
| --- | --- | --- | --- | --- |
| Flow Right  | a | b  | 0 | c  |
| Flow Top  | d | e | f | 0 |
| Flow Left | 0 | c | a | b |
| Flow Bottom  | f | 0 | d | e |