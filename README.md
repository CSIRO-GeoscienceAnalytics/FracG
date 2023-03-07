# FracG - Fault and fracture analysis and meshing software

FracG is a command line Debian GNU/Linux based software that performs analysis of discontinuity data (e.g.
faults and fractures). Currently, the framework is only tested for Ubunutu 20.04 LTS.
The input need to be provided as line-vector data in shape-file format. Prior to every analysis FracG will try
to clean up the given vector data removing duplicate entries and 
aws in the topology. The analysis can be
performed with or without obtaining the best-fit models for length distribution, meshing or including raster
data. The generic analysis are:
* Models of the length distribution (exponential, log-normal, or power-law) and of the principal orientations
(Gaussians fitted to the Kernel density estimation) are obtained by automated fitting
* Statistical parameters including scan line analysis are exported (.csv)
* Box-counting to obtain fractal dimension
* Density maps are generated (frequency, intersections, distances)
* Georeferenced graph is created
* Graph algorithms:
    * Betweenness Centrality
    * Minimum Spanning tree
    * Shortest path
    * Maximum Flow
* Classify discontinues (orientation, intersection numbers)
* Analysis can be extended by including raster data (e.g. DEM)
    * Extract raster values for line geometries
    * Build georeferenced graph including raster data as edge and vertex weights
* 2D Finite Element(FE) conforming meshes (discontinuities as side sets)
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
Note that gmsh's open cascade engine is used for the meshing..

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