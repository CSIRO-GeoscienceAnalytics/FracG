#****************************************************************/
#*				DO NOT MODIFY THIS HEADER							*/
#*					FRACG - FRACture Graph							*/
#*				Network analysis and meshing software					*/
#*																		*/
#*						(c) 2021 CSIRO									*/
#*			GNU General Public Licence version 3 (GPLv3)				*/
#*																		*/
#*						Prepared by CSIRO								*/
#*																		*/
#*					See license for full restrictions 						*/
#****************************************************************/
cmake_minimum_required(VERSION 3.9)
set(CMAKE_BUILD_TYPE Release)
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++17 -fopenmp ")
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

#Name of the Project and version
project(FracG VERSION 1.0.0 DESCRIPTION "FracG")

#Check fro all required external libaries-------------------------------
find_package(GDAL REQUIRED)
include_directories( include ${GDAL_INCLUDE_DIRS})

find_package(Boost REQUIRED COMPONENTS random system filesystem program_options timer)
include_directories( include ${BOOST_INCLUDE_DIRS})

find_package(Armadillo REQUIRED)
include_directories(${ARMADILLO_INCLUDE_DIRS})

find_package(GSL REQUIRED)
include_directories(${GSL_INCLUDE_DIRS})

find_package(OpenMP REQUIRED)
include_directories(${OPENMP_INCLUDE_DIRS})

#add gmsh API-----------------------------------------------------------
find_library(GMSH_LIB gmsh)
if(NOT GMSH_LIB)
  message(FATAL_ERROR "Could not find libgmsh")
endif()

find_path(GMSH_INC gmsh.h)
if(NOT GMSH_INC)
  message(FATAL_ERROR "Could not find gmsh.h")
endif()
include_directories(${GMSH_INC})

#include cpp and h-files------------------------------------------------
include_directories(include)
file(GLOB SOURCES "src/*.cpp")

#excetutable
add_executable(FracG ${SOURCES}) 

#import the external libaries we need 
target_link_libraries(FracG PRIVATE ${GDAL_LIBRARIES})
target_link_libraries(FracG PRIVATE ${BOOST_LIBRARIES})
target_link_libraries(FracG PRIVATE Boost::random)
target_link_libraries(FracG PRIVATE Boost::system)
target_link_libraries(FracG PRIVATE Boost::timer)
target_link_libraries(FracG PRIVATE Boost::filesystem)
target_link_libraries(FracG PRIVATE Boost::program_options)
target_link_libraries(FracG PRIVATE ${ARMADILLO_LIBRARIES})
target_link_libraries(FracG PRIVATE ${GSL_LIBRARIES})
target_link_libraries(FracG PRIVATE OpenMP::OpenMP_CXX)
target_link_libraries(FracG PRIVATE ${GMSH_LIB})

target_compile_options(FracG PRIVATE "${OpenMP_CXX_FLAGS}")

#install to global env (use sudo make install)
install(TARGETS FracG DESTINATION /usr/local/bin)
