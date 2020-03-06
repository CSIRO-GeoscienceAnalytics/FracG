#check for cmaek version adn set some compiler flags
cmake_minimum_required(VERSION 3.0)
set(CMAKE_BUILD_TYPE Debug)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14 -fopenmp ")

#Name of the Project and version
project(FracG VERSION 1.0.1 DESCRIPTION "Lineament Network Analysis")

#Check fro all required external libaries-------------------------------
find_package(GDAL REQUIRED)
include_directories( include ${GDAL_INCLUDE_DIRS})

find_package(Boost REQUIRED COMPONENTS random system filesystem)
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
target_link_libraries(FracG ${GDAL_LIBRARIES})
target_link_libraries(FracG ${BOOST_LIBRARIES})
target_link_libraries(FracG Boost::random)
target_link_libraries(FracG Boost::system)
target_link_libraries(FracG Boost::filesystem)
target_link_libraries(FracG ${ARMADILLO_LIBRARIES})
target_link_libraries(FracG ${GSL_LIBRARIES})
target_link_libraries(FracG ${OPENMP_LIBRARIES} )
target_link_libraries(FracG ${GMSH_LIB})



#install to global env (use sudo make install)
install(TARGETS FracG DESTINATION /usr/local/bin)
