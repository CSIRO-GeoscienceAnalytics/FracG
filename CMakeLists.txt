#check for cmaek version adn set some compiler flags
cmake_minimum_required(VERSION 3.0)
set(CMAKE_BUILD_TYPE Release)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14 -fopenmp")

#Name of the Project and version
project(FracG VERSION 1.0.1 DESCRIPTION "Lineament Network Analysis")

#Check fro all required external libaries
find_package(GDAL REQUIRED)
include_directories( include ${GDAL_INCLUDE_DIRS})

find_package(Boost REQUIRED COMPONENTS random)
include_directories( include ${BOOST_INCLUDE_DIRS})

find_package(Armadillo REQUIRED)
include_directories(${ARMADILLO_INCLUDE_DIRS})

find_package(OpenMP REQUIRED)
include_directories(${OPENMP_INCLUDE_DIRS})

#include cpp and h-files
include_directories(include)
file(GLOB SOURCES "src/*.cpp")

#excetutable
add_executable(FracG ${SOURCES}) 

#import the external libaries we need 
target_link_libraries(FracG ${GDAL_LIBRARIES})
target_link_libraries(FracG ${BOOST_LIBRARIES})
target_link_libraries(FracG Boost::random)
target_link_libraries(FracG ${ARMADILLO_LIBRARIES})
target_link_libraries(FracG ${OPENMP_LIBRARIES} )

#install to global env (use sudo make install)
install(TARGETS FracG DESTINATION /usr/local/bin)

