/****************************************************************/
/*				DO NOT MODIFY THIS HEADER							*/
/*					FRACG - FRACture Graph							*/
/*				Network analysis and meshing software					*/
/*																		*/
/*						(c) 2021 CSIRO									*/
/*			GNU General Public Licence version 3 (GPLv3)				*/
/*																		*/
/*						Prepared by CSIRO								*/
/*																		*/
/*					See license for full restrictions 						*/
/****************************************************************/
#ifndef fracg_util
#define fracg_util

#include <string>
#include <initializer_list>
#include <boost/filesystem.hpp>
#include "fracg.h"

namespace FracG
{
    std::string AddPrefixSuffix(boost::filesystem::path path, std::string prefix="", std::string suffix="", bool remove_extension=false);
    std::string AddPrefixSuffixSubdirs(boost::filesystem::path path, std::initializer_list<std::string> subdirs, std::string prefix="", std::string suffix="", bool remove_extension=false);
    
    std::string AddPrefixSuffix(std::string filename, std::string prefix="", std::string suffix="", bool remove_extension=false);
    std::string AddPrefixSuffixSubdirs(std::string path, std::initializer_list<std::string> subdirs, std::string prefix="", std::string suffix="", bool remove_extension=false);
    
    void CreateDir(boost::filesystem::path file);
    void CreateDir(std::string filename);
    
    std::ofstream CreateFileStream(boost::filesystem::path file);
    std::ofstream CreateFileStream(std::string filename);
    
    Direction ReadDirection(std::string direction_string);
    double DirectionAngleDegrees(Direction d);
}
#endif
