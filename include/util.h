#ifndef fracg_util
#define fracg_util

#include <string>
#include <initializer_list>
#include <boost/filesystem.hpp>

namespace FGraph
{

    std::string add_prefix_suffix(boost::filesystem::path path, std::string prefix="", std::string suffix="", bool remove_extension=false);
    std::string add_prefix_suffix_subdirs(boost::filesystem::path path, std::initializer_list<std::string> subdirs, std::string prefix="", std::string suffix="", bool remove_extension=false);
    
    std::string add_prefix_suffix(std::string filename, std::string prefix="", std::string suffix="", bool remove_extension=false);
    std::string add_prefix_suffix_subdirs(std::string path, std::initializer_list<std::string> subdirs, std::string prefix="", std::string suffix="", bool remove_extension=false);
    
    void CreateDir(boost::filesystem::path file);
    void CreateDir(std::string filename);
    
    std::ofstream CreateFileStream(boost::filesystem::path file);
    std::ofstream CreateFileStream(std::string filename);
    
}

#endif
