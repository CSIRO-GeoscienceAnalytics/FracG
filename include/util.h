#ifndef fracg_util
#define fracg_util

#include <string>
#include <initializer_list>
#include <boost/filesystem.hpp>

namespace FGraph
{

    std::string add_prefix_suffix(std::string path, std::string prefix="", std::string suffix="", bool remove_extension=false);
    std::string add_prefix_suffix(std::string path, std::initializer_list<std::string> subdirs, std::string prefix="", std::string suffix="", bool remove_extension=false);
    
}

#endif
