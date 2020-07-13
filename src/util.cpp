
#include <iostream>

#include <boost/filesystem.hpp>
#include <boost/algorithm/string.hpp>
#include "../include/util.h"

namespace FGraph
{
    namespace fs = boost::filesystem;
    namespace ba = boost::algorithm;
    
    //prepend a prefix string to a path name, while preserving the parent directory
    //append a suffix to a filename string, if that string does not already end with that suffix
    //optionally remove the existing file extension
    std::string add_prefix_suffix(std::string path, std::string prefix, std::string suffix, bool remove_extension)
    {
        fs::path base_path(path);
        if (remove_extension)
        {
            base_path = base_path.parent_path() / base_path.stem();
        }
        std::string new_filename = prefix + base_path.filename().string();
        if (!ba::ends_with(path, suffix))
        {
            new_filename += suffix;
        }
        fs::path parent_dir = base_path.parent_path();
        std::string output = (parent_dir / new_filename).string(); 
        return output;
    }
    
    //as above, but also with an optional list of subdirectories to add
    std::string add_prefix_suffix(std::string path, std::initializer_list<std::string> subdirs, std::string prefix, std::string suffix, bool remove_extension)
    {
        fs::path base(path);
        fs::path parent_dir = base.parent_path();
        for (auto subdir : subdirs) parent_dir = parent_dir / subdir;
        std::string new_path = (parent_dir / base.filename()).string();
        return add_prefix_suffix(new_path, prefix, suffix, remove_extension);
    }
}
