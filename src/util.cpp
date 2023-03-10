/****************************************************************/
/*				DO NOT MODIFY THIS HEADER						*/
/*					FRACG - FRACture Graph						*/
/*				Network analysis and meshing software			*/
/*																*/
/*						(c) 2021 CSIRO							*/
/*			GNU General Public Licence version 3 (GPLv3)		*/
/*																*/
/*						Prepared by CSIRO						*/
/*																*/
/*					See license for full restrictions 			*/
/****************************************************************/
#include <iostream>
#include <boost/filesystem.hpp>
#include <boost/algorithm/string.hpp>
#include "../include/util.h"

namespace FracG
{
    namespace fs = boost::filesystem;
    namespace ba = boost::algorithm;
    
    //prepend a prefix string to a path name, while preserving the parent directory
    //append a suffix to a filename string, if that string does not already end with that suffix
    //optionally remove the existing file extension
    std::string AddPrefixSuffix(fs::path path, std::string prefix, std::string suffix, bool remove_extension)
    {
//         std::cout << "add_ps: path " << path << ", pfx: " << prefix << ", sfx " << suffix << ", re? " << remove_extension;
        if (remove_extension)
        {
            path = path.parent_path() / path.stem();
        }
        std::string fn = path.filename().string();
        std::string new_filename = prefix + (fn == "." ? "" : fn);
//         std::cout << ",pfx fn="<<new_filename;
        if (!ba::ends_with(new_filename, suffix))
        {
            new_filename += suffix;
        }
        fs::path parent_dir = path.parent_path();
        std::string output = (parent_dir / new_filename).string(); 
//         std::cout << ",ps fn="<< new_filename << ", result = " << output << std::endl;
        return output;
    }
    
    //as above, but also with an optional list of subdirectories to add
    std::string AddPrefixSuffixSubdirs(fs::path path, std::initializer_list<std::string> subdirs, std::string prefix, std::string suffix, bool remove_extension)
    {
        fs::path parent_dir = fs::is_directory(path) ? path : path.parent_path();
        for (auto subdir : subdirs) parent_dir = parent_dir / subdir;
        fs::path new_path = parent_dir / path.filename();
        return AddPrefixSuffix(new_path, prefix, suffix, remove_extension);
    }
    
    //as above, but with a string input
    std::string AddPrefixSuffix(std::string filename, std::string prefix, std::string suffix, bool remove_extension)
    {
        fs::path path(filename);
        return AddPrefixSuffix(path, prefix, suffix, remove_extension);
    }
    
    std::string AddPrefixSuffixSubdirs(std::string filename, std::initializer_list<std::string> subdirs, std::string prefix, std::string suffix, bool remove_extension)
    {
        fs::path path(filename);
        return AddPrefixSuffixSubdirs(path, subdirs, prefix, suffix, remove_extension);
    }
    
    
    //Create a directory (the file object itsefl, if it is a directory, or othe file's parent directory
    void CreateDir(fs::path file)
    {
        fs::path dir_path = fs::is_directory(file) ? file : file.parent_path();
        fs::create_directories(dir_path);
    }
    
    void CreateDir(std::string filename)
    {
        CreateDir(fs::path(filename));
    }
    
    //Create a filestream for writing to, and its parent directory(s) if necessary
    std::ofstream CreateFileStream(fs::path file)
    {
//         std::cout << "Creating file stream for " << file << std::endl;
        CreateDir(file);
        std::ofstream out_stream;
        out_stream.open(file.string(), std::ios::out);
        return out_stream;
    }
    
    std::ofstream CreateFileStream(std::string filename)
    {
        fs::path path(filename);
        return CreateFileStream(path);
    }
    
    //convert a string into a direction
    Direction ReadDirection(std::string direction_string)
    {
        if (direction_string.length() < 1)
        {
            std::cerr << "Error: direction string is empty" << std::endl;
            return Direction::NONE;
        }
        char c = std::tolower(direction_string[0]);
        Direction direction;
        switch(c)
        {
            case 'l': direction = Direction::LEFT; break;
            case 'r': direction = Direction::RIGHT; break;
            case 't': direction = Direction::TOP; break;
            case 'b': direction = Direction::BOTTOM; break;
            case 'n': direction = Direction::NONE; break;
            default: std::cerr << "Error: Invalid direction string \""<<direction_string<<"\" given" << std::endl;
            direction = Direction::NONE; break;
        }
        return direction;
    }
    
    double DirectionAngleDegrees(Direction d)
    {
        switch (d)
        {
            case 'l': return 180;
            case 'r': return 0;
            case 't': return 90;
            case 'b': return -90;
            default:
                std::cerr << "Error: None direction has no angle" << std::endl;
                return std::nan("");
        }
    }
    
    
	//convenience function to get the opposite of a given direction
	Direction OppositeDirection(const Direction d)
	{
		switch (d)
		{
			case LEFT: return Direction::RIGHT;
			case RIGHT: return Direction::LEFT;
			case TOP: return Direction::BOTTOM;
			case BOTTOM: return Direction::TOP;
			default: return Direction::NONE;
		}
	}
    

    
}
