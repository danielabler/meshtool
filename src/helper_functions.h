#ifndef HELPERFCTS_H_
#define HELPERFCTS_H_

#include <string>
#include <iostream>
#include <sys/types.h>
#include <sys/stat.h>


//! @brief hf : Namespace for HELPERFCTS
namespace hf {


  // File access and management  
  
  std::string removeWhiteSpace(std::string mystring);
    bool hasSuffix(const std::string& str, const std::string& suffix);
  bool existsDir(const std::string& path_to_directory);
  void ensureDirExists(const std::string& path_to_directory);
  std::string getDirFromPath(const std::string path_to_directory);
  bool existsFile(const std::string& path_to_file);


}

#endif 
