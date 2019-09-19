/**
 * @file   helper_functions.cpp
 * @author Daniel Abler <daniel.abler@istb.unibe.ch>
 * @date   Feburay 2017
 * @version 1.0
 * @brief  Helper Functions
 */


#include "helper_functions.h"
// LOGGING
#include "easylogging++.h"


namespace hf {
  
// File access and management

std::string removeWhiteSpace(std::string mystring)
{
  mystring.erase(std::remove_if(mystring.begin(), mystring.end(), ::isspace), mystring.end());
  return mystring;
}


bool hasSuffix(const std::string& str, const std::string& suffix)
{
  return str.size() >= suffix.size() &&
         str.compare(str.size() - suffix.size(), suffix.size(), suffix) == 0;
}



bool existsDir(const std::string& path_to_directory)
{
struct stat info;
bool exists;

if( stat( path_to_directory.c_str(), &info ) != 0 )
{
  LOG(DEBUG) << "cannot access " << path_to_directory;
  exists = false;
}
else if( info.st_mode & S_IFDIR )  // S_ISDIR() doesn't exist on my windows
{
  LOG(DEBUG) << path_to_directory << " is a directory";
  exists = true;
}
else
{
  LOG(DEBUG) << path_to_directory << " is no directory";
  exists = false;
}

return exists;
}

std::string getDirFromPath(const std::string filename)
{
  std::string directory;
  const size_t last_slash_idx = filename.rfind('//');
  if (std::string::npos != last_slash_idx)
    {
      directory = filename.substr(0, last_slash_idx);
    }
  return directory;
}


void ensureDirExists(const std::string& path_to_directory)
{
  if (not existsDir(path_to_directory))
    {
      LOG(INFO) << "Creating directory '"<< path_to_directory;
      mkdir(path_to_directory.c_str(), ACCESSPERMS);
    }

}




bool existsFile(const std::string& path_to_file) {
// see http://techoverflow.net/blog/2013/08/21/how-to-check-if-file-exists-using-stat/
// and for a solution without 'sys/stat'
// see http://techoverflow.net/blog/2013/01/11/cpp-check-if-file-exists/
    struct stat buf;
    return (stat(path_to_file.c_str(), &buf) == 0);
}



}
