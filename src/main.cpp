/**
 * @file main.cpp
 * @author Daniel Abler <daniel.abler@istb.unibe.ch>
 * @date   Feburary May 2017
 * @brief  Commandline / config file handling for cornea meshing
 */

// STANDARD
#include <iostream>
#include <unistd.h>
// LOGGING
#include "easylogging++.h"
INITIALIZE_EASYLOGGINGPP

// MESHTOOL CLASS
#include "cornea_meshtool.h"
#include "image_meshtool.h"


void print_usage() {
  std::cout << "Usage (default settings): 	      meshtool -i path_to_input_file -o path_to_output_file -m mm|im|mi" << std::endl;
  std::cout << "Usage (with configuration file):  meshtool -c path_to_config_file \n" << std::endl;
}



/// GET PATH TO EXECUTABLE
std::string get_selfpath() {
    char buff[PATH_MAX];
    ssize_t len = ::readlink("/proc/self/exe", buff, sizeof(buff)-1);
    if (len != -1) {
        buff[len] = '\0';
        return std::string(buff);
    }
    /* handle error condition */
}


int main(int argc, char **argv) {
    // LOGGING SETTINGS
    // Load configuration from file
    el::Configurations conf("log_config.conf");
    // Reconfigure all loggers -- Now all the loggers will use configuration from file
    el::Loggers::reconfigureAllLoggers(conf);

    /// CONFIGURATION
    /// Parsing of command line arguments
    int option;
    bool c_flag = false;
    bool x_flag = false;
    bool m_flag = false;
    bool i_flag = false;
    bool o_flag = false;
    std::string path_to_config_file= "";
    std::string path_to_xsd= "";
    std::string mesh_tool_mode= "";
    std::string path_to_input_file= "";
    std::string path_to_output_file= "";
    /// Specifying the expected options
    while ((option = getopt(argc, argv,"c:x:m:i:o")) != -1) {
        switch (option) {
            case 'c' :
                c_flag = true;
                path_to_config_file = optarg;
                break;
            case 'x' :
                x_flag = true;
                path_to_xsd = optarg;
                break;
            case 'm' :
                m_flag = true;
                mesh_tool_mode = optarg;
                break;
            case 'i' :
                i_flag = true;
                path_to_input_file = optarg;
                break;
            case 'o' :
                o_flag = true;
                path_to_output_file = optarg;
                break;
            default :
                print_usage();
                LOG(FATAL) << "No parameters given -- Insufficient information to start programme. Exiting programme...";
                    exit(EXIT_FAILURE);
        }
    }

    LOG(DEBUG) << " Command line arguments: '-c' = " << c_flag << ", '-x' = " << x_flag<< ", '-m' = " << m_flag;



    /// CATCH CURRENT EXECUTION DIRECTORY (IF NOT SPECIFIED)
    LOG(DEBUG) << "Getting current execution directory...";
    std::string selfpath  = get_selfpath();
    std::string base_path = "";

    if (not selfpath.empty()){
        std::size_t found = selfpath.find_last_of("/\\");
        std::string exec_file = selfpath.substr(found+1);
        base_path = selfpath.substr(0, found+1);
        LOG(INFO) <<  "-- Current execution directory: "<< base_path;
        LOG(INFO) <<  "-- Current execution process  : "<< exec_file;
    }
    else {

        LOG(FATAL)  << "Cannot retrieve current working directory, please specify: flag -b <path-to-executable-directory>";
        exit(EXIT_FAILURE);
    }


    if (mesh_tool_mode == "image")
    {
        LOG(INFO) << "Running MeshTool in 'image' mode";

        LOG(INFO) << "Instantiating MeshTool instance";
        imt::MeshTool imageMeshTool;
        if (x_flag == true)
        {
            LOG(DEBUG) << "XSD schema location given as parameter ... updating settings";
            imageMeshTool.setPathsToXSDs("ImageMeshTool", path_to_xsd);
        }
        imageMeshTool.initConfigFile(path_to_config_file);


        /// create MeshTool parameter map (used to update settings dynamically)
        LOG(DEBUG) << "Updating MeshTool parameters...";
        std::map<std::string, std::string> MeshTool_params;
        if (i_flag == true)
        {
            LOG(DEBUG) << "Updating input file settings...";
            MeshTool_params["path_to_input_file"] = path_to_input_file;
        }
        if (o_flag == true)
        {
            LOG(DEBUG) << "Updating output file settings...";
            MeshTool_params["path_to_output_file"] = path_to_output_file;
        }
        imageMeshTool.updateConfigSettings(MeshTool_params);

        imageMeshTool.performOperation();

    }
    else if (mesh_tool_mode == "cornea")
    {
        LOG(INFO) << "Running MeshTool in 'cornea' mode";

        // INSTANTIATION
        LOG(DEBUG) << "Instantiating CorneaMeshTool";
        cmt::CorneaMeshTool corneaMeshTool;
        if (x_flag == true)
        {
            LOG(INFO) << "XSD schema location given as parameter ... updating settings";
            corneaMeshTool.setPathsToXSDs("CorneaMeshTool", path_to_xsd);
        }

        // INITIALISATION OF CONFIG FILE
        corneaMeshTool.initConfigFile(path_to_config_file);

        // PARAMETERS TO BE OVERWRITTEN

        /*
        if (o_flag == true)
          {
            LOG(INFO) << "Output path given as parameter ... updating settings";
            cornea_generation_criteria_global.path_to_output_file = path_to_output_file;
          }
        */

        // START MESHING PROCESS
        corneaMeshTool.meshCornea();
    }
    else{
        LOG(FATAL) << "Please specify MeshTool model '-m', as either 'image' or 'cornea'.";
    }


    return 0;
}
