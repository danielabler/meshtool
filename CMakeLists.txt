cmake_minimum_required(VERSION 2.8)

# about this project
# -----------------------
PROJECT(MeshTool)

SET(${PROJECT_NAME}_VERSION_MAJOR 0)
SET(${PROJECT_NAME}_VERSION_MINOR 0)
SET(${PROJECT_NAME}_VERSION_PATCH 0)
SET(${PROJECT_NAME}_VERSION "${${PROJECT_NAME}_VERSION_MAJOR}.${${PROJECT_NAME}_VERSION_MINOR}.${${PROJECT_NAME}_VERSION_PATCH}")
SET(CHIC_BMS_COPYRIGHT_YEAR "2015")

MESSAGE("Project: ${PROJECT_NAME}, Parent Directory path: ${CMAKE_SOURCE_DIR}")


# -----------------------
# CHECK PLATFORM
#------------------------
if (CMAKE_HOST_WIN32)
    MESSAGE(STATUS "Compiling on Windows")
    set(WINDOWS 1)
elseif(CMAKE_HOST_APPLE)
    MESSAGE(STATUS "Compiling on MacOS")
    set(MACOS 1)
elseif(CMAKE_HOST_UNIX)
    MESSAGE(STATUS "Compiling on Linux")
    set(LINUX 1)
endif()


# SETTINGS
# -----------------------
# Installation Dirs (relative to source dir)
#------------------------
set(INSTALL_LIB_DIR lib CACHE PATH "Installation directory for libraries")
set(INSTALL_BIN_DIR bin CACHE PATH "Installation directory for executables")
set(INSTALL_INCLUDE_DIR include CACHE PATH "Installation directory for header files")
#------------------------
# DEPENDENCIES
set(CGAL_DIR cgal_build CACHE PATH "Build dir CGAL")
set(VTK_DIR vtk_build CACHE PATH "Build dir VTK")
set(XERCESC_ROOT_DIR xerces_root CACHE PATH "Install dir Xerces")
set(XSD_DIR xsd_lib CACHE PATH "libxsd dir XSD")
MESSAGE(STATUS "CGAL_DIR:  " ${CGAL_DIR})
MESSAGE(STATUS "VTK_DIR:  " ${VTK_DIR})
MESSAGE(STATUS "XERCESC_ROOT_DIR:  " ${XERCESC_ROOT_DIR})
MESSAGE(STATUS "XSD_DIR:  " ${XSD_DIR})
#------------------------
# PATH SETTINGS
SET(BIN_INSTALL_DIR ${CMAKE_SOURCE_DIR}/${INSTALL_BIN_DIR})
SET(LIB_INSTALL_DIR ${CMAKE_SOURCE_DIR}/${INSTALL_LIB_DIR})
SET(INCLUDE_INSTALL_DIR ${CMAKE_SOURCE_DIR}/${INSTALL_INCLUDE_DIR})


# set env
# -----------------------
include_directories(${INCLUDE_INSTALL_DIR} ${PROJECT_SOURCE_DIR}/external/include)
link_directories(${LIB_INSTALL_DIR} ${PROJECT_SOURCE_DIR}/external/lib)


# Appends the cmake/modules path inside the MAKE_MODULE_PATH variable which stores the
# directories of additional CMake modules 
list(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake/modules)


# add subfolders
# -----------------------
add_subdirectory(src) 
# add_subdirectory(src-doc) 
# add_subdirectory(src-test) 


# configured flag
# -----------------------
set(CONFIGURED_ONCE TRUE CACHE INTERNAL
    "A flag showing that CMake has configured at least once.")
