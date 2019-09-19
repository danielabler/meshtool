# Include the directory itself as a path to include directories
set(CMAKE_INCLUDE_CURRENT_DIR ON)


# SETTINGS
# -----------------------
# Use Relative paths for linking to libraries
# ------------------------
SET(CMAKE_SKIP_BUILD_RPATH  FALSE)
SET(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE) 
# link libraries relative to bin location
SET(CMAKE_INSTALL_RPATH "$ORIGIN/")  
SET(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
#-----------------------	
# EXTERNAL LIBRARIES & INCLUDES

#IF(MACOS)
#    MESSAGE(STATUS "Path defintions for MacOS")
#    SET(DEPENDENCY_DIR /Users/dabler/usr_libs)
#    SET(CGAL_DIR ${DEPENDENCY_DIR}/CGAL-4.8.1_build_new)
#    SET(VTK_DIR ${DEPENDENCY_DIR}/VTK-7.1.1_build_new)
#    SET(XERCESC_ROOT_DIR ${DEPENDENCY_DIR}/xerces-c-3.2.0_install)
#    SET(XSD_DIR ${DEPENDENCY_DIR}/xsd-4.0.0-x86_64-linux-gnu/libxsd)
#ELSEIF(LINUX)
#    MESSAGE(STATUS "Path defintions for Linux")
#    SET(DEPENDENCY_DIR /home/dabler/software/binaries/)
#    SET(CGAL_DIR ${DEPENDENCY_DIR}/CGAL_build)
#    SET(VTK_DIR ${DEPENDENCY_DIR}/VTK_build)
#    SET(XERCESC_ROOT_DIR ${DEPENDENCY_DIR}/XERCES_install)
#    SET(XSD_DIR ${DEPENDENCY_DIR}/xsd-4.0.0-x86_64-linux-gnu/libxsd)
#    SET(BOOST_ROOT ${DEPENDENCY_DIR}/BOOST_source)
#ENDIF()

# avoid in-source builds
# -----------------------
set(CMAKE_DISABLE_IN_SOURCE_BUILD ON)
set(CMAKE_DISABLE_SOURCE_CHANGES  ON)
if ("${CMAKE_SOURCE_DIR}" STREQUAL "${CMAKE_BINARY_DIR}")
  message(SEND_ERROR "In-source builds are not allowed.")
endif ()


# COMPILER SETTINGS
# -----------------------
if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU" OR
    "${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
    set(CMAKE_CXX_COMPILER "g++")
    set(warnings "-Wall") # 2)
    add_definitions(-std=c++11)
endif()


# Store flags in cache but ensure that they are added only once.
if (NOT CONFIGURED_ONCE)
    set(CMAKE_CXX_FLAGS "${warnings}"
        CACHE STRING "Flags used by the compiler during all build types." FORCE)
    set(CMAKE_C_FLAGS   "${warnings}"
        CACHE STRING "Flags used by the compiler during all build types." FORCE)
endif()



# EXTERNAL LIBRARIES
# -----------------------

# --- CGAL 
find_package(CGAL COMPONENTS ImageIO)
IF(CGAL_FOUND)
    INCLUDE(${CGAL_USE_FILE})
    MESSAGE(STATUS "Found CGAL: ${CGAL_LIBRARIES}")
    MESSAGE(STATUS "Found CGAL 3RD PARTY LIBRARIES: ${CGAL_3RD_PARTY_LIBRARIES}")
    MESSAGE(STATUS "Found CGAL_ImageIO_3RD_PARTY_DEFINITIONS: ${CGAL_ImageIO_3RD_PARTY_DEFINITIONS}")
ELSE(CGAL_FOUND)
    MESSAGE(FATAL_ERROR "CGAL not found. Please set CGAL_DIR.")
ENDIF(CGAL_FOUND)


#=== TBB - for CGAL parallelisation support, see http://doc.cgal.org/latest/Mesh_3/index.html#title12
# find_package(TBB)
# if(NOT TBB_FOUND)
#         MESSAGE (STATUS "TBB not found.")
# endif(NOT TBB_FOUND)

# --- VTK 
find_package(VTK REQUIRED NO_MODULE)
IF(VTK_FOUND)
    include(${VTK_USE_FILE})
    MESSAGE(STATUS "Found VTK:")
ELSE(VTK_FOUND)
    MESSAGE(FATAL_ERROR "VTK not found..")
ENDIF(VTK_FOUND)

# --- XERCESC
find_library(XercesLibrary NAMES xerces-c PATHS ${XERCESC_ROOT_DIR}/lib)
if(NOT XercesLibrary)
  message(FATAL_ERROR "Failed to find the Xerces library.")
endif()


# SOURCE FILES
# -----------------------

# --- HEADERS

set(${PROJECT_NAME}_HDRS
    helper_functions.h
    cornea_meshtool.h
    implicit_functions.h
    zernike.h
    c3t3_subdomain_to_vtk_unstructured_grid.h
    image_meshtool.h
    VTKPolyDataToCgalPolyhedron.h
    helper_functions_vtk.h
)

# --- SRC

set(${PROJECT_NAME}_SRC
    main.cpp
    helper_functions.cpp
    cornea_meshtool.cpp
    xml-io/cornea_meshing_schema.cxx
    xml-io/imaging_meshing_schema.cxx
    #meshing_implicit_single.cpp
    zernike.cpp
    implicit_functions.cpp
    image_meshtool.cpp
    helper_functions_vtk.cpp
)


# LIBRARY TARGETS
# -----------------------
add_library(${PROJECT_NAME}_lib ${${PROJECT_NAME}_SRC})


# EXECUTABLE TARGETS
# -----------------------

add_executable(${PROJECT_NAME} main.cpp)


# INCLUDES
# -----------------------

target_include_directories(
      ${PROJECT_NAME}_lib PUBLIC  ${CMAKE_CURRENT_SOURCE_DIR}/../external/include
        ${XERCESC_ROOT_DIR}/include
      ${XSD_DIR}
      )

#include_directories ("/home/abler/bms_dependencies/CGAL-4.9_install/include")

# LINKING
# -----------------------

target_link_libraries(
 	${PROJECT_NAME}_lib
       	${VTK_LIBRARIES}
 	${CGAL_LIBRARIES} 
	${CGAL_3RD_PARTY_LIBRARIES} 
	${XercesLibrary}
	#${TBB_LIBRARY}		 # for CGAL parallelisation support, see http://doc.cgal.org/latest/Mesh_3/index.html#title12
 	)	

target_link_libraries(
 	${PROJECT_NAME}
 	${PROJECT_NAME}_lib
 	)	



# POST INSTALLATION PATH FIXING ON MAC OS
# ----------------------

if (MACOS)
    get_filename_component(LIBNAME ${XercesLibrary} NAME)
    get_filename_component(LIBPATH ${XercesLibrary} PATH)
    message(${LIBNAME} ${XercesLibrary} ${LIBPATH})
    #execute_process(COMMAND /usr/bin/install_name_tool -change "@rpath/${LIBNAME}" "${LIBPATH}" "${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}")
    #execute_process(COMMAND /usr/bin/install_name_tool -change "@rpath/libxerces-c-3.2.dylib" "${LIBPATH}/libxerces-c-3.2.dylib" "${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}")

    add_custom_command(TARGET ${PROJECT_NAME}
            POST_BUILD COMMAND
            ${CMAKE_INSTALL_NAME_TOOL} -change "@rpath/libxerces-c-3.2.dylib" "${LIBPATH}/libxerces-c-3.2.dylib"
                                                                "${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}")
endif()






 	
# INSTALLATION
# -----------------------

install(TARGETS ${PROJECT_NAME}_lib DESTINATION ${LIB_INSTALL_DIR})
install(TARGETS ${PROJECT_NAME} DESTINATION ${BIN_INSTALL_DIR})
install(FILES ${${PROJECT_NAME}_HDRS} DESTINATION ${INCLUDE_INSTALL_DIR})

#-----------------------------------------------------------------------	

set(CONFIGURED_ONCE TRUE CACHE INTERNAL
    "A flag showing that CMake has configured at least once.")


# GATHERING SHARED LIBRARIES, RELINK
# -----------------------

#if (LINUX)
#    set(APPS ${BIN_INSTALL_DIR}/${PROJECT_NAME})  # paths to executables
#    list(APPEND DIRS ${VTK_USE_FILE}
#            ${CGAL_USE_FILE}
#            ${XercesLibrary}
#            ${CGAL_LIBRARIES}
#            ${CGAL_3RD_PARTY_LIBRARIES} )
#
#    INSTALL(CODE "
#    include(BundleUtilities)
#    fixup_bundle(\"${APPS}\"   \"\"   \"${DIRS}\")
#    verify_app( \"\${APPS}\" )
#    " COMPONENT Runtime)
#endif()


