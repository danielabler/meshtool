#ifndef IMAGEMESHTOOL_H_
#define IMAGEMESHTOOL_H_

#define BOOST_PARAMETER_MAX_ARITY 12
// CONFIG HANDLING
#include "xml-io/imaging_meshing_schema.hxx"
// HELPER FCTS
#include "helper_functions.h"
#include "helper_functions_vtk.h"
#include <map>


// VTK INCLUDES
#include <vtkSmartPointer.h>
#include <vtkUnstructuredGrid.h>
#include <vtkTransform.h>
#include <vtkTransformFilter.h>
#include <vtkThreshold.h>
#include <vtkDataSetSurfaceFilter.h>
#include <vtkImageData.h>

// CGAL INCLUDES
// --- cgal primitives
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
// --- geom structures
#include <CGAL/Polyhedron_3.h>
// --- cgal image handling
#include <CGAL/Image_3.h>
// --- cgal mesh handling
#include <CGAL/Labeled_image_mesh_domain_3.h>
#include <CGAL/Mesh_triangulation_3.h>
#include <CGAL/Mesh_complex_3_in_triangulation_3.h>
#include <CGAL/Mesh_criteria_3.h>
#include <CGAL/make_mesh_3.h>
// ------ CGAL -> VTK
#include "c3t3_subdomain_to_vtk_unstructured_grid.h"
// ------ VTK -> CGAL
#include "VTKPolyDataToCgalPolyhedron.h"


// CGAL TYPE DEFS

// --- Kernel & Basic Components
//typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
struct K: public CGAL::Exact_predicates_inexact_constructions_kernel {};
typedef CGAL::Polyhedron_3<K>         Polyhedron;
// --- Labeled Image: Domain, Triangulation, Criteria
typedef CGAL::Labeled_image_mesh_domain_3<CGAL::Image_3,K> Mesh_domain;
// for CGAL parallelisation support, see http://doc.cgal.org/latest/Mesh_3/index.html#title12
// seems to have no performance effect
typedef CGAL::Mesh_triangulation_3<
				Mesh_domain,
				CGAL::Kernel_traits<Mesh_domain>::Kernel, // Same as sequential
				CGAL::Parallel_tag                        // Tag to activate parallelism
				>::type Tr;
typedef CGAL::Mesh_complex_3_in_triangulation_3<Tr> C3t3;
typedef CGAL::Mesh_criteria_3<Tr> Mesh_criteria;
typedef CGAL::Mesh_constant_domain_field_3<Mesh_domain::R, Mesh_domain::Index> Sizing_field;

// To avoid verbose function and named parameters call
using namespace CGAL::parameters;


//! @brief imt : Namespace for MESHTOOL CLASS
namespace imt {


// DATA STRUCTURES
  
struct mesh_criteria_global_t {
  double cell_radius_edge_ratio;
  double cell_size;
  double facet_angle;
  double facet_distance;
  double facet_size;
} ;

struct mesh_criteria_subdomain_t {
  double cell_size;
  int dimension;
} ;  
  
  
typedef std::map<int, mesh_criteria_subdomain_t> mesh_criteria_subdomain_map_t;

/*!
* @class MeshTool meshtool.h 
* @brief Main simulator class that encapsulates all functions 
*/
class MeshTool 
{  
  
public:
  
  //! @brief Constructor
  MeshTool();
  
  //! @brief Destructor
  ~MeshTool();
    
  
  void updateConfigSettings(std::map< std::string, std::string > config_update);
  void setPathsToXSDs(std::string library, std::string path_to_xsd);
  
  /** 
    * @brief Reads configuration file following xml schema in xml-io/MeshTool_config/MeshTool_config.xsd.
    * 	Syntactically incorrect or incomplete xml will not be accepeted.
    * @param[in] path_to_config_file Path to filesystem where configuration file is loacted.
    */ 
  void initConfigFile(std::string path_to_config_file);

  void operationCreateMeshFromImage();
  
  void operationCreateMeshFromMesh();
  
  void performOperation();

private:
  void loadImage(std::string path_to_image_file);
  void createMeshFromLabeledImage(CGAL::Image_3 cgal_image);
  void createImageFromMesh(vtkSmartPointer<vtkUnstructuredGrid> vtk_ugrid);
  std::list<vtkSmartPointer<vtkPolyData> > vtkPolyDataListFromVtkUnstructuredGrid(vtkSmartPointer<vtkUnstructuredGrid> vtk_unstructured_grid);
  vtkSmartPointer<vtkPolyData> extractSurfaceFromVtkUnstructuredGrid(vtkSmartPointer<vtkUnstructuredGrid> vtk_unstructured_grid);
protected:  
  
  // configuration settings
  std::map<std::string, std::string> configuration_params;
  mesh_criteria_global_t mesh_criteria_global;
  mesh_criteria_subdomain_map_t mesh_criteria_subdomain_map;
  double spacing[3];
  // internal domain representations
  // Internal Variables
    vtkSmartPointer<vtkUnstructuredGrid> unstructured_grid_original;
    vtkSmartPointer<vtkUnstructuredGrid> unstructured_grid_processed;  
    C3t3 cgal_c3t3_from_image_domain;
    CGAL::Image_3 cgal_image_original;
    vtkSmartPointer<vtkImageData> vtk_image_original;
    
  
};
  
} // namespace mt

#endif  // IMAGEMESHTOOL_H_