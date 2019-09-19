#ifndef CORNEAMESHTOOL_H_
#define CORNEAMESHTOOL_H_

// CONFIG HANDLING
#include "xml-io/cornea_meshing_schema.hxx"
// HELPER FCTS
#include "helper_functions.h"
#include "helper_functions_vtk.h"
#include "zernike.h"
// CGAL
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Mesh_triangulation_3.h>
#include <CGAL/Mesh_complex_3_in_triangulation_3.h>
#include <CGAL/Mesh_criteria_3.h>
#include <CGAL/Implicit_to_labeling_function_wrapper.h>
#include <CGAL/Implicit_mesh_domain_3.h>
#include <CGAL/Labeled_mesh_domain_3.h>
#include <CGAL/make_mesh_3.h>
#include <CGAL/Mesh_domain_with_polyline_features_3.h>
#include <CGAL/IO/File_medit.h>
#include "implicit_functions.h"

#include <vtkUnstructuredGrid.h>
#include <vtkSmartPointer.h>
#include <vtkIdFilter.h>
#include <vtkDataSetSurfaceFilter.h>
#include <vtkPointData.h>
#include <vtkDataArray.h>
#include <vtkUnstructuredGrid.h>
#include <vtkPolyData.h>
// ------ CGAL -> VTK
#include "c3t3_subdomain_to_vtk_unstructured_grid.h"


//! @brief hf : Namespace for HELPERFCTS
namespace cmt {

  enum boundary_type {cone_boundary, cylinder_boundary};

  struct cornea_boundary_t{
    boundary_type type;
    double angle;
    double radius;
  };

  struct boundary_cone_t{
    double z_apex;
    double alpha_rad;
  };

  struct mesh_criteria_global_t {
    double cell_radius_edge_ratio;
    double cell_size;
    double facet_angle;
    double facet_distance;
    double facet_size;
  } ;

  struct cornea_generation_criteria_global_t {
    double max_number_zernike_coeffs;
    double surface_thickness;
    double offset_along_z;
    std::string path_to_output;
  } ;


    struct lenticule_generation_criteria_global_t {
        double max_number_zernike_coeffs;
        double lenticule_radius;
        double cap_thickness;
        double lenticule_thickness;
        double lenticule_surface_distance;

    } ;

  struct surface_group {
    int surface_first;
    int rim_first;
  }; 


  struct RetValR {
    float r;
    float r_max;
  };  


  //zernike::cornea cornea_global;
 double create_cornea(double, double, double);

    double create_lenticule(double, double, double);
 /*
 surface_group is_point_on_surface(double , double , double , double );
 vtkSmartPointer<vtkUnstructuredGrid> add_nodesets(vtkSmartPointer<vtkDataSet> );
 RetValR compute_cylinder(double , double , double , double );
 double cylinder(double , double , double );
 RetValR compute_cone(double , double , double , double , double);
 double cone(double , double , double );
 */
  /*!
   * @class CorneaMeshTool cornea_meshtool.h 
   * @brief Main class that encapsulates all functions 
   */
  class CorneaMeshTool 
  {  
  
  public:
  
    //! @brief Constructor
    CorneaMeshTool();
  
    //! @brief Destructor
    ~CorneaMeshTool();
    
    void updateConfigSettings(std::map< std::string, std::string > config_update);
    void setPathsToXSDs(std::string library, std::string path_to_xsd);
    void meshCornea();
    void initConfigFile(std::string path_to_config_file);
  
  protected:
  // configuration settings
    std::map<std::string, std::string> configuration_params;
    zernike::cornea cornea;
    mesh_criteria_global_t mesh_criteria_global;
    
  };
    
 
}      // namespace cmt


#endif // CORNEAMESHTOOL_H_
