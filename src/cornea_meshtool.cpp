/**
 * @file   cornea_meshtool.cpp
 * @author Daniel Abler <daniel.abler@istb.unibe.ch>
 * @date   2017
 * @version 1.0
 * @brief  
 */

#include "cornea_meshtool.h"

// LOGGING
#include "easylogging++.h"




using namespace CGAL::parameters;
// Domain
typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef FT_to_point_function_wrapper<K::FT, K::Point_3> Function;
typedef CGAL::Implicit_multi_domain_to_labeling_function_wrapper<Function>
Function_wrapper;
typedef Function_wrapper::Function_vector Function_vector;
//typedef CGAL::Mesh_domain_with_polyline_features_3<CGAL::Labeled_mesh_domain_3<Function_wrapper, K> Mesh_domain;
//typedef CGAL::Mesh_domain_with_polyline_features_3<CGAL::Implicit_mesh_domain_3<Function,K> >  Mesh_domain;
typedef CGAL::Labeled_mesh_domain_3<Function_wrapper, K> Mesh_domain;
// Triangulation
typedef CGAL::Mesh_triangulation_3<Mesh_domain>::type Tr;
typedef CGAL::Mesh_complex_3_in_triangulation_3<Tr> C3t3;
// Mesh Criteria
typedef CGAL::Mesh_criteria_3<Tr> Mesh_criteria;
typedef Mesh_criteria::Facet_criteria    Facet_criteria;
typedef Mesh_criteria::Cell_criteria     Cell_criteria;
// Polyline
typedef K::Point_3 Point;
typedef std::vector<Point>        Polyline_3;
typedef std::list<Polyline_3>       Polylines;


// PUBLIC MEMBER FUNCTIONS 


namespace cmt {

  zernike::cornea cornea_global;
  cornea_generation_criteria_global_t cornea_generation_criteria_global;
    lenticule_generation_criteria_global_t lenticule_generation_criteria_global;
  cornea_boundary_t cornea_boundary_global;
  boundary_cone_t boundary_cone_global;

  bool create_lenticule_global = false;

  RetValR compute_cylinder(double x, double y, double z, double radius){
    RetValR result;
    result.r 	= sqrt(x*x + y*y);
    result.r_max 	= radius; 
    return result;
  }


  double create_cylinder(double x, double y, double z, double radius){
    RetValR result = compute_cylinder(x, y, z, radius); 
    if (result.r <= result.r_max && abs(z) <= 2.0){
      //std::cout << "-----------> inside \n";
      return -1;
    }
    else {
      //std::cout << "----------> outside \n";
      return 1;
    }
  }

  
  RetValR compute_cone(double x, double y, double z, double alpha_rad, double z_apex=0.0){
  RetValR result;
  result.r 	= sqrt(x*x + y*y);
  result.r_max 	= std::abs(z_apex - z) * std::tan(alpha_rad); 
  return result;
}
  
  double create_cone(double x, double y, double z, double alpha_rad, double z_apex){
  RetValR result = compute_cone(x, y, z, alpha_rad, z_apex); 
  LOG(DEBUG) << " alpha = " << alpha_rad  << " pupil_radius = " << cornea_global.pupil_radius
            << " z_apex = " << z_apex << " r_max = " << result.r_max << " r = " << result.r << " z = " << z;
  if (result.r <= result.r_max && z <= z_apex && z >= cornea_generation_criteria_global.offset_along_z - 0.2){ // small offset at lower boundary to avoid mesh distortions
    //std::cout << "-----------> inside \n";
    return -1;
  }
  else {
    //std::cout << "----------> outside \n";
    return 1;
  }
}

  

 double create_cornea(double x, double y, double z){
    // Anterior Surface
   zernike::zernike_surface_return anterior_surface  = zernike::create_zernike_surface(x, y, z, cornea_global.surface_anterior,
                                                                                       cornea_generation_criteria_global.max_number_zernike_coeffs,
                                                                                       cornea_generation_criteria_global.offset_along_z,
                                                                                       0.0,
                                                                                       cornea_generation_criteria_global.surface_thickness);
   //LOG(INFO) << "Cornea anterior surface z at 0,0 : " << anterior_surface.z_at_0_0;
    // Posterior Surface
   double shift_posterior = cornea_global.surface_distance + anterior_surface.z_at_0_0;
   //LOG(INFO) << "Cornea shift : " << shift_posterior;

     zernike::zernike_surface_return posterior_surface = zernike::create_zernike_surface(x, y, z, cornea_global.surface_posterior,
                                                                                       cornea_generation_criteria_global.max_number_zernike_coeffs,
                                                                                       cornea_generation_criteria_global.offset_along_z,
                                                                                       shift_posterior,
                                                                                       cornea_generation_criteria_global.surface_thickness);
   //LOG(INFO) << "Cornea posterior surface z at 0,0 : " << posterior_surface.z_at_0_0;

     // Additional boundary criterusm
   double boundary_criterium = -1;
   // CONE
   if (cornea_boundary_global.type == cone_boundary) {
     double alpha_rad = CGAL_PI / 180.0 * cornea_boundary_global.angle;
     double z_apex = cornea_global.pupil_radius /  std::tan(alpha_rad) +  anterior_surface.z_max;
     boundary_criterium = create_cone(x, y, z, alpha_rad, z_apex);
     LOG(DEBUG) << "--- Detected Cone boundary with angle = " << cornea_boundary_global.angle << " -> " << boundary_criterium;
     //LOG(INFO) << "--- Cone vertex Z = " << z_apex;
   }
   else if (cornea_boundary_global.type == cylinder_boundary) {
     boundary_criterium = create_cylinder( x, y, z, cornea_boundary_global.radius);
     LOG(DEBUG) << "--- Detected Cylinder boundary with radius = " << cornea_boundary_global.radius << " -> " << boundary_criterium;
   }

   if (anterior_surface.inout <=0 && posterior_surface.inout >0  && boundary_criterium <= 0){
   //  if (anterior_surface.inout <=0 && boundary_criterium <= 0){
     //LOG(INFO) << "inside"; 
      return -1;
    }
    else{
      //LOG(INFO) << "outside"; 
      return 1;
    }
  }



    double create_lenticule(double x, double y, double z){
        // Anterior Surface
        zernike::zernike_surface_return anterior_surface  = zernike::create_zernike_surface(x, y, z, cornea_global.surface_anterior,
                                                                                            cornea_generation_criteria_global.max_number_zernike_coeffs,
                                                                                            cornea_generation_criteria_global.offset_along_z,
                                                                                            lenticule_generation_criteria_global.cap_thickness,
                                                                                            lenticule_generation_criteria_global.lenticule_thickness);
        //LOG(INFO) << "Cornea anterior surface z at 0,0 : " << anterior_surface.z_at_0_0;
        // Posterior Surface
        double shift_posterior = lenticule_generation_criteria_global.lenticule_surface_distance + anterior_surface.z_at_0_0 + lenticule_generation_criteria_global.cap_thickness;
        //LOG(INFO) << "Cornea shift : " << shift_posterior;

        zernike::zernike_surface_return posterior_surface = zernike::create_zernike_surface(x, y, z, cornea_global.surface_posterior,
                                                                                            cornea_generation_criteria_global.max_number_zernike_coeffs,
                                                                                            cornea_generation_criteria_global.offset_along_z,
                                                                                            shift_posterior,
                                                                                            lenticule_generation_criteria_global.lenticule_thickness);
        //LOG(INFO) << "Cornea posterior surface z at 0,0 : " << posterior_surface.z_at_0_0;

        // Additional boundary criterusm
        double boundary_criterium = -1;
        // CONE
        if (cornea_boundary_global.type == cone_boundary) {
            double alpha_rad = CGAL_PI / 180.0 * cornea_boundary_global.angle;
            double z_apex = lenticule_generation_criteria_global.lenticule_radius /  std::tan(alpha_rad) +  anterior_surface.z_max;
            boundary_criterium = create_cone(x, y, z, alpha_rad, z_apex);
            LOG(DEBUG) << "--- Detected Cone boundary with angle = " << cornea_boundary_global.angle << " -> " << boundary_criterium;
            //LOG(INFO) << "--- Cone vertex Z = " << z_apex;
        }
        else if (cornea_boundary_global.type == cylinder_boundary) {
            boundary_criterium = create_cylinder( x, y, z, cornea_boundary_global.radius);
            LOG(DEBUG) << "--- Detected Cylinder boundary with radius = " << cornea_boundary_global.radius << " -> " << boundary_criterium;
        }

        if (anterior_surface.inout <=0 && posterior_surface.inout >0  && boundary_criterium <= 0){
        //if (anterior_surface.inout <=0 && boundary_criterium <= 0){
            //LOG(INFO) << "inside";
            return -1;
        }
        else{
            //LOG(INFO) << "outside";
            return 1;
        }
    }



//double create_lenticule(double x, double y, double z){
//
//  // Anterior Surface == anterior cornea surface, but shifted by cap thickness
//  zernike::zernike_surface_return anterior_surface  = zernike::create_zernike_surface(x, y, z, cornea_global.surface_anterior,
//                                                                                      lenticule_generation_criteria_global.max_number_zernike_coeffs,
//                                                                                      cornea_generation_criteria_global.offset_along_z,
//                                                                                      0,
//                                                                                      lenticule_generation_criteria_global.lenticule_thickness);
//  //LOG(INFO) << "lenticule anterior surface z at 0,0 : " << anterior_surface.z_at_0_0;
//  // Posterior Surface
//  //double shift_posterior = lenticule_generation_criteria_global.cap_thickness + lenticule_generation_criteria_global.lenticule_thickness;
//    double shift_posterior = cornea_global.surface_distance + anterior_surface.z_at_0_0; // true thickness
//  //LOG(INFO) << "lenticule surface shift : " << shift_posterior;
//  zernike::zernike_surface_return posterior_surface = zernike::create_zernike_surface(x, y, z, cornea_global.lenticule_surface_posterior,
//                                                                                      lenticule_generation_criteria_global.max_number_zernike_coeffs,
//                                                                                      cornea_generation_criteria_global.offset_along_z,
//                                                                                      shift_posterior,
//                                                                                      lenticule_generation_criteria_global.lenticule_thickness);
//  //LOG(INFO) << "lenticule posterior surface z at 0,0 : " << posterior_surface.z_at_0_0;
//
//
//  // Additional boundary criterusm
//  double boundary_criterium = -1;
//  // CONE
//  if (cornea_boundary_global.type == cone_boundary) {
//    double alpha_rad = CGAL_PI / 180.0 * cornea_boundary_global.angle;
//    double z_apex = lenticule_generation_criteria_global.lenticule_radius /  std::tan(alpha_rad) +  anterior_surface.z_max;
//    boundary_criterium = create_cone(x, y, z, alpha_rad, z_apex);
//    LOG(DEBUG) << "--- Detected Cone boundary with angle = " << cornea_boundary_global.angle << " -> " << boundary_criterium;
//    //LOG(INFO) << "--- Cone vertex Z = " << z_apex;
//  }
//  else if (cornea_boundary_global.type == cylinder_boundary) {
//    boundary_criterium = create_cylinder( x, y, z, cornea_boundary_global.radius);
//    LOG(DEBUG) << "--- Detected Cylinder boundary with radius = " << cornea_boundary_global.radius << " -> " << boundary_criterium;
//  }
//
//  if (anterior_surface.inout <=0 && posterior_surface.inout >0  && boundary_criterium <= 0){
// // if (anterior_surface.inout <=0  && boundary_criterium <= 0){
//    //LOG(INFO) << "inside";
//    return -1;
//  }
//  else{
//    //LOG(INFO) << "outside";
//    return 1;
//  }
//}
//


    surface_group is_point_on_surface(double x, double y, double z, double epsilon)
{
  
  // ANTERIOR SURFACE 
  zernike::RetValZ surface_ant = zernike::compute_surface_shifted(x, y, z, cornea_global.surface_anterior.zernike_coefficients,
                                                                  cornea_generation_criteria_global.max_number_zernike_coeffs,
                                                                  cornea_global.pupil_radius,
                                                                  cornea_generation_criteria_global.offset_along_z,
                                                                  0.0);
  // POSTERIOR SURFACE
  double shift_posterior = cornea_global.surface_distance + surface_ant.z_at_0_0;
  zernike::RetValZ surface_post =  zernike::compute_surface_shifted(x, y, z, cornea_global.surface_posterior.zernike_coefficients,
                                                                    cornea_generation_criteria_global.max_number_zernike_coeffs,
                                                                    cornea_global.pupil_radius,
                                                                    cornea_generation_criteria_global.offset_along_z,
                                                                    shift_posterior);  
  double radius = std::sqrt(x*x+y*y);
  RetValR boundary;
  if (cornea_boundary_global.type == cone_boundary) {
    boundary = compute_cone(x, y, z, boundary_cone_global.alpha_rad, boundary_cone_global.z_apex);
  }
  else if (cornea_boundary_global.type == cylinder_boundary) {
    boundary = compute_cylinder(x, y, z, cornea_boundary_global.radius);
      }


  /// --- 1st iteration: Edge points belong to surface
  int surface_first = 4;
  //if (std::abs(cone.r - cone.r_max) < epsilon){
  //    surface_first = 1;
  //}
  if (std::abs(surface_ant.z_zernike - z) < epsilon) {         // anterior surface
      surface_first = 1;
  }
  else if (std::abs( surface_post.z_zernike  - z) < epsilon) { // posterior surface
      surface_first = 2;
  }
  //else if (std::abs(cone.r - cone.r_max) < epsilon){
  //    surface_first = 1;
  //}
  else if (std::abs(boundary.r - boundary.r_max) < epsilon){    // rim
      surface_first = 3;
  }
  else{
    surface_first = 4; // unclassified
  }
  
  /// --- 2nd iteration: Edge points belong to rim
  int rim_first = 4;
  if (std::abs(boundary.r - boundary.r_max) < epsilon){
      rim_first = 3;
  }
  else if (std::abs(surface_ant.z_zernike - z) < epsilon) {
      rim_first = 1;
  }
  else if (std::abs( surface_post.z_zernike - z) < epsilon) {
      rim_first = 2;
  }
  else{
      rim_first = 4;
  }
  surface_group surface;
  surface.rim_first = rim_first;
  surface.surface_first = surface_first;
  LOG(DEBUG) << "surface first " << surface_first;
  LOG(DEBUG) << "rim first " << rim_first;
  return surface;
}

  
vtkSmartPointer<vtkUnstructuredGrid> add_nodesets(vtkSmartPointer<vtkDataSet> vtkdataset)
{
  LOG(INFO) <<  "======== Adding nodesets to mesh...";
  // label mesh points so that they can be reidentified later on
  vtkSmartPointer<vtkIdFilter> point_id_filter = vtkSmartPointer<vtkIdFilter>::New();
  point_id_filter->SetInputData(vtkdataset);
  point_id_filter->SetPointIds(true);
  point_id_filter->SetCellIds(false);
  point_id_filter->SetIdsArrayName("origPointIDs");
  point_id_filter->Update();
  vtkSmartPointer<vtkUnstructuredGrid> tmp_mesh = point_id_filter->GetUnstructuredGridOutput();

  // Extract Surface
  vtkSmartPointer<vtkDataSetSurfaceFilter> surfaceFilter =  vtkSmartPointer<vtkDataSetSurfaceFilter>::New();
  surfaceFilter->SetInputData(tmp_mesh);
  surfaceFilter->Update();
  vtkSmartPointer<vtkPolyData> vtkpolydata = surfaceFilter->GetOutput();
  LOG(INFO) << "--- Surface mesh has " << vtkpolydata->GetNumberOfPoints() << " points, "
            << vtkpolydata->GetNumberOfVerts() << " vertices, "
            << vtkpolydata->GetNumberOfPolys() << " polygons, "
            << vtkpolydata->GetNumberOfCells() << " cells." ;  
  // create array for holding group association on polydata
  vtkSmartPointer<vtkIntArray> surface_group_info_vtp = vtkSmartPointer<vtkIntArray>::New();
  surface_group_info_vtp->SetName("surface_group");
  surface_group_info_vtp->SetNumberOfValues(vtkpolydata->GetNumberOfPoints());
  vtkSmartPointer<vtkIntArray> rim_group_info_vtp = vtkSmartPointer<vtkIntArray>::New();
  rim_group_info_vtp->SetName("rim_group");
  rim_group_info_vtp->SetNumberOfValues(vtkpolydata->GetNumberOfPoints());
  // create array for holding group association on unstructured grid
  vtkSmartPointer<vtkIntArray> surface_group_info_vtu = vtkSmartPointer<vtkIntArray>::New();
  surface_group_info_vtu->SetName("surface_group");
  surface_group_info_vtu->SetNumberOfValues(tmp_mesh->GetNumberOfPoints());
  vtkSmartPointer<vtkIntArray> rim_group_info_vtu = vtkSmartPointer<vtkIntArray>::New();
  rim_group_info_vtu->SetName("rim_group");
  rim_group_info_vtu->SetNumberOfValues(tmp_mesh->GetNumberOfPoints());
  // initialise arrays for vtu
  for(vtkIdType i = 0; i < vtkdataset->GetNumberOfPoints(); i++)
    { 
      surface_group_info_vtu->SetValue(i, 4);
      rim_group_info_vtu->SetValue(i, 4);
    }
  // get orig point id array from polydata
  //vtkSmartPointer<vtkIntArray> orig_point_id_array = vtkpolydata->GetPointData()->GetArray("origPointIDs");
  // Iterate through all points in polydata surface
  LOG(INFO) << "--- Classifying surface points ...";
  surface_group surface;
  int orig_point_id;
  for(vtkIdType i = 0; i < vtkpolydata->GetNumberOfPoints(); i++)
    {
    double p[3];
    vtkpolydata->GetPoint(i,p);
    LOG(DEBUG) << "Point " << i << " : (" << p[0] << " " << p[1] << " " << p[2] << ")";
    // Test whether point on one of the surfaces
    surface = is_point_on_surface(p[0], p[1], p[2], 0.05);
    surface_group_info_vtp->SetValue(i, surface.surface_first);
    rim_group_info_vtp->SetValue(i, surface.rim_first);	
    // Get original point id for vtkUnstructuredGrid
    LOG(DEBUG)<< "vtkpolydata has 'origPointIDs' array:  " << vtkpolydata->GetPointData()->HasArray("origPointIDs");
    orig_point_id = vtkpolydata->GetPointData()->GetArray("origPointIDs")->GetTuple1(i);
    LOG(DEBUG)<< "Point id : " << i << " ---- orig point id : " << orig_point_id;
    // Assign to vtu array
    surface_group_info_vtu->SetValue(orig_point_id, surface.surface_first);
    rim_group_info_vtu->SetValue(orig_point_id, surface.rim_first);
    }
  LOG(INFO) << "--- ... finished";
  vtkpolydata->GetPointData()->AddArray(surface_group_info_vtp);
  vtkpolydata->GetPointData()->AddArray(rim_group_info_vtp);
  hfvtk::writeVtkData(vtkpolydata, "surface.vtp");
  // Map back to unstructuredGrid
  tmp_mesh->GetPointData()->AddArray(surface_group_info_vtu);
  tmp_mesh->GetPointData()->AddArray(rim_group_info_vtu);
  return tmp_mesh;
}
  






CorneaMeshTool::CorneaMeshTool(){
  // default settings
  this->configuration_params["path_to_config_xsd"] = "../src/xml-io/cornea_meshing_schema.xsd";
  this->configuration_params["intermediate_files_directory"] = ".";
  //this->anterior_surface
}

CorneaMeshTool::~CorneaMeshTool(){}


void CorneaMeshTool::updateConfigSettings(std::map< std::string, std::string > config_update)
{
  for ( std::map< std::string, std::string >::iterator iter = config_update.begin() ; iter != config_update.end(); iter++ )
  {
    std::string key   = iter->first;
    std::string value = iter->second;
    
    LOG(DEBUG) << "Updating setting '" << key << "' from value '" << this->configuration_params[key] << "' to '" << value << "'.";
    this->configuration_params[key] = value;
  }
}


void CorneaMeshTool::setPathsToXSDs(std::string library, std::string path_to_xsd)
{
  if (library == "CorneaMeshTool")
  {
    LOG(DEBUG) << "Path to xsd: " << hf::removeWhiteSpace(path_to_xsd);;
    this->configuration_params["path_to_config_xsd"] = hf::removeWhiteSpace(path_to_xsd);
  }
  else
  {
    LOG(DEBUG) << "There is no library binding for '" << library << "'."; 
  }
}


void CorneaMeshTool::initConfigFile(std::string path_to_config_file)
{
  LOG(INFO) << "======== Initialising MesTool config settings ... ";
  //! 1) Storing configurration file path to internal variable
  this->configuration_params["path_to_config_file"] = hf::removeWhiteSpace(path_to_config_file);
  LOG(INFO) << "--  Using config file: " << path_to_config_file;
  //! 2) Defining xml schema file (overriding schema specification in xml file)
  LOG(INFO) << "--  Using schema file: " << this->configuration_params["path_to_config_xsd"];
  xml_schema::properties props;
  props.no_namespace_schema_location ( this->configuration_params["path_to_config_xsd"] );
  //! 3) Attempt to read config file
  try
  {
    std::auto_ptr<CorneaMeshingParametersType> meshing_params (CorneaMeshingParameters(path_to_config_file, 0, props));
    LOG(INFO) << "--  Reading Config File from " << this->configuration_params["path_to_config_file"];
    // READING SURFACE INFORMATION
    zernike::zernike_coeff zernike_coeff;
    float coeff;
    zernike::j_coeff_map_type j_coeff_map_ant;
    zernike::j_coeff_map_type j_coeff_map_post;
    //---- ANTERIOR SURFACE
    LOG(INFO) << "=== reading anterior surface ...";
    CorneaMeshingParametersType::AnteriorSurface_type anterior_surface_xml = meshing_params->AnteriorSurface();
    // double index notation n, m
    ZernikeSurfaceType::ZernikeCoefficient_sequence& zernike_coefficient_seq_anterior (anterior_surface_xml.ZernikeCoefficient());
    for (ZernikeSurfaceType::ZernikeCoefficient_iterator i (zernike_coefficient_seq_anterior.begin()); i !=zernike_coefficient_seq_anterior.end(); i++)
    {
      ZernikeCoefficientType& zernike_coeff_xml (*i);
      zernike_coeff.n = zernike_coeff_xml.n();
      zernike_coeff.m = zernike_coeff_xml.m();
      coeff = zernike_coeff_xml;
      int zernike_j = zernike::zernike_nm_to_j(zernike_coeff);
      //LOG(INFO) << "m = " << zernike_coeff.m() << ", n = " << zernike_coeff.n() << ", coeff = " <<  zernike_coeff;
      LOG(INFO) << "--  n = " << zernike_coeff.n<< ", m = " << zernike_coeff.m << ", j = " << zernike_j << ", coeff = " <<  coeff;
      j_coeff_map_ant[zernike_j] = coeff;
    }
    // single index notation j
    ZernikeSurfaceType::ZernikeCoefficientSingleIndex_sequence& zernike_coefficient_single_index_seq_anterior (anterior_surface_xml.ZernikeCoefficientSingleIndex());
    for (ZernikeSurfaceType::ZernikeCoefficientSingleIndex_iterator i (zernike_coefficient_single_index_seq_anterior.begin()); i !=zernike_coefficient_single_index_seq_anterior.end(); i++)
      {
        ZernikeCoefficientType_single_index& zernike_coeff_xml (*i);
        int zernike_j = zernike_coeff_xml.j();
        j_coeff_map_ant[zernike_j] =  zernike_coeff_xml;
        LOG(INFO) << "--  j = " << zernike_j << ", coeff = " <<  j_coeff_map_ant[zernike_j];
      }
    //---- POSTERIOR SURFACE
    LOG(INFO) << "=== reading anterior surface ...";
    CorneaMeshingParametersType::PosteriorSurface_type posterior_surface_xml = meshing_params->PosteriorSurface();
    ZernikeSurfaceType::ZernikeCoefficient_sequence& zernike_coefficient_seq_posterior (posterior_surface_xml.ZernikeCoefficient());
    for (ZernikeSurfaceType::ZernikeCoefficient_iterator i (zernike_coefficient_seq_posterior.begin()); i !=zernike_coefficient_seq_posterior.end(); i++)
      {
        ZernikeCoefficientType& zernike_coeff_xml (*i);
        zernike_coeff.n = zernike_coeff_xml.n();
        zernike_coeff.m = zernike_coeff_xml.m();
        coeff = zernike_coeff_xml;
        int zernike_j = zernike::zernike_nm_to_j(zernike_coeff);
        //LOG(INFO) << "m = " << zernike_coeff.m() << ", n = " << zernike_coeff.n() << ", coeff = " <<  zernike_coeff;
        LOG(INFO) << "--  n = " << zernike_coeff.n<< ", m = " << zernike_coeff.m << ", j = " << zernike_j << ", coeff = " <<  coeff;
        j_coeff_map_post[zernike_j] = coeff;
      }
    // single index notation j
    ZernikeSurfaceType::ZernikeCoefficientSingleIndex_sequence& zernike_coefficient_single_index_seq_posterior (posterior_surface_xml.ZernikeCoefficientSingleIndex());
    for (ZernikeSurfaceType::ZernikeCoefficientSingleIndex_iterator i (zernike_coefficient_single_index_seq_posterior.begin()); i !=zernike_coefficient_single_index_seq_posterior.end(); i++)
      {
        ZernikeCoefficientType_single_index& zernike_coeff_xml (*i);
        int zernike_j = zernike_coeff_xml.j();
        j_coeff_map_post[zernike_j] =  zernike_coeff_xml;
        LOG(INFO) << "--  j = " << zernike_j << ", coeff = " <<  j_coeff_map_post[zernike_j];
      }
    //---- PUPIL RADIUS
    LOG(INFO) << "=== reading pupil radius ...";
    CorneaMeshingParametersType::PupilRadius_type pupil_radius = meshing_params->PupilRadius();
    LOG(INFO) << "--  pupil radius = " << pupil_radius;
    //---- SURFACE DISTANCE
    LOG(INFO) << "=== reading surface distance ...";
    CorneaMeshingParametersType::SurfaceDistance_type surface_distance = meshing_params->SurfaceDistance();
    LOG(INFO) << "--  surface distance = " << surface_distance;
    //---- DEFINE ZERNIKE SURFACE
    zernike::zernike_surface anterior_surface;
    anterior_surface.zernike_coefficients = j_coeff_map_ant;
    anterior_surface.pupil_radius = pupil_radius;
    zernike::zernike_surface posterior_surface;
    posterior_surface.zernike_coefficients = j_coeff_map_post;
    posterior_surface.pupil_radius = pupil_radius;
    //---- DEFINE CORNEA 
    this->cornea.surface_anterior = anterior_surface;
    this->cornea.surface_posterior= posterior_surface;
    this->cornea.pupil_radius     = pupil_radius;
    this->cornea.surface_distance = surface_distance;
    //---- DEFINE CORNEA GLOBAL
    cornea_global.surface_anterior = anterior_surface;
    cornea_global.surface_posterior= posterior_surface;
    cornea_global.pupil_radius     = pupil_radius;
    cornea_global.surface_distance = surface_distance;



      //---- DEFINE ZERNIKE SURFACE
      zernike::zernike_surface posterior_surface_lenticule;

    CorneaMeshingParametersType::PosteriorSurfaceLenticule_optional posterior_surface_lenticule_xml = meshing_params->PosteriorSurfaceLenticule();
    if (posterior_surface_lenticule_xml.present()) {
        //---- LENTICULE
        LOG(INFO) << "=== reading lenticule posterior surface ... ";
        // READING SURFACE INFORMATION
        zernike::zernike_coeff zernike_coeff_lenticule;
        float coeff_lenticule;
        zernike::j_coeff_map_type j_coeff_map_post_lenticule;
        //CorneaMeshingParametersType::PosteriorSurfaceLenticule_type posterior_surface_lenticule_xml = meshing_params->PosteriorSurfaceLenticule();
        // double index notation n, m
        ZernikeSurfaceType::ZernikeCoefficient_sequence &zernike_coefficient_seq_posterior_lenticule(
                posterior_surface_lenticule_xml->ZernikeCoefficient());
        for (ZernikeSurfaceType::ZernikeCoefficient_iterator i(zernike_coefficient_seq_posterior_lenticule.begin());
             i != zernike_coefficient_seq_posterior_lenticule.end(); i++) {
            ZernikeCoefficientType &zernike_coeff_xml(*i);
            zernike_coeff_lenticule.n = zernike_coeff_xml.n();
            zernike_coeff_lenticule.m = zernike_coeff_xml.m();
            coeff_lenticule = zernike_coeff_xml;
            int zernike_j_lenticule = zernike::zernike_nm_to_j(zernike_coeff);
            //LOG(INFO) << "m = " << zernike_coeff.m() << ", n = " << zernike_coeff.n() << ", coeff = " <<  zernike_coeff;
            LOG(INFO) << "--  n = " << zernike_coeff_lenticule.n << ", m = " << zernike_coeff_lenticule.m << ", j = "
                      << zernike_j_lenticule << ", coeff = " << coeff_lenticule;
            j_coeff_map_post_lenticule[zernike_j_lenticule] = coeff_lenticule;
        }
        // single index notation j
        ZernikeSurfaceType::ZernikeCoefficientSingleIndex_sequence &zernike_coefficient_single_index_seq_posterior_lenticule(
                posterior_surface_lenticule_xml->ZernikeCoefficientSingleIndex());
        for (ZernikeSurfaceType::ZernikeCoefficientSingleIndex_iterator i(
                zernike_coefficient_single_index_seq_posterior_lenticule.begin());
             i != zernike_coefficient_single_index_seq_posterior_lenticule.end(); i++) {
            ZernikeCoefficientType_single_index &zernike_coeff_xml(*i);
            int zernike_j_lenticule = zernike_coeff_xml.j();
            j_coeff_map_post_lenticule[zernike_j_lenticule] = zernike_coeff_xml;
            LOG(INFO) << "--  j = " << zernike_j_lenticule << ", coeff = " << j_coeff_map_ant[zernike_j_lenticule];
        }

        posterior_surface_lenticule.zernike_coefficients = j_coeff_map_post_lenticule;


        CorneaMeshingParametersType::LenticuleGenerationCriteria_optional lenticulecriteria_xml = meshing_params->LenticuleGenerationCriteria();
        if (lenticulecriteria_xml.present()) {
            //---- LENTICULE GENERATION  CRITERIA GLOBAL
            LOG(INFO) << "=== reading lenticule generation parameters ... ";

            lenticule_generation_criteria_global.max_number_zernike_coeffs = lenticulecriteria_xml->max_number_zernike_coeffs();
            LOG(INFO) << "   - max number of zernike coeffs to be considered : "
                      << lenticule_generation_criteria_global.max_number_zernike_coeffs;
            lenticule_generation_criteria_global.lenticule_thickness = lenticulecriteria_xml->lenticule_thickness();
            LOG(INFO) << "   - lenticule thickness :                             "
                      << lenticule_generation_criteria_global.lenticule_thickness;
            lenticule_generation_criteria_global.lenticule_surface_distance = lenticulecriteria_xml->lenticule_surface_distance();
            LOG(INFO) << "   - lenticule surface distance :                             "
                      << lenticule_generation_criteria_global.lenticule_thickness;
            lenticule_generation_criteria_global.lenticule_radius = lenticulecriteria_xml->lenticule_radius();
            LOG(INFO) << "   - lenticule radius:                     "
                      << lenticule_generation_criteria_global.lenticule_radius;
            lenticule_generation_criteria_global.cap_thickness = lenticulecriteria_xml->cap_thickness();
            LOG(INFO) << "   - cap thickness:                                 "
                      << lenticule_generation_criteria_global.cap_thickness;

            posterior_surface_lenticule.pupil_radius = lenticule_generation_criteria_global.lenticule_radius;

        }

        this->cornea.lenticule_surface_posterior = posterior_surface_lenticule;

        //LOG(INFO) << "=== Checking lenticule parameter consistency:";
        if (lenticule_generation_criteria_global.lenticule_surface_distance < lenticule_generation_criteria_global.lenticule_thickness)
        {
            //LOG(INFO) << "    - lenticule surface distance < lenticule thickness ";
        }
        else{
            LOG(INFO) << "    - WARNING: lenticule surface distance > lenticule thickness ";
        }




        create_lenticule_global = true;
    }







    //---- MESH CRITERIA GLOBAL
	  LOG(INFO) << "=== reading mesh parameters ... ";
    CorneaMeshingParametersType::MeshCriteria_type meshcriteria_xml = meshing_params->MeshCriteria();
	 	// Cell Size
	 	// Cell Size
	  this->mesh_criteria_global.cell_size = meshcriteria_xml.cell_size();
	  LOG(INFO) << "   - cell_size:                      " <<  this->mesh_criteria_global.cell_size;
	  // cell_radius_edge_ratio
	  this->mesh_criteria_global.cell_radius_edge_ratio = meshcriteria_xml.cell_radius_edge_ratio();
	  LOG(INFO) << "   - cell_radius_edge_ratio:         " <<  this->mesh_criteria_global.cell_radius_edge_ratio;
	  // facet_angle
	  this->mesh_criteria_global.facet_angle = meshcriteria_xml.facet_angle();
	  LOG(INFO) << "   - facet_angle:                    " <<  this->mesh_criteria_global.facet_angle;
	  // facet_distance
	  this->mesh_criteria_global.facet_distance = meshcriteria_xml.facet_distance();
	  LOG(INFO) << "   - facet_distance:                " <<  this->mesh_criteria_global.facet_distance;
	  // facet_size
	  this->mesh_criteria_global.facet_size = meshcriteria_xml.facet_size();
	  LOG(INFO) << "   - facet_size:                     " <<  this->mesh_criteria_global.facet_size;

     //---- CORNEA GENERATION  CRITERIA GLOBAL
	  LOG(INFO) << "=== reading cornea generation parameters ... ";
    CorneaMeshingParametersType::CorneaGenerationCriteria_type corneacriteria_xml = meshing_params->CorneaGenerationCriteria();
	 	cornea_generation_criteria_global.max_number_zernike_coeffs = corneacriteria_xml.max_number_zernike_coeffs();
	  LOG(INFO) << "   - max number of zernike coeffs to be considered : " <<  cornea_generation_criteria_global.max_number_zernike_coeffs;
    cornea_generation_criteria_global.surface_thickness         = corneacriteria_xml.surface_thickness();
    LOG(INFO) << "   - surface thickness :                             " <<  cornea_generation_criteria_global.surface_thickness;
    cornea_generation_criteria_global.offset_along_z            = corneacriteria_xml.offset_along_z();
    LOG(INFO) << "   - cornea offset along z axis:                     " <<  cornea_generation_criteria_global.offset_along_z;
    cornea_generation_criteria_global.path_to_output            = corneacriteria_xml.path_to_output();
    LOG(INFO) << "   - path to output:                                 " <<  cornea_generation_criteria_global.path_to_output;
    /*
    // TODO not working: 'cannot convert str to const string*'
    if ( not hf::hasSuffix(cornea_generation_criteria_global.path_to_output, "vtu" )){
      LOG(FATAL) << "!!! output path must point to '.vtu' file !!! ";
    }
    */
        
    
    //----- BOUNDARY
    LOG(INFO) << "=== reading cornea boundary parameters ...";
    CorneaMeshingParametersType::CorneaBoundary_type cornea_boundary = meshing_params->CorneaBoundary();
    CorneaBoundaryType::Cone_optional cone = cornea_boundary.Cone();
    if (cone.present()){
      ConeBoundaryType::angle_type angle = cone->angle();
      cornea_boundary_global.type = cone_boundary;
      cornea_boundary_global.angle = angle;
      LOG(INFO) << "--- Cone angle : " << cornea_boundary_global.angle;
      // Compute z_apex here, stays the same throughout
      // dummy surface probe
      zernike::RetValZ surface_ant = zernike::compute_surface_shifted(0, 0, 0, cornea_global.surface_anterior.zernike_coefficients,
                                                                      cornea_generation_criteria_global.max_number_zernike_coeffs,
                                                                      cornea_global.pupil_radius,
                                                                      cornea_generation_criteria_global.offset_along_z,
                                                                      0.0);
      double alpha_rad = CGAL_PI / 180.0 * cornea_boundary_global.angle;
      double z_apex = cornea_global.pupil_radius /  std::tan(alpha_rad) +  surface_ant.z_max;
      boundary_cone_global.z_apex = z_apex;
      boundary_cone_global.alpha_rad = alpha_rad;
      LOG(INFO) << "--- Cone Apex at z = " << boundary_cone_global.z_apex;
    }
    
    CorneaBoundaryType::Cylinder_optional cylinder = cornea_boundary.Cylinder();
    if (cylinder.present()){
      CylinderBoundaryType::radius_type radius = cylinder->radius();
      cornea_boundary_global.type = cylinder_boundary;
      cornea_boundary_global.radius = radius;
      LOG(INFO) << "--- Cylinder radius : " << cornea_boundary_global.radius;
    }

    // WRITING Computed Config Params as python dict
    std::string path = hf::getDirFromPath(cornea_generation_criteria_global.path_to_output);
    hf::ensureDirExists(path);
    std::string path_to_python_dict_file = path + "/params.py";
    ofstream myfile;
    myfile.open (path_to_python_dict_file);
    myfile << "params = {} \n";
    myfile << "params['z_apex'] = " << boundary_cone_global.z_apex << "\n";
    myfile.close();

    // CorneaBoundaryType::ZernikeCoefficient_sequence& zernike_coefficient_seq_posterior (posterior_surface_xml.ZernikeCoefficient());
    // for (ZernikeSurfaceType::ZernikeCoefficient_iterator i (zernike_coefficient_seq_posterior.begin()); i !=zernike_coefficient_seq_posterior.end(); i++)
    //   {
    //     ZernikeCoefficientType& zernike_coeff_xml (*i);
    //     zernike_coeff.n = zernike_coeff_xml.n();
    //     zernike_coeff.m = zernike_coeff_xml.m();
    //     coeff = zernike_coeff_xml;
    //     int zernike_j = zernike::zernike_nm_to_j(zernike_coeff);
    //     //LOG(INFO) << "m = " << zernike_coeff.m() << ", n = " << zernike_coeff.n() << ", coeff = " <<  zernike_coeff;
    //     LOG(INFO) << "--  n = " << zernike_coeff.n<< ", m = " << zernike_coeff.m << ", j = " << zernike_j << ", coeff = " <<  coeff;
    //     j_coeff_map_post[zernike_j] = coeff;
    //   }

  }

  catch (const xml_schema::exception& e)
  {
    LOG(FATAL) << ">>>>> Problem reading config file: " << e;
  }
  LOG(DEBUG) << "Successfully read config settings. ";     
}

void CorneaMeshTool::meshCornea()
{
  LOG(INFO) << "======== Starting to mesh cornea ...";
  // Define functions
  Function f1(&create_cornea);
  Function_vector v;
  v.push_back(f1);

    if (create_lenticule_global)
    {
        Function f2(&create_lenticule);
        v.push_back(f2);
    }
  //std::vector<std::string> vps;
  //vps.push_back("+-");
  // Domain (Warning: Sphere_3 constructor uses square radius !)
  //Mesh_domain domain(Function_wrapper(v, vps), K::Sphere_3(CGAL::ORIGIN, 5.*5.));
  Mesh_domain domain(Function_wrapper(v), K::Sphere_3(CGAL::ORIGIN, 30.*30.));
  //Mesh_domain domain(v, K::Sphere_3(CGAL::ORIGIN, 30.*30.));
  //Mesh_domain domain(f1, K::Sphere_3(CGAL::ORIGIN, 30.*30.));
  // Set mesh criteria
  /*
  - facet_angle. This parameter controls the shape of surface facets. Actually, it is a lower bound for the angle (in degree) of surface facets. When boundary surfaces are smooth, the termination of the meshing process is guaranteed if the angular bound is at most 30 degrees.
  - facet_size. Thi sparameter controls the size of surface facets. Actually, each surface facet has a surface Delaunay ball which is a ball circumscribing the surface facet and centered on the surface patch. The parameter facet_size is either a constant or a spatially variable scalar field, providing an upper bound for the radii of surface Delaunay balls.
  - facet_distance. This parameter controls the approximation error of boundary and subdivision surfaces. Actually, it is either a constant or a spatially variable scalar field. It provides an upper bound for the distance between the circumcenter of a surface facet and the center of a surface Delaunay ball of this facet.
  - facet_topology. This parameters controls the set of topological constraints which have to be verified by each surface facet. By default, each vertex of a surface facet has to be located on a surface patch, on a curve segment, or on a corner. It can also be set to check whether the three vertices of a surface facet belongs to the same surface patch. This has to be done cautiously, as such a criteria needs that each surface patches intersection is an input 1-dimensional feature.
  The criteria for mesh cells are governed by two parameters:
  - cell_radius_edge_ratio. This parameter controls the shape of mesh cells (but can't filter slivers, as we discussed earlier). Actually, it is an upper bound for the ratio between the circumradius of a mesh tetrahedron and its shortest edge. There is a theoretical bound for this parameter: the Delaunay refinement process is guaranteed to terminate for values of cell_radius_edge_ratio bigger than 2.
  - cell_size. This parameter controls the size of mesh tetrahedra. It is either a scalar or a spatially variable scalar field. It provides an upper bound on the circumradii of the mesh tetrahedra.
  */
  Facet_criteria facet_criteria( this->mesh_criteria_global.facet_angle, // angle
                                 this->mesh_criteria_global.facet_size,  // size
                                 this->mesh_criteria_global.facet_distance); // approximation
  Cell_criteria cell_criteria( this->mesh_criteria_global.cell_radius_edge_ratio, // radius-edge ratio,
                               this->mesh_criteria_global.cell_size); // size
  Mesh_criteria criteria(facet_criteria, cell_criteria);

  /*
  // TODO make features work: include boundary mode here (cone, cylinder)!
  LOG(INFO) << "--- Adding edges as additional features to cornea mesh ...";
  // Create edge for surface 1
  zernike::coord_2D_spher coords_sphere;
  coords_sphere.r = cornea_global.pupil_radius;
  Polylines polylines_ant (1);
  Polyline_3& polyline_ant = polylines_ant.front();
  Polylines polylines_post (1);
  Polyline_3& polyline_post = polylines_post.front();
  // sample zernike surface along edges
  for(int phi = 0; phi < 3600; ++phi)
  {
    coords_sphere.phi = phi/100.0*CGAL_PI/180.0;
    zernike::coord_2D_cart coords_cart = zernike::spher_to_cart_2D(coords_sphere);
    // ANTERIOR
    zernike::RetValZ surface_ant = zernike::compute_surface_shifted( coords_cart.x, coords_cart.y, 0,
                                                           cornea_global.surface_anterior.zernike_coefficients,
                                                           cornea_generation_criteria_global.max_number_zernike_coeffs,
                                                           cornea_global.pupil_radius,
                                                           cornea_generation_criteria_global.offset_along_z,
                                                          0.0);
    Point p_ant (coords_cart.x, coords_cart.y, surface_ant.z_zernike);
    polyline_ant.push_back(p_ant);
    LOG(DEBUG) << "-- anterior surface point: x = " << coords_cart.x << " ,y = " << coords_cart.y  << " ,z = " << surface_ant.z_zernike
              << "(r = " << coords_sphere.r << ", phi = " << phi << " )";
    // POSTERIOR
    double shift_posterior = cornea_global.surface_distance + surface_ant.z_at_0_0;
    zernike::RetValZ surface_post = zernike::compute_surface_shifted( coords_cart.x, coords_cart.y, 0,
                                                           cornea_global.surface_posterior.zernike_coefficients,
                                                           cornea_generation_criteria_global.max_number_zernike_coeffs,
                                                           cornea_global.pupil_radius,
                                                           cornea_generation_criteria_global.offset_along_z,
                                                           shift_posterior);
    Point p_post (coords_cart.x, coords_cart.y, surface_post.z_zernike);
    polyline_post.push_back(p_post);
    LOG(DEBUG) << "-- posterio surface point: x = " << coords_cart.x << " ,y = " << coords_cart.y  << " ,z = " << surface_post.z_zernike
              << "(r = " << coords_sphere.r << ", phi = " << phi << " )";
  }
  polyline_ant.push_back(polyline_ant.front()); // close the line
  polyline_post.push_back(polyline_post.front()); // close the line
  // Insert edge in domain
  domain.add_features(polylines_ant.begin(), polylines_ant.end());
  domain.add_features(polylines_post.begin(), polylines_post.end());
 */

  LOG(INFO) << "--- Generating mesh ...";
  // Mesh generation with feature preservation
  C3t3 c3t3 = CGAL::make_mesh_3<C3t3>(domain, criteria);

  // Output
  //std::ofstream medit_file("out.mesh");
  //CGAL::output_to_medit(medit_file, c3t3);

  LOG(INFO) << "--- Converting cgal mesh to vtkUnstructuredGrid...";
  vtkSmartPointer<vtkUnstructuredGrid> unstructuredGrid = c3t3_subdomain_to_vtk_unstructured_grid(c3t3);
  vtkSmartPointer<vtkUnstructuredGrid> unstructuredGrid_nodesets = add_nodesets(unstructuredGrid);
  LOG(INFO) << "--- Writing mesh to " << cornea_generation_criteria_global.path_to_output;
  std::string path = hf::getDirFromPath(cornea_generation_criteria_global.path_to_output);
  hf::ensureDirExists(path);
  hfvtk::writeVtkUnstructuredGrid(unstructuredGrid_nodesets, cornea_generation_criteria_global.path_to_output);

}

} // namespace cmt



