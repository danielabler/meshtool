/**
 * @file   meshtool.cpp
 * @author Daniel Abler <daniel.abler@istb.unibe.ch>
 * @date    2015
 * @version 1.0
 * @brief  main file for meshtool
 */

#include "image_meshtool.h"
#include <CGAL/read_vtk_image_data.h>
#include <CGAL/Image_3_vtk_interface.h>

// LOGGING
#include "easylogging++.h"

// This allows to output intermediate files for debugging
// files are stored in his->configuration_params["intermediate_files_directory"]
#define __WITH_INT_OUTPUT__


namespace imt {
   
// PUBLIC MEMBER FUNCTIONS 

// PUBLIC MEMBER FUNCTIONS 
  
MeshTool::MeshTool(){
  // default settings
  this->configuration_params["path_to_MeshTool_config_xsd"] = "../src/xml-io/imaging_meshing_schema.xsd";
  this->configuration_params["intermediate_files_directory"] = ".";
}

MeshTool::~MeshTool(){}


void MeshTool::updateConfigSettings(std::map< std::string, std::string > config_update)
{
  for ( std::map< std::string, std::string >::iterator iter = config_update.begin() ; iter != config_update.end(); iter++ )
  {
    std::string key   = iter->first;
    std::string value = iter->second;
    
    LOG(DEBUG) << "Updating setting '" << key << "' from value '" << this->configuration_params[key] << "' to '" << value << "'.";
    this->configuration_params[key] = value;
  }
}



void MeshTool::setPathsToXSDs(std::string library, std::string path_to_xsd)
{
  if (library == "ImageMeshTool")
  {
    this->configuration_params["path_to_MeshTool_config_xsd"] = hf::removeWhiteSpace(path_to_xsd);
  }
  else
  {
    LOG(DEBUG) << "There is no library binding for '" << library << "'."; 
  }
}

void MeshTool::initConfigFile(std::string path_to_config_file)
{
  LOG(INFO) << "Initialising MesTool config settings ... ";
    //! 1) Storing configuration file path to internal variable
    this->configuration_params["path_to_config_file"] = hf::removeWhiteSpace(path_to_config_file);
    LOG(INFO) << "Using config file: " << path_to_config_file;
    //! 2) Defining xml schema file (overriding schema specification in xml file)
    xml_schema::properties props;
    props.no_namespace_schema_location ( this->configuration_params["path_to_MeshTool_config_xsd"] );
    LOG(INFO) << "Using schema file: " << this->configuration_params["path_to_MeshTool_config_xsd"];
    //! 3) Attempt to read config file
    try
    {
      std::auto_ptr<MeshTool_Config_t> config (MeshTool_Config(path_to_config_file, 0, props));
      LOG(INFO) << "--------------------------------------------------";
      LOG(INFO) << "Reading Config File from " << this->configuration_params["path_to_config_file"]; 	     
      //! Processed Variables:
      //! - _MAIN CONFIG_   
      
      // IF image to mesh
      MeshTool_Config_t::Operation_Image2Mesh_optional& i2m = config->Operation_Image2Mesh();
      if (i2m.present())
      {
	this->configuration_params["operation_mode"] = "image_to_mesh";
	//! OPERATION IMAGE_TO_MESH
	MeshTool_Config_t::Operation_Image2Mesh_type& image_to_mesh = config->Operation_Image2Mesh().get();
	//! path_to_input_file
	this->configuration_params["path_to_input_file"] = image_to_mesh.path_to_input_file();
	LOG(INFO) << "PARAMETER path_to_input_file: 		" << this->configuration_params["path_to_input_file"];
	//! path_to_output_file
	this->configuration_params["path_to_output_file"] = image_to_mesh.path_to_output_file();
	LOG(INFO) << "PARAMETER path_to_output_file: 		" << this->configuration_params["path_to_output_file"];
	
	//! Mesh Criteria Global
	LOG(INFO) << "PARAMETER Mesh Criteria Global: ";
	Operation_Image2Mesh_t::MeshCriteria_global_type& mesh_criteria_global = image_to_mesh.MeshCriteria_global();
	// Cell Size
	this->mesh_criteria_global.cell_size = mesh_criteria_global.cell_size();
	LOG(INFO) << "   - cell_size:                      " <<  this->mesh_criteria_global.cell_size;
	// cell_radius_edge_ratio
	this->mesh_criteria_global.cell_radius_edge_ratio = mesh_criteria_global.cell_radius_edge_ratio();
	LOG(INFO) << "   - cell_radius_edge_ratio:         " <<  this->mesh_criteria_global.cell_radius_edge_ratio;
	// facet_angle
	this->mesh_criteria_global.facet_angle = mesh_criteria_global.facet_angle();
	LOG(INFO) << "   - facet_angle:                    " <<  this->mesh_criteria_global.facet_angle;
	// facet_distance
	this->mesh_criteria_global.facet_distance = mesh_criteria_global.facet_distance();
	LOG(INFO) << "   - facet_distance:                " <<  this->mesh_criteria_global.facet_distance;
	// facet_size
	this->mesh_criteria_global.facet_size = mesh_criteria_global.facet_size();
	LOG(INFO) << "   - facet_size:                     " <<  this->mesh_criteria_global.facet_size;
	
      /**
       * CGAL Mesh Criteria:
       * [- edge_size: 			a scalar field (resp. a constant) providing a space varying (resp. a uniform) upper bound for the lengths of curve segment edges. This parameter has to be set to a positive value when 1-dimensional features protection is used.]
       * 
       * - facet_angle: 		a lower bound for the angles (in degrees) of the surface mesh facets.
       * - facet_size:  		a scalar field (resp. a constant) describing a space varying (resp. a uniform) upper-bound or for the radii of the surface Delaunay balls.
       * - facet_distance: 		a scalar field (resp. a constant) describing a space varying (resp. a uniform) upper bound for the same distance.
       * [- facet_topology: 		the set of topological constraints which have to be verified by each surface facet. The default value is CGAL::FACET_VERTICES_ON_SURFACE. See Mesh_facet_topology manual page to get all possible values.]
       * 
       * - cell_radius_edge_ratio: 	an upper bound for the radius-edge ratio of the mesh tetrahedra.
       * - cell_size: 			a scalar field (resp. a constant) describing a space varying (resp. a uniform) upper-bound for the circumradii of the mesh tetrahedra.
       * 
       * see http://doc.cgal.org/latest/Mesh_3/classCGAL_1_1Mesh__criteria__3.html
       */
      
	//! Mesh Criteria Subdomain
	Operation_Image2Mesh_t::MeshCriteria_SubDomain_sequence& mesh_criteria_subdom_seq = image_to_mesh.MeshCriteria_SubDomain();
	for (Operation_Image2Mesh_t::MeshCriteria_SubDomain_iterator i (mesh_criteria_subdom_seq.begin()); i !=mesh_criteria_subdom_seq.end(); i++)
	{
	  MeshCriteria_SubDomain_t& mesh_criteria_subdomain (*i);
	  mesh_criteria_subdomain_t mesh_criteria;
	  mesh_criteria.cell_size = mesh_criteria_subdomain.cell_size();
	  mesh_criteria.dimension = mesh_criteria_subdomain.dimension();
	  int subdomain_id	= mesh_criteria_subdomain.domain_id();
	  LOG(INFO) << "PARAMETER Mesh Criteria Subdomain: " << subdomain_id;
	  LOG(INFO) << "   - cell_size:                     " <<  mesh_criteria.cell_size;
	  LOG(INFO) << "   - dimension:                     " <<  mesh_criteria.dimension;
	  this->mesh_criteria_subdomain_map[subdomain_id] = mesh_criteria;
	  
	}
      }
      
      
           // IF image to mesh
      MeshTool_Config_t::Operation_Remesh_optional& m2m = config->Operation_Remesh();
      if (m2m.present())
      {
	this->configuration_params["operation_mode"] = "mesh_to_mesh";
	//! OPERATION MESH_TO_MESH
	MeshTool_Config_t::Operation_Remesh_type& remesh = config->Operation_Remesh().get();
	//! path_to_input_file
	this->configuration_params["path_to_input_file"] = remesh.path_to_input_file();
	LOG(INFO) << "PARAMETER path_to_input_file: 		" << this->configuration_params["path_to_input_file"];
	//! path_to_output_file
	this->configuration_params["path_to_output_file"] = remesh.path_to_output_file();
	LOG(INFO) << "PARAMETER path_to_output_file: 		" << this->configuration_params["path_to_output_file"];
	
	// Spacing
	this->spacing[0] = remesh.spacing_xyz();
	this->spacing[1] = remesh.spacing_xyz();
	this->spacing[2] = remesh.spacing_xyz();
	//! Mesh Criteria Global
	LOG(INFO) << "PARAMETER Mesh Criteria Global: ";
	Operation_Remesh_t::MeshCriteria_global_type& mesh_criteria_global = remesh.MeshCriteria_global();
	// Cell Size
	this->mesh_criteria_global.cell_size = mesh_criteria_global.cell_size();
	LOG(INFO) << "   - cell_size:                      " <<  this->mesh_criteria_global.cell_size;
	// cell_radius_edge_ratio
	this->mesh_criteria_global.cell_radius_edge_ratio = mesh_criteria_global.cell_radius_edge_ratio();
	LOG(INFO) << "   - cell_radius_edge_ratio:         " <<  this->mesh_criteria_global.cell_radius_edge_ratio;
	// facet_angle
	this->mesh_criteria_global.facet_angle = mesh_criteria_global.facet_angle();
	LOG(INFO) << "   - facet_angle:                    " <<  this->mesh_criteria_global.facet_angle;
	// facet_distance
	this->mesh_criteria_global.facet_distance = mesh_criteria_global.facet_distance();
	LOG(INFO) << "   - facet_distance:                " <<  this->mesh_criteria_global.facet_distance;
	// facet_size
	this->mesh_criteria_global.facet_size = mesh_criteria_global.facet_size();
	LOG(INFO) << "   - facet_size:                     " <<  this->mesh_criteria_global.facet_size;
	
	//! Mesh Criteria Subdomain
	Operation_Remesh_t::MeshCriteria_SubDomain_sequence& mesh_criteria_subdom_seq = remesh.MeshCriteria_SubDomain();
	for (Operation_Remesh_t::MeshCriteria_SubDomain_iterator i (mesh_criteria_subdom_seq.begin()); i !=mesh_criteria_subdom_seq.end(); i++)
	{
	  MeshCriteria_SubDomain_t& mesh_criteria_subdomain (*i);
	  mesh_criteria_subdomain_t mesh_criteria;
	  mesh_criteria.cell_size = mesh_criteria_subdomain.cell_size();
	  mesh_criteria.dimension = mesh_criteria_subdomain.dimension();
	  int subdomain_id	= mesh_criteria_subdomain.domain_id();
	  LOG(INFO) << "PARAMTER Mesh Criteria Subdomain: " << subdomain_id;
	  LOG(INFO) << "   - cell_size:                     " <<  mesh_criteria.cell_size;
	  LOG(INFO) << "   - dimension:                     " <<  mesh_criteria.dimension;
	  this->mesh_criteria_subdomain_map[subdomain_id] = mesh_criteria;
	  
	}
      }
      
      
      
      LOG(INFO) << "--------------------------------------------------";
    }
    catch (const xml_schema::exception& e)
    {
      LOG(FATAL) << "Problem reading config file: " << e;
    }
    
    LOG(DEBUG) << "Successfully read config settings. ";     
} 



void MeshTool::loadImage(std::string path_to_image_file)
{
    // initialise vtk, cgal images
    vtkSmartPointer<vtkImageData> vtk_image = vtkSmartPointer<vtkImageData>::New();
    CGAL::Image_3 cgal_image;

    std::string inria_suffix    = "inr";
    std::string inria_suffix_gz = "inr.gz";
    if ( hf::hasSuffix(path_to_image_file , inria_suffix) or hf::hasSuffix(path_to_image_file, inria_suffix_gz))
    {
        LOG(INFO) << "Input is INRIA file -- loading directly with CGAL.";
        const char* char_input_file = path_to_image_file.c_str();
        cgal_image.read(char_input_file);
    }
    else
    {
        LOG(INFO) << "Input is not INRIA file -- loading with VTK.";

        vtk_image = hfvtk::loadImage(path_to_image_file);
	
	// Check if vtk_image has scalar data of type VTK_UNISGNED_CHAR -- otherwise CGAL behaves strangely...
	LOG(DEBUG) << "VTK_IMAGE has Scalar data of type : " << vtk_image->GetScalarTypeAsString();
	if (not (vtk_image->GetScalarType()==3))
	{
	  LOG(INFO) << "Label information must be provided as 'unsigned char' -- trying to convert";

	  vtkSmartPointer<vtkImageCast> cast = vtkSmartPointer<vtkImageCast>::New();
	  cast->SetInputData(vtk_image);
	  cast->SetOutputScalarTypeToUnsignedChar();
	  cast->Update();

	  vtk_image = cast->GetOutput();
	}
	
	LOG(DEBUG) << "Storing original ctk image to internal variable";
	this->vtk_image_original = vtk_image;
    
	
	
        LOG(DEBUG) << "Converting vtk image to cgal image";
        //cgal_image.read_vtk_image_data(vtk_image);
	cgal_image = CGAL::read_vtk_image_data(vtk_image);
    }

    LOG(DEBUG) << "Storing CGAL image to internal variable";
    this->cgal_image_original = cgal_image;
}


void MeshTool::operationCreateMeshFromImage()
{
 //! 1) Load Image as cgal:image_3
      loadImage(this->configuration_params["path_to_input_file"]);  
 
 //! 2) create mesh from labeled cgal=image_3
     createMeshFromLabeledImage(this->cgal_image_original);
 
 //! 3) Align output mesh with original input image
     double origin[3];
     this->vtk_image_original->GetOrigin(origin);
     
     vtkSmartPointer<vtkTransform> transform = vtkSmartPointer<vtkTransform>::New();
     transform->Translate(origin);
     vtkSmartPointer<vtkTransformFilter> transform_filter = vtkSmartPointer<vtkTransformFilter>::New();
     transform_filter->SetInputData(this->unstructured_grid_processed);
     transform_filter->SetTransform(transform);
     transform_filter->Update();
     this->unstructured_grid_processed = vtkUnstructuredGrid::SafeDownCast(transform_filter->GetOutput());
     
 //! 3) Output mesh
     hfvtk::writeVtkUnstructuredGrid(this->unstructured_grid_processed, this->configuration_params["path_to_output_file"]);
}

void MeshTool::createMeshFromLabeledImage(CGAL::Image_3 cgal_image)
{
    //! @todo normalisation of cgal meshing parameters to voxel sizes
  
//     double vox_x = cgal_image.vx();
//     double vox_y = cgal_image.vy();
//     double vox_z = cgal_image.vz();

//     double N;
//     if (Configuration::Instance()->CgalMeshNormalised())
//     {
//         N = ( vox_x + vox_y + vox_z ) / 3.0 ;
//         LOG(DEBUG) << "Using normalised mesh criteria, N = " << N;
//     }
//     else
//     {
//         N = 1.0;
//         LOG(DEBUG) << "Using mesh criteria as given, N = " << N;
//     }
    LOG(INFO) << "Image to be meshed has:";
    LOG(INFO) << "-- Voxel sizes: " << "vx = " << cgal_image.vx() << ", vy = " << cgal_image.vy() << ", vz = " << cgal_image.vz();
    LOG(INFO) << "-- dimensions:  " << " x = " << cgal_image.xdim() << ",  y = " << cgal_image.ydim() << ",  z = " << cgal_image.zdim();

    // Defining mesh criteria
    Mesh_domain domain(cgal_image);

    // Mesh criteria
    LOG(INFO) << "Mesh Criteria: ";
    mesh_criteria_global_t global_mesh_crit = this->mesh_criteria_global;
    mesh_criteria_subdomain_map_t subdomain_map = this->mesh_criteria_subdomain_map;
    
    Sizing_field size(global_mesh_crit.cell_size);
    
    // meshcriteria for subdomains
    for (mesh_criteria_subdomain_map_t::iterator it = subdomain_map.begin(); it != subdomain_map.end(); it++)
    {
      int subddomain_id = it->first;
      mesh_criteria_subdomain_t  subdomain_mesh_crit = it->second;
      size.set_size(
		    subdomain_mesh_crit.cell_size,
                    subdomain_mesh_crit.dimension,
                    subddomain_id
                     );
      LOG(INFO) << "-- Cell Size subdomain "  << subddomain_id << " : " << subdomain_mesh_crit.cell_size;
    }
   
    LOG(INFO) <<   "-- Facet angle:            " << global_mesh_crit.facet_angle;
    LOG(INFO) <<   "-- Facet size:             " << global_mesh_crit.facet_size;
    LOG(INFO) <<   "-- Facet distance:         " << global_mesh_crit.facet_distance;
    LOG(INFO) <<   "-- Cell-radius-edge-ratio: " << global_mesh_crit.cell_radius_edge_ratio;
    LOG(INFO) <<   "-- Cell size:              " << global_mesh_crit.cell_size;

    // Mesh criteria
    Mesh_criteria criteria(
        facet_angle    		= global_mesh_crit.facet_angle,
        facet_size     		= global_mesh_crit.facet_size,
        facet_distance 		= global_mesh_crit.facet_distance,
        cell_radius_edge_ratio 	= global_mesh_crit.cell_radius_edge_ratio,
        cell_size 		= size
    );
    
    
    LOG(INFO) << "Creating mesh...";
    C3t3 c3t3 = CGAL::make_mesh_3<C3t3>(domain, criteria);

    LOG(INFO) << "Storing Mesh to local cgal c3t3 mesh instance...";
    this->cgal_c3t3_from_image_domain = c3t3;

    LOG(INFO) << "Converting cgal mesh to vtkUnstructuredGrid...";
    vtkSmartPointer<vtkUnstructuredGrid> unstructuredGrid = c3t3_subdomain_to_vtk_unstructured_grid(c3t3);

    this->unstructured_grid_processed = unstructuredGrid;
}


void MeshTool::operationCreateMeshFromMesh()
{ 
  LOG(DEBUG)<< "trying to load grid from : " <<  this->configuration_params["path_to_input_file"];
 //! 1) Load mesh
 vtkSmartPointer<vtkUnstructuredGrid> vtk_ugrid = hfvtk::loadVtkUnstructuredGrid(this->configuration_params["path_to_input_file"]);
 //! 2) Convert mesh to image
 double spacing[3];
 spacing[0] = 1;
 spacing[1] = 1;
 spacing[2] = 1;
 vtkSmartPointer<vtkImageData> vtk_image = hfvtk::createVtkImageFromVtkUnstructuredGrid(vtk_ugrid, spacing);
 std::string path_to_temp_image = "temp_image.vti";
 hfvtk::writeVtkImageData(vtk_image, path_to_temp_image);
 //! 3) Run operationCreateMeshFromImage with modified input path
 this->configuration_params["path_to_input_file"] = path_to_temp_image;
 operationCreateMeshFromImage();
}



void MeshTool::performOperation()
{
  if (this->configuration_params["operation_mode"] == "image_to_mesh")
  {
    operationCreateMeshFromImage();
  }
  else if (this->configuration_params["operation_mode"] == "mesh_to_mesh")
  {
    operationCreateMeshFromMesh();
  }
  else
  {
    LOG(DEBUG) << "operation mode '"<<this->configuration_params["operation_mode"]<< "' invalid."; 
  }
}





// PRIVATE MEMBER FUNCTIONS
  
} // end namespace mt