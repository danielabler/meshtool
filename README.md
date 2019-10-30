# README #

## Introduction ##

* *MeshTool* is a CGAL-based c++ tool for creating FEM meshes. The tool produces VTK 3D meshes (vtkUnstructuredGrid) as output.
 Currently, two applications are supported:
  * Creation of cornea meshes from Zernike-parameters of anterior and posterior surface.
  * Creation of meshes from image segmentations. 
* Folder  <base_dir>\pre-post-processing\ contains scripts for generating Abaqus input files from these VTK meshes.


## Usage ##

### Installation

Quickest way to install is to use the docker or singularity container as explained [here](https://github.com/danielabler/dockerfiles/tree/master/meshtool).
### Commandline Parameters ###

From <base_dir>\bin:

```#!shell

./bin/MeshTool -c <path_to_configuration_file> -m <mode>
```

Where 

* **base_dir** is the installation path of the MeshTool project.
* **path_to_configuration_file** is the path (relative or absolute) to a xml configuration file.
* **mode** is the operational model of the tool, currently one of `'cornea'`, or `'image'`.


### Image Meshing Mode ###

In this mode, `MeshTool` generates a VTK mesh from a segmented image.
An example configuration is provided in `test-data/image/config/image-to-mesh_from-mhd.xml`.
Required parameters are:
* File settings:
  * path_to_input_file
  * path_to_output_file
* Global meshing parameters:
  These are applicable to all subdomains of the mesh and based on
  [CGAL mesh criteria](http://doc.cgal.org/latest/Mesh_3/classCGAL_1_1Mesh__criteria__3.html).    
  * cell_radius_edge_ratio
  * cell_size
  * facet_angle
  * facet_distance
  * facet_size

If multiple labels (image colors) exist in the source image, a subdomain will be created for each label.
Each subdomain can be meshed according to a specific set of parameters. 
These domain-specific meshing parameters are optional:

* cell_size: size of the cell
* dimension: 3
* domain_id: color value of the respective subdomain
  
    
The image meshing command has the following form: 

```#!shell
./bin/MeshTool -c test-data/test.xml -x <base_dir>/src/xml-io/imaging_meshing_schema.xsd -m 'image'
```

where 'base_dir' is the **absolute** path to the installation path of the MeshTool project. 

### Cornea Meshing Mode ###

Cornea and meshing parameters are specified in xml format. The corresponding xml schema and an example xml are in <base_dir>\src\xml-io. 
Paths to the configuration xml-file and the xml-schema file need to be passed to the executable, using commandline options '-c' and '-x', respectively:

From <base_dir>:

```#!shell
./bin/MeshTool -c test-data/test.xml -x <base_dir>/src/xml-io/cornea_meshing_schema.xsd -m 'cornea'
```

Main Configuration options:

* Cornea

    * Anterior/Posterior Surface: List of Zernike coefficients in single index (*<ZernikeCoefficientSingleIndex j='0'>*) or n-m notation (*<ZernikeCoefficient n='1' m='2'>*)
    * Pupil Radius
    * Distance between posterior and anterior surface
    * Information for cornea generation, such as max number of zernike coefficients to be considered for generation, output file, etc
    * Cornea boundary, currently: cone with custom angle, or cylinder with custom radius

* Lenticule:

    * Anterior surface for lenticule
    * Posterior surface for lenticule
    * Lenticule configuration, in particular
        * lenticule radius
        * lenticule cap thickness: distance between anterior cornea surface and anterior lenticule surface
        * lenticule surface distance: distance between anterior and posterior lenticule surface
        * lenticule thickness: only relevant for generation of 3D volume by intersection of 2 volumes; 
          value should greater than  lenticule surface distance
    * Meshing criteria for lenticule

* Mesh Criteria, based on [CGAL mesh criteria](http://doc.cgal.org/latest/Mesh_3/classCGAL_1_1Mesh__criteria__3.html)

### Conversion to Abaqus .inp ###

* <base_dir>\pre-post-processing\create_abaqus.py provides an example configuration for generating Abaqus *.inp input files from the resulting VTK *.vtu meshes.
* <base_dir>\pre-post-processing\create_abq.py provides a command line interface for this conversion:
  ```shell script
  python  <base_dir>\pre-post-processing\create_abq.py -i <path_to_vtu_input> -o <path_to_inp_output> -p <base_dir>\pre-post-processing
  ```


## Development ##

### Dependencies ###

* VTK
* CGAL
* [Xerces library](https://xerces.apache.org/xerces-c/) (needed by [CodeSynthesis XSD](http://www.codesynthesis.com/products/xsd/))
* [CodeSynthesis XSD](http://www.codesynthesis.com/products/xsd/)

==> Adapt DEPENDENCY_DIR and CGAL_DIR, VTK_DIR, XERCESC_ROOT_DIR, XSD_DIR in <base_dir>/src/CMakeLists.txt for your system

### Compilation ###

* go to main directory <base_dir>

```#!shell
mkdir build
cd build

cmake ..
make
make install
```

#### Compilation MacOS ####

* I had difficulties including the full path to boost libraries in the compiled meshtool executable;
  instead boost libraries were linked as if they were present in same directory as the meshtool binary.
* Compiling CGAL against system boost libraries (instead of custom boost installations) solved this issue...
  There should be more robust ways of fixing library paths using 'rpath', but it seems that 'rpath' related
  settings in MacOS are somewhat problematic:
  https://gitlab.kitware.com/cmake/cmake/issues/16589 .
 
* To check libraries linked into executable do:
  ```#!shell

  otool -L path/to/binary
  ```
  
* Installation of package with 'fixup_bundle' does not work. The problem seems to be related to the way apps are
  organized in MacOS. The bundling helper appears to require a certain folder structure to work correctly.
  http://cmake.3232098.n2.nabble.com/fixup-bundles-on-macosx-for-non-app-target-td7591503.html
  
  The executable can still be run from <build-dir>/src/MeshTool 

### Modifying XML schema ###

MeshTool uses the [codesynthesis binding xml binding libraries](http://www.codesynthesis.com/products/xsd/) for xml serialisation and reading. In case xml schemas are modified, binding libraries have to be re-generated. 
See the [CodeSynthesis user guide](http://www.codesynthesis.com/projects/xsd/documentation/cxx/tree/guide/#3) for further information.

* For configuration file:

```#!shell
cd <base_dir>/src/xml-io
path_to_xsd_library/bin/xsd cxx-tree cornea_meshing_schema.xsd
```



## Usage on SKULL ##

MeshTool and dependencies are installed on skull.istb.unibe.ch on user account *abler*:

* dependencies: /home/abler/bms_dependencies/
* repository:   /home/abler/bms_dependencies/cornea_meshing

### Recompile ###

```#!shell
cd /home/abler/cornea_meshing/build
cmake ..
make install
```

### Example Meshing ###

From <base_dir> = '/home/abler/cornea_meshing'.

* For cornea meshing:

  ```#!shell
  ./bin/MeshTool -x /home/abler/cornea_meshing/src/xml-io/cornea_meshing_schema.xsd -c test-data/cornea/config/test_lenticule.xml -m 'cornea'
  ```

* For image meshing:

  ```#!shell
  ./bin/MeshTool -x /home/abler/cornea_meshing/src/xml-io/imaging_meshing_schema.xsd -c test-data/image/config/image-to-mesh_from-mhd.xml -m 'image'
  ```

### Conversion to Abaqus inp ###

If not done already, add vtk to python path, adding the following to your .profile should work:

```#!shell
export PYTHONPATH=/home/abler/bms_dependencies/VTK-6.2.0_build/Wrapping/Python/:/home/abler/bms_dependencies/VTK-6.2.0_build/bin:/home/aler/bms_dependencies/VTK-6.2.0_build/lib:$PYTHONPATH
export LD_LIBRARY_PATH=/home/abler/bms_dependencies/VTK-6.2.0_build/bin:/home/abler/bms_dependencies/VTK-6.2.0_build/lib:$LD_LIBRARY_PATH
```

Activate python virtual-environment in shell:

```#!shell
source /home/abler/cornea_meshing/pyenv/bin/activate  
```

Then, run analysis script:

```#!shell
python /home/abler/cornea_meshing/pre-post-processing/create_abaqus_test.py
```