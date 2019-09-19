"""
Created in September 2016

@author: djs.abler@gmail.com

"""


#=====================================================================
# General Imports
#=====================================================================

import vtk
from vtk.numpy_interface import dataset_adapter as dsa
from vtk.numpy_interface import dataset_adapter as dsa
from vtk.numpy_interface import algorithms as algs
from vtk.util import numpy_support

import os
import numpy as np
import shutil

#=====================================================================
# Import of own modules, classes, functions
#=====================================================================

from CommonHelperFunctionsBasic import *



# ==============================================================================
# COMMAND LINE
# ==============================================================================

def create_cmd_line_from_dict(path_to_executable, cmd_arg_dict, nameless_args_list=None):
    cmd_string = path_to_executable
    if not nameless_args_list == None:
        for item in nameless_args_list:
            cmd_string = "%s %s" % (cmd_string, item)
    for cla, clavalue in cmd_arg_dict.items():
        cmd_string = "%s %s %s" % (cmd_string, cla, clavalue)
    return cmd_string

# ==============================================================================
# FILE MANAGEMENT
# ==============================================================================

def copy_from_server(base_path_remote, base_path_local, sim_group, model_name, server_settings_dict):
    remote_dir = os.path.join(base_path_remote, sim_group, model_name)
    local_dir = os.path.join(base_path_local, sim_group, model_name)
    endings_file_dict = {  # "field_vars"    : "_output-dict.pkl",
        # "connectivity"  : "_connectivity.txt",
        # "nodal_coords"  : "_nodal_coordinates.txt",
        # "nodes"         : "_nodes.txt",
        # "elements"      : "_elements.txt",
        "volumes": "_job_volumes.txt"}
    user_name = server_settings_dict["user_name"]
    server_name = server_settings_dict["server_name"]
    certificate = server_settings_dict["certificate"]
    # make dir on local
    ensure_dir_exists(local_dir)
    # copy all files from remote to local
    for item, ending in endings_file_dict.items():
        file_name = model_name + ending
        path_to_file_remote = os.path.join(remote_dir, file_name)
        path_to_file_local = os.path.join(local_dir, file_name)
        print("Copying file %s from remote:%s to local:%s" % (file_name, path_to_file_remote, path_to_file_local))
        scp_string = "scp" + " -i " + certificate + " " + user_name + "@" + server_name + ":" + path_to_file_remote + " " + path_to_file_local
        print(scp_string)
        os.system(scp_string)


def copy_file_from_server(path_to_file_remote, path_to_file_local, server_settings_dict):
    user_name = server_settings_dict["user_name"]
    server_name = server_settings_dict["server_name"]
    print(
    "Trying to download: '%s' from '%s@%s' to '%s'" % (path_to_file_remote, user_name, server_name, path_to_file_local))
    if server_settings_dict.has_key("certificate"):
        certificate = server_settings_dict["certificate"]
        scp_string = "scp" + " -i " + certificate + " " + user_name + "@" + server_name + ":" + path_to_file_remote + " " + path_to_file_local
        # print scp_string
        os.system(scp_string)
    else:
        print("-- No certificate ... try password authentification")
        if server_settings_dict.has_key("password"):
            password = server_settings_dict["password"]
            scp_string = "sshpass" + " -p \"" + password + "\" scp " + user_name + "@" + server_name + ":" + path_to_file_remote + " " + path_to_file_local
            # print scp_string
            os.system(scp_string)
        else:
            print("-- Password authentication not working -- exit.")
            scp_string = None
    print("--------------------------------------------------------------")
    return scp_string


def copyFile(src, dest):
    try:
        shutil.copy(src, dest)
    # eg. src and dest are the same file
    except shutil.Error as e:
        print('Error: %s' % e)
    # eg. source or destination doesn't exist
    except IOError as e:
        print('Error: %s' % e.strerror)


def copy_files_to_dir(source_dir, target_dir):
    for root, dirs, files in os.walk(source_dir):
        for file in files:
            path_file = os.path.join(root,file)
            shutil.copy2(path_file,target_dir)


# ==============================================================================
# TEXT FILE MANIPULATION
# ==============================================================================

def replaceTextInFile(replacement_dict, path_to_input_file, path_to_output_file):
    with open(path_to_input_file) as infile, open(path_to_output_file, 'w') as outfile:
        for line in infile:
            for src, target in replacement_dict.iteritems():
                line = line.replace(src, target)
            outfile.write(line)
            
            
# ==============================================================================
# VTK
# ==============================================================================

def threshold_vtk_data(data_vtu_in, array_type, array_name, lower_thr=None, upper_thr=None):
    threshold_filter = vtk.vtkThreshold()
    threshold_filter.SetInputData(data_vtu_in)
    if array_type == "cell":
        threshold_filter.SetInputArrayToProcess(0, 0, 0, vtk.vtkDataObject.FIELD_ASSOCIATION_CELLS, array_name)
    elif array_type == "point":
        threshold_filter.SetInputArrayToProcess(0, 0, 0, vtk.vtkDataObject.FIELD_ASSOCIATION_POINTS, array_name)
    if (lower_thr == None) and (upper_thr != None):
        threshold_filter.ThresholdByUpper(upper_thr)
        print("-- select %ss where %s > %s" % (array_type, array_name, upper_thr))
    elif (lower_thr != None) and (upper_thr == None):
        print("-- select %ss where %s < %s" % (array_type, array_name, lower_thr))
        threshold_filter.ThresholdByLower(lower_thr)
    elif (lower_thr != None) and (upper_thr != None):
        print("-- select %ss where %s between %s and %s" % (array_type, array_name, lower_thr, upper_thr))
        threshold_filter.ThresholdBetween(lower_thr, upper_thr)
    threshold_filter.Update()
    data_vtu_out = threshold_filter.GetOutput()
    return data_vtu_out


def compute_vtk_volume(data_vtu_in, method="vtk"):
    volume = 0
    if method == "vtk":
        quality = vtk.vtkMeshQuality()
        for cell_id in range(data_vtu_in.GetNumberOfCells()):
            cell_curr = data_vtu_in.GetCell(cell_id)
            vol_curr = quality.TetVolume(cell_curr)
            volume = volume + vol_curr
    elif method == "abaqus":
        if data_vtu_in.GetCellData().HasArray("EVOL"):
            cell_volume_array = data_vtu_in.GetCellData().GetArray("EVOL")
            for cell_id in range(cell_volume_array.GetNumberOfTuples()):
                vol_curr = cell_volume_array.GetValue(cell_id)
                volume = volume + vol_curr
    print("Volume from '%s': %f" % (method, volume))
    return volume


def read_vtk_data(_path_to_file):
    if os.path.exists(_path_to_file):
        extension = get_file_extension(_path_to_file)
        if extension == 'vtk':
            reader = vtk.vtkDataSetReader()
        elif extension == 'vtp':
            reader = vtk.vtkXMLPolyDataReader()
        elif extension == 'stl':
            reader = vtk.vtkSTLReader()
        elif extension == 'vtu':
            reader = vtk.vtkXMLUnstructuredGridReader()
        reader.SetFileName(_path_to_file)
        reader.Update()
        return reader.GetOutput()
    else:
        print("Path does not exist")


def write_vtk_data(_data, _path_to_file):
    ensure_dir_exists(_path_to_file)
    writer = vtk.vtkXMLDataSetWriter()
    if vtk.VTK_MAJOR_VERSION <= 5:
        writer.SetInput(_data)
    else:
        writer.SetInputData(_data)
    writer.SetFileName(_path_to_file)
    writer.Update()
    writer.Write()


def load_image(path_to_image):
    file_ext = get_file_extension(path_to_image)
    if (file_ext=="vti"):
        print("Opening as  '.vti' file.")
        image_reader = vtk.vtkXMLImageDataReader()
    elif (file_ext=="nii"):
        print("Opening as  '.nii' file.")
        image_reader = vtk.vtkNIFTIImageReader()
    else:
        print("Attempting to load as other vtk image.")
        reader_factory = vtk.vtkImageReader2Factory()
        image_reader = reader_factory.CreateImageReader2(path_to_image)
    image_reader.SetFileName(path_to_image)
    image_reader.Update()
    image = image_reader.GetOutput()
    return image

def write_image(_data, path_to_image):
    #== make sure data type is unsigned char
    cast = vtk.vtkImageCast()
    if vtk.VTK_MAJOR_VERSION <= 5:
        cast.SetInput(_data)
    else:
        cast.SetInputData(_data)
    cast.SetOutputScalarTypeToUnsignedChar()
    cast.Update()
    _data = cast.GetOutput()
    #== write
    file_ext = get_file_extension(path_to_image)
    if (file_ext=="vti"):
        print("Writing as  '.vti' file.")
        image_writer = vtk.vtkXMLImageDataWriter()
    elif (file_ext=="nii"):
        print("Writing as  '.nii' file.")
        image_writer = vtk.vtkNIFTIImageWriter()
        print("VTK seems to change image orientation of NIFTI. Make sure to check image orientation relative to original image")
    elif (file_ext=="mhd" or file_ext=="mha"):
        print("Writing as .mhd/.raw or .mha image.")
        image_writer = vtk.vtkMetaImageWriter()
        image_writer.SetCompression(False)
    else:
        print("No valid image file extension specified!")
    if vtk.VTK_MAJOR_VERSION <= 5:
        image_writer.SetInput(_data)
    else:
        image_writer.SetInputData(_data)
    image_writer.SetFileName(path_to_image)
    #image_writer.Update()
    image_writer.Write()
    print("Image has been saved as %s"%(path_to_image))


def make_vtk_id_list(it):
    vil = vtk.vtkIdList()
    for i in it:
        vil.InsertNextId(int(i))
    return vil



def vtu_has_arrays(mesh_vtu, array_dict):
    mesh_vtu_wrapped = dsa.WrapDataObject(mesh_vtu)
    cell_array_names  = mesh_vtu_wrapped.CellData.keys()
    point_array_names = mesh_vtu_wrapped.PointData.keys()
    missing_arrays_dict = {}
    for name, type in array_dict.items():
        if type == "cell_array":
            if not name in cell_array_names:
                missing_arrays_dict[name] = type
        elif type == "point_array":
            if not name in point_array_names:
                missing_arrays_dict[name] = type
    if len(missing_arrays_dict) != 0:
        status = False
    else:
        status = True
        
    if not status:
        print("Missing arrays: ")
        for name, type  in missing_arrays_dict.items():
            print("-- %s : %s"%(name, type))
    return status, missing_arrays_dict
    
    

def get_node_ids_from_domain_surface(mesh_vtu, domain_id_list):
    subgrid = get_domain_with_node_ids(mesh_vtu, domain_id_list)
    # -- extract surface
    surfaceFilter = vtk.vtkDataSetSurfaceFilter()
    surfaceFilter.SetInputData(subgrid)
    surfaceFilter.Update()
    surface = surfaceFilter.GetOutput()
    # -- extract node ids on surface
    id_mapping_surface_to_orig = {}
    id_mapping_orig_to_surface = {}
    if surface.GetNumberOfPoints() > 0:
        surface_wrapped = dsa.WrapDataObject(surface)
        extracted_point_ids_orig    = surface_wrapped.PointData['originalPointIDs']
        extracted_point_ids_surface = range(algs.shape(extracted_point_ids_orig)[0])
        id_mapping_surface_to_orig = dict(zip(extracted_point_ids_surface, extracted_point_ids_orig))
        id_mapping_orig_to_surface = dict(zip(extracted_point_ids_orig, extracted_point_ids_surface))
    else:
        print("No selection")
    #return id_mapping_surface_to_orig, id_mapping_orig_to_surface
    return list(extracted_point_ids_orig)
    
    
def get_node_ids_from_domain(mesh_vtu, domain_id_list):
    subgrid = get_domain_with_node_ids(mesh_vtu, domain_id_list)
    subgrid_wrapped = dsa.WrapDataObject(subgrid)
    extracted_point_ids_orig = subgrid_wrapped.PointData['originalPointIDs']
    return list(extracted_point_ids_orig)
    

def get_domain_with_node_ids(mesh_vtu, domain_id_list):
    # -- add copy of element point ids to mesh
    pointIdFilter = vtk.vtkIdFilter()
    pointIdFilter.SetInputData(mesh_vtu)
    pointIdFilter.SetPointIds(True)
    pointIdFilter.SetCellIds(False)
    pointIdFilter.SetIdsArrayName("originalPointIDs")
    pointIdFilter.Update()
    mesh = pointIdFilter.GetUnstructuredGridOutput()
    # -- NUMPY - VTK
    mesh_wrapped        = dsa.WrapDataObject(mesh)
    # -- extract element block ids
    material_array_name = "ElementBlockIds"
    element_block_id_array   = mesh_wrapped.CellData[material_array_name]
    # -- create empty array for selection
    element_selection_array_name = "ElementSelection" 
    element_selection_array = np.zeros(algs.shape(element_block_id_array)[0])    
    # -- assign selection id = 1 for elements in domain_id_list
    for domain_id in domain_id_list:
        element_selection_array[ element_block_id_array == domain_id ] = 1
    # -- add selection array to mesh    
    mesh_wrapped.CellData.append(element_selection_array, element_selection_array_name)
    # -- extract parts of mesh that have been selected
    selector = vtk.vtkThreshold()
    selector.SetInputData(mesh)
    selector.SetInputArrayToProcess(0, 0, 0, 
                                 vtk.vtkDataObject.FIELD_ASSOCIATION_CELLS,
                                                element_selection_array_name)
    selector.ThresholdBetween(1, 1)
    selector.Update()
    subgrid = selector.GetOutput()
    return subgrid

def get_unique_integer_list_from_vtu_array(mesh_vtu, array_type, array_name):
    list_unique = []
    mesh_vtu_wrapped = dsa.WrapDataObject(mesh_vtu)
    if array_type == "cell_array":
        data = mesh_vtu_wrapped.CellData
    elif array_type == "point_array":
        data = mesh_vtu_wrapped.PointData
    if array_name in data.keys():
        array = data[array_name]
        list_unique = list(set(array))
    return list_unique
    

    
# ==============================================================================
# NAMING FUNCTIONS
# ==============================================================================
