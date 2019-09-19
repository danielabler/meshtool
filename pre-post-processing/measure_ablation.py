 
# -*- coding: utf-8 -*-
"""
Created on Fri Nov  6 10:01:03 2015

@author: danabl
"""

#==============================================================================
# General Imports
#==============================================================================
import os
import vtk



def get_file_extension(path_to_file):
    # need to test for existence of '.'
    # return None of component has no file extension
    file_name = os.path.split(path_to_file)[-1]
    file_name_split = file_name.split(".")
    if file_name_split[-1] == file_name:
        # there is no '.' in file_name
        return None
    else:
        return file_name_split[-1]
        
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
        print "Path does not exist"

    
def get_surface(mesh):
    surfaceFilter = vtk.vtkDataSetSurfaceFilter()
    surfaceFilter.SetInputData(mesh)
    surfaceFilter.Update()
    surface = surfaceFilter.GetOutput()    
    return surface
    
def get_center_z_positions(mesh):
    surface_mesh = get_surface(mesh)
    # build tree    
    tree = vtk.vtkOBBTree()
    tree.SetDataSet(surface_mesh)
    tree.BuildLocator()
    # line along z through x=0, y=0
    line_p0 = [0.0, 0.0, -10.0] 
    line_p1 = [0.0, 0.0, 10.0]
    # intersect surface mesh with line
    intersect_points = vtk.vtkPoints()
    tree.IntersectWithLine(line_p0, line_p1, intersect_points, None)
    # extract intersection points
    pointsVTKIntersectionData = intersect_points.GetData()
    noPointsVTKIntersection = pointsVTKIntersectionData.GetNumberOfTuples()
    pointsIntersection = []
    for idx in range(noPointsVTKIntersection):
        _tup = pointsVTKIntersectionData.GetTuple3(idx)
        pointsIntersection.append(_tup)
    return pointsIntersection
    
#=====================================================================
# Import of own modules, classes, functions
#=====================================================================
base_path = "/data/ownCloud/MY_PROJECTS/BrainTumourModeling/brain_tumour_modeling_software"
os.sys.path.append(base_path)


test_data_dir = os.path.join("/data/ownCloud/MY_PROJECTS/CorneaMeshing")
test_data_input  = os.path.join(test_data_dir)
test_data_output = os.path.join(test_data_dir)

path_to_vtu_reference = os.path.join(test_data_input, "mesh_cornea_50k.vtu")
mesh_reference  = read_vtk_data(path_to_vtu_reference)
intersection_points_reference = get_center_z_positions(mesh_reference)
    
for i in range(1,8):    
    path_to_vtu_ablated   = os.path.join(test_data_input, "mesh_ablated_%i_50k.vtu"%(i))
    mesh_ablated    = read_vtk_data(path_to_vtu_ablated)
    intersection_points_ablated   = get_center_z_positions(mesh_ablated)
    # difference 
    diff_anterior  = - (intersection_points_reference[0][2] - intersection_points_ablated[0][2])
    diff_posterior = - (intersection_points_reference[1][2] - intersection_points_ablated[1][2])

    print("Mesh ID=%i"%(i))
    print("Ablation depth = %f"%(diff_anterior))
    print("Difference posterior = depth = %f"%(diff_posterior))
    print("------------")