 
# -*- coding: utf-8 -*-
"""
Created on Fri Nov  6 10:01:03 2015

@author: danabl
"""

#==============================================================================
# General Imports
#==============================================================================
import os
import socket
#==============================================================================
# PATH SETTINGS
#==============================================================================
hostname = socket.gethostname()
if hostname == "danabl-ISTB":
    repository_base_path = "/home/danabl/SoftwareDevelopment/cornea_meshing"
elif hostname == "istb-skull.unibe.ch":
    repository_base_path = "/home/abler/cornea_meshing"

test_data_dir = os.path.join(repository_base_path, "test-data")
test_data_input  = os.path.join(test_data_dir)
test_data_output = os.path.join(test_data_dir)

file_name        = "test"
path_to_vtu_in   = os.path.join(test_data_input, file_name + ".vtu")
path_to_vtu_out  = os.path.join(test_data_input, file_name + "_fem.vtu")

#=====================================================================
# Import of own modules, classes, functions
#=====================================================================
python_base_path = os.path.join(repository_base_path, "pre-post-processing")
os.sys.path.append(python_base_path)
import CommonHelperFunctions as chf
import FemModel as fem

#=====================================================================
# PROGRAMME
#=====================================================================

mesh = chf.read_vtk_data(path_to_vtu_in)
model = fem.FemModel("test_model", mesh)
#== ABAQUS
path_to_abq_file = os.path.join(test_data_output, file_name + "_fem.inp")
model.export_to_abaqus(path_to_abq_file)

