 
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
#=====================================================================
# Import of own modules, classes, functions
#=====================================================================

hostname = socket.gethostname()
if hostname == "danabl-ISTB":
    repository_base_path = "/home/danabl/SoftwareDevelopment/cornea_meshing"
elif hostname == "istb-skull.unibe.ch":
    repository_base_path = "/home/abler/cornea_meshing"

python_base_path = os.path.join(repository_base_path, "pre-post-processing")
os.sys.path.append(python_base_path)
import CommonHelperFunctions as chf
import FemModel as fem

#==============================================================================
# PATH SETTINGnS
#==============================================================================
if hostname == "danabl-ISTB":
    data_base_path = "/data/ownCloud/Projects/Cornea/CorneaAblation/corneas/"
elif hostname == "istb-skull.unibe.ch":
    data_base_path = "/home/abler/CorneaAblation/corneas/"

model_name     = "pre_OD_cuenca"

data_dir    = os.path.join(data_base_path, "dir_" + model_name)
data_input  = os.path.join(data_dir)
data_output = os.path.join(data_dir)

file_name        = model_name
path_to_vtu_in   = os.path.join(data_input, file_name + ".vtu")
path_to_vtu_out  = os.path.join(data_input, file_name + "_fem.vtu")


#=====================================================================
# PROGRAMME
#=====================================================================

#-- Try to get z_apex from 'params.py' file
os.sys.path.append(data_dir)
import params

mesh = chf.read_vtk_data(path_to_vtu_in)
model = fem.FemModel("test_model", mesh)
model.cone_apex_z = params.params["z_apex"]

#== ABAQUS
path_to_abq_file = os.path.join(data_output, file_name + ".inp")
model.export_to_abaqus(path_to_abq_file)

