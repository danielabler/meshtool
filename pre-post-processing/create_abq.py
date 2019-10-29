 
# -*- coding: utf-8 -*-
"""
Created on Fri Nov  6 10:01:03 2015

@author: danabl
"""

#==============================================================================
# General Imports
#==============================================================================
import os
import sys, getopt
import socket

#==============================================================================
# PATH SETTINGS
#==============================================================================

#=====================================================================
# PROGRAMME
#=====================================================================


def main(argv):

    try:
        opts, args = getopt.getopt(argv, "hi:o:p:",["in=", "out=", "pythonpath="])
    except getopt.GetoptError:
        print("create_abq.py -i <path_to_vtu_file> -o <path_to_abq_file>")
        print("or")
        print("create_abq.py -i <path_to_vtu_file> -o <path_to_abq_file> -p <path_to_python_functions>")
        sys.exit(2)

    path_file_vtu = ""
    path_file_abq = ""
    #path_python = "/home/dabler/Documents/SoftwareDevelopment/meshtool_github/pre-post-processing"
    path_python = "/home/mesher/software/MESHTOOL_source/pre-post-processing"

    for opt, arg in opts:
        if opt == "-h":
            print("create_abq.py -i <path_to_vtu_file> -o <path_to_abq_file>")
            print("(assumes that python functions located at '%s')"%path_python)
            print("or")
            print("create_abq.py -i <path_to_vtu_file> -o <path_to_abq_file> -p <path_to_python_functions>")
            sys.exit()
        elif opt in ("-i", "--in"):
            path_file_vtu = arg.strip()
        elif opt in ("-o", "--out"):
            path_file_abq = arg.strip()
        elif opt in ("-p", "--pythonpath"):
            path_python = arg.strip()

    os.sys.path.append(path_python)
    import CommonHelperFunctions as chf
    import FemModel as fem

    if path_file_vtu == "" and path_file_abq == "":
        print("create_abq.py -i <path_to_vtu_file> -o <path_to_abq_file>")
        print("or")
        print("create_abq.py -i <path_to_vtu_file> -o <path_to_abq_file> -p <path_to_python_functions>")
    elif path_file_vtu == "":
        print("Path to input file missing: '-i <path_to_vtu_file>'")
    elif path_file_abq == "":
        print("Path to output file missing: '-o <path_to_abq_file>'")
    else:
        mesh = chf.read_vtk_data(path_file_vtu)
        model = fem.FemModel("cornea", mesh)
        model.export_to_abaqus(path_file_abq)



if __name__ == "__main__":
    main(sys.argv[1:])
