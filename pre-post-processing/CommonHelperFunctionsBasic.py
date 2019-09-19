"""
Created in September 2016

@author: djs.abler@gmail.com


Collection of functions used by other scripts.
The functions included here should only use standard python libraries,
i.e. no external dependencies, such as vtk, pandas etc
"""

# ==============================================================================
# IMPORTS
# ==============================================================================

import os


# ==============================================================================
# FILE MANAGEMENT
# ==============================================================================


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


def ensure_dir_exists(path):
    if get_file_extension(path) == None:
        # assume that 'path' is directory, add trailing '/'
        path = path + '/'
    if os.path.exists(os.path.dirname(path)):
        return True
    else:
        os.makedirs(os.path.dirname(path))
        return False
