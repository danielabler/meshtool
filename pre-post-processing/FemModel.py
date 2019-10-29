#=====================================================================
# General Imports
#=====================================================================

import os
import numpy as np
import vtk
import pandas as pd
from vtk.numpy_interface import dataset_adapter as dsa
from vtk.util import numpy_support

#=====================================================================
# Import of own modules, classes, functions
#=====================================================================

import CommonHelperFunctions as chf

#=====================================================================
# FemModel
#=====================================================================

class FemModel():

    def __init__(self, model_name, mesh):
        self.model_name = model_name
        self.mesh = vtk.vtkUnstructuredGrid()
        self.mesh.DeepCopy(mesh)
        self.mesh_wrapped = dsa.WrapDataObject(self.mesh)
        self.material_prop_key      = "MAT"
        self.boundary_condition_key = "BC"
        self.cell_type_vtk_to_abq_map = {10 : "C3D4H" }  
        self.cone_apex_z = 0.0

    def export_to_abaqus(self, path_to_file=None, ):
        self.path_to_abq_file = path_to_file
        chf.ensure_dir_exists(path_to_file)
        # -- id maps
        # -- FILE SETTINGS
        abq_file = open(path_to_file, 'w')
        self.__write_abq_comment(abq_file, ["Abaqus input file", "exported from VKT"])
        # -- 1) Write Header
        self.__write_abq_header(abq_file)
        # -- INPUT statement
        self.__write_abq_include(abq_file, "./parameters.inp")
        # -- 2) Write Nodes
        self.__write_abq_comment(abq_file, ["NODE DEFINITIONS"])
        self.__write_abq_nodes(abq_file, "./sets/nodes.inp")
        # -- 3) write Elements
        # -- 3.1) ALL ELEMENTS
        self.__write_abq_comment(abq_file, ["ELEMENT DEFINITIONS"])
        all_element_vtk_ids = range(self.mesh.GetNumberOfCells())
        self.__write_abq_element_set_2(abq_file, "ELSET_ALL", all_element_vtk_ids, "./sets/elements.inp")
        # -- ELEMENTSETs, NODESETs, SURFACES go into sets.inp
        rel_path_to_sets_file = "./sets/sets.inp"
        abs_path_to_sets_file = self.__create_abs_from_rel_path( self.path_to_abq_file, rel_path_to_sets_file)
        chf.ensure_dir_exists(abs_path_to_sets_file)
        sets_file = open(abs_path_to_sets_file, 'w')
        self.__write_abq_include(abq_file, rel_path_to_sets_file)
        #-- 3.) Elsets for cornes substructures
        elsets = chf.get_unique_integer_list_from_vtu_array(self.mesh, "cell_array", "ElementBlockIds")
        elset_dict = {1 : "cornea", 3: "lenticule"}
        for elsetid in elsets:
            elset_name = "ELSET_%s"%(elset_dict[elsetid])
            print("Creating ElementSets for '%s'"%elset_name)
            selected_cell_vtk_ids = np.where(self.mesh_wrapped.CellData["ElementBlockIds"]==elsetid)[0]
            self.__write_abq_set(sets_file, elset_name, selected_cell_vtk_ids, set_type="element")
            #self.__write_abq_element_set_2(abq_file, elset_name, selected_cell_vtk_ids, "./sets/elements.inp")
        self.__write_abq_comment(abq_file, ["SET DEFINITIONS"])
        #-- 3.3) Nodesets for groups
        # rim_group, id==3 -> rim with all nodes on rim belonging to rim
        # surface group, id == 1 -> anterior surface, with rim nodes belonging to surface
        # surface group, id == 2 -> posterior surface, with rim nodes belonging to surface
        for group_name in ["rim_group", "surface_group"]:
            groups = chf.get_unique_integer_list_from_vtu_array(self.mesh, "point_array", group_name)
            print("Creating Nodesets for goups: %s"%groups)
            group_elset_dict = {}
            for group_id in groups:
                group_id_full = "%s_%03i"%(group_name, group_id)
                node_set_name = "NODES_%s_%03i"%(group_name, group_id)
                groups_name    = "GROUP_%s_%03i"%(group_name, group_id)
                selected_point_vtk_ids = np.where(self.mesh_wrapped.PointData[group_name]==group_id)[0]
                self.__write_abq_set(sets_file, node_set_name, selected_point_vtk_ids, set_type="node")
                group_temp_dict = {}
                group_temp_dict["node_set_name"] = node_set_name
                group_temp_dict["group_name"] = groups_name
                group_elset_dict[group_id_full] = group_temp_dict
        print(group_elset_dict)
        # -- 3.4) Write Surfaces for
        #    - "surface_group" -> 1, 2
        #    - "rim_group" -> 3
        # 1) Surface_group == 1
        self.__write_abq_comment(sets_file, ["SURFACES"])
        print("Creating surface groups ... (may take some time)")
        group_name = "surface_group"
        group_id = 1
        print(" ---- group name: %s, group id: %d -> anterior surface"%(group_name, group_id))
        surface_node_ids = np.where(self.mesh_wrapped.PointData[group_name]==group_id)[0]
        print(surface_node_ids)
        surface = self.__get_surfaces(surface_node_ids)
        self.__write_abq_surface(sets_file, surface, "anterior_surface")
        group_id = 2
        print(" ---- group name: %s, group id: %d -> posterior surface"%(group_name, group_id))
        surface_node_ids = np.where(self.mesh_wrapped.PointData[group_name]==group_id)[0]
        surface = self.__get_surfaces(surface_node_ids)
        self.__write_abq_surface(sets_file, surface, "posterior_surface")
        group_name = "rim_group"
        group_id = 3
        print(" ---- group name: %s, group id: %d -> rim surface"%(group_name, group_id))
        surface_node_ids = np.where(self.mesh_wrapped.PointData[group_name]==group_id)[0]
        surface = self.__get_surfaces(surface_node_ids)
        self.__write_abq_surface(sets_file, surface, "rim_surface")
        # -- close sets file
        sets_file.close()

        # -- 4) write Material Properties
        self.__write_abq_comment(abq_file, ["MATERIAL"])
        abq_file.write("**\n")
        self.__write_abq_include(abq_file, "./sets/orientation_cornea.inp")
        abq_file.write("**\n")
        abq_file.write("*SOLID SECTION, ELSET=ELSET_ALL, MATERIAL=CORNEA_MAT, ORIENTATION=OR_CORNEA_SUP\n")
        abq_file.write("**\n")
        abq_file.write("*MATERIAL, NAME=CORNEA_MAT\n")
        abq_file.write("**\n")
        abq_file.write("*Anisotropic Hyperelastic, user, local direction=2, formulation=invariant, type=incompressible, properties = 16\n")
        abq_file.write("<C10>,<C01>,<C20>,<C11>,<C02>,<D1>,<D2>,<K1>,\n")
        abq_file.write("<K2>,<K3>,<K4>,<KAPPA1>,<KAPPA2>,<DVOL>,<I40>,<I60>\n")
        abq_file.write("**\n")
        self.__write_abq_comment(abq_file, ["*Anisotropic Hyperelastic, holzapfel, local direction=2", "<C10>,  <DVOL>,  <K1>,  <K2>, <KAPPA>"])
        # -- DENSITY
        self.__write_abq_comment(abq_file, ["*Density", "<DensityMat>"])
        # -- MATERIAL DAMPING
        self.__write_abq_comment(abq_file, ["*Damping, alpha=<AlphaMat>"])
        # -- BOUDNARY DEFINITION
        surface_set_rim_name = "NODES_rim_group_003" # nodes rim, including those of ant/post surface
        self.__write_abq_comment(abq_file, ["BOUNDARY DEFINTION"])
        abq_file.write("*Map solution, unbalanced stress=ramp\n")
        abq_file.write("**\n")
        abq_file.write("*TRANSFORM, NSET=%s, TYPE=S\n"%(surface_set_rim_name))
        abq_file.write("0,0,%f,0,0,0\n"%(self.cone_apex_z))
        abq_file.write("**\n")
        abq_file.write("*Boundary\n")
        abq_file.write("%s, 2, 3, 0.0\n"%(surface_set_rim_name))

        # -- STEP
        self.__write_abq_comment(abq_file, ["STEP"])
        abq_file.write("*Step, name=Equilibrium, nlgeom=YES, inc=1000\n")
        abq_file.write("*Static\n")
        abq_file.write("<StepMin>, 1., <StepTol>, <StepMax>\n")
        abq_file.write("**\n")
        # -- LOADS
        self.__write_abq_comment(abq_file, ["LOADS"])
        abq_file.write("**\n")
        abq_file.write("*Dsload\n")
        abq_file.write("POSTERIOR_SURFACE, P, <IOP>\n")
        abq_file.write("** Name: IOP   Type: Pressure\n")
        abq_file.write("**\n")
        # -- OUTPUTS
        self.__write_abq_comment(abq_file, ["OUTPUT REQUESTS"])
        abq_file.write("*Restart, write, frequency=0\n")
        self.__write_abq_comment(abq_file, ["FIELD OUTPUT: F-Output-1"])
        abq_file.write("*Output, field, variable=PRESELECT\n")
        self.__write_abq_comment(abq_file, ["HISTORY OUTPUT: H-Output-1"])
        abq_file.write("*Output, history, variable=PRESELECT\n")
        # -- End Step
        abq_file.write("*End Step\n")
        # -- Close File
        abq_file.close()                 
        

    def __add_node_bc(self, bc_type_name, bc_selection):
        node_id_list = chf.get_node_ids_from_domain_surface(self.mesh, bc_selection)
        point_array = np.zeros(self.mesh.GetNumberOfPoints())
        point_array[node_id_list] = 1
        name = self.boundary_condition_key + "_" + bc_type_name
        self.mesh_wrapped.PointData.append(point_array, name)
            
    def __match_surface(self, n_1, n_2, n_3, n_4):
        """
        identifies facet in surface from pattern
        e.g. (True, False, True, True) means that nodes 1, 3, 4 are part of surface, node 2 is not
        """
        surface = None
        if (n_1 and n_2 and n_3 and not n_4):
            surface = "S1"
        elif (n_1 and n_2 and n_4 and not n_3):
            surface = "S2"
        elif (n_2 and n_3 and n_4 and not n_1):
            surface = "S3"
        elif (n_3 and n_4 and n_1 and not n_2):
            surface = "S4"
        else:
            #print("None applies")
            pass
        return surface

           

    def save_fem_model(self, path_to_file):
        chf.write_vtk_data(self.mesh, path_to_file)
        
        
    def __create_connectivity_table(self):
        print("-- creating connectivity table")
        connect_pd = pd.DataFrame(columns=("cell_id", "node_1", "node_2", "node_3", "node_4"))
        cnt = 0
        n = self.mesh.GetNumberOfCells()
        for _cell_id in range(n):
            cell = self.mesh.GetCell(_cell_id)
            id_list = cell.GetPointIds()
            point_id_list = []
            for _point_id in range(id_list.GetNumberOfIds()):
                point_id = id_list.GetId(_point_id)
                point_id_list.append(point_id)
            connect_pd.loc[_cell_id] = [_cell_id] + point_id_list
            cnt = cnt + 1
            if cnt%5000==0:
                print("- %i / %i"%(cnt, n))
        connect_pd = connect_pd.astype(int)
        self.connectivity_table = connect_pd

    def __get_surface_elements(self, id_list):
        if not hasattr(self, "connectivity_table"):
            self.__create_connectivity_table()
        # Identify cell ids that contain at least one surface node
        cell_id_list = []
        for surface_node_id in id_list:
            query_exp = "node_1==%i or node_2==%i or node_3==%i or node_4==%i"%(surface_node_id, surface_node_id, surface_node_id, surface_node_id)
            cell_id_list = cell_id_list + list(self.connectivity_table.query(query_exp)["cell_id"].values)
            #print(query_exp, cell_id_list)
        unique_cell_ids = list(set(cell_id_list))
        return unique_cell_ids

    def __get_surfaces(self, id_list):
        # reduce to those elements that contain at least 1 node
        unique_cell_ids = self.__get_surface_elements(id_list)
        # identify those cells that contain 3 surface nodes
        surface_pd = pd.DataFrame(columns=("cell_id", "surface"))
        sel_connect_pd = self.connectivity_table.loc[unique_cell_ids]
        print(sel_connect_pd)
        for index, cell_id, node_1, node_2, node_3, node_4 in sel_connect_pd.itertuples():
            node_1_in = node_1 in id_list
            node_2_in = node_2 in id_list
            node_3_in = node_3 in id_list
            node_4_in = node_4 in id_list
            #print("info: ", index, cell_id, node_1, node_2, node_3, node_4 )
            surface = self.__match_surface(node_1_in, node_2_in, node_3_in, node_4_in)
            if not surface==None:
                surface_pd.loc[cell_id] = [cell_id, surface]
            #print("surface: ",surface)
        surface_pd.cell_id = surface_pd.cell_id.astype(int)
        return surface_pd
        
    def __write_abq_surface(self, abq_file, surface_pd, name):
        # -- write element sets (1 for each suface type)
        surface_types = surface_pd.surface.unique()
        element_set_name = "elset_" + name + "_all"
        self.__write_abq_element_set(abq_file, element_set_name, surface_pd.cell_id.values)   
        for surface_type in surface_types:
            element_set_name = "elset_" + name + "_" + surface_type
            selected_elements = surface_pd[surface_pd.surface == surface_type].cell_id.values
            self.__write_abq_element_set(abq_file, element_set_name, selected_elements)            
        # -- write surfaces
        abq_file.write("**\n")
        abq_file.write("*SURFACE, NAME=%s, TYPE=ELEMENT\n"%(name))
        for surface_type in surface_types:
            element_set_name = "elset_" + name + "_" + surface_type
            abq_file.write("%s, %s\n"%(element_set_name, surface_type))
        abq_file.write("**\n")
      
        
    def __write_abq_output_request(self, abq_file):
        abq_file.write("**\n")
        abq_file.write("*EL FILE, FREQ=1 \n")
        abq_file.write("TEMP, IVOL\n")
        abq_file.write("**\n")
        abq_file.write("*OUTPUT, FIELD, FREQ=1 \n")
        abq_file.write("*ELEMENT OUTPUT\n")
        abq_file.write("TEMP, NFLUX, EVOL, THE, ELEN, ELSE, S, LE, HFL\n")
        abq_file.write("*NODE OUTPUT\n")
        abq_file.write("COORD, U, RF, CF, NT, RFL\n")
        abq_file.write("**\n")
        
            
        
    def __convert_vtk_id_to_abq_id(self, vtk_id):
        return vtk_id + 1
        
    def __write_abq_header(self, abq_file):
        abq_file.write("**\n")
        abq_file.write('*HEADING \n')
        abq_file.write("**\n")
         

    def __write_abq_comment(self, abq_file, text_list):
        abq_file.write("**\n")
        for line in text_list:
            abq_file.write("** %s\n"%(line))
        abq_file.write("**\n")
        
    def __write_abq_nodes(self, abq_file, rel_path_to_target_file = None):
        print("-- Writing  nodes ")
        abq_file.write("**\n")
        abq_file.write('*NODE \n')
        if rel_path_to_target_file == None:
            target_file = abq_file
        else:
            print("                ... to %s "%(rel_path_to_target_file))
            abs_path_to_target_file = self.__create_abs_from_rel_path( self.path_to_abq_file, rel_path_to_target_file)
            chf.ensure_dir_exists(abs_path_to_target_file)
            target_file = open(abs_path_to_target_file, 'w')
            self.__write_abq_include(abq_file, rel_path_to_target_file)
        n_points = self.mesh.GetNumberOfPoints()
        for vtk_node_id in range(n_points):
            abq_node_id = self.__convert_vtk_id_to_abq_id(vtk_node_id)
            point = self.mesh.GetPoint(vtk_node_id)
            x, y, z = point
            target_file.write("%i, %f, %f, %f\n"%(abq_node_id, x, y, z))
        abq_file.write("**\n")
            
    def __write_abq_include(self, abq_file, input_string):
        abq_file.write("**\n")
        abq_file.write('*INCLUDE, Input=%s \n'%(input_string))
        abq_file.write("**\n")            

    def __write_abq_element_set_2(self, abq_file, element_set_name, selected_element_vtk_ids, rel_path_to_target_file = None):
        print("-- Writing element set '%s' "%(element_set_name))
        cell_types = set(self.mesh_wrapped.CellTypes[selected_element_vtk_ids])
        if len(cell_types) == 1:
            cell_type = cell_types.pop()
            if cell_type in self.cell_type_vtk_to_abq_map.keys():                    
                element_set_type = self.cell_type_vtk_to_abq_map[cell_type]
                self.__write_abq_element_set_with_nodes(abq_file, element_set_name, element_set_type, selected_element_vtk_ids, rel_path_to_target_file)                    
            else:
                element_set_type = cell_types
                print("No mapping for cell type %s"%(cell_type))
        else:
            print("None or mixed element types (n=%i) in element set"%(len(cell_types)))
        return element_set_type


    def __write_abq_element_set_with_nodes(self, abq_file, element_set_name, element_set_type, selected_element_vtk_ids, rel_path_to_target_file = None):
        element_type =  element_set_type
        abq_file.write("**\n")
        abq_file.write("*ELEMENT, ELSET=%s, TYPE=%s \n"%(element_set_name, element_type))
        if rel_path_to_target_file == None:
            target_file = abq_file
        else:
            print("                ... to %s "%(rel_path_to_target_file))
            abs_path_to_target_file = self.__create_abs_from_rel_path( self.path_to_abq_file, rel_path_to_target_file)
            chf.ensure_dir_exists(abs_path_to_target_file)
            target_file = open(abs_path_to_target_file, 'w')
            self.__write_abq_include(abq_file, rel_path_to_target_file)

        for vtk_element_id in selected_element_vtk_ids:
            abq_element_id = self.__convert_vtk_id_to_abq_id(vtk_element_id)
            cell  = self.mesh.GetCell(vtk_element_id)
            cell_pts_vtkIdlist = cell.GetPointIds()
            n_pts = cell_pts_vtkIdlist.GetNumberOfIds()
            target_file.write("%i"%(abq_element_id))
            for i in range(n_pts):
                vtk_node_id = cell_pts_vtkIdlist.GetId(i)
                abq_node_id = self.__convert_vtk_id_to_abq_id(vtk_node_id)
                target_file.write(", %i"%(abq_node_id))
            target_file.write("\n")
        abq_file.write("**\n")
                
                
    def __write_abq_element_set(self, abq_file, element_set_name, selected_element_vtk_ids):
        abq_file.write("**\n")
        abq_file.write("*ELSET, ELSET=%s \n"%(element_set_name))
        item_in_line_count = 0
        for vtk_element_id in selected_element_vtk_ids:
            abq_element_id = self.__convert_vtk_id_to_abq_id(vtk_element_id)
            item_in_line_count = item_in_line_count +1;
            if (item_in_line_count == 5):
              abq_file.write("\n")
              item_in_line_count = 0;
            abq_file.write("%i,"%(abq_element_id))
        abq_file.write("\n")
        abq_file.write("**\n")
        
    def __write_abq_set(self, abq_file, set_name, selected_vtk_ids, set_type="element"):
        abq_file.write("**\n")
        if set_type == "element":
            abq_file.write("*ELSET, ELSET=%s \n"%(set_name))
        elif set_type == "node":
            abq_file.write("*NSET, NSET=%s \n"%(set_name))
        item_in_line_count = 0
        for vtk_id in selected_vtk_ids:
            abq_id = self.__convert_vtk_id_to_abq_id(vtk_id)
            item_in_line_count = item_in_line_count +1;
            if (item_in_line_count == 5):
              abq_file.write("\n")
              item_in_line_count = 0;
            abq_file.write("%i,"%(abq_id))
        abq_file.write("\n")
        abq_file.write("**\n")
        
    def __write_abq_node_sets(self, abq_file, node_set_name, selected_nodes_vtk_ids):
        abq_file.write("**\n")
        abq_file.write("*Nset, nset=%s \n"%(node_set_name))
        item_in_line_count = 0
        for vtk_node_id in selected_nodes_vtk_ids:
            abq_node_id = self.__convert_vtk_id_to_abq_id(vtk_node_id)
            item_in_line_count = item_in_line_count +1;
            if (item_in_line_count == 5):
              abq_file.write("\n")
              item_in_line_count = 0;
            abq_file.write("%i,"%(abq_node_id))
        abq_file.write("\n")
        abq_file.write("**\n")
                
                
    def __write_abq_material_as_elset(self, abq_file, elset_name, material_name):
        abq_file.write("**\n")
        abq_file.write("*SOLID SECTION, ELSET=%s, MATERIAL=%s \n"%(elset_name, material_name))
        abq_file.write("1 \n")
        abq_file.write("**\n")
            
    def __write_abq_material_properties(self, abq_file, material_name, property_dict):
        abq_file.write("**\n")
        abq_file.write("*MATERIAL,NAME=%s \n"%(material_name))
        # -- MASS
        if property_dict.has_key("Density"):            
            abq_file.write("*Density \n")         
            abq_file.write("%f \n"%(property_dict["Density"]))
        # -- HEAT TRANSFER
        if property_dict.has_key("Thermal_Conductivity"):
            abq_file.write("*Conductivity \n")           
            abq_file.write("%f, \n"%(property_dict["Thermal_Conductivity"]))
        if property_dict.has_key("Specific_Heat"):
            abq_file.write("*Specific Heat \n")           
            abq_file.write("%f, \n"%(property_dict["Specific_Heat"]))
        # -- HETVAL SUBROUTINE
        if property_dict.has_key("HETVAL"):
            abq_file.write("*Heat Generation \n")           
        # -- MECHANICS
        if property_dict.has_key("Temp_dep_Material_Props"):
            abq_file.write("*Elastic \n")           
            for temp, prop_dict in property_dict["Temp_dep_Material_Props"].items():
                poisson_factor_ref_mat = prop_dict["Poisson_Factor"]
                e_module_ref_mat       = prop_dict["E_Module"]
                poisson_ratio = self.material_properties["Poisson_Factor"][poisson_factor_ref_mat]
                e_module      = self.material_properties["E_Module"][e_module_ref_mat]
                abq_file.write("%f, %f, %f \n"%(e_module, poisson_ratio, temp))
        else:
            if property_dict.has_key("E_Module") & property_dict.has_key("Poisson_Factor"):
                abq_file.write("*Elastic \n")           
                abq_file.write("%f, %f \n"%(property_dict["E_Module"], property_dict["Poisson_Factor"]))

        abq_file.write("**\n")
 

    def __create_abs_from_rel_path(self, main_file_path, sub_file_rel_path):
        sub_file_abs_path = os.path.join( os.path.dirname(main_file_path) , os.path.normpath(sub_file_rel_path))
        return sub_file_abs_path
