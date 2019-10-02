
# ----------------------------------------
# USAGE:
# ----------------------------------------
# python3 visualize_Porcupine_plot.py pdb_file node_definition_file eigenvector_file output_vis_state_file_name output_vis_state_structural_file_name 

# ----------------------------------------
# PREAMBLE:
# ----------------------------------------

import sys
import MDAnalysis
import numpy as np
from scipy import interpolate
import matplotlib as mpl
import matplotlib.pyplot as plt

pdb_file = sys.argv[1]
node_definition_file = sys.argv[2]
eigenvector_file = sys.argv[3]
vis_state_file_name = sys.argv[4]

# ----------------------------------------
# SETTING COLORS
# ----------------------------------------
color_defs = {}
# colorid 33 is to be used to color eigenvector reps
color_defs[33] = "color change rgb 33 1.0 0.0 0.0\n"        # red

# ----------------------------------------
# VECTOR RADIUSES
# ----------------------------------------
vector_max_magnitude = 5.0  # units of Angstroms
vector_max_radius = 0.5
cone_radius_scale_factor = 5/3.

# ----------------------------------------
# CREATING MDANALYSIS SELECTIONS BY READING IN THE NODE DEFINITION FILE
# ----------------------------------------
u = MDAnalysis.Universe(pdb_file)
u_all = u.select_atoms('all')

pos_list = []
vmd_list = 'serial '
count = 0
with open(node_definition_file,'r') as f:
    for line in f:
        line_list = line.split()
        if line_list[0][0] == '#':
            continue
        temp = u.select_atoms('bynum %s'%(int(line_list[2]))).atoms[0]
        if temp.name != line_list[1]:
            print('Atom selection does not return the same atom as the node_definition file states. Node Index: %d'%(line_list[0]))
            sys.exit()
        pos_list.append(temp.position)
        vmd_list += '%s '%(int(line_list[2]))
        count += 1

pos_list = np.array(pos_list)
nNodes_range = list(range(count))

# ----------------------------------------
# LOADING IN AND ANALYZING THE EIGENVECTOR FILE
# ----------------------------------------
eigenvector = np.loadtxt(eigenvector_file)
node_components = eigenvector.reshape((count,3))
magnitude_scale_factor = vector_max_magnitude/np.max(np.sum(np.square(node_components),axis=1)**0.5)
node_components *= magnitude_scale_factor

component_magnitudes = np.sum(np.square(node_components),axis=1)**0.5   # not a bug since these values should never be negative
magnitudes_min_max = (np.min(component_magnitudes),np.max(component_magnitudes))
print(magnitudes_min_max)

start_middle_end_coords = np.zeros((count,3,3))
start_middle_end_coords[:,0] = pos_list     # setting start of drawn vector to be node position
start_middle_end_coords[:,1] = pos_list + 0.9 * node_components     # setting "middle" of drawn vector
start_middle_end_coords[:,2] = pos_list + node_components   # setting end of drawn vector

cylinder_radius = [vector_max_radius * (i/vector_max_magnitude) for i in component_magnitudes]
cone_radius = [i * cone_radius_scale_factor for i in cylinder_radius]
print(np.min(cylinder_radius),np.max(cylinder_radius),np.min(cone_radius),np.max(cone_radius))

# ----------------------------------------
# CREATING VIS STATE FILE
# ----------------------------------------
with open(vis_state_file_name,'w') as W:
    ### starting lines
    W.write('#!/usr/local/bin/vmd\nset viewplist {}\nset fixedlist {}\n\n')

    ### setting colorids
    W.write('# setting colorid rgb values\n')
    W.write(color_defs[33])

    ### drawing porcupine vectors
    W.write('\n# Draw vectors at the nodes\nmol new\n')
    W.write('draw material AOEdgy\n')
    W.write('graphics top color 34 \n')
    for i in nNodes_range:
        W.write('draw cylinder {' + str(start_middle_end_coords[i][0][0]) + ' ' + str(start_middle_end_coords[i][0][1]) + ' ' + str(start_middle_end_coords[i][0][2]) + '} {' + str(start_middle_end_coords[i][1][0]) + ' ' + str(start_middle_end_coords[i][1][1]) + ' ' + str(start_middle_end_coords[i][1][2]) + '} radius ' + str(cylinder_radius[i]) + ' resolution 85 filled 0\n')
        W.write('draw cone {' + str(start_middle_end_coords[i][1][0]) + ' ' + str(start_middle_end_coords[i][1][1]) + ' ' + str(start_middle_end_coords[i][1][2]) + '} {' + str(start_middle_end_coords[i][2][0]) + ' ' + str(start_middle_end_coords[i][2][1]) + ' ' + str(start_middle_end_coords[i][2][2]) + '} radius ' + str(cone_radius[i]) + ' resolution 85\n')

    W.write('mol delrep 0 top\nmol rename top vectors\n### setting viewpoints\nset viewpoints([molinfo top]) {{{1 0 0 -0.345217} {0 1 0 -0.726605} {0 0 1 -0.474629} {0 0 0 1}} {{0.224004 -0.148068 -0.963274 0} {-0.783484 -0.615181 -0.08763 0} {-0.579614 0.774343 -0.253812 0} {0 0 0 1}} {{0.100633 0 0 0} {0 0.100633 0 0} {0 0 0.100633 0} {0 0 0 1}} {{1 0 0 0} {0 1 0 0} {0 0 1 0} {0 0 0 1}}}\nlappend viewplist [molinfo top]\n\n')
    
    ### prepping the molecule and reps
    W.write('mol new ' + pdb_file + ' type pdb first 0 last -1 step 1 filebonds 1 autobonds 1 waitfor all\n')
    W.write('mol delrep 0 top\n')
    W.write('mol representation NewCartoon 0.160000 50.000000 4.100000 0\n')
    W.write('mol color Name\n')
    W.write('mol selection {all}\n')
    W.write('mol material AOEdgy\n')
    W.write('mol addrep top\n')
    W.write('mol representation Licorice 0.10000 85.000000 85.000000\n')
    W.write('mol color Name\n')
    W.write('mol selection {same residue as ' + vmd_list + '}\n')
    W.write('mol material AOEdgy\n')
    W.write('mol addrep top\n\n')
    W.write('mol representation CPK 1.000000 0.300000 85.000000 85.000000\n')
    W.write('mol color Name\n')
    W.write('mol selection {' + vmd_list + '}\n')
    W.write('mol material AOEdgy\n')
    W.write('mol addrep top\n\n')
    
    W.write('### setting viewpoints\nset viewpoints([molinfo top]) {{{1 0 0 -0.345217} {0 1 0 -0.726605} {0 0 1 -0.474629} {0 0 0 1}} {{0.224004 -0.148068 -0.963274 0} {-0.783484 -0.615181 -0.08763 0} {-0.579614 0.774343 -0.253812 0} {0 0 0 1}} {{0.100633 0 0 0} {0 0.100633 0 0} {0 0 0.100633 0} {0 0 0 1}} {{1 0 0 0} {0 1 0 0} {0 0 1 0} {0 0 0 1}}}\nlappend viewplist [molinfo top]\nset topmol [molinfo top]\n\nforeach v $viewplist { \n  molinfo $v set {center_matrix rotate_matrix scale_matrix global_matrix} $viewpoints($v)\n}\nforeach v $fixedlist {\n  molinfo $v set fixed 1\n}\nunset viewplist\nunset fixedlist\n')

print('Finished creaing vis-state file', vis_state_file_name)

