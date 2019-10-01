
# ----------------------------------------
# USAGE:
# ----------------------------------------
# python3 visualize_P_edge.py pdb_file [COM][ATOMIC] edge_metric_file output_vis_state_file_name output_vis_state_structural_file_name average_structure_file distance_cutoff value_cutoff

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
node_definition = str(sys.argv[2])
edge_metric_file = sys.argv[3]          # used to color and size edge cylinders
vis_state_file_name = sys.argv[4]
vis_structure_file_name = sys.argv[5]
average_structure = np.loadtxt(sys.argv[6])
distance_cutoff = float(sys.argv[7])
distance2_cutoff = distance_cutoff**2
value_threshold = float(sys.argv[8])       # think about how to implement this...

# ----------------------------------------
# SETTING COLORS
# ----------------------------------------
color_defs = {}
# colorid 33 is to be used to color nodes
color_defs[32] = "color change rgb 32 0.0 0.0 0.0\n"        # black
####    edges to be colored from lightgray to blue
edge_high_centrality_rgb = np.array([0.0,0.0,1.0])         # blue
edge_low_centrality_rgb = np.array([0.827,0.827,0.827])    # lightgray

# ----------------------------------------
# CREATING/FILLING COLOR DICTIONARIES
# ----------------------------------------
edge_possible_colorids = list(range(33,1057))
nEdge_colorids = len(edge_possible_colorids)
cmap_positions = np.linspace(0,1,nEdge_colorids)

cdict = {'red':[], 'green':[], 'blue':[]}
for colorid in edge_possible_colorids:
    index = colorid - edge_possible_colorids[0]
    thiscolor = np.abs(edge_low_centrality_rgb + index * (edge_high_centrality_rgb - edge_low_centrality_rgb)/(nEdge_colorids-1))
    ### VMD Stuff
    color_defs[colorid] = "color change rgb " + str(colorid) + " " + str(thiscolor[0]) + " " + str(thiscolor[1]) + " " + str(thiscolor[2]) + '\n'
    ### MATPLOTLIB Stuff
    cdict['red'].append((cmap_positions[index],thiscolor[0],thiscolor[0]))
    cdict['green'].append((cmap_positions[index],thiscolor[1],thiscolor[1]))
    cdict['blue'].append((cmap_positions[index],thiscolor[2],thiscolor[2]))

# ----------------------------------------
# EDGE RADIUSES
# ----------------------------------------
edge_large_radius = 0.5
edge_small_radius = 0.05
edge_radius_difference = edge_large_radius - edge_small_radius

# ----------------------------------------
# CREATING MDANALYSIS SELECTIONS, COLLECTING COM POSITIONS
# ----------------------------------------
u = MDAnalysis.Universe(pdb_file)
u_all = u.select_atoms('all')

selection_list = []
com_list = []
count = 0
with open('node_selection.txt','w') as f:
    # ----------------------------------------
    # SUBSTRATE SELECTION - CENTER OF MASS OF RESIDUES
    # ----------------------------------------
    if node_definition.upper() == 'COM':
        substrate_selection = u.select_atoms('protein')
        u_all.translate(-substrate_selection.center_of_mass())
        nResidues_range = list(range(substrate_selection.n_residues))
        for i in nResidues_range:
            temp = substrate_selection.residues[i].atoms
            selection_list.append(temp)
            temp_resname = temp.resnames[0]
            temp_resid = temp.resids[0]
            f.write("%02d   %s   %s\n" %(count,temp_resname,temp_resid))
            com_list.append(temp.center_of_mass())
            count += 1

    # ----------------------------------------
    # SUBSTRATE SELECTION - ATOMS
    # ----------------------------------------
    elif node_definition.upper() == 'ATOMIC':
        substrate_selection = u.select_atoms('protein')
        nAtoms_range = list(range(substrate_selection.n_atoms))
        for i in nAtoms_range:
            temp = substrate_selection.atoms[i]
            selection_list.append(temp)
            temp_resname = temp.resname
            temp_resid = temp.resid
            f.write("%02d   %s   %s\n" %(count,temp_resname,temp_resid))
            com_list.append(temp.pos)
            count += 1

nNodes_range = list(range(count))

# ----------------------------------------
# PAIR-WISE DISTANCE ANALYSIS; APPLYING DISTANCE CUTOFF
# ----------------------------------------
# THIS CAN BE SOME MUCH MORE PYTHONIC... NUMPY MASKED ARRAY
distance_cutoff_matrix = np.zeros((count,count),dtype=np.int)
for i in nNodes_range[:-1]:
    for j in nNodes_range[i+1:]:
        dist2 = np.sum((average_structure[i] - average_structure[j])**2)
        if dist2 < distance2_cutoff:
            distance_cutoff_matrix[i,j] = 1

# ----------------------------------------
# LOADING IN AND ANALYZING THE EDGE METRIC FILE
# ----------------------------------------
edge_matrix = np.zeros((count,count),dtype=np.float64)
with open(edge_metric_file,'r') as W:
    for line in W:
        temp = line.split()
        node1 = int(temp[0]) - 1 
        node2 = int(temp[1]) - 1 
        if temp[2] != 0. and distance_cutoff_matrix[node1,node2] == 1:
            edge_matrix[node1,node2] = edge_matrix[node2,node1] = float(temp[2])

edge_metric_min_max = (0.0,1.0)
#edge_metric_min_max = (np.min(edge_matrix),np.max(edge_matrix))
range_ = edge_metric_min_max[1] - edge_metric_min_max[0]
delta = (range_)/(nEdge_colorids-1)
print(edge_metric_min_max, range_, delta)

# prepping the edge_radius_colorid_data array
edge_radius_colorid_data = []
for i in nNodes_range[:-1]:
    for j in nNodes_range[i+1:]:
        if edge_matrix[i,j] > value_threshold:     #applying a second cutoff, specifically on the centrality value
            edge_radius_colorid_data.append([i,j,edge_matrix[i,j],0,0]) # org: node1,node2,edge_metric,radius,color

edge_radius_colorid_data = np.array(edge_radius_colorid_data)
for i in edge_radius_colorid_data:
    # radius
    i[-2] = ((i[2] - edge_metric_min_max[0])/range_)*edge_radius_difference + edge_small_radius
    # color
    i[-1] = edge_possible_colorids[0] + int((i[2] - edge_metric_min_max[0])/delta)

print(len(edge_radius_colorid_data))
print(np.min(edge_radius_colorid_data[:,-2]), np.max(edge_radius_colorid_data[:,-2]), np.min(edge_radius_colorid_data[:,-1]), np.max(edge_radius_colorid_data[:,-1]))

# sorting the edge_radius_colorid_data array by colorid to prevent repetitive graphics calls in VMD
unsorted_colorids = edge_radius_colorid_data[:,-1]
idx = unsorted_colorids.argsort()
edge_radius_colorid_data = edge_radius_colorid_data[idx]

# ----------------------------------------
# CREATING VIS STATE FILE
# ----------------------------------------
with open(vis_state_file_name,'w') as W:
    ### starting lines
    W.write('#!/usr/local/bin/vmd\nset viewplist {}\nset fixedlist {}\n\n')

    ### setting colorids
    W.write('# setting colorid rgb values\n')
    W.write(color_defs[32])
    for i in edge_possible_colorids:
        W.write(color_defs[i])

    ### drawing node spheres
    W.write('\n# Draw spheres at the nodes\nmol new\n')
    W.write('draw material AOEdgy\n')
    W.write('graphics top color 32 \n')
    for i in nNodes_range:
        W.write('draw sphere {' + str(com_list[i][0]) + ' ' + str(com_list[i][1]) + ' ' + str(com_list[i][2]) + '} resolution 85 radius 0.25\n')

    W.write('mol delrep 0 top\nmol rename top node_spheres\n### setting viewpoints\nset viewpoints([molinfo top]) {{{1 0 0 -0.0931692} {0 1 0 0.209481} {0 0 1 0.0970525} {0 0 0 1}} {{-0.937397 -0.311889 -0.154964 0} {-0.319652 0.947151 0.0273406 0} {0.13825 0.0751647 -0.98754 0} {0 0 0 1}} {{0.0429779 0 0 0} {0 0.0429779 0 0} {0 0 0.0429779 0} {0 0 0 1}} {{1 0 0 0.03} {0 1 0 -0.0100001} {0 0 1 0} {0 0 0 1}}}\nlappend viewplist [molinfo top]\n\n')

    W.write('\n# Draw cylinders between the nodes\nmol new\n')
    ### drawing edge cylinders
    colorid = 99999
    for i in edge_radius_colorid_data:
        if i[-1] != colorid:
            W.write('graphics top color ' + str(int(i[-1])) + '\n')
            colorid = i[-1]
        W.write('draw cylinder {' + str(com_list[int(i[0])][0]) + ' ' + str(com_list[int(i[0])][1]) + ' ' + str(com_list[int(i[0])][2]) + '} {' + str(com_list[int(i[1])][0]) + ' ' + str(com_list[int(i[1])][1]) + ' ' + str(com_list[int(i[1])][2]) + '} radius ' + str(i[-2]) + ' resolution 85 filled 0\n')

    W.write('mol delrep 0 top\nmol rename top edges\n### setting viewpoints\nset viewpoints([molinfo top]) {{{1 0 0 -0.0931692} {0 1 0 0.209481} {0 0 1 0.0970525} {0 0 0 1}} {{-0.937397 -0.311889 -0.154964 0} {-0.319652 0.947151 0.0273406 0} {0.13825 0.0751647 -0.98754 0} {0 0 0 1}} {{0.0429779 0 0 0} {0 0.0429779 0 0} {0 0 0.0429779 0} {0 0 0 1}} {{1 0 0 0.03} {0 1 0 -0.0100001} {0 0 1 0} {0 0 0 1}}}\nlappend viewplist [molinfo top]\n\n')

    ### prepping the molecule and reps
    W.write('mol new ' + vis_structure_file_name + ' type pdb first 0 last -1 step 1 filebonds 1 autobonds 1 waitfor all\n')
    W.write('mol delrep 0 top\n')
    W.write('mol representation NewCartoon 0.160000 50.000000 4.100000 0\n')
    W.write('mol color Name\n')
    W.write('mol selection {all}\n')
    W.write('mol material AOEdgy\n')
    W.write('mol addrep top\n')

    W.write('### setting viewpoints\nset viewpoints([molinfo top]) {{{1 0 0 -0.0931692} {0 1 0 0.209481} {0 0 1 0.0970525} {0 0 0 1}} {{-0.937397 -0.311889 -0.154964 0} {-0.319652 0.947151 0.0273406 0} {0.13825 0.0751647 -0.98754 0} {0 0 0 1}} {{0.0429779 0 0 0} {0 0.0429779 0 0} {0 0 0.0429779 0} {0 0 0 1}} {{1 0 0 0.03} {0 1 0 -0.0100001} {0 0 1 0} {0 0 0 1}}}\nlappend viewplist [molinfo top]\nset topmol [molinfo top]\n\nforeach v $viewplist { \n  molinfo $v set {center_matrix rotate_matrix scale_matrix global_matrix} $viewpoints($v)\n}\nforeach v $fixedlist {\n  molinfo $v set fixed 1\n}\nunset viewplist\nunset fixedlist\n')

print('Finished creaing vis-state file', vis_state_file_name)

u_all.write(vis_structure_file_name)

# ----------------------------------------
# CREATING COLORBAR IMAGE 
# ----------------------------------------
cmap = mpl.colors.LinearSegmentedColormap('my_cmap',cdict,nEdge_colorids)
fig, ax = plt.subplots(figsize=(2,8))
fig.subplots_adjust(right=0.5)
norm = mpl.colors.Normalize(vmin=edge_metric_min_max[0],vmax=edge_metric_min_max[1])
cb = mpl.colorbar.ColorbarBase(ax,cmap=cmap,spacing='uniform',orientation='vertical',norm=norm,ticks=[0,0.25,0.50,0.75,1])
cb.set_label(r'P$_{Edge}$',size=16)
cb.set_ticklabels([sys.argv[8],'0.25','0.50','0.75','1.0'])
plt.savefig('color_bar.%s_%s.png'%(sys.argv[7],sys.argv[8]),dpi=600,transparent=True)
plt.close()

