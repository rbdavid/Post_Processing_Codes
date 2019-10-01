
# ----------------------------------------
# USAGE:
# ----------------------------------------
# python3 visualize_node_dC_dk.py pdb_file [COM][ATOMIC] node_metric_file output_vis_state_file_name output_vis_state_structural_file_name output_color_bar_figure_name

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
node_definition = str(sys.argv[2])      # currently only allows for single-residue coarse graining...
node_metric_file = sys.argv[3]          # used to color node spheres
vis_state_file_name = sys.argv[4]
vis_structure_file_name = sys.argv[5]
color_bar_figure_name = sys.argv[6]

ss_vmd_selection = 'resid 50 104 130 225 304 337 431 433'

# ----------------------------------------
# SETTING COLORS
# ----------------------------------------
color_defs = {}
node_high_centrality_rgb = np.array([0.25098039215686274, 0.0, 0.29411764705882354])                # #40004b; purple
node_mid_centrality_rgb  = np.array([0.8862745098039215, 0.8862745098039215, 0.8862745098039215])   # #e2e2e2; gray
node_low_centrality_rgb  = np.array([0.0, 0.26666666666666666, 0.10588235294117647])                # #00441b; green

# ----------------------------------------
# CREATING/FILLING COLOR DICTIONARIES
# ----------------------------------------
node_possible_colorids = list(range(33,1057))
nNode_colorids = len(node_possible_colorids)
low_to_mid_colorids = list(range(33,545))
mid_to_high_colorids = list(range(545,1057))
cmap_positions = np.linspace(0,1,nNode_colorids)
print(cmap_positions, cmap_positions[0], cmap_positions[-1])

cdict = {'red':[], 'green':[], 'blue':[]}
for colorid in low_to_mid_colorids:
    index = colorid - low_to_mid_colorids[0]
    thiscolor = np.abs(node_low_centrality_rgb + index * (node_mid_centrality_rgb - node_low_centrality_rgb)/(len(low_to_mid_colorids)-1))
    ### VMD Stuff
    color_defs[colorid] = "color change rgb " + str(colorid) + " " + str(thiscolor[0]) + " " + str(thiscolor[1]) + " " + str(thiscolor[2]) + '\n'
    ### MATPLOTLIB Stuff
    cdict['red'].append((cmap_positions[index],thiscolor[0],thiscolor[0]))
    cdict['green'].append((cmap_positions[index],thiscolor[1],thiscolor[1]))
    cdict['blue'].append((cmap_positions[index],thiscolor[2],thiscolor[2]))

for colorid in mid_to_high_colorids:
    index = colorid - mid_to_high_colorids[0]
    thiscolor = np.abs(node_mid_centrality_rgb + index * (node_high_centrality_rgb - node_mid_centrality_rgb)/(len(mid_to_high_colorids)-1))
    ### VMD Stuff
    color_defs[colorid] = "color change rgb " + str(colorid) + " " + str(thiscolor[0]) + " " + str(thiscolor[1]) + " " + str(thiscolor[2]) + '\n'
    ### MATPLOTLIB Stuff
    cdict['red'].append((cmap_positions[colorid-node_possible_colorids[0]],thiscolor[0],thiscolor[0]))
    cdict['green'].append((cmap_positions[colorid-node_possible_colorids[0]],thiscolor[1],thiscolor[1]))
    cdict['blue'].append((cmap_positions[colorid-node_possible_colorids[0]],thiscolor[2],thiscolor[2]))

# ----------------------------------------
# NODE RADIUSES
# ----------------------------------------
node_large_radius = 1.5
node_small_radius = 0.1
node_radius_difference = node_large_radius - node_small_radius

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
# LOADING IN AND ANALYZING THE NODE METRIC FILE
# ----------------------------------------
node_metric_data = np.loadtxt(node_metric_file)[:,1]  # data organization: node id (1 indexed), dC/dk values; only grabbing 2nd column

node_metric_min_max = (np.min(node_metric_data),np.max(node_metric_data))
node_metric_abs_max = np.max(np.abs(node_metric_min_max))
print node_metric_min_max, node_metric_abs_max
node_metric_min_max = (-node_metric_abs_max,node_metric_abs_max)
range_ = node_metric_min_max[1] - node_metric_min_max[0]
delta = (range_)/(nNode_colorids-1)
print node_metric_min_max, delta, range_

scale_factor = 99.99/(node_metric_abs_max)
print 'scale factor = ', scale_factor
node_vmd_selection = 'resid'
node_colorids = np.zeros(count,dtype=np.int)
node_radius = np.zeros(count,dtype=np.float)
for i in nNodes_range:
    selection_list[i].tempfactors = node_metric_data[i]*scale_factor
    node_colorids[i] = node_possible_colorids[0] + int((node_metric_data[i]-node_metric_min_max[0])/delta)
    node_radius[i] = node_small_radius + (2*np.abs(node_metric_data[i])/range_)*node_radius_difference
    if np.abs(node_metric_data[i]) > 0.1*node_metric_abs_max:  #XXX... just some arb. cutoff; 10% of max value
        node_vmd_selection += ' ' + str(i+1)

print(np.min(node_colorids), np.max(node_colorids), np.min(node_radius), np.max(node_radius))

####
with open(vis_state_file_name,'w') as W:
    ### starting lines
    W.write('#!/usr/local/bin/vmd\nset viewplist {}\nset fixedlist {}\n\n')

    ### setting colorids
    W.write('# setting colorid rgb values\n')
    for i in node_possible_colorids:
        W.write(color_defs[i])

    ### drawing node spheres
    W.write('\n# Draw spheres at the nodes\nmol new\n')
    W.write('draw material AOEdgy\n')
    for i in nNodes_range:
        if np.abs(node_metric_data[i]) > 0.25*node_metric_abs_max:  #XXX... just some arb. cutoff; 25% of max value
            W.write('# node ' + str(i) + 'P_node value ' + str(node_metric_data[i]) + '\n')
            W.write('graphics top color ' + str(node_colorids[i]) + '\n')
            W.write('draw sphere {' + str(com_list[i][0]) + ' ' + str(com_list[i][1]) + ' ' + str(com_list[i][2]) + '} resolution 85 radius ' + str(node_radius[i]) + ' \n')

    W.write('mol delrep 0 top\nmol rename top graphics\n### setting viewpoints\nset viewpoints([molinfo top]) {{{1 0 0 -0.0931692} {0 1 0 0.209481} {0 0 1 0.0970525} {0 0 0 1}} {{-0.937397 -0.311889 -0.154964 0} {-0.319652 0.947151 0.0273406 0} {0.13825 0.0751647 -0.98754 0} {0 0 0 1}} {{0.0429779 0 0 0} {0 0.0429779 0 0} {0 0 0.0429779 0} {0 0 0 1}} {{1 0 0 0.03} {0 1 0 -0.0100001} {0 0 1 0} {0 0 0 1}}}\nlappend viewplist [molinfo top]\n\n')

    ### prepping the molecule and reps
    W.write('mol new ' + vis_structure_file_name + ' type pdb first 0 last -1 step 1 filebonds 1 autobonds 1 waitfor all\n')
    W.write('mol delrep 0 top\n')
    W.write('mol representation NewCartoon 0.160000 50.000000 4.100000 0\n')
    W.write('mol color Beta\n')
    W.write('mol selection {all}\n')
    W.write('mol material AOEdgy\n')
    W.write('mol addrep top\n')
    W.write('mol scaleminmax top 0 -99.99 99.99\n')
    W.write('mol representation Licorice 0.150000 85.000000 85.000000\n')
    W.write('mol color Name\n')
    W.write('mol selection {' + node_vmd_selection + '}\n')
    W.write('mol material AOEdgy\n')
    W.write('mol addrep top\n\n')
    W.write('mol representation Licorice 0.150000 85.000000 85.000000\n')
    W.write('mol color Name\n')
    W.write('mol selection {' + ss_vmd_selection + '}\n')
    W.write('mol material AOEdgy\n')
    W.write('mol addrep top\n\n')
    
    W.write('### setting viewpoints\nset viewpoints([molinfo top]) {{{1 0 0 -0.0931692} {0 1 0 0.209481} {0 0 1 0.0970525} {0 0 0 1}} {{-0.937397 -0.311889 -0.154964 0} {-0.319652 0.947151 0.0273406 0} {0.13825 0.0751647 -0.98754 0} {0 0 0 1}} {{0.0429779 0 0 0} {0 0.0429779 0 0} {0 0 0.0429779 0} {0 0 0 1}} {{1 0 0 0.03} {0 1 0 -0.0100001} {0 0 1 0} {0 0 0 1}}}\nlappend viewplist [molinfo top]\nset topmol [molinfo top]\n\nforeach v $viewplist { \n  molinfo $v set {center_matrix rotate_matrix scale_matrix global_matrix} $viewpoints($v)\n}\nforeach v $fixedlist {\n  molinfo $v set fixed 1\n}\nunset viewplist\nunset fixedlist\n')

print('Finished creating vis-state file', vis_state_file_name)

u_all.write(vis_structure_file_name)

### MATPLOTLIB Stuff
cmap = mpl.colors.LinearSegmentedColormap('my_cmap',cdict,nNode_colorids)
fig, ax = plt.subplots(figsize=(3,8))
fig.subplots_adjust(right=0.5)
norm = mpl.colors.Normalize(vmin=node_metric_min_max[0]*10**3,vmax=node_metric_min_max[1]*10**3)
cb = mpl.colorbar.ColorbarBase(ax,cmap=cmap,spacing='uniform',orientation='vertical',norm=norm) #,ticks=[0,lower_bound,0.25,0.50,0.75,1]
cb.set_label(r'$\sum_{j} k_{i,j} \frac{d||C_{m,n}||^2}{dk_{i,j}}$ (10$^{-3}$ $\AA^{4}$)',size=16)
#cb.set_ticklabels([lower_bound_str,'0.25','0.50','0.75','1.0'])
plt.savefig(color_bar_figure_name,dpi=600,transparent=True)
plt.close()

