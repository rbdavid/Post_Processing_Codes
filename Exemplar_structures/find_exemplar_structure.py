
# ----------------------------------------
# USAGE:
# ----------------------------------------
# python3 find_exemplar_structure.py average_structure trajectory_location start end 

# ----------------------------------------
# PREAMBLE:
# ----------------------------------------

import sys
import MDAnalysis
from MDAnalysis.analysis.align import rotation_matrix
import numpy as np

average_structure = np.loadtxt(sys.argv[1])
pdb = sys.argv[2]
traj_loc = sys.argv[3]
start = int(sys.argv[4])
end = int(sys.argv[5])

trajectory_list = [traj_loc%(i) for i in list(range(start,end+1))]

atom_selection = '(resid 55 56 59 66 67 72:77 79:81 83 84 86 and name CA) or (resid 52 68 69 and name CA CZ) or (resid 62 70 78 85 and name CA CG) or (resid 71 and name CA OG1) or (resid 82 and name CA CD) or (resid 87 and name CA CE) or (resid 88 and name CA SG)'

# ----------------------------------------
# ANALYZE TRAJECTORIES:
# ----------------------------------------

u = MDAnalysis.Universe(pdb)
u_all = u.select_atoms('all')
u_selection = u.select_atoms(atom_selection)
nAtoms = u_selection.n_atoms
if nAtoms != average_structure.shape[0]:
    print('Average structure and atom selections do not have the same length.')
    sys.exit()

minRMSD = 9999.
traj_step = ['',-1] # traj_string and step number 

for traj in trajectory_list:
    u.load_new(traj)
    print('Analyzing ',traj)
    for ts in u.trajectory:
        u_all.translate(-u_selection.center_of_mass())
        R,d = rotation_matrix(u_selection.positions,average_structure)
        if d < minRMSD:
            minRMSD = d
            traj_step = [traj,ts.frame]

print(minRMSD,traj_step[0],traj_step[1])

u.load_new(traj_step[0])
u.trajectory[traj_step[1]]
u_all.translate(-u_selection.center_of_mass())
R,d = rotation_matrix(u_selection.positions,average_structure)
u_all.rotate(R)
u_all.write('exemplar_structure.pdb')

