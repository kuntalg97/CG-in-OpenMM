import mdtraj as md

# Load the trajectory
traj = md.load('output.dcd', top='water512_CG.pdb')

# Visualize trajectory (Example using nglview)
import nglview as nv
view = nv.show_mdtraj(traj)
view
