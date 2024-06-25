import mdtraj as md
import matplotlib.pyplot as plt

# Load the trajectory
traj = md.load('output.dcd', top='water512_CG.pdb')

# Define atom pairs for which to compute the RDF (e.g., all pairs of oxygen atoms in water)
oxygen_atoms = traj.topology.select('name O')

# Compute RDF
rdf, r = md.compute_rdf(traj, pairs=traj.topology.select_pairs(oxygen_atoms, oxygen_atoms), r_range=(0.0, 1.5), bin_width=0.01)

# Plot RDF
plt.plot(r, rdf)
plt.xlabel('Distance (nm)')
plt.ylabel('g(r)')
plt.title('Radial Distribution Function')
plt.show()
