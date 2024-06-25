from openmm.app import *
from openmm import *
from openmm.unit import *
from sys import stdout

pdb = PDBFile('water512_CG.pdb')
table_file = 'bw_sp2b.table'

coordinates_in_angstroms = pdb.positions
coordinates_in_nanometers = coordinates_in_angstroms / 10.0

energies = []  # energy values
distances = []  # distance values

with open(table_file, 'r') as file:
    for _ in range(5):
        next(file)

    for line in file:
        columns = line.split()

        distance = float(columns[1])
        energy = float(columns[2])
        distance = distance*0.1   # converts to nm
        energy = energy*4.184     # converts to kJ/mol

        distances.append(distance)
        energies.append(energy)

# Create an OpenMM System
system = openmm.System()

forcefield = ForceField('field.xml')
system = forcefield.createSystem(pdb.topology)

box_size = 2.5022 * nanometers
system.setDefaultPeriodicBoxVectors([box_size, 0, 0], [0, box_size, 0], [0, 0, box_size])

# Add particles to the system 
#for atom in pdb.topology.atoms():
#    mass = 18.015400
#    system.addParticle(mass)

N = system.getNumParticles()
print ("Number of particles =", N)

# Create a CustomCompoundBondForce with a tabulated potential
nonbond_force = openmm.CustomNonbondedForce('energy(r)')

for atom in pdb.topology.atoms():
    nonbond_force.addParticle()

# Define the tabulated function
tabulated_function = openmm.Continuous1DFunction(energies, 0.0, 15, False)
nonbond_force.addTabulatedFunction('energy', tabulated_function)

# Add the bond term to the force
#for i in range(0,N-2):
#    for j in range(i+1,N-1): 
#        bond_force.addBond([i, j], [])  # assuming bond between particles 0 and 1

# Add the CustomCompoundBondForce to the System
system.addForce(nonbond_force)

integrator = LangevinMiddleIntegrator(300*kelvin, 1/picosecond, 0.005*picoseconds)
#integrator = NoseHooverIntegrator(300*kelvin, 1/picosecond, 0.005*picoseconds, 4, 1, 1) 
simulation = Simulation(pdb.topology, system, integrator)
simulation.context.setPositions(coordinates_in_nanometers)
simulation.context.setPeriodicBoxVectors([box_size, 0, 0], [0, box_size, 0], [0, 0, box_size])
simulation.minimizeEnergy()
simulation.reporters.append(PDBReporter('output.pdb', 1000))
simulation.reporters.append(StateDataReporter(stdout, 1000, step=True, potentialEnergy=True, temperature=True))
#simulation.reporters.append(ForceReporter('forces.txt', 1000))
simulation.step(2500000)  
