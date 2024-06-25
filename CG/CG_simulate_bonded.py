from openmm.app import *
from openmm import *
from openmm.unit import *
from sys import stdout

pdb = PDBFile('water512_CG.pdb')
table_file = 'bw_sp2b.table'

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

# Add particles to the system (assuming you have two particles)
for atom in pdb.topology.atoms():
    mass = 18.015400
    system.addParticle(mass)

N = system.getNumParticles()
print ("Number of particles =", N)

# Create a CustomCompoundBondForce with a tabulated potential
bond_force = openmm.CustomCompoundBondForce(2, 'energy(distance(p1,p2))')

# Define the tabulated function
tabulated_function = openmm.Continuous1DFunction(energies, 0.0, 10.0, False)
bond_force.addTabulatedFunction('energy', tabulated_function)

# Add the bond term to the force
for i in range(0,N-2):
    for j in range(i+1,N-1): 
        bond_force.addBond([i, j], [])  # assuming bond between particles 0 and 1

# Add the CustomCompoundBondForce to the System
system.addForce(bond_force)

integrator = LangevinMiddleIntegrator(300*kelvin, 1/picosecond, 0.005*picoseconds)
#integrator = NoseHooverIntegrator(300*kelvin, 1/picosecond, 0.005*picoseconds, 4, 1, 1) 
simulation = Simulation(pdb.topology, system, integrator)
simulation.context.setPositions(pdb.positions)
#simulation.minimizeEnergy()
simulation.reporters.append(PDBReporter('output.pdb', 1000))
simulation.reporters.append(StateDataReporter(stdout, 1000, step=True, potentialEnergy=True, temperature=True))
#simulation.reporters.append(ForceReporter('forces.txt', 1000))
simulation.step(2500000)  
