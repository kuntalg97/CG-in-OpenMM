def read_lammps_tabulated_pair_potential(file_path):
    distances = []
    energies = []

    with open(file_path, 'r') as file:
        for _ in range(5):
            next(file)

        for line in file:
            # Split the line into columns
            columns = line.split()

            distance = float(columns[1])
            energy = float(columns[2])

            distances.append(distance)
            energies.append(energy)

    return distances, energies

# Example usage
file_path = 'bw_sp2b.table'
distances, energies = read_lammps_tabulated_pair_potential(file_path)

# Print the results or perform further processing
#for i in range(len(distances)):
#    print(f"Distance: {distances[i]}, Energy: {energies[i]}")

n_lines = len(open(file_path).readlines())
n_lines -= 5

for i in range(1,n_lines+1):
    print("%16.6f%16.6f"%(distances[i-1],energies[i-1]))
