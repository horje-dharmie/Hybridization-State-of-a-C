#!/home/damiojedeji/anaconda3/bin/python3

import numpy as np
import mdtraj as md

# Load the trajectory file
# Make sure to replace 'path_to_your_file' with the actual file path
traj = md.load('a-C.gro')

# Select carbon atoms
carbon_indices = traj.top.select('element C')

# Compute pairwise distances between carbon atoms
atom_pairs = np.array([(i, j) for i in carbon_indices for j in carbon_indices if i < j])
distances = md.compute_distances(traj, atom_pairs)

# Identify covalent bonds using a cutoff of 1.9 Å
cutoff = 0.19
bond_indices = np.where(distances < cutoff)

# Create a dictionary to store bonds for each atom
atom_bonds = {i: [] for i in carbon_indices}
for idx in range(len(atom_pairs)):
    if distances[0, idx] < cutoff:
        i, j = atom_pairs[idx]
        atom_bonds[i].append(j)
        atom_bonds[j].append(i)

# Calculate coordination number for each atom
coordination_numbers = {atom: len(bonds) for atom, bonds in atom_bonds.items()}

# Determine hybridization state based on coordination number
hybridization_states = {}
for atom, coord_number in coordination_numbers.items():
    if coord_number == 4:
        hybridization_states[atom] = 'sp3'
    elif coord_number == 3:
        hybridization_states[atom] = 'sp2'
    elif coord_number == 2:
        hybridization_states[atom] = 'sp'
    else:
        hybridization_states[atom] = 'unknown'  # Handle unexpected cases

# Count the number of sp2 and sp3 atoms
sp2_count = sum(1 for state in hybridization_states.values() if state == 'sp2')
sp3_count = sum(1 for state in hybridization_states.values() if state == 'sp3')

# Calculate the ratio of sp2 to sp3
if sp3_count > 0:
    ratio_sp2_to_sp3 = sp2_count / sp3_count
else:
    ratio_sp2_to_sp3 = float('inf')  # Handle case where there are no sp3 atoms

# Print the results
print(f"Number of sp2 atoms: {sp2_count}")
print(f"Number of sp3 atoms: {sp3_count}")
print(f"Ratio of sp2 to sp3: {ratio_sp2_to_sp3}")

# Assign beta factors based on hybridization states
# Pronounced differences in beta factors
hybridization_to_beta = {'sp3': 10, 'sp2': 50, 'sp': 90, 'unknown': 0}
bfactors = np.zeros(traj.n_atoms)
for atom in carbon_indices:
    bfactors[atom] = hybridization_to_beta[hybridization_states[atom]]

# Save the PDB file with beta factors
traj.save_pdb('amorphous_carbon_hybridization.pdb', bfactors=bfactors)

