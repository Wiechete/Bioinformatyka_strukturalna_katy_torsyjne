from Bio import PDB
import numpy as np
import sys
import os

# funkcja obliczajaca katy torsyjne
#def calculate_rna_torsion_angles(structure):
#    torsion_angles = []
#    # petla sprawdzajaca czy residue jest nukleotydem na podstawie fosforu i atomow 05
#    for model in structure:
#        for chain in model:
#            for residue in chain:
#                if residue.id[0] == " " and residue.has_id('P') and residue.has_id('O5\''):
#                    atoms = [residue['P'], residue['O5\''], residue['C5\''], residue['C4\''],
#                             residue['C3\''], residue['O3\''], residue['O4\''], residue['C1\''],
#                             residue['C2\''], residue['C4\''], residue['C3\''], residue['C1\'']]
                    
#                    angles = []
#                    for i in range(len(atoms) - 3):# obliczanie katow
#                        # Zamien numpy arrays na Bio.PDB.Vectors
#                        v1 = PDB.Vector(list(atoms[i].get_vector()))
#                        v2 = PDB.Vector(list(atoms[i + 1].get_vector()))
#                        v3 = PDB.Vector(list(atoms[i + 2].get_vector()))
#                        v4 = PDB.Vector(list(atoms[i + 3].get_vector()))

#                        angles.append(PDB.vectors.calc_dihedral(v1, v2, v3, v4))# oblicz katy miedzy vectorami
                    
#                    torsion_angles.append(angles)
    
#    return np.array(torsion_angles)

def calculate_rna_torsion_angles(structure):
    torsion_angles = []
    # petla sprawdzajaca, czy residue jest nukleotydem
    for model in structure:
        for chain in model:
            for residue in chain:
                if residue.id[0] == " ":
                    # Atomowa definicja katow
                    atoms_alpha = [residue['O3\''], residue['P'], residue['O5\''], residue['C5\'']]
                    atoms_beta = [residue['P'], residue['O5\''], residue['C5\''], residue['C4\'']]
                    atoms_gamma = [residue['O5\''], residue['C5\''], residue['C4\''], residue['C3\'']]
                    atoms_delta = [residue['C5\''], residue['C4\''], residue['C3\''], residue['O3\'']]
                    atoms_epsilon = [residue['C4\''], residue['C3\''], residue['O3\''], residue['P']]
                    atoms_zeta = [residue['C3\''], residue['O3\''], residue['P'], residue['O5\'']]
                    
                    # Sprawdzanie dostepnosci atomow chi
                    try:
                        atoms_chi = [residue['N1'], residue['C2'], residue['N3'], residue['C4'],
                                     residue['C5'], residue['C6'], residue['O2'], residue['O4']]
                    except KeyError:
                        continue  # Pomijamy nukleotydy, ktore nie maja wszystkich atomow potrzebnych do obliczenia chi

                    # Obliczanie katow
                    angles_alpha = calculate_dihedral_angle(atoms_alpha)
                    angles_beta = calculate_dihedral_angle(atoms_beta)
                    angles_gamma = calculate_dihedral_angle(atoms_gamma)
                    angles_delta = calculate_dihedral_angle(atoms_delta)
                    angles_epsilon = calculate_dihedral_angle(atoms_epsilon)
                    angles_zeta = calculate_dihedral_angle(atoms_zeta)
                    angles_chi = calculate_dihedral_angle(atoms_chi)

                    torsion_angles.append([angles_alpha, angles_beta, angles_gamma, angles_delta, angles_epsilon, angles_zeta, angles_chi])

    return np.array(torsion_angles)

# obliczanie katow
def calculate_dihedral_angle(atoms):
    vectors = [PDB.Vector(list(atom.get_vector())) for atom in atoms]
    return PDB.vectors.calc_dihedral(vectors[0], vectors[1], vectors[2], vectors[3])

# zapisywanie katow torsyjnych
def save_torsion_angles_to_file(angles, filename):
    np.savetxt(filename, angles, fmt='%f', delimiter='\t')

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python script.py <pdb_file>")
        sys.exit(1)

    pdb_file = sys.argv[1]
    filename = os.path.basename(pdb_file)
    pdb_id = filename[0:4]

structure = PDB.PDBParser(QUIET=True).get_structure(pdb_id, filename)
torsion_angles = calculate_rna_torsion_angles(structure)
save_torsion_angles_to_file(torsion_angles, f"{pdb_id}_torsion_angles.txt")

