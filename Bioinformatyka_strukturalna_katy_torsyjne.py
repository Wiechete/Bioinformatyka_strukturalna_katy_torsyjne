from Bio import PDB
import numpy as np
import sys
import os

def calculate_rna_torsion_angles(structure):
    torsion_angles = []
    
    for model in structure:
        for chain in model:
            for residue in chain:
                if residue.id[0] == " " and residue.has_id('P') and residue.has_id('O5\''):
                    atoms = [residue['P'], residue['O5\''], residue['C5\''], residue['C4\''],
                             residue['C3\''], residue['O3\''], residue['O4\''], residue['C1\''],
                             residue['C2\''], residue['C4\''], residue['C3\''], residue['C1\'']]
                    
                    angles = []
                    for i in range(len(atoms) - 3):
                        # Convert numpy arrays to Bio.PDB.Vectors
                        v1 = PDB.Vector(list(atoms[i].get_vector()))
                        v2 = PDB.Vector(list(atoms[i + 1].get_vector()))
                        v3 = PDB.Vector(list(atoms[i + 2].get_vector()))
                        v4 = PDB.Vector(list(atoms[i + 3].get_vector()))

                        angles.append(PDB.vectors.calc_dihedral(v1, v2, v3, v4))
                    
                    torsion_angles.append(angles)
    
    return np.array(torsion_angles)

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

