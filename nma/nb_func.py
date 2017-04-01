import numpy as np

__all__ = ["hydrophobic"]

# Connect all pairs of non-bonded atoms from hydrophobic residues with springs 
def hydrophobic(structure, non_bonded):
  output = []
  atom_resid = [str(atom.residue)[:3] for atom in structure.top.atoms]
  hydrophobic = ['GLY', 'ALA', 'VAL', 'LEU', 'ILE', 'PRO', 'PHE', 'MET', 'TRP']
  for pair in non_bonded:
    if np.all(np.array([atom_resid[pair[0]] in hydrophobic, atom_resid[pair[1]] in hydrophobic])):
      output.append(pair)
  return set(output)
