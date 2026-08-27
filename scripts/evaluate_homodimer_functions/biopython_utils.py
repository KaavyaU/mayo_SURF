"""
biopython_utils.py — BioPython-based structural analysis for homodimer evaluation.

Provides clash scoring, interface residue detection, and sequence extraction.
No hotspot logic — homodimers don't use predefined hotspots.
"""

import numpy as np
from scipy.spatial import cKDTree
from Bio.PDB import PDBParser, Selection
from Bio.PDB.Polypeptide import is_aa

THREE_TO_ONE = {
    'ALA': 'A', 'CYS': 'C', 'ASP': 'D', 'GLU': 'E', 'PHE': 'F',
    'GLY': 'G', 'HIS': 'H', 'ILE': 'I', 'LYS': 'K', 'LEU': 'L',
    'MET': 'M', 'ASN': 'N', 'PRO': 'P', 'GLN': 'Q', 'ARG': 'R',
    'SER': 'S', 'THR': 'T', 'VAL': 'V', 'TRP': 'W', 'TYR': 'Y',
}


def calculate_clash_score(pdb_file, threshold=2.4, only_ca=False):
    """
    Count interchain clashes (atom pairs closer than `threshold` Å).
    Skips intra-residue and bonded-neighbor pairs.
    """
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure('protein', pdb_file)

    atoms = []
    atom_info = []

    for model in structure:
        for chain in model:
            for residue in chain:
                for atom in residue:
                    if atom.element == 'H':
                        continue
                    if only_ca and atom.get_name() != 'CA':
                        continue
                    atoms.append(atom.coord)
                    atom_info.append((chain.id, residue.id[1], atom.get_name()))

    tree = cKDTree(atoms)
    pairs = tree.query_pairs(threshold)

    valid_clashes = set()
    for (i, j) in pairs:
        chain_i, res_i, _ = atom_info[i]
        chain_j, res_j, _ = atom_info[j]

        # Same residue
        if chain_i == chain_j and res_i == res_j:
            continue
        # Bonded neighbors
        if chain_i == chain_j and abs(res_i - res_j) == 1:
            continue
        # For all-atom mode, only count interchain clashes
        if not only_ca and chain_i == chain_j:
            continue

        valid_clashes.add((i, j))

    return len(valid_clashes)


def get_interface_residues(pdb_file, binder_chain='A', target_chain='B', distance_cutoff=4.0):
    """
    Return {resnum: one_letter_aa} for binder residues within `distance_cutoff`
    of any target atom.
    """
    parser = PDBParser(QUIET=True)
    model = parser.get_structure('complex', pdb_file)[0]

    chain_ids = [c.id for c in model]
    if binder_chain not in chain_ids:
        raise ValueError(f"Binder chain '{binder_chain}' not found. Available: {chain_ids}")
    if target_chain not in chain_ids:
        raise ValueError(f"Target chain '{target_chain}' not found. Available: {chain_ids}")

    binder_atoms = list(Selection.unfold_entities(model[binder_chain], 'A'))
    target_atoms = list(Selection.unfold_entities(model[target_chain], 'A'))

    binder_coords = np.array([a.coord for a in binder_atoms])
    target_coords = np.array([a.coord for a in target_atoms])

    binder_tree = cKDTree(binder_coords)
    target_tree = cKDTree(target_coords)

    contact_pairs = binder_tree.query_ball_tree(target_tree, distance_cutoff)

    interface_residues = {}
    for binder_idx, close_target_indices in enumerate(contact_pairs):
        if not close_target_indices:
            continue
        residue = binder_atoms[binder_idx].get_parent()
        resname = residue.get_resname()
        if resname in THREE_TO_ONE:
            interface_residues[residue.id[1]] = THREE_TO_ONE[resname]

    return interface_residues


def get_binder_sequence(pdb_file, binder_chain='A'):
    """Extract the one-letter amino acid sequence of the binder chain."""
    parser = PDBParser(QUIET=True)
    model = parser.get_structure('protein', pdb_file)[0]

    if binder_chain not in model:
        return ''

    seq = ''
    for residue in model[binder_chain]:
        if is_aa(residue, standard=True):
            seq += THREE_TO_ONE.get(residue.get_resname(), 'X')

    return seq
