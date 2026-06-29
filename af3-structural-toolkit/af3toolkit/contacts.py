"""
Structural contact detection from an AF3-output PDB file.

Two kinds of contacts:
  - find_contacts: residue-residue contacts between two chains (e.g. the
    receptor fragment and the Reelin fragment), within a distance cutoff.
  - find_chain_ion_contacts: residues of a chain coordinating an ion
    (e.g. calcium, relevant for LDL-receptor-family LA repeats).

Both return LOCAL PDB residue numbering -- use sequence_mapping.build_position_maps
separately to translate these into full-sequence positions and exon/repeat IDs.
"""

import pandas as pd
from Bio.PDB import PDBParser, NeighborSearch
from Bio.PDB.Polypeptide import is_aa


def find_contacts(pdb_file, chain1_id="A", chain2_id="B", distance_cutoff=3.5):
    """
    Find residue-residue contacts between two chains within distance_cutoff (Angstroms).

    Returns
    -------
    (chain1_contacts, chain2_contacts) : tuple of lists
        Each is a sorted list of (residue_name, local_pdb_position) tuples.
    """
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("complex", pdb_file)
    model = structure[0]

    chain1_atoms = [atom for res in model[chain1_id] if is_aa(res) for atom in res]
    chain2_atoms = [atom for res in model[chain2_id] if is_aa(res) for atom in res]

    ns = NeighborSearch(chain2_atoms)

    chain1_contacts = set()
    chain2_contacts = set()

    for atom in chain1_atoms:
        neighbors = ns.search(atom.coord, distance_cutoff)
        for neighbor in neighbors:
            res1 = atom.get_parent()
            res2 = neighbor.get_parent()
            if res1 != res2:
                chain1_contacts.add((res1.get_resname(), res1.get_id()[1]))
                chain2_contacts.add((res2.get_resname(), res2.get_id()[1]))

    return (
        sorted(chain1_contacts, key=lambda x: x[1]),
        sorted(chain2_contacts, key=lambda x: x[1]),
    )


def find_chain_ion_contacts(pdb_file, chain_id="A", ion_resname="CA", distance_cutoff=3.5):
    """
    Find residues of a chain coordinating a given ion (default calcium, "CA").

    Returns a sorted list of (residue_name, local_pdb_position) tuples.
    """
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("complex", pdb_file)
    model = structure[0]

    chain_atoms = [atom for res in model[chain_id] if is_aa(res) for atom in res]

    ion_atoms = []
    for chain in model:
        for res in chain:
            if res.get_resname() == ion_resname:
                ion_atoms.extend(res.get_atoms())

    ns = NeighborSearch(ion_atoms)

    contacting_residues = set()
    for atom in chain_atoms:
        neighbors = ns.search(atom.coord, distance_cutoff)
        for neighbor in neighbors:
            res = atom.get_parent()
            contacting_residues.add((res.get_resname(), res.get_id()[1]))

    return sorted(contacting_residues, key=lambda x: x[1])


def format_contacts_as_dataframe(chain1_contacts, chain2_contacts, chain1_id="A", chain2_id="B"):
    """
    Combine the two contact lists from find_contacts() into a single dataframe,
    padding the shorter list so both columns are the same length.
    """
    sorted_a = sorted(chain1_contacts, key=lambda x: x[1])
    sorted_b = sorted(chain2_contacts, key=lambda x: x[1])

    max_len = max(len(sorted_a), len(sorted_b))
    sorted_a += [("", "")] * (max_len - len(sorted_a))
    sorted_b += [("", "")] * (max_len - len(sorted_b))

    return pd.DataFrame({
        f"Chain {chain1_id} Residue": [res[0] for res in sorted_a],
        f"Chain {chain1_id} PDB Pos": [res[1] for res in sorted_a],
        f"Chain {chain2_id} Residue": [res[0] for res in sorted_b],
        f"Chain {chain2_id} PDB Pos": [res[1] for res in sorted_b],
    })
