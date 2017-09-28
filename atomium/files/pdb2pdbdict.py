"""This module handles the conversion of Pdb objects to PDB data
dictionaries."""

from .pdbstring2pdbdict import atoms_to_chains, atoms_to_residues

def pdb_to_pdb_dict(pdb):
    """Converts a :py:class:`.Pdb` to a data ``dict``

    :param Pdb pdb: The Pdb to save..
    :rtype: ``dict``"""

    return {
     "deposition_date": pdb._deposition_date,
     "code": pdb._code,
     "title": pdb._title,
     "models": [structure_to_model_dict(model) for model in pdb._models]
    }


def structure_to_model_dict(structure):
    """Converts an :py:class:`.AtomicStructure` to a model ``dict``

    :param AtomicStructure structure: the structure to convert.
    :rtype: ``dict``"""

    atoms, heteroatoms = [], []
    structure_atoms = sorted(
     structure.atoms(), key=lambda a: a.atom_id() if a.atom_id() else 10000000
    )
    for atom in structure_atoms:
        list_ = heteroatoms if atom.chain() is None else atoms
        list_.append(atom_to_atom_dict(atom))
    chains = atoms_to_chains(atoms)
    molecules = atoms_to_residues(heteroatoms)
    return {"chains": chains, "molecules": molecules}


def atom_to_atom_dict(atom):
    """Converts an :py:class:`.Atom` to an atom ``dict``

    :param Atom atom: the atom to convert.
    :rtype: ``dict``"""

    id_, residue_name, chain_id, residue_id, insert_code = "", "", "", "", ""
    if atom.residue():
        id_ = atom.residue().residue_id()
        residue_name = atom.residue().name()
        chain_id = atom.chain().chain_id() if atom.chain() is not None else ""
        residue_id = int("".join([c for c in id_ if c.isdigit()]))
        insert_code = id_[-1] if id_ and id_[-1].isalpha() else ""
    elif atom.molecule():
        id_ = atom.molecule().molecule_id()
        residue_name = atom.molecule().name()
        chain_id = id_[0] if id_ and id_[0].isalpha() else None
        residue_id = int("".join([c for c in id_ if c.isdigit()]))
    return {
     "atom_id": atom.atom_id(), "atom_name": atom.name(), "alt_loc": None,
     "residue_name": residue_name, "full_id": id_,
     "chain_id": chain_id, "residue_id": residue_id, "insert_code": insert_code,
     "x": atom.x(), "y": atom.y(), "z": atom.z(),
     "occupancy": 1.0,
     "element": atom.element(), "charge": atom.charge(),
     "temp_factor": atom.bfactor() if atom.bfactor() else None,
    }
