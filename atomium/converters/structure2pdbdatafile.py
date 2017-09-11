"""Contains functions for converting PdbDataFiles from AtomicStructures."""

from ..files.pdbdatafile import PdbDataFile

def structure_to_pdb_data_file(structure, model=1):
    """Converts an :py:class:`.AtomicStructure` to a :py:class:`.PdbDataFile`.

    :param AtomicStructure atom: The structure to convert.
    :param int model: The model number of the structure.
    :rtype: ``PdbDataFile``"""

    data_file = PdbDataFile()
    data_file.atoms, data_file.heteroatoms = [], []
    for atom in structure.atoms():
        if atom.residue():
            data_file.atoms.append(atom_to_atom_dict(atom, model=model))
        else:
            data_file.heteroatoms.append(atom_to_atom_dict(atom, model=model))
    data_file.atoms = sorted(data_file.atoms, key=lambda k: k["atom_id"])
    data_file.heteroatoms = sorted(data_file.heteroatoms, key=lambda k: k["atom_id"])
    data_file.connections = atoms_to_connections(structure.atoms())
    return data_file


def atom_to_atom_dict(atom, model=1):
    """Converts an :py:class:`.Atom` to an atom ``dict``.

    :param Atom atom: The atom to convert.
    :param int model: The model number that the atom is in.
    :rtype: ``dict``"""

    residue_name, chain_id, residue_id, insert_code = None, None, None, None
    if atom.residue():
        id_ = atom.residue().residue_id()
        residue_name = atom.residue().name()
        chain_id = atom.chain().chain_id() if atom.chain() else None
        residue_id = int("".join([c for c in id_ if c.isdigit()]))
        insert_code = id_[-1] if id_ and id_[-1].isalpha() else None
    elif atom.molecule():
        id_ = atom.molecule().molecule_id()
        residue_name = atom.molecule().name()
        chain_id = id_[0] if id_ and id_[0].isalpha() else None
        residue_id = int("".join([c for c in id_ if c.isdigit()]))
    return {
     "atom_id": atom.atom_id(), "atom_name": atom.name(), "alt_loc": None,
     "residue_name": residue_name,
     "chain_id": chain_id, "residue_id": residue_id, "insert_code": insert_code,
     "x": atom.x(), "y": atom.y(), "z": atom.z(),
     "occupancy": 1.0,
     "element": atom.element(), "charge": atom.charge(),
     "temperature_factor": atom.bfactor() if atom.bfactor() else None,
     "model": model
    }


def atoms_to_connections(atoms):
    """Gets a connections ``list`` from a collection of atoms. Only atoms that
    are not part of residues will have their bonds used.

    :param atoms: The :py:class:`.Atom` objects to use.
    :rtype: ``list``"""

    connections = []
    for atom in atoms:
        if not atom.residue():
            connections.append({
             "atom": atom.atom_id(),
             "bond_to": sorted([a.atom_id() for a in atom.bonded_atoms()])
            })
    return sorted(connections, key=lambda k: k["atom"])
