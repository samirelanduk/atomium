"""Contains functions for converting Models to PdbDataFiles."""

def atom_to_atom_dict(atom):
    """Converts an :py:class:`.Atom` to an atom ``dict``.

    :param Atom atom: The atom to convert."""
    
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
     "atom_id": atom.atom_id(), "atom_name": atom.atom_name(), "alt_loc": None,
     "residue_name": residue_name,
     "chain_id": chain_id, "residue_id": residue_id, "insert_code": insert_code,
     "x": atom.x(), "y": atom.y(), "z": atom.z(),
     "occupancy": 1.0, "temperature_factor": None,
     "element": atom.element(), "charge": atom.charge()
    }
