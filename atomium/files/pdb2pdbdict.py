"""This module handles the conversion of Pdb objects to PDB data
dictionaries."""

from .pdbstring2pdbdict import atoms_to_chains, atoms_to_residues
from ..models.data import CODES

def pdb_to_pdb_dict(pdb):
    """Converts a :py:class:`.Pdb` to a data ``dict``.

    :param Pdb pdb: The Pdb to save.
    :rtype: ``dict``"""

    pdb_dicts = [structure_to_pdb_dict(model) for model in pdb._models]
    pdb_dict = pdb_dicts[0]
    for d in pdb_dicts[1:]:
        pdb_dict["models"].append(d["models"][0])
    pdb_dict["sequences"] = sequences_from_model(pdb._models[0])
    pdb_dict["deposition_date"] = pdb._deposition_date
    pdb_dict["code"] = pdb._code
    pdb_dict["title"] = pdb._title
    pdb_dict["resolution"] = pdb._resolution
    pdb_dict["rfactor"] = pdb._rfactor
    pdb_dict["rfree"] = pdb._rfree
    pdb_dict["rcount"] = pdb._rcount
    pdb_dict["organism"] = pdb._organism
    pdb_dict["expression_system"] = pdb._expression_system
    pdb_dict["technique"] = pdb._technique
    pdb_dict["classification"] = pdb._classification
    pdb_dict["keywords"] = pdb._keywords
    pdb_dict["biomolecules"] = pdb._biomolecules
    return pdb_dict


def structure_to_pdb_dict(structure):
    """Converts an :py:class:`.AtomStructure` to a model ``dict``.

    :param AtomStructure structure: the structure to convert.
    :rtype: ``dict``"""

    atoms, heteroatoms = [], []
    structure_atoms = sorted(
     structure.atoms(), key=lambda a: a.id
    )
    for atom in structure_atoms:
        list_ = heteroatoms if atom.residue is None else atoms
        list_.append(atom_to_atom_dict(atom))
    chains = atoms_to_chains(atoms, heteroatoms)
    connections = structure_to_connections(structure)
    return {
     "models": [chains], "connections": connections,
     "deposition_date": None, "code": None, "title": None, "resolution": None,
     "organism": None, "expression_system": None, "technique": None,
     "classification": None, "rfactor": None, "rfree": None, "rcount": None,
     "keywords": [], "biomolecules": [], "sequences": {}
    }


def atom_to_atom_dict(atom):
    """Converts an :py:class:`.Atom` to an atom ``dict``.

    :param Atom atom: the atom to convert.
    :rtype: ``dict``"""

    id_, residue_name, chain_id, residue_id, insert_code = "", "", "", "", ""
    if atom.residue:
        id_ = atom.residue.id
        residue_name = atom.residue.name
        chain_id = atom.chain.id if atom.chain is not None else ""
        residue_id = int("".join([c for c in id_ if c.isdigit()]))
        insert_code = id_[-1] if id_ and id_[-1].isalpha() else ""
    elif atom.ligand:
        id_ = atom.ligand.id
        residue_name = atom.ligand.name
        chain_id = atom.chain.id if atom.chain is not None else ""
        residue_id = int("".join([c for c in id_ if c.isdigit()]))
    return {
     "atom_id": atom.id, "atom_name": atom.name, "alt_loc": None,
     "residue_name": residue_name, "full_id": id_,
     "chain_id": chain_id, "residue_id": residue_id, "insert_code": insert_code,
     "x": atom.x, "y": atom.y, "z": atom.z,
     "occupancy": 1.0, "anisotropy": atom.anisotropy,
     "element": atom.element, "charge": atom.charge,
     "temp_factor": atom.bfactor if atom.bfactor else None,
    }


def structure_to_connections(structure):
    """Gets a connections ``list`` from an :py:class:`.AtomStructure`. Only
    atoms that are not part of residues will have their bonds used.

    :param AtomStructure structure: The structure to use.
    :rtype: ``list``"""

    connections = []
    for atom in structure.atoms():
        if not atom.residue:
            connections.append({
             "atom": atom.id,
             "bond_to": sorted([a.id for a in atom.bonded_atoms])
            })
    return sorted(connections, key=lambda k: k["atom"])


def sequences_from_model(model):
    """Takes a model and creates a sequences ``dict`` from its chains.

    :param Model model: the model to parse.
    :rtype: ``dict``"""

    sequences = {}
    lookup = {v: k for k, v in CODES.items() if k not in ["HIP", "HIE"]}
    for chain in model.chains():
        if chain.rep_sequence:
            sequences[chain.id] = [
             lookup.get(c, "???") for c in chain.rep_sequence
            ]
    return sequences
