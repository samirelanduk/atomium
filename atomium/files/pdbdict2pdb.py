"""This module handles the conversion of PDB data dictionaries to Pdb
objects."""

from .pdb import Pdb
from ..structures import Model, Chain, Residue, Molecule, Atom

def pdb_dict_to_pdb(pdb_dict):
    """Converts a data ``dict`` to a :py:class:`.Pdb`

    :param dict pdb_dict: The data dictionary to load.
    :rtype: :py:class:`.Pdb`"""

    pdb = Pdb()
    pdb._deposition_date = pdb_dict["deposition_date"]
    pdb._code = pdb_dict["code"]
    pdb._title = pdb_dict["title"]
    pdb._models = [model_dict_to_model(d) for d in pdb_dict["models"]]
    return pdb


def model_dict_to_model(model_dict):
    """Converts a model ``dict`` to a :py:class:`.Model`

    :param dict model_dict: The model dictionary to load.
    :rtype: :py:class:`.Model`"""

    model = Model()
    for chain in model_dict["chains"]:
        model.add_chain(chain_dict_to_chain(chain))
    for molecule in model_dict["molecules"]:
        model.add_molecule(residue_dict_to_residue(molecule, molecule=True))
    return model


def chain_dict_to_chain(chain_dict):
    """Converts a chain ``dict`` to a :py:class:`.Chain`

    :param dict chain_dict: The chain dictionary to load.
    :rtype: :py:class:`.Chain`"""

    residues = [residue_dict_to_residue(res) for res in chain_dict["residues"]]
    return Chain(*residues, chain_id=chain_dict["chain_id"])


def residue_dict_to_residue(residue_dict, molecule=False):
    """Converts a residue ``dict`` to a :py:class:`.Residue` or
    :py:class:`.Molecule`.

    :param dict residue_dict: The residue dictionary to load.
    :param bool molecule: if ``True``, a :py:class:`.Molecule` will be returned\
    instead of a :py:class:`.Residue`.
    :rtype: :py:class:`.Residue` or :py:class:`.Molecule`"""

    atoms = [atom_dict_to_atom(atom) for atom in residue_dict["atoms"]]
    MolClass = Molecule if molecule else Residue
    residue = MolClass(*atoms, name=residue_dict["name"])
    residue._id = residue_dict["id"]
    return residue


def atom_dict_to_atom(atom_dict):
    """Converts an atom ``dict`` to a :py:class:`.Atom`

    :param dict atom_dict: The atom dictionary to load.
    :rtype: :py:class:`.Atom`"""

    return Atom(
     atom_dict["element"], atom_dict["x"], atom_dict["y"], atom_dict["z"],
     atom_id=atom_dict["atom_id"], name=atom_dict["atom_name"],
     charge=atom_dict["charge"], bfactor=atom_dict["temp_factor"]
    )
