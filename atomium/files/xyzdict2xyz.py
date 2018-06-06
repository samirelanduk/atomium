"""This module handles the conversion of XYZ data dictionaries to Xyz
objects."""

from .xyz import Xyz
from ..models import Model, Atom

def xyz_dict_to_xyz(xyz_dict):
    """Converts a data ``dict`` to a :py:class:`.Xyz`

    :param dict xyz_dict: The data dictionary to load.
    :rtype: :py:class:`.Xiz`"""

    model = xyz_dict_to_model(xyz_dict)
    xyz = Xyz()
    xyz._model = model
    xyz._title = xyz_dict["title"]
    return xyz


def xyz_dict_to_model(xyz_dict):
    """Converts a data ``dict`` to a :py:class:`.Model`

    :param dict xyz_dict: The model dictionary to load.
    :rtype: :py:class:`.Model`"""

    atoms = [atom_dict_to_atom(atom) for atom in xyz_dict["atoms"]]
    return Model(*atoms)


def atom_dict_to_atom(atom_dict):
    """Converts an atom ``dict`` to a :py:class:`.Atom`

    :param dict atom_dict: The atom dictionary to load.
    :rtype: :py:class:`.Atom`"""

    return Atom(
     atom_dict["element"], atom_dict["x"], atom_dict["y"], atom_dict["z"]
    )
