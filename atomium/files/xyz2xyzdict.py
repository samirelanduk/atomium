"""This module handles the conversion of Xyz objects to XYZ data
dictionaries."""

def xyz_to_xyz_dict(xyz):
    """Converts a :py:class:`.Xyz` to a data ``dict``

    :param Xyz xyz: The Xyz to save..
    :rtype: ``dict``"""

    xyz_dict = structure_to_xyz_dict(xyz._model)
    xyz_dict["title"] = xyz._title
    return xyz_dict


def structure_to_xyz_dict(structure):
    """Converts an :py:class:`.AtomStructure` to a model ``dict``.

    :param AtomStructure structure: the structure to convert.
    :rtype: ``dict``"""

    return {
     "title": None,
     "atoms": [atom_to_atom_dict(atom) for atom in structure.atoms()]
    }


def atom_to_atom_dict(atom):
    """Converts an :py:class:`.Atom` to an atom ``dict``.

    :param Atom atom: the atom to convert.
    :rtype: ``dict``"""

    return {
     "element": atom.element, "x": atom.x, "y": atom.y, "z": atom.z,
    }
