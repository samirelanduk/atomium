"""This module handles the conversion of XYZ data dictionaries to .xyz
filestrings."""

def xyz_dict_to_xyz_string(xyz_dict):
    """Converts a data ``dict`` to a .xyz filestring.

    :param dict xyz_dict: The data dictionary to pack.
    :rtype: ``str``"""

    from .utilities import lines_to_string
    lines = []
    pack_header(lines, xyz_dict)
    pack_structure(lines, xyz_dict)
    return lines_to_string(lines)


def pack_header(lines, xyz_dict):
    """Adds a .xyz title to a list of lines.

    :param list lines: The record lines to add to.
    :param dict xyz_dict: The data dictionary to pack."""

    lines.append(str(len(xyz_dict["atoms"])))
    if xyz_dict["title"] is not None:
        lines.append(xyz_dict["title"])


def pack_structure(lines, xyz_dict):
    """Adds .xyz atom lines to a list of lines.

    :param list lines: The record lines to add to.
    :param dict xyz_dict: The data dictionary to pack."""

    for atom in xyz_dict["atoms"]:
        lines.append("{}{:11.3f}{:11.3f}{:11.3f}".format(
         atom["element"], atom["x"], atom["y"], atom["z"]
        ))
