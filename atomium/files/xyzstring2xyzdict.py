"""This module handles the conversion of .xyz filestrings to XYZ data
dictionaries."""

def xyz_string_to_xyz_dict(filestring):
    """Converts the string of a .xyz file to a parsed data ``dict``.

    :param str filestring: The filestring to parse.
    :rtype: ``dict``"""

    from .utilities import string_to_lines
    xyz_dict = {}
    lines = string_to_lines(filestring)
    extract_header(xyz_dict, lines)
    extract_structure(xyz_dict, lines)
    return xyz_dict


def extract_header(xyz_dict, lines):
    """Takes a ``dict`` and adds header information to it by parsing file lines.

    :param dict xyz_dict: the ``dict`` to update.
    :param list lines: the file lines to read from."""

    try:
        int(lines[0])
        lines.pop(0)
    except: pass
    try:
        float(lines[0].split()[1])
        xyz_dict["title"] = None
    except:
        xyz_dict["title"] = lines.pop(0).strip()


def extract_structure(xyz_dict, lines):
    """Takes a ``dict`` and adds structure information to it by parsing file
    lines.

    :param dict xyz_dict: the ``dict`` to update.
    :param list lines: the file lines to read from."""

    xyz_dict["atoms"] = [{
     "element": line.split()[0],
     "x": float(line.split()[1]),
     "y": float(line.split()[2]),
     "z": float(line.split()[3])
    } for line in lines]
