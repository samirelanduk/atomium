"""Contains functions for dealing with the .xyz file format."""

import re
from copy import deepcopy
from .data import DATA_DICT, ATOM_DICT, MODEL_DICT

def xyz_string_to_xyz_dict(filestring):
    """Takes the filecontents of a .xyz file and produces an atomium data
    dictionary from them.

    :param str filestring: The contents of a .xyz file.
    :rtype: ``dict``"""

    lines = list(filter(lambda l: bool(l.strip()), filestring.split("\n")))
    xyz_dict = {"header_lines": [], "atom_lines": []}
    element_pattern = r"[A-Z]"
    float_pattern = r"\s{1,}\-?\d*\.?\d*"
    pattern = re.compile(element_pattern + (3 * float_pattern))
    while not pattern.match(lines[0]):
        xyz_dict["header_lines"].append(lines.pop(0))
    for line in lines:
        xyz_dict["atom_lines"].append(line)
    return xyz_dict


def xyz_dict_to_data_dict(xyz_dict):
    """Takes a basic .xyz dict and turns it into a standard atomium data
    dictionary.

    :param dict xyz_dict: The .xyz dictionary.
    :rtype: ``dict``"""

    d = deepcopy(DATA_DICT)
    d["description"]["title"] = xyz_dict["header_lines"][-1]
    d["models"].append(deepcopy(MODEL_DICT))
    for atom_line in xyz_dict["atom_lines"]:
        a = deepcopy(ATOM_DICT)
        chunks = atom_line.split()
        a["element"] = chunks[0]
        a["x"], a["y"], a["z"] = [float(n) for n in chunks[1:]]
        d["models"][0]["atoms"].append(a)
    return d


def file_to_xyz_string(file):
    """Takes a :py:class:`.File` and turns it into a .xyz filestring that
    represents it.

    :param File file: the File to convert.
    :rtype: ``str``"""

    lines = []
    lines.append(str(len(file.model.atoms())))
    if file.title is not None:
        lines.append(file.title)
    for atom in file.model.atoms():
        lines.append("{}{:11.3f}{:11.3f}{:11.3f}".format(
         atom.element, atom.x, atom.y, atom.z
        ))
    return "\n".join(lines)
