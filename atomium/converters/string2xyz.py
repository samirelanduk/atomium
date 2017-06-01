from ..xyz.xyz import Xyz
from ..structures.models import Model
from ..structures.atoms import Atom
from .strings import string2lines

def string_to_xyz(s):
    """Converts a string taken from a .xyz file and turns it into a
    :py:class:`.Xyz` object.

    :param str s: The string to convert."""

    lines = string2lines(s)
    remove_atom_num(lines)
    comment = extract_comment(lines)
    atoms = [parse_atom(line) for line in lines]
    atoms = filter(None, atoms)
    xyz = Xyz(comment)
    xyz._model = Model(*atoms)
    return xyz


def remove_atom_num(lines):
    """Takes a list of strings representing file lines and removes the first
    line if it is an integer.

    :param list lines: The file lines.
    :rtype: ``list``"""

    try:
        int(lines[0])
        lines.pop(0)
    except (ValueError, IndexError) as e: pass


def parse_atom(line):
    """Takes a line from an .xyz file and creates an :py:class:`.Atom` from it.
    If the line canot be parsed it returns ``None``.

    :param str line: The line to parse.
    :rtype: ``Atom`` or ``None``"""

    try:
        element, x, y, z = line.split()
        return Atom(element, float(x), float(y), float(z))
    except: pass


def extract_comment(lines):
    """Checks the first line in a list of lines, and if it isn't an atom line,
    removes and returns it. If it is an atom line, an empty string is returned.

    :param list lines: The file lines.
    :rtype: ``str``"""

    if not lines or parse_atom(lines[0]):
        return ""
    else:
        return lines.pop(0)
