"""Contains tools for turning an AtomicStructure into a .xyz file."""

from ..structures.models import AtomicStructure

def structure_to_xyz_string(structure, comment=""):
    """Takes an :py:class:`.AtomicStructure` and turns it into a string in
    .xyz format.

    :param AtomicStructure structure: the structure to convert.
    :param str comment: A comment to add to the file (default is empty string).
    :rtype: ``str``"""

    if not isinstance(structure, AtomicStructure):
        raise TypeError("{} is not a AtomicStructure".format(str(structure)))
    if not isinstance(comment, str):
        raise TypeError("comment {} is not a string".format(str(comment)))
    atom_num = str(len(structure.atoms()))
    atoms = ["{}{:11.3f}{:11.3f}{:11.3f}".format(
     atom.element(), atom.x(), atom.y(), atom.z()
    ) for atom in sorted(structure.atoms(), key=lambda k: k.element())]
    return "\n".join([atom_num, comment, *atoms])
