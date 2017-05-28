"""Contains tools for turning a Model into a .xyz file."""

from ..structures.models import Model

def model_to_xyz_string(model, comment=""):
    """Takes a :py:class:`.Model` and turns it into a string in .xyz format.

    :param Model model: the Model to convert.
    :param str comment: A comment to add to the file (default is empty string).
    :rtype: ``str``"""

    if not isinstance(model, Model):
        raise TypeError("{} is not a Model".format(str(model)))
    if not isinstance(comment, str):
        raise TypeError("comment {} is not a string".format(str(comment)))
    atom_num = str(len(model.atoms()))
    atoms = ["{}{:11.3f}{:11.3f}{:11.3f}".format(
     atom.element(), atom.x(), atom.y(), atom.z()
    ) for atom in sorted(model.atoms(), key=lambda k: k.element())]
    return "\n".join([atom_num, comment, *atoms])
