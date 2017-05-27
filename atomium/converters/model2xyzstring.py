"""Contains tools for turning a Model into a .xyz file."""

def model_to_xyz_string(model):
    """Takes a :py:class:`.Model` and turns it into a string in .xyz format.

    :param Model model: the Model to convert.
    :rtype: ``str``"""

    atom_num = str(len(model.atoms()))
    atoms = ["{}{:11.3f}{:11.3f}{:11.3f}".format(
     atom.element(), atom.x(), atom.y(), atom.z()
    ) for atom in sorted(model.atoms(), key=lambda k: k.element())]
    return "\n".join([atom_num, "", *atoms])
