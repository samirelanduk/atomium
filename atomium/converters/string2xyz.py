from ..xyz.xyz import Xyz
from ..structures.models import Model

def string_to_xyz(s):
    """Converts a string taken from a .xyz file and turns it into a
    :py:class:`.Xyz` object.

    :param str s: The string to convert."""
    
    xyz = Xyz()
    xyz._model = Model()
    return xyz
