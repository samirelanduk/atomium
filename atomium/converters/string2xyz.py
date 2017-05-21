from ..xyz.xyz import Xyz
from ..structures.models import Model
from ..structures.atoms import Atom
from .strings import string2lines

def string_to_xyz(s):
    """Converts a string taken from a .xyz file and turns it into a
    :py:class:`.Xyz` object.

    :param str s: The string to convert."""

    lines = string2lines(s)
    xyz = Xyz()
    if len(lines) > 2:
        try:
            int(lines[0])
            lines = lines[1:]
        except:
            pass
        try:
            element, x, y, z = lines[0].split()
            float(x), float(y), float(z)
        except:
            xyz._comment = lines[0]
            lines = lines[1:]
    xyz._model = Model()
    for line in lines:
        element, x, y, z = line.split()
        xyz._model._atoms.add(Atom(element, float(x), float(y), float(z)))
    return xyz
