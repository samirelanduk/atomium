from matrices.checks import is_numeric

"""Contains the Atom class and its associated classes."""

class Atom:
    """Represents atoms in three dimensional space. Every atom has an element
    and a set of Cartesian coordinates.

    :param str element: The atom's element symbol. It doesn't need to be on the\
    Periodic Table, but it does need to be 1 or 2 characters.
    :param number x: The atom's x coordinate.
    :param number y: The atom's y coordinate.
    :param number z: The atom's z coordinate.
    :raises TypeError: if the element is not str.
    :raises ValueError: if the element is not 1 or 2 characters.
    :raises TypeError: if the coordinates are not numeric."""

    def __init__(self, element, x, y, z):
        if not isinstance(element, str):
            raise TypeError("Element '{}' is not str".format(element))
        if not 0 < len(element) < 3:
            raise ValueError("Element {} is not 1 or 2 chars".format(element))
        if not is_numeric(x):
            raise TypeError("x coordinate '{}' is not numeric".format(x))
        if not is_numeric(y):
            raise TypeError("y coordinate '{}' is not numeric".format(y))
        if not is_numeric(z):
            raise TypeError("z coordinate '{}' is not numeric".format(z))
        self._element = element
        self._x = x
        self._y = y
        self._z = z
