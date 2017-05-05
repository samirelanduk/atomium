"""Contains the Atom class and its associated classes."""

class Atom:
    """Represents atoms in three dimensional space. Every atom has an element
    and a set of Cartesian coordinates.

    :param str element: The atom's element symbol.
    :param number x: The atom's x coordinate.
    :param number y: The atom's y coordinate.
    :param number z: The atom's z coordinate."""

    def __init__(self, element, x, y, z):
        self._element = element
        self._x = x
        self._y = y
        self._z = z
