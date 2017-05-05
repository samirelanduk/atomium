from matrices.checks import is_numeric

"""Contains the Atom class and its associated classes."""

class Atom:
    """Represents an atom in three dimensional space. Every atom has an element
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


    def __repr__(self):
        return "<{} Atom at ({}, {}, {})>".format(
         self._element, self._x, self._y, self._z
        )


    def element(self, element=None):
        """Returns the atom's element symbol. If a value is given, the element
        symbol will be updated, but it must be ``str`` and it must be 1 or 2
        characters.

        :param str element: If given, the atom's element will be set to this.
        :raises TypeError: if the element is not str.
        :raises ValueError: if the element given is not 1 or 2 characters.
        :rtype: ``str``"""

        if element is None:
            return self._element
        else:
            if not isinstance(element, str):
                raise TypeError("Element '{}' is not str".format(element))
            if not 0 < len(element) < 3:
                raise ValueError("Element {} isn't 1 - 2 chars".format(element))
            self._element = element


    def x(self, x=None):
        """Returns the atom's x coordinate. If a value is given, the x
        coordinate woll be updated, but it must be numeric.

        :param number x: If given, the atom's x coordinate will be set to this.
        :raises TypeError: if the x coordinate given is not numeric.
        :rtype: ``int`` or ``float``"""

        if x is None:
            return self._x
        else:
            if not is_numeric(x):
                raise TypeError("x coordinate '{}' is not numeric".format(x))
            self._x = x


    def y(self, y=None):
        """Returns the atom's y coordinate. If a value is given, the y
        coordinate woll be updated, but it must be numeric.

        :param number y: If given, the atom's y coordinate will be set to this.
        :raises TypeError: if the y coordinate given is not numeric.
        :rtype: ``int`` or ``float``"""

        if y is None:
            return self._y
        else:
            if not is_numeric(y):
                raise TypeError("y coordinate '{}' is not numeric".format(y))
            self._y = y


    def z(self, z=None):
        """Returns the atom's z coordinate. If a value is given, the z
        coordinate woll be updated, but it must be numeric.

        :param number z: If given, the atom's z coordinate will be set to this.
        :raises TypeError: if the z coordinate given is not numeric.
        :rtype: ``int`` or ``float``"""

        if z is None:
            return self._z
        else:
            if not is_numeric(z):
                raise TypeError("z coordinate '{}' is not numeric".format(z))
            self._z = z
