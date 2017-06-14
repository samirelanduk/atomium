from math import sqrt
from matrices.checks import is_numeric

"""Contains the classes for atoms and their bonds."""

class Atom:
    """Represents an atom in three dimensional space. Every atom has an element
    and a set of Cartesian coordinates.

    :param str element: The atom's element symbol. It doesn't need to be on the\
    Periodic Table, but it does need to be 1 or 2 characters.
    :param number x: The atom's x coordinate.
    :param number y: The atom's y coordinate.
    :param number z: The atom's z coordinate.
    :param int atom_id: A unique integer ID for the atom. The class keeps track\
    of IDs that have already been used, though you can free up the ID by\
    changing or garbage collecting the atom that has the ID you want.
    :raises TypeError: if the element is not str.
    :raises ValueError: if the element is not 1 or 2 characters.
    :raises TypeError: if the coordinates are not numeric.
    :raises TypeError: if the atom_id is not int.
    :raises ValueError: if you give an atom_id that has already been used."""

    known_ids = set()

    def __init__(self, element, x, y, z, atom_id=None):
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
        if atom_id is not None and not isinstance(atom_id, int):
            raise TypeError("ID {} is not an integer".format(atom_id))
        if atom_id in Atom.known_ids:
            raise ValueError("There's already an atom of ID {}".format(atom_id))
        self._element = element
        self._x = x
        self._y = y
        self._z = z
        self._id = atom_id
        if atom_id is not None: Atom.known_ids.add(atom_id)
        self._bonds = set()


    def __repr__(self):
        return "<{} Atom at ({}, {}, {})>".format(
         self._element, self._x, self._y, self._z
        )


    def __setattr__(self, attr, value):
        if attr == "_id" and "_id" in self.__dict__ and value is not None:
            Atom.known_ids.remove(self._id)
            Atom.known_ids.add(value)
        self.__dict__[attr] = value


    def __del__(self):
        if "_id" in self.__dict__ and self._id in Atom.known_ids:
            Atom.known_ids.remove(self._id)


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
        coordinate will be updated, but it must be numeric.

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
        coordinate will be updated, but it must be numeric.

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
        coordinate will be updated, but it must be numeric.

        :param number z: If given, the atom's z coordinate will be set to this.
        :raises TypeError: if the z coordinate given is not numeric.
        :rtype: ``int`` or ``float``"""

        if z is None:
            return self._z
        else:
            if not is_numeric(z):
                raise TypeError("z coordinate '{}' is not numeric".format(z))
            self._z = z


    def atom_id(self, atom_id=None):
        """Returns the atom's unique integer ID. If a value is given, the ID
        will be updated, provided it is a unique integer.

        :param int atom_id: If given, the ID will be set to this.
        :raises TypeError: if the ID given is not numeric.
        :raises ValueError: if the ID given is already in use."""

        if atom_id is None:
            return self._id
        else:
            if not isinstance(atom_id, int):
                raise TypeError("Atom ID '{}' is not int".format(atom_id))
            self._id = atom_id


    def bonds(self):
        """Returns the :py:class:`.Bond` objects that the atom is connected to.

        :rtype: ``set``"""

        return set(self._bonds)


    def bonded_atoms(self):
        """Returns all the atoms that are bonded to this atom.

        :rtype: ``set``"""

        atoms = set()
        [atoms.update(bond.atoms()) for bond in self.bonds()]
        if atoms: atoms.remove(self)
        return atoms


    def bond(self, other):
        """Bonds the atom to some other atom by creating a :py:class:`.Bond`
        between them."""

        Bond(self, other)


    def unbond(self, other):
        """Breaks the bond between this atom and another.

        :param Atom other: The atom to unbond from.
        :raises TypeError: if something other than an :py:class:`Atom` is given.
        :raises ValueError: if the atom given isn't bonded to begin with."""

        if not isinstance(other, Atom):
            raise TypeError("Cannot unbond non-atom {}".format(other))
        for bond in self.bonds():
            if other in bond.atoms():
                bond.destroy()
                return
        raise ValueError("{} cannot unbond non-bonded {}".format(self, other))


    def mass(self):
        """Returns the atom's mass according to the Periodic Table, based on the
        atom's :py:meth:`element`. If the element doesn't match any symbol on
        the Periodic Table, a mass of 0 will be returned.

        The element lookup is case-insensitive.

        :rtype: ``float``"""

        return PERIODIC_TABLE.get(self._element.upper(), 0)


    def distance_to(self, other):
        """Returns the distance (in whatever units the coordinates are defined
        in) between this atom and another. You can also give a (x, y, z) tuple
        instead of another atom if you so wish.

        :param Atom other: The other atom.
        :raises TypeError: if something other than an :py:class:`Atom` is\
        given and that object isn't in the form (x, y, z).
        :rtype: ``float``"""

        x, y, z = None, None, None
        if not isinstance(other, Atom):
            try:
                assert len(other) == 3
                x, y, z = other
            except:
                raise TypeError("'{}' is not an Atom".format(other))
        else:
            x, y, z = other._x, other._y, other._z
        x_sum = pow((x - self._x), 2)
        y_sum = pow((y - self._y), 2)
        z_sum = pow((z - self._z), 2)
        return sqrt(x_sum + y_sum + z_sum)



class Bond:
    """Represents a chemical bond between an :py:class:`.Atom` and another. It
    doesn't matter what order the atoms are given, as they are just stored
    unordered in a set anyway.

    :param Atom atom1: The first atom.
    :param Atom atom2: The second atom.
    :raises TypeError: if non :py:class:`.Atom` objects are given.
    :raises ValueError: if the two atoms are the same atom."""

    def __init__(self, atom1, atom2):
        if not isinstance(atom1, Atom):
            raise TypeError("bond atom {} is not an atom".format(atom1))
        if not isinstance(atom2, Atom):
            raise TypeError("bond atom {} is not an atom".format(atom2))
        if atom1 is atom2:
            raise ValueError("Cannot bond atom {} to itself".format(atom1))
        self._atoms = set((atom1, atom2))
        atom1._bonds.add(self), atom2._bonds.add(self)


    def __repr__(self):
        atom1, atom2 = self._atoms
        return "<{}-{} Bond>".format(atom1.element(), atom2.element())


    def atoms(self):
        """Returns the two :py:class:`.Atom` objects that the bond connects.
        They are given as an unordered set.

        :rtype: ``set``"""

        return set(self._atoms)


    def length(self):
        """Returns the length of the bond, defined as the distance between its
        two atoms.

        :rtype: ``float``"""

        atom1, atom2 = self._atoms
        return atom1.distance_to(atom2)


    def destroy(self):
        """Destroys the bond and removes it from its atoms' ``set`` of bonds.

        Usually those are the only pointers to the bond and so calling this
        method will cause the bond to be removed by garbage collection. If you
        do have other variables pointing to the bond though, the object will
        remain in memory."""

        atom1, atom2 = self._atoms
        atom1._bonds.remove(self)
        atom2._bonds.remove(self)




PERIODIC_TABLE = {
 "H": 1.0079, "HE": 4.0026, "LI": 6.941, "BE": 9.0122, "B": 10.811,
 "C": 12.0107, "N": 14.0067, "O": 15.9994, "F": 18.9984, "NE": 20.1797,
 "NA": 22.9897, "MG": 24.305, "AL": 26.9815, "SI": 28.0855, "P": 30.9738,
 "S": 32.065, "CL": 35.453, "K": 39.0983, "AR": 39.948, "CA": 40.078,
 "SC": 44.9559, "TI": 47.867, "V": 50.9415, "CR": 51.9961, "MN": 54.938,
 "FE": 55.845, "NI": 58.6934, "CO": 58.9332, "CU": 63.546, "ZN": 65.39,
 "GA": 69.723, "GE": 72.64, "AS": 74.9216, "SE": 78.96, "BR": 79.904,
 "KR": 83.8, "RB": 85.4678, "SR": 87.62, "Y": 88.9059, "ZR": 91.224,
 "NB": 92.9064, "MO": 95.94, "TC": 98, "RU": 101.07, "RH": 102.9055,
 "PD": 106.42, "AG": 107.8682, "CD": 112.411, "IN": 114.818, "SN": 118.71,
 "SB": 121.76, "I": 126.9045, "TE": 127.6, "XE": 131.293, "CS": 132.9055,
 "BA": 137.327, "LA": 138.9055, "CE": 140.116, "PR": 140.9077, "ND": 144.24,
 "PM": 145, "SM": 150.36, "EU": 151.964, "GD": 157.25, "TB": 158.9253,
 "DY": 162.5, "HO": 164.9303, "ER": 167.259, "TM": 168.9342, "YB": 173.04,
 "LU": 174.967, "HF": 178.49, "TA": 180.9479, "W": 183.84, "RE": 186.207,
 "OS": 190.23, "IR": 192.217, "PT": 195.078, "AU": 196.9665, "HG": 200.59,
 "TL": 204.3833, "PB": 207.2, "BI": 208.9804, "PO": 209, "AT": 210, "RN": 222,
 "FR": 223, "RA": 226, "AC": 227, "PA": 231.0359, "TH": 232.0381, "NP": 237,
 "U": 238.0289, "AM": 243, "PU": 244, "CM": 247, "BK": 247, "CF": 251,
 "ES": 252, "FM": 257, "MD": 258, "NO": 259, "RF": 261, "LR": 262, "DB": 262,
 "BH": 264, "SG": 266, "MT": 268, "RG": 272, "HS": 277
}
