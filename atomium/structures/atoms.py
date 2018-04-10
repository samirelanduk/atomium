"""This module contains the classes for atoms and their bonds."""

import math
import numpy as np

def atom_query(func):
    """Decorator which can be applied to any function which returns atoms. It
    lets you query the output.

    :param function func: The function to enhance.
    :rtype: ``function``"""

    def new(*args, id=None, name=None,
     element=None, hydrogen=True, het=True, metal=True, **kwargs):
        atoms = func(*args, **kwargs)
        if id:
            atoms = set(filter(lambda a: a._id == id, atoms))
        if name:
            atoms = set(filter(lambda a: a._name == name, atoms))
        if element:
            atoms = set(filter(
             lambda a: a._element.lower() == element.lower(), atoms
            ))
        if not hydrogen:
            atoms = set(filter(lambda a: a._element.lower() != "h", atoms))
        if not het:
            atoms = set(filter(lambda a: a._residue is not None, atoms))
        if not metal:
            atoms = set(filter(lambda a: a._element.upper() not in METALS, atoms))
        return atoms
    new.__name__ = func.__name__
    new.__doc__ = func.__doc__
    return new



class Atom:
    """Represents an atom in three dimensional space. Every atom has an element
    and a set of Cartesian coordinates.

    :param str element: The atom's element symbol. It doesn't need to be on the\
    Periodic Table.
    :param number x: The atom's x coordinate.
    :param number y: The atom's y coordinate.
    :param number z: The atom's z coordinate.
    :param int id: A unique integer ID for the atom. This is supposed to\
    be unique. If you do not assign one, one will be generated.
    :param str name: The atom's name.
    :param number charge: The charge of the atom.
    :param number bfactor: The B-factor of the atom (its uncertainty).
    :raises TypeError: if the element is not str.
    :raises TypeError: if the coordinates are not numeric.
    :raises TypeError: if the id is not int.
    :raises TypeError: if the name is not str.
    :raises TypeError: if the charge is not numeric.
    :raises TypeError: if the bfactor is not numeric."""

    def __init__(self, element, x=0, y=0, z=0, id=0, name=None, charge=0,
                 bfactor=0):
        if not isinstance(element, str):
            raise TypeError("Element '{}' is not str".format(element))
        if any(not isinstance(coord, (int, float)) for coord in (x, y, z)):
            raise TypeError("Coordinates {} not numeric".format((x, y, z)))
        if not isinstance(id, int):
            raise TypeError("ID {} is not an integer".format(id))
        if name is not None and not isinstance(name, str):
            raise TypeError("name {} is not a string".format(name))
        if not isinstance(charge, (float, int)):
            raise TypeError("charge '{}' is not numeric".format(charge))
        if not isinstance(bfactor, (float, int)):
            raise TypeError("bfactor '{}' is not numeric".format(bfactor))
        self._element = element
        self._x = x
        self._y = y
        self._z = z
        self._id = id
        self._name = name
        self._charge = charge
        self._bfactor = bfactor
        self._bonds = set()
        self._residue, self._chain, self._molecule = None, None, None
        self._model, self._complex = None, None


    def __repr__(self):
        return "<{} Atom {} at ({}, {}, {})>".format(
         self._element, self._id ,
         self._x, self._y, self._z
        )


    @property
    def element(self):
        """The atom's element symbol. This is used to calculate its mass using a
        Periodic Table.

        :raises TypeError: if the element is set to non-str.
        :rtype: ``str``"""

        return self._element


    @element.setter
    def element(self, element):
        if not isinstance(element, str):
            raise TypeError("Element '{}' is not str".format(element))
        self._element = element


    @property
    def x(self):
        """The atom's x-coordinate.

        :raises TypeError: if the x coordinate given is not numeric.
        :rtype: ``float``"""

        return self._x


    @x.setter
    def x(self, x):
        if not isinstance(x, (float, int)):
            raise TypeError("x coordinate '{}' is not numeric".format(x))
        self._x = x


    @property
    def y(self):
        """The atom's y-coordinate.

        :raises TypeError: if the y coordinate given is not numeric.
        :rtype: ``float``"""

        return self._y


    @y.setter
    def y(self, y):
        if not isinstance(y, (float, int)):
            raise TypeError("y coordinate '{}' is not numeric".format(y))
        self._y = y


    @property
    def z(self):
        """The atom's z-coordinate.

        :raises TypeError: if the z coordinate given is not numeric.
        :rtype: ``float``"""

        return self._z


    @z.setter
    def z(self, z):
        if not isinstance(z, (float, int)):
            raise TypeError("z coordinate '{}' is not numeric".format(z))
        self._z = z


    @property
    def id(self):
        """Yhe atom's unique integer ID.

        :raises TypeError: if the ID given is not an integer.
        :rtype: ``int``"""

        return self._id


    @property
    def name(self):
        """The atom's name. This is often used to determine what 'kind' of atom
        it is.

        :raises TypeError: if the name is set to non-str.
        :rtype: ``str``"""

        return self._name


    @name.setter
    def name(self, name):
        if not isinstance(name, str):
            raise TypeError("Name '{}' is not str".format(name))
        self._name = name


    @property
    def bfactor(self):
        """The atom's B-factor - the uncertainty in its position.

        :raises TypeError: if the bfactor is set to be non-numeric.
        :rtype: ``float``"""

        return self._bfactor


    @bfactor.setter
    def bfactor(self, bfactor=None):
        if not isinstance(bfactor, (float, int)):
            raise TypeError("bfactor '{}' is not numeric".format(bfactor))
        self._bfactor = bfactor


    @property
    def charge(self):
        """The atom's charge.

        :raises TypeError: if the charge is set to be non-numeric.
        :rtype: ``float``"""

        return self._charge


    @charge.setter
    def charge(self, charge=None):
        if not isinstance(charge, (float, int)):
            raise TypeError("charge '{}' is not numeric".format(charge))
        self._charge = charge


    @property
    def location(self):
        """The atom's Cartesian coordinates.

        :rtype: ``tuple``"""

        return (self._x, self._y, self._z)


    def trim(self, places):
        """Rounds the coordinate values to a given number of decimal places.
        Useful for removing floating point rounding errors after transformation.

        :param int places: The number of places to round the coordinates to. If\
        ``None``, no rounding will be done."""

        if places is not None:
            self._x = round(self._x, places)
            self._y = round(self._y, places)
            self._z = round(self._z, places)


    def translate(self, dx=0, dy=0, dz=0, trim=12):
        """Translates an atom in 3D space.

        :param float dx: The distance to move in the x direction.
        :param float dy: The distance to move in the y direction.
        :param float dz: The distance to move in the z direction.
        :param int trim: The amount of rounding to do to the atom's coordinates\
        after translating - the default is 12 decimal places but this can be\
        set to ``None`` if no rounding is to be done."""

        try:
            dx, dy, dz = dx
        except TypeError: pass
        self._x += dx
        self._y += dy
        self._z += dz
        self.trim(trim)


    def move_to(self, x, y, z):
        """Moves the atom to the coordinates given.

        :param number x: The atom's new x coordinate.
        :param number y: The atom's new y coordinate.
        :param number z: The atom's new z coordinate.
        :raises TypeError: if non-numeric coordinates are given."""

        if any(not isinstance(coord, (int, float)) for coord in (x, y, z)):
            raise TypeError("Coordinates {} not numeric".format((x, y, z)))
        self._x, self._y, self._z = x, y, z


    def generate_rotation_matrix(self, angle, axis):
        """Generates a rotation matrix that would rotate this atom by a given
        angle around a given axis.

        :param float angle: Angle in radians.
        :param str axis: The axis to rotate around. Can be `x`, `y` or `z`."""

        axis = [1 if i == ["x", "y", "z"].index(axis) else 0 for i in range(3)]
        axis = np.asarray(axis)
        axis = axis / math.sqrt(np.dot(axis, axis))
        a = math.cos(angle / 2)
        b, c, d = -axis * math.sin(angle / 2)
        aa, bb, cc, dd = a * a, b * b, c * c, d * d
        bc, ad, ac, ab, bd, cd = b * c, a * d, a * c, a * b, b * d, c * d
        return np.array([
         [aa + bb - cc - dd, 2 * (bc + ad), 2 * (bd - ac)],
         [2 * (bc - ad), aa + cc - bb - dd, 2 * (cd + ab)],
         [2 * (bd + ac), 2 * (cd - ab), aa + dd - bb - cc]
        ])


    def rotate(self, angle, axis, degrees=False, trim=12):
        """Rotates an atom in 3D space.

        :param float angle: Angle in radians.
        :param str axis: The axis to rotate around. Can be `x`, `y` or `z`.
        :param bool degrees: if ``True`` the angle will be interpreted as\
        degrees.
        :param int trim: The amount of rounding to do to the atom's coordinates\
        after translating - the default is 12 decimal places but this can be\
        set to ``None`` if no rounding is to be done."""

        if axis not in ("x", "y", "z"):
            raise ValueError("{} is not a valid axis".format(axis))
        angle = math.radians(angle) if degrees else angle
        matrix = self.generate_rotation_matrix(angle, axis)
        vector = matrix.dot(self.location)
        self._x, self._y, self._z = vector
        self.trim(trim)


    @property
    def mass(self):
        """The atom's molar mass according to the Periodic Table, based on the
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

        try:
            x, y, z = other
        except:
            if not isinstance(other, Atom):
                raise TypeError("'{}' is not an Atom".format(other))
            x, y, z = other.location
        x_sum = pow((x - self._x), 2)
        y_sum = pow((y - self._y), 2)
        z_sum = pow((z - self._z), 2)
        return math.sqrt(x_sum + y_sum + z_sum)


    @property
    def residue(self):
        """Returns the :py:class:`.Residue` the atom is part of, or ``None`` if
        it is not part of one.

        :rtype: ``Residue``"""

        return self._residue


    @property
    def chain(self):
        """Returns the :py:class:`.Chain` the atom is part of, or ``None`` if
        it is not part of one.

        :rtype: ``Chain``"""

        return self._chain


    @property
    def molecule(self):
        """Returns the :py:class:`.Molecule` the atom is part of, or ``None`` if
        it is not part of one.

        :rtype: ``Molecule``"""

        return self._molecule


    @property
    def model(self):
        """Returns the :py:class:`.Model` the atom is part of, or ``None`` if
        it is not part of one.

        :rtype: ``Model``"""

        return self._model


    @property
    def complex(self):
        """Returns the :py:class:`.Complex` the atom is part of, or ``None`` if
        it is not part of one.

        :rtype: ``Complex``"""

        return self._complex


    @property
    def bonds(self):
        """The atomic :py:class:`.Bond` objects that the atom is associated
        with.

        :rtype: ``set``"""

        return set(self._bonds)


    @atom_query
    def bonded_atoms(self):
        """Returns all the atoms that are bonded to this atom.

        :param int id: if given, only atoms whose ID matches this will be\
        returned.
        :param str name: if given, only atoms whose name matches this will be\
        returned.
        :param str element: if given, only atoms whose element matches this\
        will be returned.
        :param bool hydrogen: If ``False``, hydrogen atoms will be excluded.
        :param bool het: If ``False``, non-chain atoms will be excluded.
        :param bool metal: If ``False``, metal atoms will be excluded.
        :rtype: ``set``"""

        atoms = set()
        [atoms.update(bond.atoms()) for bond in self.bonds]
        if atoms: atoms.remove(self)
        return atoms


    def bond_to(self, other):
        """Bonds the atom to some other atom by creating a :py:class:`.Bond`
        between them."""

        if other not in self.bonded_atoms():
            Bond(self, other)


    def unbond_from(self, other):
        """Breaks the bond between this atom and another.

        :param Atom other: The atom to unbond from.
        :raises TypeError: if something other than an :py:class:`Atom` is given.
        :raises ValueError: if the atom given isn't bonded to begin with."""

        if not isinstance(other, Atom):
            raise TypeError("Cannot unbond non-atom {}".format(other))
        for bond in self.bonds:
            if other in bond.atoms() and other is not self:
                bond.destroy()
                return
        raise ValueError("{} cannot unbond non-bonded {}".format(self, other))


    def bond_with(self, other):
        """Returns the :py:class:`.Bond` between this atom and another.

        :param Atom other: The atom to get the bond with.
        :raises TypeError: if something other than an :py:class:`Atom` is given.
        :rtype: ``Bond``"""

        if not isinstance(other, Atom):
            raise TypeError("Cannot get bond with non-atom {}".format(other))
        if other is self:
            return None
        for bond in self.bonds:
            if other in bond.atoms():
                return bond


    def copy(self):
        """Returns a copy of the atom. The new atom will have the same element,
        location, name, charge, ID and bfactor as the original, but will not be
        part of any model or other molecule, and will not have the any bonds.

        :rtype: ``Atom``"""

        return Atom(self._element, self._x, self._y, self._z, id=self._id,
         name=self._name, charge=self._charge, bfactor=self._bfactor)


    @atom_query
    def nearby_atoms(self, cutoff, *args, **kwargs):
        """Returns all atoms in the associated :py:class:`.Model` that are
        within a given distance (in the units of the atom coordinates) of this
        atom.

        :param cutoff: The distance cutoff to use.
        :param int id: if given, only atoms whose ID matches this will be\
        returned.
        :param str name: if given, only atoms whose name matches this will be\
        returned.
        :param str element: if given, only atoms whose element matches this\
        will be returned.
        :param bool hydrogen: If ``False``, hydrogen atoms will be excluded.
        :param bool het: If ``False``, non-chain atoms will be excluded.
        :param bool metal: If ``False``, metal atoms will be excluded.
        :rtype: ``set``"""

        if self._model:
            atoms =  self._model.atoms_in_sphere(
             *self.location, cutoff, *args, **kwargs
            )
            try:
                atoms.remove(self)
            except: pass
            return atoms
        return set()



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
        return "<{}-{} Bond>".format(atom1.element, atom2.element)


    @atom_query
    def atoms(self):
        """Returns the two :py:class:`.Atom` objects that the bond connects.
        They are given as an unordered set.

        :param int id: if given, only atoms whose ID matches this will be\
        returned.
        :param str name: if given, only atoms whose name matches this will be\
        returned.
        :param str element: if given, only atoms whose element matches this\
        will be returned.
        :param bool hydrogen: If ``False``, hydrogen atoms will be excluded.
        :param bool het: If ``False``, non-chain atoms will be excluded.
        :param bool metal: If ``False``, metal atoms will be excluded.
        :rtype: ``set``"""

        return set(self._atoms)


    @property
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

METALS = [
 "LI", "BE", "NA", "MG", "AL", "K", "CA", "SC", "TI", "V", "CR", "MN", "FE",
 "CO", "NI", "CU", "ZN", "HA", "RB", "SR", "Y", "ZR", "NB", "MO", "TC", "RU",
 "RH", "PD", "AG", "CD", "IN", "SN", "CS", "BA", "LA", "CE", "PR", "ND", "PM",
 "SM", "EU", "GD", "TB", "DY", "HO", "ER", "TM", "YB", "LU", "HF", "TA", "W",
 "RE", "OS", "IR", "PT", "AU", "HG", "TL", "PB", "BI", "PO", "FR", "RA", "AC",
 "TH", "PA", "U", "NP", "PU", "AM", "CM", "BK", "CF", "ES", "FM", "MD", "NO",
 "LR", "RF", "DB", "SG", "BH", "HS", "MT", "DS", "RG", "CN", "UUT", "FL", "LV"
]
