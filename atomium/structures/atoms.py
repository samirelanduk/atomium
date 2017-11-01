"""This module contains the classes for atoms and their bonds."""

import weakref
from math import sqrt, acos, degrees
from points import Vector

class Atom:
    """Represents an atom in three dimensional space. Every atom has an element
    and a set of Cartesian coordinates.

    :param str element: The atom's element symbol. It doesn't need to be on the\
    Periodic Table, but it does need to be 1 or 2 characters.
    :param number x: The atom's x coordinate.
    :param number y: The atom's y coordinate.
    :param number z: The atom's z coordinate.
    :param int atom_id: A unique integer ID for the atom. This is supposed to\
    be unique.
    :param str name: The atom's name.
    :param number charge: The charge of the atom.
    :raises TypeError: if the element is not str.
    :raises ValueError: if the element is not 1 or 2 characters.
    :raises TypeError: if the coordinates are not numeric.
    :raises TypeError: if the atom_id is not int.
    :raises TypeError: if the name is not str.
    :raises TypeError: if the charge is not numeric.
    :raises TypeError: if the bfactor is not numeric."""

    def __init__(self, element, x, y, z, atom_id=None, name=None, charge=0,
                 bfactor=0):
        if not isinstance(element, str):
            raise TypeError("Element '{}' is not str".format(element))
        if not 0 < len(element) < 3:
            raise ValueError("Element {} is not 1 or 2 chars".format(element))
        if not isinstance(x, (float, int)):
            raise TypeError("x coordinate '{}' is not numeric".format(x))
        if not isinstance(y, (float, int)):
            raise TypeError("y coordinate '{}' is not numeric".format(y))
        if not isinstance(z, (float, int)):
            raise TypeError("z coordinate '{}' is not numeric".format(z))
        if atom_id is not None and not isinstance(atom_id, int):
            raise TypeError("ID {} is not an integer".format(atom_id))
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
        self._id = atom_id
        self._name = name
        self._charge = charge
        self._bfactor = bfactor
        self._bonds = set()
        self._residue, self._chain, self._molecule = None, None, None
        self._model = None


    def __repr__(self):
        return "<{} Atom {}{}at ({}, {}, {})>".format(
         self._element, self._id if self._id else "", " " if self._id else "",
         self._x, self._y, self._z
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
        coordinate will be updated, but it must be numeric.

        :param number x: If given, the atom's x coordinate will be set to this.
        :raises TypeError: if the x coordinate given is not numeric.
        :rtype: ``int`` or ``float``"""

        if x is None:
            return self._x
        else:
            if not isinstance(x, (float, int)):
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
            if not isinstance(y, (float, int)):
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
            if not isinstance(z, (float, int)):
                raise TypeError("z coordinate '{}' is not numeric".format(z))
            self._z = z


    def location(self):
        """Returns the atom's Cartesian coordinates.

        :rtype: ``tuple``"""

        return (self._x, self._y, self._z)


    def translate(self, dx, dy, dz):
        """Translates an atom in 3D space.

        :param float dx: The distance to move in the x direction.
        :param float dy: The distance to move in the y direction.
        :param float dz: The distance to move in the z direction."""

        self._x += dx
        self._y += dy
        self._z += dz


    def rotate(self, angle, axis):
        """Rotates an atom in 3D space.

        :param float angle: Angle in radians.
        :param str axis: The axis to rotate around."""
        
        vector = Vector(self._x, self._y, self._z)
        vector.rotate(angle, axis)
        self._x, self._y, self._z = vector.values()


    def atom_id(self):
        """Returns the atom's unique integer ID.

        :rtype: ``int``"""

        return self._id


    def name(self, name=None):
        """Returns the atom's name. If a value is given, the name will be
        updated, provided it is a string.

        :param str name: If given, the name will be set to this.
        :raises TypeError: if the name given is not str."""

        if name is None:
            return self._name
        else:
            if not isinstance(name, str):
                raise TypeError("Atom name '{}' is not str".format(name))
            self._name = name


    def residue(self):
        """Returns the :py:class:`.Residue` the atom is part of, or ``None`` if
        it is not part of one.

        :rtype: ``Residue``"""

        return self._residue


    def chain(self):
        """Returns the :py:class:`.Chain` the atom is part of, or ``None`` if
        it is not part of one.

        :rtype: ``Chain``"""

        return self._chain


    def molecule(self):
        """Returns the :py:class:`.Molecule` the atom is part of, or ``None`` if
        it is not part of one.

        :rtype: ``Molecule``"""

        return self._molecule


    def model(self):
        """Returns the :py:class:`.Model` the atom is part of, or ``None`` if
        it is not part of one.

        :rtype: ``Model``"""

        return self._model


    def charge(self, charge=None):
        """Returns the atom's charge. If a value is given, the charge will be
        updated, but it must be numeric.

        :param number charge: If given, the atom's charge will be set to this.
        :raises TypeError: if the charge given is not numeric.
        :rtype: ``int`` or ``float``"""

        if charge is None:
            return self._charge
        else:
            if not isinstance(charge, (float, int)):
                raise TypeError("charge '{}' is not numeric".format(charge))
            self._charge = charge


    def bfactor(self, bfactor=None):
        """Returns the atom's B-factor - the uncertainty in its position. If a
        value is given, the bfactor will be updated, but it must be numeric.

        :param number bfactor: If given, the atom's bfactor will be set to this.
        :raises TypeError: if the bfactor given is not numeric.
        :rtype: ``int`` or ``float``"""

        if bfactor is None:
            return self._bfactor
        else:
            if not isinstance(bfactor, (float, int)):
                raise TypeError("bfactor '{}' is not numeric".format(bfactor))
            self._bfactor = bfactor


    def bonds(self):
        """Returns the :py:class:`.Bond` objects that the atom is associated
        with.

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

        if other not in self.bonded_atoms():
            Bond(self, other)


    def unbond(self, other):
        """Breaks the bond between this atom and another.

        :param Atom other: The atom to unbond from.
        :raises TypeError: if something other than an :py:class:`Atom` is given.
        :raises ValueError: if the atom given isn't bonded to begin with."""

        if not isinstance(other, Atom):
            raise TypeError("Cannot unbond non-atom {}".format(other))
        for bond in self.bonds():
            if other in bond.atoms() and other is not self:
                bond.destroy()
                return
        raise ValueError("{} cannot unbond non-bonded {}".format(self, other))


    def bond_with(self, other):
        """Returns the :py:class:`.Bond` between this atom and another.

        :param Atom other: The atom to get the bond with.
        :raises TypeError: if something other than an :py:class:`Atom` is given.
        :rtype: ``Bond`` (or ``None`` if no bond exists)."""

        if not isinstance(other, Atom):
            raise TypeError("Cannot get bond with non-atom {}".format(other))
        if other is self:
            return None
        for bond in self.bonds():
            if other in bond.atoms():
                return bond


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
            x, y, z = other.location()
        x_sum = pow((x - self._x), 2)
        y_sum = pow((y - self._y), 2)
        z_sum = pow((z - self._z), 2)
        return sqrt(x_sum + y_sum + z_sum)


    def nearby(self, cutoff, *args, **kwargs):
        """Returns all atoms in the associated :py:class:`.Model` that are
        within a given distance (in the units of the atom coordinates) of this
        atom.

        :param int cutoff: The distance cutoff to use.
        :param str element: If given, only atoms whose element matches this\
        will be returned.
        :param str exclude: If given, only atoms whose element doesn't match\
        this will be returned.
        :param int atom_id: If given, only atoms whose atom ID matches this\
        will be returned (this will only return one atom).
        :param str name: If given, only atoms whose name matches this will be\
        returned."""

        if self._model:
            atoms = self._model.atoms(*args, **kwargs)
            try:
                atoms.remove(self)
            except: pass
            return set(filter(lambda a: a.distance_to(self) <= cutoff, atoms))
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


    def vector(self, target):
        """Returns the ``Vector`` which represents the bond. You must specify
        which atom you want the Vector to point towards.

        The Vector comes from the ``points`` library.

        :param Atom target: The atom the Vector will point towards.
        :raises TypeError: if a non-atom is given.
        :raises ValueError: if the atom given isn't in the Bond.
        :rtype: ``Vector``"""

        if not isinstance(target, Atom):
            raise TypeError("bond atom {} is not an atom".format(target))
        if target not in self._atoms:
            raise ValueError("{} is not in bond {}".format(target, self))
        base = [atom for atom in self._atoms if atom is not target][0]
        values = [c1 - c2 for c1, c2 in zip(target.location(), base.location())]
        return Vector(*values)


    def angle_with(self, other, degrees=False):
        """Returns the angle between this Bond and another. If the Bonds share
        an atom, the atom will be used as the source of the two vectors. If
        they don't, the closest atoms will be used as the base.

        :param Bond other: The other Bond.
        :param bool degrees: If ``True``, the angle will be returned in degrees.
        :raises TypeError: if a non-Bond is given.
        :rtype: ``float``"""

        if not isinstance(other, Bond):
            raise TypeError("{} is not a Bond".format(other))
        pairs = [
         [set([a1, a2]), a1.distance_to(a2)]
        for a1 in self._atoms for a2 in other._atoms]
        pair = min(pairs, key=lambda k: k[1])[0]
        v1 = self.vector([a for a in self._atoms if a not in pair][0])
        v2 = other.vector([a for a in other._atoms if a not in pair][0])
        return v1.angle_with(v2, degrees=degrees)


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
