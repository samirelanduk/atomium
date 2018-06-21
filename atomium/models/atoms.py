"""Contains the atom class."""

import re
import numpy as np
from .data import PERIODIC_TABLE, METALS

class Atom:
    """An atom in space - a point particle with a location, element, charge etc.

    Atoms are the building blocks of all structures in atomium.

    :param str element: The atom's elemental symbol.
    :param number x: The atom's x coordinate.
    :param number y: The atom's y coordinate.
    :param number z: The atom's z coordinate.
    :param int id: An integer ID for the atom.
    :param str name: The atom's name.
    :param number charge: The charge of the atom.
    :param number bfactor: The B-factor of the atom (its uncertainty).
    :param list anisotropy: The directional uncertainty of the atom."""

    def __init__(self, element, x=0, y=0, z=0, id=0, name=None, charge=0,
                 bfactor=0, anisotropy=(0, 0, 0, 0, 0, 0)):
        self._element = str(element)
        self._x, self._y, self._z = x, y, z
        self._id, self._name = int(id), str(name) if name else None
        self._charge = float(charge)
        self._bfactor = float(bfactor)
        self._anisotropy = list(anisotropy)
        self._bonded_atoms = set()
        for attr in ("ligand", "residue", "chain", "model"):
            self.__dict__["_" + attr] = None


    def __repr__(self):
        return "<{}{} Atom{} at ({}, {}, {})>".format(
         self._element,
         " ({})".format(self._name) if self._name else "",
         " {}".format(self._id) if self._id else "",
         self._x, self._y, self._z
        )


    def __str__(self):
        return "<Atom{} ({})>".format(
         " {}".format(self._id) if self._id else "",
         self._name if self._name else self._element
        )


    @property
    def element(self):
        """The atom's element symbol. This is used to calculate its mass using a
        Periodic Table.

        :rtype: ``str``"""

        return self._element


    @element.setter
    def element(self, element):
        self._element = element


    @property
    def x(self):
        """The atom's x-coordinate.

        :rtype: ``float``"""

        return self._x


    @x.setter
    def x(self, x):
        self._x = x


    @property
    def y(self):
        """The atom's y-coordinate.

        :rtype: ``float``"""

        return self._y


    @y.setter
    def y(self, y):
        self._y = y


    @property
    def z(self):
        """The atom's z-coordinate.

        :rtype: ``float``"""

        return self._z


    @z.setter
    def z(self, z):
        self._z = z


    @property
    def id(self):
        """The atom's unique integer ID. It cannot be updated - the ID the atom
        is created with is its ID forever.

        :rtype: ``int``"""

        return self._id


    @property
    def name(self):
        """The atom's name. This is often used to determine what 'kind' of atom
        it is.

        :rtype: ``str``"""

        return self._name


    @name.setter
    def name(self, name):
        self._name = str(name)


    @property
    def charge(self):
        """The atom's charge - usually just zero, or 'neutral'.

        :rtype: ``float``"""

        return self._charge


    @charge.setter
    def charge(self, charge):
        self._charge = float(charge)


    @property
    def bfactor(self):
        """The atom's B-factor - the uncertainty in its position in all
        directions.

        :rtype: ``float``"""

        return self._bfactor


    @bfactor.setter
    def bfactor(self, bfactor):
        self._bfactor = float(bfactor)


    @property
    def anisotropy(self):
        """The atom's directional uncertainty, represented by a list of six
        numbers.

        :rtype: ``list``"""

        return self._anisotropy


    @property
    def location(self):
        """The atom's Cartesian coordinates in ``(x, y, z)`` format.

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
        """Translates an atom in 3D space. You can provide three values, or a
        single vector.

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
        :param number z: The atom's new z coordinate."""

        self._x, self._y, self._z = x, y, z


    def transform(self, matrix, trim=12):
        """Transforms the atom using a 3x3 matrix supplied. This is useful if
        the :py:meth:`.rotate` method isn't powerful enough for your needs.

        :param array matrix: A NumPy matrix representing the transformation.\
        You can supply a list of lists if you like and it will be converted to\
        a NumPy matrix.
        :param int trim: The amount of rounding to do to the atom's coordinates\
        after transforming - the default is 12 decimal places but this can be\
        set to ``None`` if no rounding is to be done."""

        vector = np.array(matrix).dot(self.location)
        self._x, self._y, self._z = vector
        self.trim(trim)


    def rotate(self, angle, axis, *args, **kwargs):
        """Rotates the atom by an angle in radians, around one of the the three
        axes.

        :param float angle: The angle to rotate by in radians.
        :param str axis: the axis to rotate around.
        :param int trim: The amount of rounding to do to the atom's coordinates\
        after rotating - the default is 12 decimal places but this can be\
        set to ``None`` if no rounding is to be done."""

        try:
            axis = [1 if i == "xyz".index(axis) else 0 for i in range(3)]
        except ValueError:
            raise ValueError("'{}' is not a valid axis".format(axis))
        axis = np.asarray(axis)
        axis = axis / np.sqrt(np.dot(axis, axis))
        a = np.cos(angle / 2)
        b, c, d = -axis * np.sin(angle / 2)
        aa, bb, cc, dd = a * a, b * b, c * c, d * d
        bc, ad, ac, ab, bd, cd = b * c, a * d, a * c, a * b, b * d, c * d
        self.transform(np.array([
         [aa + bb - cc - dd, 2 * (bc + ad), 2 * (bd - ac)],
         [2 * (bc - ad), aa + cc - bb - dd, 2 * (cd + ab)],
         [2 * (bd + ac), 2 * (cd - ab), aa + dd - bb - cc]
        ]), *args, **kwargs)


    @property
    def mass(self):
        """The atom's molar mass according to the Periodic Table, based on the
        atom's :py:meth:`element`. If the element doesn't match any symbol on
        the Periodic Table, a mass of 0 will be returned.

        The element lookup is case-insensitive.

        :rtype: ``float``"""

        return PERIODIC_TABLE.get(self._element.upper(), 0)


    @property
    def is_metal(self):
        """Checks whether the atom's element matches a metal element.

        The element lookup is case-insensitive.

        :rtype: ``bool``"""

        return self._element.upper() in METALS


    def distance_to(self, other):
        """Returns the distance (in whatever units the coordinates are defined
        in) between this atom and another. You can also give a (x, y, z) tuple
        instead of another atom if you so wish.

        :param Atom other: The other atom (or location tuple).
        :rtype: ``float``"""

        try:
            x, y, z = other
        except:
            x, y, z = other.location
        x_sum = pow((x - self._x), 2)
        y_sum = pow((y - self._y), 2)
        z_sum = pow((z - self._z), 2)
        return np.sqrt(x_sum + y_sum + z_sum)


    @property
    def ligand(self):
        """Returns the :py:class:`.Ligand` the atom is part of, or ``None`` if
        it is not part of one.

        :rtype: ``Ligand``"""

        return self._ligand


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
    def model(self):
        """Returns the :py:class:`.Model` the atom is part of, or ``None`` if
        it is not part of one.

        :rtype: ``Model``"""

        return self._model


    @property
    def bonded_atoms(self):
        """Returns the atoms bonded to this one.

        :rtype: ``set``"""

        return set(self._bonded_atoms)


    def bond_to(self, other):
        """Bonds the atom to some other atom. The two atoms will be placed
        inside each other's :py:meth:`.bonded_atoms`.

        :param Atom other: The atom to bond to."""

        self._bonded_atoms.add(other)
        other._bonded_atoms.add(self)


    def unbond_from(self, other):
        """Unbonds the atom from some other atom. If they aren't bonded to
        begin with, nothing bad will happen.

        :param Atom other: The atom to unbond from."""

        try:
            self._bonded_atoms.remove(other)
        except KeyError: pass
        try:
            other._bonded_atoms.remove(self)
        except KeyError: pass


    def nearby_atoms(self, cutoff, *args, **kwargs):
        """Returns all atoms in the associated :py:class:`.Model` that are
        within a given distance (in the units of the atom coordinates) of this
        atom. If the atom is not part of a model, no atoms will be returned.

        {}

        :param float cutoff: The radius to search within.
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


    def nearby_residues(self, *args, ligands=False, **kwargs):
        """Returns all residues in the associated :py:class:`.Model` that are
        within a given distance (in the units of the atom coordinates) of this
        atom. If the atom is not part of a model, no residues will be returned.

        {}

        :param float cutoff: the distance cutoff to use.
        :param bool ligands: if ``True``, ligands will be returned too.
        :rtype: ``set``"""

        nearby_atoms = self.nearby_atoms(*args, **kwargs)
        residues = set()
        for atom in nearby_atoms:
            residues.add(atom.residue)
            if ligands: residues.add(atom.ligand)
        try:
            residues.remove(None)
        except: pass
        return residues


    def copy(self):
        """Returns a copy of the atom. The new atom will have the same element,
        location, name, charge, ID, bfactor etc. as the original, but will not
        be part of any model or other molecule, and will not have any bonds.

        :rtype: ``Atom``"""

        return Atom(
         element=self._element, x=self._x, y=self._y, z=self._z, id=self._id,
         name=self._name, charge=self._charge, bfactor=self._bfactor,
         anisotropy=self.anisotropy
        )



QUERY_DOCSTRING = """You can specify which atoms should be searched in this
        function. Any atom property can be specified such as ``name='CA'``.
        String properties can be searched by regex, as in ``element='[^C]'``.
        Numeric properties can be searched by threshold, as in
        ``mass__gt=20``."""


def atom_query(func):
    """Decorator which can be applied to any function which returns atoms. It
    lets you query the output.

    The new function looks for keyword arguments which match atom attributes, or
    which are atom attributes with ``_regex`` or ``__`` in them. It then gets
    the atoms that the unmodified function would return, with the remaining
    arguments, and then filters them with each specification in the query.

    It will also update the new function's docstring with the explanatory text
    above.

    :param function func: The function to enhance.
    :rtype: ``function``"""

    def new(*args, **kwargs):
        query = {}
        for k, v in kwargs.items():
            if k in Atom.__dict__ or k.endswith("_regex") or "__" in k:
                query[k] = v
        for key in query: del kwargs[key]
        atoms = func(*args, **kwargs)
        if atoms:
            atom_attributes = list(list(atoms)[0].__dict__.keys())
            for k, v in query.items():
                attr = k.split("__")[0].split("_regex")[0]
                if "_" + attr in atom_attributes: attr = "_" + attr
                if k.endswith("_regex"):
                    atom_filter = lambda a: re.match(v, a.__dict__[attr])
                else:
                    comp = "__eq__"
                    if "__" in k:
                        comp = "__{}__".format(k.split("__")[1])
                    def atom_filter(a):
                        return a.__getattribute__(attr).__getattribute__(comp)(v)
                atoms = [a for a in atoms if atom_filter(a)]
            return set(atoms)
        return set()
    new.__name__ = func.__name__
    new.__doc__ = func.__doc__.format(QUERY_DOCSTRING)
    return new

Atom.nearby_atoms.__doc__ = Atom.nearby_atoms.__doc__.format(QUERY_DOCSTRING)
Atom.nearby_residues.__doc__ = Atom.nearby_residues.__doc__.format(QUERY_DOCSTRING)
