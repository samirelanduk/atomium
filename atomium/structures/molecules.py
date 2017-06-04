"""Contains classes for structures made of atoms."""

from collections import Counter
from math import sqrt
from geometrica import translate, rotate
from .atoms import Atom

class AtomicStructure:
    """Represents structures made of :py:class:`.Atom` objects, which tends to
    be rather a lot of things in practice. This class would not generally be
    instantiated directly, and is here to be a parent class to other, more
    specific entities.

    AtomicStructures are containers of their atoms.

    :param \*atoms: The :py:class:`.Atom` objects that make up the structure.
    :raises TypeError: if non-atoms are given."""

    def __init__(self, *atoms):
        if not all(isinstance(atom, Atom) for atom in atoms):
            non_atoms = [atom for atom in atoms if not isinstance(atom, Atom)]
            raise TypeError(
             "AtomicStructures need atoms, not '{}'".format(non_atoms[0])
            )
        self._atoms = set(atoms)


    def __repr__(self):
        return "<{} ({} atoms)>".format(self.__class__.__name__, len(self._atoms))


    def __contains__(self, member):
        return member in self._atoms


    def atoms(self, element=None, atom_id=None):
        """Returns the :py:class:`.Atom` objects in the structure. You can
        filter these by element if you wish.

        :param str element: If given, only atoms whose element matches this\
        will be returned.
        :param int atom_id: If given, only atoms whose atom ID matches this\
        will be returned (this will only return one atom).
        :rtype: ``set``"""

        atoms = set(self._atoms)
        if element:
            atoms = set(filter(lambda a: a.element() == element, atoms))
        if atom_id:
            atoms = set(filter(lambda a: a.atom_id() == atom_id, atoms))
        return atoms


    def atom(self, *args, **kwargs):
        """Returns the first :py:class:`.Atom` that matches the criteria given.

        Note that atoms are stored unordered in a ``set`` so if more than one
        atom matches the criteria you give, it might not be the same atom that
        is returned each time you call this method.

        :param str element: If given, only atoms whose element matches this\
        will be searched.
        :param int atom_id: If given, only atoms whose atom ID matches this\
        will be searched.
        :rtype: ``Atom``"""

        atoms = self.atoms(*args, **kwargs)
        for atom in atoms: return atom


    def add_atom(self, atom):
        """Adds an :py:class:`.Atom` to the structure.

        :param Atom atom: The atom to add.
        :raises TypeError: if the atom given is not an Atom."""

        if not isinstance(atom, Atom):
            raise TypeError("Can only add atoms, not '{}'".format(atom))
        self._atoms.add(atom)


    def remove_atom(self, atom):
        """Removes an :py:class:`.Atom` from the structure.

        :param Atom atom: The atom to remove."""

        self._atoms.remove(atom)


    def mass(self):
        """Returns the mass of the structure in Daltons, based on the masses of
        its atoms.

        :rtype: ``float``"""

        return sum([atom.mass() for atom in self._atoms])


    def formula(self):
        """Returns the formula (count of each atom) of the structure.

        :rtype: ``Counter``"""

        return Counter([atom.element() for atom in self._atoms])


    def translate(self, dx, dy, dz):
        """Translates the structure through space, updating all atom
        coordinates accordingly.

        :param Number dx: The distance to move in the x direction.
        :param Number dy: The distance to move in the y direction.
        :param Number dz: The distance to move in the z direction."""

        atoms = list(self._atoms)
        points = translate(atoms, dx, dy, dz)
        for index, atom in enumerate(atoms):
            atom._x, atom._y, atom._z = points[index]


    def rotate(self, axis, angle):
        """Rotates the structure about an axis, updating all atom coordinates
        accordingly.

        :param str axis: The axis to rotate around. Can only be 'x', 'y' or 'z'.
        :param Number angle: The angle in degrees. Rotation is right handed."""

        atoms = list(self._atoms)
        points = rotate(atoms, axis, angle)
        for index, atom in enumerate(atoms):
            atom._x, atom._y, atom._z = points[index]


    def center_of_mass(self):
        """Returns the center of mass of the structure. This is the average of
        all the atom coordinates, weighted the mass of each atom.

        :returns: (x, y, z) ``tuple``"""

        mass = self.mass()
        average_x = sum([atom._x * atom.mass() for atom in self._atoms]) / mass
        average_y = sum([atom._y * atom.mass() for atom in self._atoms]) / mass
        average_z = sum([atom._z * atom.mass() for atom in self._atoms]) / mass
        return (average_x, average_y, average_z)


    def radius_of_gyration(self):
        """The radius of gyration of a structure is a measure of how extended it
        is. It is the root mean square deviation of the atoms' distance from the
        structure's :py:meth:`.center_of_mass`.

        :rtype: ``float``"""

        center_of_mass = self.center_of_mass()
        square_deviation = sum(
         [atom.distance_to(center_of_mass) ** 2 for atom in self._atoms]
        )
        mean_square_deviation = square_deviation / len(self._atoms)
        return sqrt(mean_square_deviation)


    def to_file_string(self, file_format, description=""):
        """Converts a structure to a filestring. Currently supported file formats
        are: .xyz.

        :param str file_format: The file format to use, in lowercase.
        :param str description: A structure description to put in the file.
        :raises ValueError: if an unsopported file format is given."""

        if file_format == "xyz":
            from ..converters.structure2xyzstring import structure_to_xyz_string
            return structure_to_xyz_string(self, description)
        else:
            raise ValueError("{} is not a valid file type".format(file_format))


    def save(self, path, *args, **kwargs):
        """Saves the structure to file, in the format implied by the extension
        of the path you provide (i.e. giving a path ``/path/to/file.xyz`` will
        save as .xyz).

        :param str path: The path to save to. The extension you provide here is\
        important as atomium will use that to determine what file format to\
        save as.
        :param str description: A structure description to put in the file."""

        file_format = path.split(".")[-1].lower()
        s = self.to_file_string(file_format, *args, **kwargs)
        from ..converters.strings import string_to_file
        string_to_file(s, path)
