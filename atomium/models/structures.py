"""Contains the base Atom Structure class."""

import re
import numpy as np
import rmsd
from collections import Counter
from .atoms import atom_query, QUERY_DOCSTRING

class AtomStructure:
    """A structure made of atoms. In practice this class usually acts as a base
    class to things with more recognisable names, but it can also be used as a
    generic container of :py:class:`.Atom` objects, or as a base class for some
    user defined custom class.

    All atom structures are containers of their atoms, and of any
    sub-structures they might contain.

    :param \*atoms: The atoms the structure is to be made of. The atoms will\
    be updated with awareness of the new structure they are part of if a\
    sub-class is used. You can also pass in other atom structures here, and\
    their atoms will be used.
    :param str id: The structure's ID.
    :param str name: The structure's name."""

    CLASS_NAMES = ("ligand", "residue", "chain", "model")

    def __init__(self, *atoms, id=None, name=None):
        self._atoms = set()
        for atom in atoms:
            try:
                self._atoms.update(atom._atoms)
            except AttributeError: self._atoms.add(atom)
        self._id_atoms = {id: set() for id in set([a.id for a in self._atoms])}
        update = self.__class__.__name__.lower() in self.CLASS_NAMES
        for atom in self._atoms:
            self._id_atoms[atom.id].add(atom)
            if update:
                atom.__dict__["_" + self.__class__.__name__.lower()] = self
        self._id = str(id) if id else None
        self._name = str(name) if name else None


    def __repr__(self):
        return "<{}{} ({}{} atom{})>".format(
         self.__class__.__name__,
         " {}".format(self._name) if self._name else "",
         "{}, ".format(self._id) if self._id else "",
         len(self._atoms),
         "s" if len(self._atoms) != 1 else ""
        )


    def __contains__(self, member):
        try:
            atoms = member._atoms
        except AttributeError:
            atoms = {member}
        return atoms.issubset(self._atoms)


    @atom_query
    def atoms(self):
        """Returns the :py:class:`.Atom` objects in the structure.

        {}

        :rtype: ``set``"""

        return set(self._atoms)


    def atom(self, *args, **kwargs):
        """Returns the first :py:class:`.Atom` that matches the criteria given.

        Note that atoms are stored unordered in a ``set`` so if more than one
        atom matches the criteria you give, it might not be the same atom that
        is returned each time you call this method.

        {}

        :rtype: ``Atom``"""

        if "id" in kwargs:
            atoms = self._id_atoms.get(kwargs["id"])
        elif len(args) == 1:
            atoms = self._id_atoms.get(args[0])
        else:
            atoms = self.atoms(*args, **kwargs)
        if not atoms: atoms = set()
        for atom in atoms: return atom


    def pairwise_atoms(self, *args, **kwargs):
        """A generator which yeilds all the pairwise atom combinations of the
        structure. There will be no duplicates in the returned generator, and
        the number of returned pairs will be a triangle number.

        :rtype: ``tuple``"""

        atoms = list(self.atoms(*args, **kwargs))
        for a_index in range(len(atoms) - 1):
            for o_index in range(a_index + 1, len(atoms)):
                yield {atoms[a_index], atoms[o_index]}


    def add(self, obj):
        """Adds an atom or other :py:class:`.AtomStructure` to this structure.

        :param obj: The atom or structure to add."""

        try:
            atoms = obj._atoms
        except AttributeError: atoms = {obj}
        for atom in atoms:
            if atom.id in self._id_atoms:
                self._id_atoms[atom.id].add(atom)
            else:
                self._id_atoms[atom.id] = {atom}
            class_name = self.__class__.__name__.lower()
            if class_name in self.CLASS_NAMES:
                atom.__dict__["_" + class_name] = self
            self._atoms.add(atom)


    def remove(self, obj):
        """Removes an atom or other :py:class:`.AtomStructure` from this
        structure.

        :param obj: The atom or structure to remove."""

        try:
            atoms = obj._atoms
        except AttributeError: atoms = {obj}
        for atom in atoms:
            try:
                self._id_atoms[atom.id].remove(atom)
                if not self._id_atoms[atom.id]: del self._id_atoms[atom.id]
                atom.__dict__["_" + self.__class__.__name__.lower()] = None
                self._atoms.remove(atom)
            except KeyError: pass


    def _get(self, object_name, id=None, name=None, id_regex=None,
             name_regex=None, **kwargs):
        objects = set()
        for atom in self._atoms:
            objects.add(atom.__dict__["_" + object_name])
        try:
            objects.remove(None)
        except: pass
        if id:
            objects = [o for o in objects if o._id == id]
        if id_regex:
            objects = [o for o in objects if re.match(id_regex, o._id)]
        if name:
            objects = [o for o in objects if o._name == name]
        if name_regex:
            objects = [o for o in objects if re.match(name_regex, o._name)]
        if object_name == "ligand":
            if "water" in kwargs and not kwargs["water"]:
                objects = [o for o in objects if o._name not in ("HOH", "WAT")]
        return set(objects)


    @property
    def id(self):
        """The structure's identifier. Once created it is not modifiable.

        :rtype: ``str``"""

        return self._id


    @property
    def name(self):
        """The structure's name.

        :rtype: ``str``"""

        return self._name


    @name.setter
    def name(self, name):
        self._name = name


    def trim(self, places):
        """Rounds the coordinate values to a given number of decimal places.
        Useful for removing floating point rounding errors after transformation.

        :param int places: The number of places to round the coordinates to. If\
        ``None``, no rounding will be done."""

        for atom in self._atoms:
            atom.trim(places)


    def translate(self, *args, **kwargs):
        """Translates the structure through space, updating all atom
        coordinates accordingly. You can provide three values, or a single
        vector.

        :param Number dx: The distance to move in the x direction.
        :param Number dy: The distance to move in the y direction.
        :param Number dz: The distance to move in the z direction.
        :param int trim: The amount of rounding to do to the atoms' coordinates\
        after translating - the default is 12 decimal places but this can be\
        set to ``None`` if no rounding is to be done."""

        for atom in self._atoms:
            atom.translate(*args, **kwargs)


    def transform(self, *args, **kwargs):
        """Transforms the structure using a 3x3 matrix supplied. This is useful
        if the :py:meth:`.rotate` method isn't powerful enough for your needs.

        :param array matrix: A NumPy matrix representing the transformation.\
        You can supply a list of lists if you like and it will be converted to\
        a NumPy matrix.
        :param int trim: The amount of rounding to do to the atoms' coordinates\
        after transforming - the default is 12 decimal places but this can be\
        set to ``None`` if no rounding is to be done."""

        for atom in self._atoms:
            atom.transform(*args, **kwargs)


    def rotate(self, *args, **kwargs):
        """Rotates the structure about an axis, updating all atom coordinates
        accordingly.

        :param Number angle: The angle in radians.
        :param str axis: The axis to rotate around. Can only be 'x', 'y' or 'z'.
        :param int trim: The amount of rounding to do to the atoms' coordinates\
        after translating - the default is 12 decimal places but this can be\
        set to ``None`` if no rounding is to be done."""

        for atom in self._atoms:
            atom.rotate(*args, **kwargs)


    @property
    def mass(self):
        """The mass of the structure in Daltons - just the sum of its atoms'
        masses.

        :rtype: ``float``"""

        return round(sum([atom.mass for atom in self._atoms]), 12)


    @property
    def charge(self):
        """Returns the charge of the structure, based on the charges of its
        atoms.

        :rtype: ``float``"""

        return round(sum([atom.charge for atom in self._atoms]), 12)


    @property
    def formula(self):
        """Returns the formula (count of each atom) of the structure.

        :rtype: ``Counter``"""

        return Counter([atom.element for atom in self._atoms])


    @property
    def center_of_mass(self):
        """Returns the center of mass of the structure. This is the average of
        all the atom coordinates, weighted by the mass of each atom.

        :returns: (x, y, z) ``tuple``"""

        mass = self.mass
        average_x = sum([atom.x * atom.mass for atom in self._atoms]) / mass
        average_y = sum([atom.y * atom.mass for atom in self._atoms]) / mass
        average_z = sum([atom.z * atom.mass for atom in self._atoms]) / mass
        return (average_x, average_y, average_z)


    @property
    def radius_of_gyration(self):
        """The radius of gyration of a structure is a measure of how extended it
        is. It is the root mean square deviation of the atoms' distance from the
        structure's :py:meth:`.center_of_mass`.

        :rtype: ``float``"""

        center_of_mass = self.center_of_mass
        atoms = self._atoms
        square_deviation = sum(
         [atom.distance_to(center_of_mass) ** 2 for atom in atoms]
        )
        mean_square_deviation = square_deviation / len(atoms)
        return np.sqrt(mean_square_deviation)


    def pairing_with(self, structure):
        """Takes another structure with the same number of atoms as this one,
        and attempts to find the nearest equivalent of every atom in this
        structure, in that structure.

        Atoms will be aligned first by element, then by name, then by number of
        bonds, then IDs, and finally by memory address - this last metric is
        used to ensure that even when allocation is essentially random, it is at
        least the same every time two structures are aligned.

        :param AtomStructure structure: the structure to pair with.
        :raises ValueError: if the other structure has a different number of\
        atoms.
        :rtype: ``dict``"""

        atoms, other_atoms = list(self._atoms), list(structure._atoms)
        if len(atoms) != len(other_atoms):
            raise ValueError("{} and {} have different numbers of atoms".format(
             self, structure
            ))
        for l in atoms, other_atoms:
            l.sort(key=lambda a: (
             a.element, a.name, len(a.bonded_atoms), a.id, id(a)
            ))
        return {a1: a2 for a1, a2 in zip(atoms, other_atoms)}


    def superimpose_onto(self, other):
        """Superimoses this structure onto another - it will be translated so
        that its center of mass matches the other structure's, then rotated so
        as to minimise the RMSD.

        The other structure must have the same number of atoms.

        :param AtomStructure other: The structure to superimpose onto. This\
        structure does not move."""

        center, other_center = self.center_of_mass, other.center_of_mass
        self.translate(-center[0], -center[1], -center[2])
        atoms, other_atoms = zip(*self.pairing_with(other).items())
        P = np.array([[*atom.location] for atom in atoms])
        Q = [(
         a.x - other_center[0], a.y - other_center[1], a.z - other_center[2]
        ) for a in other_atoms]
        P = rmsd.kabsch_rotate(P, Q)
        for atom, row in zip(atoms, P):
            atom.move_to(*row)
        self.translate(other_center[0], other_center[1], other_center[2])


    def rmsd_with(self, structure, superimpose=False):
        """Calculates the Root Mean Square Deviation between this structure and
        another.

        You can get the RMSD either of the coordinates as they are, or of
        superimposed coordinates.

        :param AtomStructure structure: the structure to check against.
        :param bool superimpose: if ``True``, the structure will be\
        superimosed first (and then moved back).
        :raises TypeError: if the other structure is not an\
        :py:class:`.AtomStructure`.
        :raises ValueError: if the other structure has a different number of\
        atoms.
        :rtype: ``float``"""

        pairing = self.pairing_with(structure)
        if superimpose:
            atoms = list(self._atoms)
            locations = [atom.location for atom in atoms]
            self.superimpose_onto(structure)
        sd = sum(a1.distance_to(a2) ** 2 for a1, a2 in pairing.items())
        if superimpose:
            for atom, location in zip(atoms, locations):
                atom.move_to(*location)
        msd = sd / len(pairing)
        return np.sqrt(msd)


    def grid(self, size=1, margin=0):
        """A generator which models a grid around the structure and returns the
        coordinates of all the points in that grid. The origin is always one of
        those points, and the grid will be a box.

        :param int size: The spacing between grid points. The default is 1.
        :param int margin: How far to extend the grid beyond the structure\
        coordinates. The default is 0.
        :rtype: ``tuple``"""

        atom_locations = [atom.location for atom in self._atoms]
        dimension_values = []
        for dimension in range(3):
            coordinates = [loc[dimension] for loc in atom_locations]
            min_, max_ = min(coordinates) - margin, max(coordinates) + margin
            values = [0]
            while values[0] > min_: values.insert(0, values[0] - size)
            while values[-1] < max_: values.append(values[-1] + size)
            dimension_values.append(values)
        for x in dimension_values[0]:
            for y in dimension_values[1]:
                for z in dimension_values[2]:
                    yield (x, y, z)


    @atom_query
    def atoms_in_sphere(self, x, y, z, radius):
        """Returns all the atoms in a given sphere within the structure.

        {}

        :rtype: ``set``"""

        atoms = filter(
         lambda a: a.distance_to((x, y, z)) <= radius, self._atoms
        )
        return set(atoms)


    def nearby_atoms(self, *args, **kwargs):
        """Gets all atoms within a given cutoff of this structure's atoms.

        {}

        :param float cutoff: the distance cutoff to use.
        :rtype: ``set``"""

        atoms = set()
        for atom in self._atoms:
            atoms.update(atom.nearby_atoms(*args, **kwargs))
        for atom in self._atoms:
            try:
                atoms.remove(atom)
            except: pass
        return atoms


    def nearby_residues(self, *args, **kwargs):
        """Gets all atoms within a given cutoff of this structure's atoms.

        {}

        :param float cutoff: the distance cutoff to use.
        :param bool ligands: if ``True``, ligands will be returned too.
        :rtype: ``set``"""

        residues = set()
        for atom in self._atoms:
            residues.update(atom.nearby_residues(*args, **kwargs))
        for het in self._get("residue").union(self._get("ligand")):
            try:
                residues.remove(het)
            except: pass
        return residues


    def copy(self):
        """Returns a copy of the structure, with its own distinct atoms.
        Its atoms will have the same ID, location, element, charge, name and
        bfactor as their counterparts in the original, but will have no bonds,
        and be a member of no model - even if their counterparts do and are.

        :rtype: ``AtomStructure``"""

        copy = self.__class__(*[a.copy() for a in self.atoms()])
        copy._id, copy._name = self._id, self._name
        return copy


    def to_file_string(self, file_format, description=None):
        """Converts a structure to a filestring. Currently supported file
        formats are: .xyz and .pdb.

        :param str file_format: The file format to use, in lowercase.
        :param str description: A structure description to put in the file.
        :raises ValueError: if an unsupported file format is given."""

        if file_format == "xyz":
            from ..files.xyz2xyzdict import structure_to_xyz_dict
            from ..files.xyzdict2xyzstring import xyz_dict_to_xyz_string
            xyz_dict = structure_to_xyz_dict(self)
            xyz_dict["title"] = description
            return xyz_dict_to_xyz_string(xyz_dict)
        elif file_format == "pdb":
            from ..files.pdb2pdbdict import structure_to_pdb_dict
            from ..files.pdbdict2pdbstring import pdb_dict_to_pdb_string
            pdb_dict = structure_to_pdb_dict(self)
            pdb_dict["title"] = description
            return pdb_dict_to_pdb_string(pdb_dict)
        else:
            raise ValueError("{} is not a valid file type".format(file_format))


    def save(self, path, *args, **kwargs):
        """Saves the structure to file, in the format implied by the extension
        of the path you provide (i.e. giving a path ``/path/to/file.xyz`` will
        save as .xyz etc.).

        :param str path: The path to save to. The extension you provide here is\
        important as atomium will use that to determine what file format to\
        save as.
        :param str description: A structure description to put in the file."""

        file_format = path.split(".")[-1].lower()
        s = self.to_file_string(file_format, *args, **kwargs)
        from ..files.utilities import string_to_file
        string_to_file(s, path)


AtomStructure.atom.__doc__ =\
 AtomStructure.atom.__doc__.format(QUERY_DOCSTRING)
AtomStructure.nearby_atoms.__doc__ =\
 AtomStructure.nearby_atoms.__doc__.format(QUERY_DOCSTRING)
AtomStructure.nearby_residues.__doc__ =\
 AtomStructure.nearby_residues.__doc__.format(QUERY_DOCSTRING)
