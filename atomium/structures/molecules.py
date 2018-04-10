"""This module contains classes for structures made of atoms."""

import math
from collections import Counter
from itertools import combinations
from functools import reduce
import operator
import numpy as np
import rmsd
from .atoms import Atom, atom_query

class AtomicStructure:
    """Represents structures made of :py:class:`.Atom` objects, which tends to
    be quite a lot of things in practice. This class would not generally be
    instantiated directly, and is here to be a parent class to other, more
    specific entities.

    AtomicStructures are containers of their atoms, and support the ``in``
    keyword.

    :param \*atoms: The :py:class:`.Atom` objects that make up the structure.\
    These can also be AtomicStructures themsevles, in which case the atoms of\
    that structure will be used in its place.
    :raises TypeError: if non-atoms or AtomicStructures are given."""

    def __init__(self, *atoms):
        self._atoms = set()
        for atom in atoms:
            if not isinstance(atom, (Atom, AtomicStructure)):
                raise TypeError(
                 "AtomicStructures need atoms, not '{}'".format(atom)
                )
        [self._atoms.add(atom) if isinstance(atom, Atom)
         else self._atoms.update(atom._atoms) for atom in atoms]
        self._id_atoms = {id: set() for id in set([a.id for a in self._atoms])}
        for atom in self._atoms:
            self._id_atoms[atom.id].add(atom)


    def __repr__(self):
        return "<{} ({} atoms)>".format(
         self.__class__.__name__, len(self._atoms)
        )


    def __contains__(self, member):
        return member in self._atoms


    @atom_query
    def atoms(self):
        """Returns the :py:class:`.Atom` objects in the structure. You can
        filter these by element if you wish.

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


    def atom(self, *args, **kwargs):
        """Returns the first :py:class:`.Atom` that matches the criteria given.

        Note that atoms are stored unordered in a ``set`` so if more than one
        atom matches the criteria you give, it might not be the same atom that
        is returned each time you call this method.

        :param int id: if given, only atoms whose ID matches this will be\
        returned.
        :param str name: if given, only atoms whose name matches this will be\
        returned.
        :param str element: if given, only atoms whose element matches this\
        will be returned.
        :param bool hydrogen: If ``False``, hydrogen atoms will be excluded.
        :param bool het: If ``False``, non-chain atoms will be excluded.
        :param bool metal: If ``False``, metal atoms will be excluded.
        :rtype: ``Atom``"""

        if "id" in kwargs:
            atoms = self._id_atoms.get(kwargs["id"])
            if not atoms: atoms = set()
        elif len(args) == 1:
            atoms = self._id_atoms.get(args[0])
        else:
            atoms = self.atoms(*args, **kwargs)
        if not atoms: atoms = set()
        for atom in atoms: return atom


    def add_atom(self, atom):
        """Adds an :py:class:`.Atom` to the structure.

        :param Atom atom: The atom to add.
        :raises TypeError: if the atom given is not an Atom."""

        if not isinstance(atom, Atom):
            raise TypeError("Can only add atoms, not '{}'".format(atom))
        if atom.id in self._id_atoms:
            self._id_atoms[atom.id].add(atom)
        else:
            self._id_atoms[atom.id] = {atom}
        atom.__dict__["_" + self.__class__.__name__.lower()] = self
        self._atoms.add(atom)


    def remove_atom(self, atom):
        """Removes an :py:class:`.Atom` from the structure.

        :param Atom atom: The atom to remove."""

        try:
            self._id_atoms[atom.id].remove(atom)
            if not self._id_atoms[atom.id]: del self._id_atoms[atom.id]
            atom.__dict__["_" + self.__class__.__name__.lower()] = None
            self._atoms.remove(atom)
        except KeyError: pass


    def pairwise_atoms(self, *args, **kwargs):
        """A generator which yeilds all the pairwise atom combinations of the
        structure. There will be no duplicates in the returned generator, and
        the number of returned pairs will be a triangle number.

        :rtype: ``tuple``"""

        atoms = list(self.atoms(*args, **kwargs))
        for a_index in range(len(atoms) - 1):
            for o_index in range(a_index + 1, len(atoms)):
                yield {atoms[a_index], atoms[o_index]}


    def add(self, structure):
        """Adds an atomic structure to this one - its atoms will be incorporated
        and so will its residues/molecules/chains if it has any.

        :param AtomicStructure structure: The structure to add.
        :raises TypeError: if the structure is not an AtomicStructure."""

        if not isinstance(structure, AtomicStructure):
            raise TypeError("{} is not an atomic structure".format(structure))
        for atom in structure._atoms:
            self.add_atom(atom)


    def remove(self, structure):
        """Removes an atomic structure from this one.

        :param AtomicStructure structure: The structure to add.
        :raises TypeError: if the structure is not an AtomicStructure."""

        if not isinstance(structure, AtomicStructure):
            raise TypeError("{} is not an atomic structure".format(structure))
        for atom in structure._atoms:
            self.remove_atom(atom)


    def residues(self, id=None, name=None):
        """Returns all the :py:class:`.Residue` objects in the structure which
        match the given criteria.

        :param str id: Filter by residue ID.
        :param str name: Filter by name.
        :rtype: ``Residue``"""

        res = set()
        for atom in self._atoms: res.add(atom.residue)
        try: res.remove(None)
        except KeyError: pass
        if id: res = set(filter(lambda r: r.id == id, res))
        if name: res = set(filter(lambda r: r.name == name, res))
        return res


    def residue(self, *args, **kwargs):
        """Returns the first :py:class:`.Residue` object in the structure which
        matches the given criteria.

        :param str id: Filter by residue ID.
        :param str name: Filter by name.
        :rtype: ``Residue``"""

        residues = self.residues(*args, **kwargs)
        for res in residues: return res


    def molecules(self, id=None, name=None, generic=False, water=True):
        """Returns all the :py:class:`.Molecule` objects in the structure which
        match the given criteria.

        :param str id: Filter by molecule ID.
        :param str name: Filter by name.
        :param bool generic: if ``True``, chains will be excluded.
        :param bool water: if ``False``, water molecules will be excluded.
        :rtype: ``Molecule``"""

        molecules = set()
        for atom in self._atoms: molecules.add(atom.molecule)
        try: molecules.remove(None)
        except KeyError: pass
        if id: molecules = set(filter(lambda r: r.id == id, molecules))
        if name: molecules = set(filter(lambda r: r.name == name, molecules))
        if generic:
            from .chains import Chain
            molecules = set(filter(
             lambda m: not isinstance(m, (Residue, Chain)), molecules
            ))
        if not water:
            molecules = set(filter(
             lambda m: m.name not in ("HOH", "WAT"), molecules
            ))
        return molecules


    def molecule(self, *args, **kwargs):
        """Returns the first :py:class:`.Molecule` object in the structure which
        matches the given criteria.

        :param str id: Filter by molecule ID.
        :param str name: Filter by name.
        :param bool generic: if ``True``, chains will be excluded.
        :param bool water: if ``False``, water molecules will be excluded.
        :rtype: ``Molecule``"""

        molecules = self.molecules(*args, **kwargs)
        for mol in molecules: return mol


    def chains(self, id=None, name=None):
        """Returns all the :py:class:`.Chain` objects in the structure which
        match the given criteria.

        :param str id: Filter by chain ID.
        :param str name: Filter by name.
        :rtype: ``Chain``"""

        chains = set()
        for atom in self._atoms:
            chains.add(atom.chain)
        try:
            chains.remove(None)
        except KeyError: pass
        if id:
            chains = set(filter(lambda r: r.id == id, chains))
        if name:
            chains = set(filter(lambda r: r.name == name, chains))
        return chains


    def chain(self, *args, **kwargs):
        """Returns the first :py:class:`.Chain` object in the structure which
        matches the given criteria.

        :param str id: Filter by chain ID.
        :param str name: Filter by name.
        :rtype: ``Chain``"""

        chains = self.chains(*args, **kwargs)
        for chain in chains: return chain


    def complexes(self, id=None, name=None):
        """Returns all the :py:class:`.Complex` objects in the structure which
        match the given criteria.

        :param str id: Filter by complex ID.
        :param str name: Filter by name.
        :rtype: ``Complex``"""

        complexes = set()
        for atom in self._atoms:
            complexes.add(atom.complex)
        try:
            complexes.remove(None)
        except KeyError: pass
        if id:
            complexes = set(filter(lambda r: r.id == id, complexes))
        if name:
            complexes = set(filter(lambda r: r.name == name, complexes))
        return complexes


    def complex(self, *args, **kwargs):
        """Returns the first :py:class:`.Complex` object in the structure which
        matches the given criteria.

        :param str id: Filter by complex ID.
        :param str name: Filter by name.
        :rtype: ``Complex``"""

        complexes = self.complexes(*args, **kwargs)
        for complex in complexes: return complex


    def trim(self, places):
        """Rounds the coordinate values to a given number of decimal places.
        Useful for removing floating point rounding errors after transformation.

        :param int places: The number of places to round the coordinates to. If\
        ``None``, no rounding will be done."""

        for atom in self.atoms():
            atom.trim(places)


    def translate(self, *args, **kwargs):
        """Translates the structure through space, updating all atom
        coordinates accordingly.

        :param Number dx: The distance to move in the x direction.
        :param Number dy: The distance to move in the y direction.
        :param Number dz: The distance to move in the z direction.
        :param int trim: The amount of rounding to do to the atom's coordinates\
        after translating - the default is 12 decimal places but this can be\
        set to ``None`` if no rounding is to be done."""

        for atom in self._atoms:
            atom.translate(*args, **kwargs)


    def rotate(self, angle, axis, degrees=False, trim=12):
        """Rotates the structure about an axis, updating all atom coordinates
        accordingly.

        :param Number angle: The angle in radians.
        :param str axis: The axis to rotate around. Can only be 'x', 'y' or 'z'.
        :param bool degrees: if ``True`` the angle will be interpreted as\
        degrees.
        :param int trim: The amount of rounding to do to the atom's coordinates\
        after translating - the default is 12 decimal places but this can be\
        set to ``None`` if no rounding is to be done."""

        if axis not in ("x", "y", "z"):
            raise ValueError("{} is not a valid axis".format(axis))
        angle = math.radians(angle) if degrees else angle
        matrix = Atom.generate_rotation_matrix(None, angle, axis)
        atoms = list(self._atoms)
        for atom, vector in zip(atoms, [atom.location for atom in atoms]):
            atom._x, atom._y, atom._z = matrix.dot(vector)
        self.trim(trim)


    @property
    def mass(self):
        """The mass of the structure in Daltons.

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
        return math.sqrt(mean_square_deviation)


    def pairing_with(self, structure):
        """Takes another structure with the same number of atoms as this one,
        and attempts to find the nearest equivalent of every atom in this
        structure, in that structure.

        Atoms will be aligned first by element, then by name, then by number of
        bonds, then IDs, and finally by memory address - this last metric is
        used to ensure that even when allocation is essentially random, it is at
        least the same every time two structures are aligned.

        :param AtomicStructure structure: the structure to pair with.
        :raises TypeError: if the other structure is not an\
        :py:class:`.AtomicStructure`.
        :raises ValueError: if the other structure has a different number of\
        atoms.
        :rtype: ``dict``"""

        if not isinstance(structure, AtomicStructure):
            raise TypeError("{} is not an AtomicStructure".format(structure))
        atoms, other_atoms = list(self._atoms), list(structure._atoms)
        if len(atoms) != len(other_atoms):
            raise ValueError("{} and {} have different numbers of atoms".format(
             self, structure
            ))
        for l in atoms, other_atoms:
            l.sort(key=lambda a: (
             a.element, a.name, len(a.bonds), a.id, id(a)
            ))
        return {a1: a2 for a1, a2 in zip(atoms, other_atoms)}


    def superimpose_onto(self, other):
        """Superimoses this structure onto another - it will be translated so
        that its center of mass matches the other structure's, then rotated so
        as to minimise the RMSD.

        The other structure must have the same number of atoms.

        :param AtomicStructure other: The structure to superimpose onto. This\
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

        :param AtomicStructure structure: the structure to check against.
        :param bool superimpose: if ``True``, the structure will be\
        superimosed first (and then moved back).
        :raises TypeError: if the other structure is not an\
        :py:class:`.AtomicStructure`.
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
        return math.sqrt(msd)


    def copy(self):
        """Returns a copy of the structure, with its own distinct atoms.

        Its atoms will have the same ID, location, element, charge, name and
        bfactor as their counterparts in the original, but will have no bonds,
        and be a member of no model - even if their counterparts do and are.

        :rtype: ``AtomicStructure``"""

        copy = self.__class__(*[a.copy() for a in self.atoms()])
        copy._id, copy._name = self._id, self._name
        return copy


    def grid(self, size=1, margin=0):
        """A generator which models a grid around the structure and returns the
        coordinates of all the points in that grid. The origin is always one of
        those points.

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

        :param x: The x-coordinate of the centre of the sphere.
        :param y: The y-coordinate of the centre of the sphere.
        :param z: The z-coordinate of the centre of the sphere.
        :param radius: The radius of the sphere.
        :param int id: if given, only atoms whose ID matches this will be\
        returned.
        :param str name: if given, only atoms whose name matches this will be\
        returned.
        :param str element: if given, only atoms whose element matches this\
        will be returned.
        :param bool hydrogen: If ``False``, hydrogen atoms will be excluded.
        :param bool het: If ``False``, non-chain atoms will be excluded.
        :param bool metal: If ``False``, metal atoms will be excluded.
        :raises TypeError: if the model is not an atomium model object.
        :raises TypeError: if the coordinates are not numeric.
        :raises TypeError: if the radius is not numeric.
        :raises ValueError: if the radius is negative.
        :rtype: ``set``"""

        if any(not isinstance(c, (int, float)) for c in (x, y, z)):
            raise TypeError("({}, {}, {}) not valid coordinate".format(x, y, z))
        if not isinstance(radius, (int, float)):
            raise TypeError("{} is not a valid radius".format(radius))
        if radius < 0:
            raise ValueError("{} is not a valid radius".format(radius))
        atoms = filter(
         lambda a: a.distance_to((x, y, z)) <= radius, self._atoms
        )
        return set(atoms)


    def to_file_string(self, file_format, description=None):
        """Converts a structure to a filestring. Currently supported file formats
        are: .xyz and .pdb.

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



class Molecule(AtomicStructure):
    """Base class: :py:class:`AtomicStructure`

    A Molecule is a collection of atoms which form a unit of some kind.

    :param \*atoms: The :py:class:`.Atom` objects that make up the structure.\
    These can also be :py:class:`.AtomicStructure` objects, in which case the\
    atoms of that structure will be used in its place.
    :param str molecule_id: A unique str ID for the molecule. Uniqueness is not\
    actually enforced.
    :param str id: The molecule ID.
    :param str name: A name for the molecule.
    :raises TypeError: if non-atoms are given.
    :raises TypeError: if the ID or name is not str."""

    def __init__(self, *atoms, id=None, name=None):
        AtomicStructure.__init__(self, *atoms)
        if id is not None and not isinstance(id, str):
            raise TypeError("ID {} is not a string".format(id))
        if name is not None and not isinstance(name, str):
            raise TypeError("Molecule name {} is not a string".format(name))
        self._id = id
        self._name = name
        for atom in self._atoms:
            atom._molecule = self


    def __repr__(self):
        id_, name = "", ""
        if self._id: id_ = self._id + " "
        if self._name: name = self._name + ", "
        return "<{} {}({}{} atoms)>".format(
         self.__class__.__name__, id_,
         name, len(self._atoms)
        )


    @property
    def id(self):
        """The molecule's unique string ID.

        :rtype: ``str``"""

        return self._id


    @property
    def name(self):
        """The molecule's name.

        :raises TypeError: if the name given is not str."""

        return self._name


    @name.setter
    def name(self, name):
        if not isinstance(name, str):
            raise TypeError("Molecule name '{}' is not str".format(name))
        self._name = name


    @property
    def model(self):
        """The :py:class:`.Model` that the Molecule is part of.

        :rtype: ``Model``"""

        for atom in self._atoms:
            return atom.model


    def site(self, cutoff=4, water=False, main_chain=False, carbon=True):
        """Returns the :py:class:`.Site` that encompasses this molecule. This is
        all the residues with a non-hydrogen atom within 4 Angstroms of a
        non-hydrogen atom in the molecule.

        :param float cutoff: determines the distance cutoff to use. The efault\
        is 4.
        :param bool water: If ``True``, water molecules will be\
        obtained as well as residues.
        :param bool main_chain: If ``True`` main chain atoms will be considered\
        when determining binding residues.
        :param bool carbon: If ``False`` carbon atoms will not be considered.

        :rtype: :py:class:`.Site`"""

        from .chains import Site
        atoms, nearby = self.atoms(hydrogen=False), set()
        for atom in atoms:
            nearby.update(atom.nearby_atoms(cutoff, hydrogen=False))
        if not main_chain:
            nearby = [a for a in nearby if
             a.name not in ["C", "CA", "O", "N"] or not a.residue]
        if not carbon:
            nearby = [a for a in nearby if a.element != "C"]
        residues = [atom.residue for atom in nearby if atom not in atoms]
        if water:
            residues += [
             atom.molecule for atom in nearby if atom not in atoms
              and atom.molecule != None and atom.molecule.name == "HOH"
             ]
        residues = set([residue for residue in residues if residue])
        return Site(*residues, ligand=self)



class Residue(Molecule):
    """Base class: :py:class:`Molecule`

    A Residue is a subunit of some sort of polymer, such as an amino acid
    residue.

    :param \*atoms: The :py:class:`.Atom` objects that make up the structure.\
    These can also be :py:class:`.AtomicStructure` objects, in which case the\
    atoms of that structure will be used in its place.
    :param str residue_id: A unique str ID for the residue. Uniqueness is not\
    actually enforced.
    :param str id: The residue ID.
    :param str name: A name for the residue.
    :raises TypeError: if non-atoms are given.
    :raises TypeError: if the ID or name is not str."""

    def __init__(self, *atoms, **kwargs):
        Molecule.__init__(self, *atoms, **kwargs)
        self._next, self._previous = None, None
        for atom in self._atoms:
            atom._residue = self


    @property
    def full_name(self):
        """Returns the full name of the reside if it is one of the 20 canonical
        amino acids - otherwise it just returns the name itself.

        :rtype: ``str``"""

        if self._name:
            return RESIDUES.get(self._name.upper(), self._name)


    @property
    def next(self):
        """Residues can be linked to each other in a linear chain. This property
        returns the :py:class:`.Residue` downstream of this one. Alternatively,
        if you supply a residue, that residue will be assigned as the 'next' one
        downstream to this, and this residue will be upstream to that.

        Note that is a separate concept from bonds. Creating a connection of
        this kind implies no, and requires no, explicit bonding.

        :param Residue residue: The residue to connect to. If ``None`` is\
        given, any existing connection downstream of this residue will be\
        broken.
        :raises TypeError: if a non-residue is given.
        :raises ValueError: if you try to connect a residue to itself.
        :rtype: ``Residue``"""

        return self._next


    @next.setter
    def next(self, residue):
        if residue is None:
            if self._next: self._next._previous = None
            self._next = None
        elif not isinstance(residue, Residue):
            raise TypeError("{} is not a Residue".format(residue))
        elif residue is self:
            raise ValueError("Cannot link {} to itself".format(self))
        else:
            self._next = residue
            residue._previous = self


    @property
    def previous(self):
        """Residues can be linked to each other in a linear chain. This property
        returns the :py:class:`.Residue` upstream of this one. Alternatively,
        if you supply a residue, that residue will be assigned as the 'previous'
        one upstream to this, and this residue will be downstream from that.

        Note that is a separate concept from bonds. Creating a connection of
        this kind implies no, and requires no, explicit bonding.

        :param Residue residue: The residue to connect to. If ``None`` is\
        given, any existing connection upstream of this residue will be\
        broken.
        :raises TypeError: if a non-residue is given.
        :raises ValueError: if you try to connect a residue to itself.
        :rtype: ``Residue``"""

        return self._previous


    @previous.setter
    def previous(self, residue):
        if residue is None:
            if self._previous: self._previous._next = None
            self._previous = None
        elif not isinstance(residue, Residue):
            raise TypeError("{} is not a Residue".format(residue))
        elif residue is self:
            raise ValueError("Cannot link {} to itself".format(self))
        else:
            self._previous = residue
            residue._next = self


    @property
    def chain(self):
        """The :py:class:`.Chain` that the Residue is part of.

        :rtype: ``Chain``"""

        for atom in self._atoms:
            return atom.chain



RESIDUES = {
 "GLY": "glycine", "ALA": "alanine", "VAL": "valine", "LEU": "leucine",
 "ILE": "isoleucine", "MET": "methionine", "PHE": "phenylalanine",
 "TRP": "tryptophan", "PRO": "proline", "SER": "serine", "THR": "threonine",
 "CYS": "cysteine", "TYR": "tyrosine", "ASN": "asparagine", "GLN": "glutamine",
 "ASP": "aspartic acid", "GLU": "glutamic acid", "LYS": "lysine",
 "ARG": "arginine", "HIS": "histidine", "HOH": "water"
}
