"""This module contains classes for structures made of atoms."""

from collections import Counter
from itertools import combinations
from points import Vector
from points import Matrix
import numpy as np
import rmsd
import points
import weakref
from math import sqrt, pi, degrees, floor, ceil
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
        atoms_ = set()
        for atom in atoms:
            if isinstance(atom, Atom):
                atoms_.add(atom)
            elif isinstance(atom, AtomicStructure):
                atoms_.update(atom.atoms())
            else:
                raise TypeError(
                 "AtomicStructures need atoms, not '{}'".format(atom)
                )
        self._id_atoms = {a.atom_id(): a for a in atoms_ if a.atom_id()}
        self._atoms = set([a for a in atoms_ if not a.atom_id()])


    def __repr__(self):
        return "<{} ({} atoms)>".format(
         self.__class__.__name__, len(self._atoms) + len(self._id_atoms)
        )


    def __contains__(self, member):
        return member in self._atoms or member in self._id_atoms.values()


    @atom_query
    def atoms(self):
        """Returns the :py:class:`.Atom` objects in the structure. You can
        filter these by element if you wish.

        :param str element: If given, only atoms whose element matches this\
        will be returned.
        :param str exclude: If given, only atoms whose element doesn't match\
        this will be returned.
        :param int atom_id: If given, only atoms whose atom ID matches this\
        will be returned (this will only return one atom).
        :param str name: If given, only atoms whose name matches this will be\
        returned.
        :rtype: ``set``"""

        atoms = self._atoms.union(set(self._id_atoms.values()))
        return atoms


    def atom(self, *args, **kwargs):
        """Returns the first :py:class:`.Atom` that matches the criteria given.

        Note that atoms are stored unordered in a ``set`` so if more than one
        atom matches the criteria you give, it might not be the same atom that
        is returned each time you call this method.

        :param int atom_id: If given, only atoms whose atom ID matches this\
        will be searched.
        :param str element: If given, only atoms whose element matches this\
        will be searched.
        :param str exclude: If given, only atoms whose element doesn't match\
        this will be returned.
        :param str name: If given, only atoms whose name matches this will be\
        searched.
        :rtype: ``Atom``"""

        if "atom_id" in kwargs:
            return self._id_atoms.get(kwargs["atom_id"])
        if len(args) == 1:
            return self._id_atoms.get(args[0])
        atoms = self.atoms(*args, **kwargs)
        for atom in atoms: return atom


    def add_atom(self, atom):
        """Adds an :py:class:`.Atom` to the structure.

        :param Atom atom: The atom to add.
        :raises TypeError: if the atom given is not an Atom."""

        if not isinstance(atom, Atom):
            raise TypeError("Can only add atoms, not '{}'".format(atom))
        if atom.atom_id():
            self._id_atoms[atom.atom_id()] = atom
        else:
            self._atoms.add(atom)


    def remove_atom(self, atom):
        """Removes an :py:class:`.Atom` from the structure.

        :param Atom atom: The atom to remove."""

        if atom.atom_id():
            del self._id_atoms[atom.atom_id()]
        else:
            self._atoms.remove(atom)


    def pairwise_atoms(self, *args, **kwargs):
        """A generator which yeilds all the pairwise atom combinations of the
        structure. There will be no duplicates in the returned generator, and
        the number of returned pairs will be a triangle number.

        :rtype: ``list``"""

        atoms = list(self.atoms(*args, **kwargs))
        for a_index in range(len(atoms) - 1):
            for o_index in range(a_index + 1, len(atoms)):
                yield [atoms[a_index], atoms[o_index]]


    def grid(self, size=1, margin=0):
        """A generator which models a grid around the structure and returns the
        coordinates of all the points in that grid. The origin is always one of
        those points.

        :param int size: The spacing between grid points. The default is 1.
        :param int margin: How far to extend the grid beyond the structure\
        coordinates. The default is 0.
        :rtype: ``tuple``"""

        atom_locations = [atom.location() for atom in self.atoms()]
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


    def mass(self):
        """Returns the mass of the structure in Daltons, based on the masses of
        its atoms.

        :rtype: ``float``"""

        return sum([atom.mass() for atom in self.atoms()])


    def charge(self):
        """Returns the charge of the structure, based on the charges of its
        atoms.

        :rtype: ``float``"""

        return sum([atom.charge() for atom in self.atoms()])


    def formula(self):
        """Returns the formula (count of each atom) of the structure.

        :rtype: ``Counter``"""

        return Counter([atom.element() for atom in self.atoms()])


    def round(self, places):
        """Rounds the coordinate values of all the atoms to a given number of
        decimal places. Useful for removing floating point rounding errors after
        rotation.

        :param int places: The number of places to round the coordinates to."""

        for atom in self.atoms():
            atom.round(places)


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
        atoms, other_atoms = list(self.atoms()), list(structure.atoms())
        if len(atoms) != len(other_atoms):
            raise ValueError("{} and {} have different numbers of atoms".format(
             self, structure
            ))
        for l in atoms, other_atoms:
            l.sort(key=lambda a: (
             a.element(), a.name(), len(a.bonds()), a.atom_id(), id(a)
            ))
        return {a1: a2 for a1, a2 in zip(atoms, other_atoms)}


    def rmsd_with(self, structure):
        """Calculates the Root Mean Square Deviation between this structure and
        another.

        No translation or rotation is performed beforehand - atom coordinates
        are compared as-is.

        :raises TypeError: if the other structure is not an\
        :py:class:`.AtomicStructure`.
        :raises ValueError: if the other structure has a different number of\
        atoms.
        :rtype: ``float``"""

        pairing = self.pairing_with(structure)
        sd = sum(a1.distance_to(a2) ** 2 for a1, a2 in pairing.items())
        msd = sd / len(pairing)
        return sqrt(msd)


    def translate(self, dx, dy, dz):
        """Translates the structure through space, updating all atom
        coordinates accordingly.

        :param Number dx: The distance to move in the x direction.
        :param Number dy: The distance to move in the y direction.
        :param Number dz: The distance to move in the z direction."""

        for atom in self.atoms():
            atom.translate(dx, dy, dz)


    def rotate(self, angle, axis):
        """Rotates the structure about an axis, updating all atom coordinates
        accordingly.

        :param Number angle: The angle in radians.
        :param str axis: The axis to rotate around. Can only be 'x', 'y' or 'z'."""

        for atom in self.atoms():
            atom.rotate(angle, axis)


    def superimpose_onto(self, structure):
        this_center, other_center = self.center_of_mass(), structure.center_of_mass()
        self.translate(-this_center[0], -this_center[1], -this_center[2])
        structure.translate(-other_center[0], -other_center[1], -other_center[2])
        atoms, other_atoms = zip(*self.pairing_with(structure).items())
        P = np.array([[*atom.location()] for atom in atoms])
        Q = np.array([[*atom.location()] for atom in other_atoms])
        P = rmsd.kabsch_rotate(P, Q)
        for atom, row in zip(atoms, P):
            atom._x, atom._y, atom._z = row
        structure.translate(other_center[0], other_center[1], other_center[2])
        self.translate(other_center[0], other_center[1], other_center[2])


    def center_of_mass(self):
        """Returns the center of mass of the structure. This is the average of
        all the atom coordinates, weighted by the mass of each atom.

        :returns: (x, y, z) ``tuple``"""

        mass = self.mass()
        atoms = self.atoms()
        average_x = sum([atom._x * atom.mass() for atom in atoms]) / mass
        average_y = sum([atom._y * atom.mass() for atom in atoms]) / mass
        average_z = sum([atom._z * atom.mass() for atom in atoms]) / mass
        return (average_x, average_y, average_z)


    def radius_of_gyration(self):
        """The radius of gyration of a structure is a measure of how extended it
        is. It is the root mean square deviation of the atoms' distance from the
        structure's :py:meth:`.center_of_mass`.

        :rtype: ``float``"""

        center_of_mass = self.center_of_mass()
        atoms = self.atoms()
        square_deviation = sum(
         [atom.distance_to(center_of_mass) ** 2 for atom in atoms]
        )
        mean_square_deviation = square_deviation / len(atoms)
        return sqrt(mean_square_deviation)


    @atom_query
    def atoms_in_sphere(self, x, y, z, radius):
        """Returns all the atoms in a given sphere within the structure.

        :param x: The x-coordinate of the centre of the sphere.
        :param y: The y-coordinate of the centre of the sphere.
        :param z: The z-coordinate of the centre of the sphere.
        :param radius: The radius of the sphere.
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
         lambda a: a.distance_to((x, y, z)) <= radius, self.atoms()
        )
        return set(atoms)


    def find_pattern(self, pattern):
        print("")
        if not isinstance(pattern, AtomicStructure):
            raise TypeError("pattern {} is not AtomicStructure".format(pattern))

        """Lesk has described an algorithm primarily for pattern
        searching in proteins but which can also be used with
        the smaller molecules present in the CCDB. The
        algorithm assigns a structure atom, S,, as a candidate
        match for a pattern atom, P,“, if S, has other structure
        atoms of the appropriate atomic types at the same distances
        from it as does P,. All of the structure atoms
        which are not matched to any pattern atom are removed
        from consideration and the candidate matches are
        checked again to ensure that the deletions have not
        affected any of the current structure atom-to-pattern
        atom matches. This process is repeated until no more
        eliminations can be made. """

        # Get respective atom sets
        pattern_atoms, structure_atoms = list(pattern.atoms()), self.atoms()

        # Start dictionary of potential matches
        matches = {p_atom: set() for p_atom in pattern_atoms}

        # For each atom in pattern, get distances to other pattern atoms
        nearby = {atom: [(a, atom.distance_to(a)) for a in atom.nearby(3)
         if a in pattern.atoms()] for atom in pattern_atoms}

        index = 1
        while True:
            print("Iteration", index)
            print("There are {} structure atoms under consideration".format(len(structure_atoms)))
            matches = {p_atom: set() for p_atom in pattern_atoms}
            # Go through each pattern atom and try and find matches
            for p_atom in list(pattern_atoms):
                # What are the other atoms in this pattern?
                neighbors = filter(lambda a: a is not p_atom, pattern_atoms)

                # Go through each atom in the structure and see if it is a match
                for atom in structure_atoms:

                    # First of all, does it have the same name as the pattern atom?
                    if atom.name() == p_atom.name():

                        # If so, go through each pattern neighbor and see if there
                        # is a corresponding structure atom
                        for p_atom_neighbor in neighbors:
                            # Get all structure atoms that are broadly at the same
                            # distance from this structure atom as the pattern atom
                            # neighbor is
                            distance = p_atom_neighbor.distance_to(p_atom)
                            atoms_at_distance = atom.nearby(
                             distance + 1, atoms=set(structure_atoms)
                            ) - atom.nearby(
                             distance - 1, atoms=set(structure_atoms)
                            )
                            # Do any of them have the same name as this pattern
                            # neighbor? If not, stop looking through the patterns as
                            # this structure atom is not suitable
                            if p_atom_neighbor.name() not in [a.name() for a in atoms_at_distance]:
                                break
                        else:
                            # If all of the neighbor atoms were ok, this is a match
                            matches[p_atom].add(atom)

            # Remove structure atoms which are not potential matches

            reduced_structure_atoms = set([a for sublist in matches.values() for a in sublist])
            for p_atom in matches:
                print("Pattern atom:", p_atom)
                for match in matches[p_atom]:
                    print("\tMatch:", match)
            if structure_atoms == reduced_structure_atoms:
                break
            structure_atoms = reduced_structure_atoms
            index += 1

        # Generate all possible combinations of the structure atoms
        possibles = combinations(structure_atoms, len(pattern_atoms))
        print(possibles)
        print(len(list(possibles)), "possible combinations of these structure atoms")

        for possible in possibles:
            pass
        print("")


    def copy(self):
        """Returns a copy of the structure, with its own distinct atoms.

        Note that the returned structure will just be a plain
        :py:class:`.AtomicStructure` and not a chain or model or whatever the
        original structure was. Its atoms will have the same ID, location,
        element, charge, name and bfactor as their counterparts in the original,
        but will have no bonds, and be a member of no molecule or model - even
        if their counterparts do and are.

        :rtype: ``AtomicStructure``"""

        return AtomicStructure(*[a.copy() for a in self.atoms()])


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
    :param str name: A name for the molecule.
    :raises TypeError: if non-atoms are given.
    :raises TypeError: if the molecule_id is not str."""

    def __init__(self, *atoms, molecule_id=None, name=None):
        AtomicStructure.__init__(self, *atoms)
        if molecule_id is not None and not isinstance(molecule_id, str):
            raise TypeError("ID {} is not a string".format(molecule_id))
        if name is not None and not isinstance(name, str):
            raise TypeError("Molecule name {} is not a string".format(name))
        self._id = molecule_id
        self._name = name
        for atom in self._atoms:
            atom._molecule = self
        for atom in self._id_atoms:
            self._id_atoms[atom]._molecule = self


    def __repr__(self):
        id_, name = "", ""
        if self._id: id_ = self._id + " "
        if self._name: name = self._name + ", "
        return "<{} {}({}{} atoms)>".format(
         self.__class__.__name__, id_,
         name, len(self._atoms) + len(self._id_atoms)
        )


    def molecule_id(self, molecule_id=None):
        """Returns the molecule's unique string ID.

        :rtype: ``str``"""

        return self._id


    def name(self, name=None, full=False):
        """Returns the molecule's name. If a value is given, the name will be
        updated, provided it is a string.

        :param int name: If given, the name will be set to this.
        :param bool full: If ``True`` the full English name will be returned.
        :raises TypeError: if the name given is not str."""

        if name is None:
            return RESIDUES.get(self._name, self._name) if full else self._name
        else:
            if not isinstance(name, str):
                raise TypeError("Molecule name '{}' is not str".format(name))
            self._name = name


    def add_atom(self, atom, *args, **kwargs):
        """Adds an :py:class:`.Atom` to the molecule.

        :param Atom atom: The atom to add.
        :raises TypeError: if the atom given is not an Atom."""

        AtomicStructure.add_atom(self, atom, *args, **kwargs)
        atom._molecule = self


    def remove_atom(self, atom, *args, **kwargs):
        """Removes an :py:class:`.Atom` from the molecule.

        :param Atom atom: The atom to remove."""

        AtomicStructure.remove_atom(self, atom, *args, **kwargs)
        atom._molecule = None


    def model(self):
        """Returns the :py:class:`.Model` that the Molecule is part of.

        :rtype: ``Model``"""

        for atom in self.atoms():
            return atom.model()


    def site(self, include_water=False, main_chain=False):
        """Returns the :py:class:`.Site` that encompasses this molecule. This is
        all the residues with a non-hydrogen atom within 4 Angstroms of a
        non-hydrogen atom in the molecule.

        :param bool include_water: If ``True``, water molecules will be\
        obtained as well as residues.
        :param bool main_chain: If ``True`` main chain atoms will be considered\
        when determining binding residues.

        :rtype: :py:class:`.Site`"""

        from .chains import Site
        atoms, nearby = self.atoms(hydrogen=False), set()
        for atom in atoms:
            nearby.update(atom.nearby(4, hydrogen=False))
        if not main_chain:
            nearby = [a for a in nearby if
             a.name() not in ["C", "CA", "O", "N"] or not a.residue()]
        residues = [atom.residue() for atom in nearby if atom not in atoms]
        if include_water:
            residues += [
             atom.molecule() for atom in nearby if atom not in atoms
              and atom.molecule() != None and atom.molecule().name() == "HOH"
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
    :param str name: A name for the residue.
    :raises TypeError: if non-atoms are given.
    :raises TypeError: if the residue_id is not str."""

    def __init__(self, *atoms, residue_id=None, **kwargs):
        if residue_id: kwargs["molecule_id"] = residue_id
        Molecule.__init__(self, *atoms, **kwargs)
        self._next, self._previous = None, None
        self._side_chains = set()
        for atom in self._atoms:
            atom._residue = self
        for atom in self._id_atoms:
            self._id_atoms[atom]._residue = self


    def __contains__(self, member):
        return Molecule.__contains__(self, member)


    def residue_id(self, residue_id=None):
        """Returns the residue's unique string ID.

        :rtype: ``str``"""

        return self._id


    def add_atom(self, atom, *args, **kwargs):
        """Adds an :py:class:`.Atom` to the residue.

        :param Atom atom: The atom to add.
        :raises TypeError: if the atom given is not an Atom."""

        Molecule.add_atom(self, atom, *args, **kwargs)
        atom._residue = self


    def remove_atom(self, atom, *args, **kwargs):
        """Removes an :py:class:`.Atom` from the residue.

        :param Atom atom: The atom to remove."""

        Molecule.remove_atom(self, atom, *args, **kwargs)
        atom._residue = None


    def next(self, residue=""):
        """Residues can be linked to each other in a linear chain. This method
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

        if residue is None:
            if self._next: self._next._previous = None
            self._next = None
        elif residue == "":
            return self._next
        elif not isinstance(residue, Residue):
            raise TypeError("{} is not a Residue".format(residue))
        elif residue is self:
            raise ValueError("Cannot link {} to itself".format(self))
        else:
            self._next = residue
            residue._previous = self


    def previous(self, residue=""):
        """Residues can be linked to each other in a linear chain. This method
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

        if residue is None:
            if self._previous: self._previous._next = None
            self._previous = None
        elif residue == "":
            return self._previous
        elif not isinstance(residue, Residue):
            raise TypeError("{} is not a Residue".format(residue))
        elif residue is self:
            raise ValueError("Cannot link {} to itself".format(self))
        else:
            self._previous = residue
            residue._next = self


    def chain(self):
        """Returns the :py:class:`.Chain` that the Residue is part of.

        :rtype: ``Chain``"""

        for atom in self.atoms():
            return atom.chain()



RESIDUES = {
 "GLY": "glycine", "ALA": "alanine", "VAL": "valine", "LEU": "leucine",
 "ILE": "isoleucine", "MET": "methionine", "PHE": "phenylalanine",
 "TRP": "tryptophan", "PRO": "proline", "SER": "serine", "THR": "threonine",
 "CYS": "cysteine", "TYR": "tyrosine", "ASN": "asparagine", "GLN": "glutamine",
 "ASP": "aspartic acid", "GLU": "glutamic acid", "LYS": "lysine",
 "ARG": "arginine", "HIS": "histidine", "HOH": "water"
}
