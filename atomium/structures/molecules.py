"""This module contains classes for structures made of atoms."""

from collections import Counter
import weakref
from math import sqrt
from .atoms import Atom

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


    def atoms(self, atom_id=None, element=None, exclude=None, name=None):
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
        if element:
            atoms = set(filter(lambda a: a.element() == element, atoms))
        if exclude:
            atoms = set(filter(lambda a: a.element() != exclude, atoms))
        if atom_id:
            atoms = set(filter(lambda a: a.atom_id() == atom_id, atoms))
        if name:
            atoms = set(filter(lambda a: a.name() == name, atoms))
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


    def pairwise_atoms(self):
        """A generator which yeilds all the pairwise atom combinations of the
        structure. There will be no duplicates in the returned generator, and
        the number of returned pairs will be a triangle number.

        :rtype: ``list``"""

        atoms = list(self.atoms())
        for a_index in range(len(atoms) - 1):
            for o_index in range(a_index + 1, len(atoms)):
                yield [atoms[a_index], atoms[o_index]]


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


    def name(self, name=None):
        """Returns the molecule's name. If a value is given, the name will be
        updated, provided it is a string.

        :param int name: If given, the name will be set to this.
        :raises TypeError: if the name given is not str."""

        if name is None:
            return self._name
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


    def site(self, include_water=False):
        """Returns the :py:class:`.Site` that encompasses this molecule. This is
        all the residues with a non-hydrogen atom within 4 Angstroms of a
        non-hydrogen atom in the molecule.

        :param bool include_water: If ``True``, water molecules will be\
        obtained as well as residues.

        :rtype: :py:class:`.Site`"""

        from .chains import Site
        atoms, nearby = self.atoms(exclude="H"), set()
        for atom in atoms:
            nearby.update(atom.nearby(4, exclude="H"))
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
