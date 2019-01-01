"""Structure classes."""

import numpy as np
import rmsd
from collections import Counter, OrderedDict
from .base import StructureClass, query, StructureSet

class AtomStructure:
    """A structure made of atoms. This contains various useful methods that rely
    on a ``atoms()`` method, which the inheriting object must supply itself.

    The class would never be instantiated directly."""

    @property
    def id(self):
        """The structure's unique ID.

        :rtype: ``str``"""

        try:
            return self._id
        except AttributeError as e:
            raise AttributeError(
             "{}s don't have IDs".format(self.__class__.__name__)
            )


    @property
    def name(self):
        """The structure's name.

        :rtype: ``str``"""

        try:
            return self._name
        except AttributeError as e:
            raise AttributeError(
             "{}s don't have names".format(self.__class__.__name__)
            )


    @name.setter
    def name(self, name):
        self._name = name


    @property
    def mass(self):
        """The structure's mass - the sum of all its atoms' masses.

        :rtype: ``float``"""

        return round(sum([atom.mass for atom in self.atoms()]), 12)


    @property
    def charge(self):
        """The structure's charge - the sum of all its atoms' charges.

        :rtype: ``float``"""

        return round(sum([atom.charge for atom in self.atoms()]), 12)


    @property
    def formula(self):
        """The structure's formula as a ``Counter`` dictionary - the count of
        all its atoms' elements.

        :rtype: ``Counter``"""

        return Counter([atom.element for atom in self.atoms()])


    @property
    def center_of_mass(self):
        """Returns the center of mass of the structure. This is the average of
        all the atom coordinates, weighted by the mass of each atom.

        :rtype: ``tuple``"""

        mass = self.mass
        average_x = sum([atom.x * atom.mass for atom in self.atoms()]) / mass
        average_y = sum([atom.y * atom.mass for atom in self.atoms()]) / mass
        average_z = sum([atom.z * atom.mass for atom in self.atoms()]) / mass
        return (average_x, average_y, average_z)


    @property
    def radius_of_gyration(self):
        """The radius of gyration of a structure is a measure of how extended it
        is. It is the root mean square deviation of the atoms' distance from the
        structure's :py:meth:`.center_of_mass`.

        :rtype: ``float``"""

        center_of_mass = self.center_of_mass
        atoms = self.atoms()
        square_deviation = sum(
         [atom.distance_to(center_of_mass) ** 2 for atom in atoms]
        )
        mean_square_deviation = square_deviation / len(atoms)
        return np.sqrt(mean_square_deviation)


    def pairwise_atoms(self, *args, **kwargs):
        """A generator which yeilds all the pairwise atom combinations of the
        structure. There will be no duplicates in the returned generator, and
        the number of returned pairs will be a triangle number.

        :rtype: ``tuple``"""

        atoms = list(self.atoms(*args, **kwargs))
        for a_index in range(len(atoms) - 1):
            for o_index in range(a_index + 1, len(atoms)):
                yield {atoms[a_index], atoms[o_index]}


    def translate(self, dx=0, dy=0, dz=0, trim=12):
        """Translates the structure through space, updating all atom
        coordinates accordingly. You can provide three values, or a single
        vector.

        :param Number dx: The distance to move in the x direction.
        :param Number dy: The distance to move in the y direction.
        :param Number dz: The distance to move in the z direction.
        :param int trim: The amount of rounding to do to the atoms' coordinates\
        after translating - the default is 12 decimal places but this can be\
        set to ``None`` if no rounding is to be done."""

        try:
            _,_,_ = dx
            vector = dx
        except TypeError: vector = (dx, dy, dz)
        Atom.translate_atoms(vector, *self.atoms())
        self.trim(trim)


    def transform(self, matrix, trim=12):
        """Transforms the structure using a 3x3 matrix supplied. This is useful
        if the :py:meth:`.rotate` method isn't powerful enough for your needs.

        :param array matrix: A NumPy matrix representing the transformation.\
        You can supply a list of lists if you like and it will be converted to\
        a NumPy matrix.
        :param int trim: The amount of rounding to do to the atoms' coordinates\
        after transforming - the default is 12 decimal places but this can be\
        set to ``None`` if no rounding is to be done."""

        Atom.transform_atoms(matrix, *self.atoms())
        self.trim(trim)


    def rotate(self, angle, axis, trim=12):
        """Rotates the structure about an axis, updating all atom coordinates
        accordingly.

        :param Number angle: The angle in radians.
        :param str axis: The axis to rotate around. Can only be 'x', 'y' or 'z'.
        :param int trim: The amount of rounding to do to the atoms' coordinates\
        after translating - the default is 12 decimal places but this can be\
        set to ``None`` if no rounding is to be done."""

        Atom.rotate_atoms(angle, axis, *self.atoms())
        self.trim(trim)


    def trim(self, places):
        """Rounds the coordinate values to a given number of decimal places.
        Useful for removing floating point rounding errors after transformation.

        :param int places: The number of places to round the coordinates to. If\
        ``None``, no rounding will be done."""

        for atom in self.atoms():
            atom.trim(places)


    def pairing_with(self, structure):
        """Takes another structure with the same number of atoms as this one,
        and attempts to find the nearest equivalent of every atom in this
        structure, in that structure.

        Atoms will be aligned first by ID (if equal), then element, then by
        name, and finally by memory address - this last metric is
        used to ensure that even when allocation is essentially random, it is at
        least the same every time two structures are aligned.

        :param AtomStructure structure: the structure to pair with.
        :raises ValueError: if the other structure has a different number of\
        atoms.
        :rtype: ``dict``"""

        atoms = self.atoms()
        other_atoms = structure.atoms()
        if len(atoms) != len(other_atoms):
            raise ValueError("{} and {} have different numbers of atoms".format(
             self, structure
            ))
        pair = {}
        common_ids = set(a._id for a in atoms) & set(a._id for a in other_atoms)
        id_atoms = {a._id: a for a in atoms}
        id_other_atoms = {a._id: a for a in other_atoms}
        for id_ in common_ids:
            pair[id_atoms[id_]] = id_other_atoms[id_]
            atoms.remove(id_atoms[id_])
            other_atoms.remove(id_other_atoms[id_])
        atoms, other_atoms = list(atoms), list(other_atoms)
        for l in atoms, other_atoms:
            l.sort(key=lambda a: (
             a._element, a._name, id(a)
            ))
        return {**pair, **{a1: a2 for a1, a2 in zip(atoms, other_atoms)}}


    def rmsd_with(self, structure):
        """Calculates the Root Mean Square Deviation between this structure and
        another.

        You can get the RMSD either of the coordinates as they are, or of
        superimposed coordinates.

        :param AtomStructure structure: the structure to check against.
        :raises ValueError: if the other structure has a different number of\
        atoms.
        :rtype: ``float``"""

        pairing = self.pairing_with(structure)
        coords1, coords2 = [[a.location for a in atoms]
         for atoms in zip(*pairing.items())]
        c1, c2 = self.center_of_mass, structure.center_of_mass
        coords1 = [[x - c1[0], y - c1[1], z - c1[2]] for x, y, z in coords1]
        coords2 = [[x - c2[0], y - c2[1], z - c2[2]] for x, y, z in coords2]
        return round(rmsd.kabsch_rmsd(coords1, coords2), 12)


    def grid(self, size=1, margin=0):
        """A generator which models a grid around the structure and returns the
        coordinates of all the points in that grid. The origin is always one of
        those points, and the grid will be a box.

        :param int size: The spacing between grid points. The default is 1.
        :param int margin: How far to extend the grid beyond the structure\
        coordinates. The default is 0.
        :rtype: ``tuple``"""

        atom_locations = [atom.location for atom in self.atoms()]
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


    def atoms_in_sphere(self, location, radius, *args, **kwargs):
        """Returns all the atoms in a given sphere within this structure.

        :param tuple location: the centre of the sphere.
        :param float radius: the radius of the sphere.
        :rtype: ``set``"""

        return {a for a in self.atoms(*args, **kwargs)
         if a.distance_to(location) <= radius}


    def nearby_atoms(self, *args, **kwargs):
        """Returns all atoms within a given distance of this structure,
        excluding the structure's own atoms.

        :param float cutoff: the distance cutoff to use.
        :rtype: ``set``"""

        atoms = set()
        for atom in self.atoms():
            atoms.update(atom.nearby_atoms(*args, **kwargs))
        return atoms - self.atoms()


    def equivalent_to(self, other):
        """Two structures are equivalent if (1) they have the same number of
        atoms and (2) their atoms can be paired using :py:meth:`.pairing_with`
        such that each atom is equivalent to its pair.

        :param AtomStructure other: the structure to compare with.
        :rtype: ``bool``"""

        try:
            mapping = self.pairing_with(other)
            for atom1, atom2 in mapping.items():
                if not atom1.equivalent_to(atom2): return False
            return True
        except: return False


    def save(self, path):
        from .utilities import save
        ext = path.split(".")[-1]
        if ext == "cif":
            from .mmcif import structure_to_mmcif_string
            string = structure_to_mmcif_string(self)
        elif ext == "mmtf":
            from .mmtf import structure_to_mmtf_string
            string = structure_to_mmtf_string(self)
        elif ext == "pdb":
            from .pdb import structure_to_pdb_string
            string = structure_to_pdb_string(self)
        else:
            raise ValueError("Unsupported file extension: " + ext)
        save(string, path)



class Molecule:
    """A molecule is a top-level constituent of a :py:class:`.Model` - a chain,
    a ligand, or a water molecule."""

    @property
    def model(self):
        """Returns the molecules :py:class:`.Model`.

        :rtype: ``Model``"""

        return self._model



class Het:
    """A direct container of atoms, such as a residue or ligand. Though never
    instantiated directly, there is an initaliser method for setting up the
    atom dictionary."""

    def __init__(self, *atoms):
        for atom in atoms: atom._structure = self
        self._atoms = StructureSet(*atoms)


    def __contains__(self, atom):
        return atom in self._atoms.structures


    @property
    def chain(self):
        """Returns the :py:class:`.Chain` the structure is part of (if a
        residue) or associated with (if a ligand).

        :rtype: ``Chain``"""

        return self._chain


    @property
    def model(self):
        """Returns the :py:class:`.Model` the structure is part of, via its
        chain.

        :rtype: ``Model``"""

        try:
            return self._chain._model
        except AttributeError: return None


    def add(self, atom):
        """Adds an :py:class:`.Atom` to the structure.

        :param Atom atom: the atom to add."""

        self._atoms.add(atom)
        atom._structure = self


    def remove(self, atom):
        """Removes an :py:class:`.Atom` from the structure.

        :param Atom atom: the atom to remove."""

        self._atoms.remove(atom)
        atom._structure = None


    def nearby_structures(self, *args, **kwargs):
        """Returns all other het structures within a given distance of this
        structure, excluding itself.

        :param float cutoff: the distance cutoff to use.
        :rtype: ``set``"""

        structures = set()
        for atom in self._atoms.structures:
            structures.update(atom.nearby_structures(*args, **kwargs))
        return structures



class Model(AtomStructure, metaclass=StructureClass):
    """The universe in which all other molecules live, interact, and generally
    exist.

    It is a cotainer of its molecules, residues, and atoms.

    :param \*molecules: The chains, ligands, and waters that will inhabit the\
    model."""

    def __init__(self, *molecules):
        self._chains = StructureSet()
        self._ligands = StructureSet()
        self._waters = StructureSet()
        for mol in molecules:
            mol._model = self
            d = (self._chains if isinstance(mol, Chain) else self._waters
             if mol._water else self._ligands)
            d.add(mol)


    def __repr__(self):
        chains = "{} chains".format(len(self._chains))
        if len(self._chains) == 1: chains = chains[:-1]
        ligands = "{} ligands".format(len(self._ligands))
        if len(self._ligands) == 1: ligands = ligands[:-1]
        return "<Model ({}, {})>".format(chains, ligands)


    def __contains__(self, obj):
        return (obj in self.molecules() or obj in self.residues()
         or obj in self.atoms())


    def chains(self):
        """Returns the model's chains.

        :rtype: ``set``"""

        return self._chains


    def ligands(self):
        """Returns the model's ligands.

        :rtype: ``set``"""

        return self._ligands


    def waters(self):
        """Returns the model's water ligands.

        :rtype: ``set``"""

        return self._waters


    def molecules(self):
        """Returns all of the model's molecules (chains, ligands, waters).

        :rtype: ``set``"""

        return self._chains + self._ligands + self._waters


    def residues(self):
        """Returns all of the model's residues in all its chains.

        :rtype: ``set``"""

        res = []
        for chain in self._chains.structures:
            res += chain.residues()
        return StructureSet(*res)


    def atoms(self):
        """Returns all of the model's atoms in all its molecules.

        :rtype: ``set``"""

        atoms = set()
        for mol in self.molecules():
            try:
                atoms.update(mol._atoms.structures)
            except:
                for res in mol._residues.structures:
                    atoms.update(res._atoms.structures)
        return StructureSet(*atoms)


    def dehydrate(self):
        """Removes all water ligands from the model."""

        self._waters = StructureSet()


    def add(self, molecule):
        """Add a molecule to the model.

        :param molecule: the chain/ligand to add."""

        d = (self._chains if isinstance(molecule, Chain) else self._waters
         if molecule._water else self._ligands)
        d.add(molecule)
        molecule._model = self


    def remove(self, molecule):
        """Removes a molecule from the model.

        :param molecule: the chain/ligand to remove."""

        d = (self._chains if isinstance(molecule, Chain) else self._waters
         if molecule._water else self._ligands)
        d.remove(molecule)
        molecule._model = None



class Chain(AtomStructure, Molecule, metaclass=StructureClass):
    """A sequence of residues. Unlike other structures, they are iterables, and
    have a length.

    Residues can also be accessed using indexing.

    :param \*residues: The residues that will make up the chain.
    :param str id: the chain's unique ID.
    :param str internal_id: the internal ID used for transformations.
    :param str sequence: the actual sequence the chain should have."""

    def __init__(self, *residues, id=None, internal_id=None, sequence=""):
        self._id, self._internal_id, self._sequence = id, internal_id, sequence
        for res in residues: res._chain = self
        self._residues = StructureSet(*residues)
        self._model = None


    def __repr__(self):
        return "<Chain {} ({} residues)>".format(self._id, len(self._residues))


    def __len__(self):
        return len(self._residues)


    def __iter__(self):
        return iter(self._residues.structures)


    def __getitem__(self, key):
        return self.residues()[key]


    def __contains__(self, obj):
        return obj in self._residues.structures or obj in self.atoms()


    @property
    def length(self):
        """Returns the number of residues in the chain.

        :rtype: ``int``"""

        return len(self)


    def residues(self):
        """Returns the residues in the chain.

        :rtype: ``tuple``"""

        return self._residues


    def ligands(self):
        """Returns all the ligands associated with the chain - but only if the
        chain is part of a model.

        :rtype: ``set``"""

        return StructureSet() if self._model is None else StructureSet(
         *[l for l in self._model._ligands.structures if l._chain is self]
        )


    def atoms(self):
        """Returns all the atoms in with the chain.

        :rtype: ``set``"""

        atoms = set()
        for res in self._residues.structures:
            atoms.update(res._atoms.structures)
        return StructureSet(*atoms)


    @property
    def sequence(self):
        """Returns the sequence associated with the chain. Note that this is the
        sequence that the molecule actually has in real life - some may be
        missing from the actual sequence of residues in the structure.

        :rtype: ``str``"""

        return self._sequence


    @sequence.setter
    def sequence(self, sequence):
        self._sequence = sequence


    def copy(self):
        """Creates a copy of the chain, with new atoms and residues.

        :rtype: ``Chain``"""

        return Chain(*[r.copy() for r in self.residues()], id=self._id,
         internal_id=self._internal_id, sequence=self._sequence)



class Ligand(AtomStructure, Molecule, Het, metaclass=StructureClass):
    """A small molecule, usually associated with a polymer chain.

    :param \*atoms: The atoms that will make up the ligand.
    :param str id: the ligand's unique ID.
    :param str name: the ligand's name.
    :param str internal_id: the internal ID used for transformations.
    :param Chain chain: the chain the ligand is associated with.
    :param bool water: if ``True``, the ligand will be treated as water."""

    def __init__(self, *atoms, id=None, name=None, internal_id=None, chain=None,
                 water=False):
        self._id, self._name, self._internal_id = id, name, internal_id
        self._chain, self._water = chain, water
        Het.__init__(self, *atoms)


    def __repr__(self):
        return "<{} {} ({})>".format(
         "Water" if self._water else "Ligand", self._name, self._id
        )


    def atoms(self):
        """Returns the atoms in the ligand.

        :rtype: ``set``"""

        return self._atoms


    @property
    def water(self):
        """Returns ``True`` if the ligand is a water ligand.

        :rtype: ``bool``"""

        return self._water


    def copy(self):
        """Creates a copy of the ligand, with new atoms, and no chain.

        :rtype: ``Ligand``"""

        return Ligand(*[a.copy() for a in self.atoms()], id=self._id,
         name=self._name, internal_id=self._internal_id, water=self._water)



class Residue(AtomStructure, Het, metaclass=StructureClass):
    """A small subunit within a chain.

    :param \*atoms: The atoms the residue is to be made of.
    :param str id: The residue's ID.
    :param str name: The residue's name."""

    from atomium import data as __data

    def __init__(self, *atoms, id=None, name=None):
        self._id, self._name = id, name
        self._next, self._previous = None, None
        self._chain = None
        Het.__init__(self, *atoms)


    def __repr__(self):
        return "<Residue {} ({})>".format(self._name, self._id)


    def atoms(self):
        """Returns the atoms in the residue.

        :rtype: ``set``"""

        return self._atoms


    @property
    def next(self):
        """Residues can be linked to each other in a linear chain. This property
        returns the :py:class:`.Residue` downstream of this one. Alternatively,
        if you supply a residue, that residue will be assigned as the 'next' one
        downstream to this, and this residue will be upstream to that.
        Note that is a separate concept from bonds.

        :raises ValueError: if you try to connect a residue to itself.
        :rtype: ``Residue``"""

        return self._next


    @next.setter
    def next(self, next):
        if next is None:
            if self._next: self._next._previous = None
            self._next = None
        elif next is self:
            raise ValueError("Cannot link {} to itself".format(self))
        else:
            self._next = next
            next._previous = self


    @property
    def previous(self):
        """Residues can be linked to each other in a linear chain. This property
        returns the :py:class:`.Residue` upstream of this one. Alternatively,
        if you supply a residue, that residue will be assigned as the 'previous'
        one upstream to this, and this residue will be downstream to that.

        :raises ValueError: if you try to connect a residue to itself.
        :rtype: ``Residue``"""

        return self._previous


    @previous.setter
    def previous(self, previous):
        if previous is None:
            if self._previous: self._previous._next = None
            self._previous = None
        elif previous is self:
            raise ValueError("Cannot link {} to itself".format(self))
        else:
            self._previous = previous
            previous._next = self


    @property
    def code(self):
        """Returns the single letter code, based on its three letter name - or
        just 'X' if it doesn't match anything.

        :rtype: ``str``"""

        return self.__data.CODES.get(self._name, "X")


    @property
    def full_name(self):
        """Returns the residue's full name, based on its three letter name - or
        just the three letter name if it doesn't match anything.

        :rtype: ``str``"""

        return self.__data.FULL_NAMES.get(self._name, self._name)


    def copy(self):
        """Creates a copy of the residue, with new atoms.

        :rtype: ``Residue``"""

        return Residue(
         *[a.copy() for a in self.atoms()], id=self._id, name=self._name
        )



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
    :param number bvalue: The B-value of the atom (its uncertainty).
    :param list anisotropy: The directional uncertainty of the atom."""

    from atomium import data as __data

    __slots__ = [
     "_element", "_x", "_y", "_z", "_id", "_name", "_charge",
     "_bvalue", "_anisotropy", "_structure", "_model", "_flag"
    ]

    def __init__(self, element, x, y, z, id, name, charge, bvalue, anisotropy):
        self._element, self._x, self._y, self._z = element, x, y, z
        self._id, self._name, self._charge = id, name, charge
        self._bvalue, self._anisotropy = bvalue, anisotropy
        self._structure, self._model = None, None


    def __repr__(self):
        return "<Atom {} ({})>".format(self._id, self._name)


    @staticmethod
    def translate_atoms(vector, *atoms):
        """Translates multiple atoms using some vector.

        :param vector: the three values representing the delta position.
        :param \*atoms: the atoms to translate."""

        for atom in atoms:
            atom._x += vector[0]
            atom._y += vector[1]
            atom._z += vector[2]


    @staticmethod
    def transform_atoms(matrix, *atoms):
        """Transforms multiple atoms using some matrix.

        :param matrix: the transformation matrix.
        :param \*atoms: the atoms to transform."""

        atoms = list(atoms)
        locations = [a.location for a in atoms]
        output = np.dot(np.array(matrix), np.array(locations).transpose())
        for atom, location in zip(atoms, output.transpose()):
            atom._x, atom._y, atom._z = location


    @staticmethod
    def rotate_atoms(angle, axis, *atoms, **kwargs):
        """Rotates multiple atoms using an axis and an angle.

        :param float angle: the angle to rotate by in radians.
        :param str axis: the axis to rotate around (x, y, or z).
        :param \*atoms: the atoms to rotate."""

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
        Atom.transform_atoms(np.array([
         [aa + bb - cc - dd, 2 * (bc + ad), 2 * (bd - ac)],
         [2 * (bc - ad), aa + cc - bb - dd, 2 * (cd + ab)],
         [2 * (bd + ac), 2 * (cd - ab), aa + dd - bb - cc]
        ]), *atoms, **kwargs)


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
        self._name = name


    @property
    def charge(self):
        """The atom's charge - usually just zero, or 'neutral'.

        :rtype: ``float``"""

        return self._charge


    @charge.setter
    def charge(self, charge):
        self._charge = charge


    @property
    def bvalue(self):
        """The atom's B-value - the uncertainty in its position in all
        directions.

        :rtype: ``float``"""

        return self._bvalue


    @bvalue.setter
    def bvalue(self, bvalue):
        self._bvalue = bvalue


    @property
    def anisotropy(self):
        """The atom's directional uncertainty, represented by a list of six
        numbers.

        :rtype: ``list``"""

        return self._anisotropy


    @property
    def location(self):
        """The atom's location in space.

        :rtype: ``tuple``"""

        return (self._x, self._y, self._z)


    def distance_to(self, other):
        """Returns the distance (in whatever units the coordinates are defined
        in) between this atom and another. You can also give a (x, y, z) tuple
        instead of another atom if you so wish.

        :param Atom other: The other atom (or location tuple).
        :rtype: ``float``"""

        try:
            x, y, z = other
        except: x, y, z = other.location
        x_sum = pow((x - self._x), 2)
        y_sum = pow((y - self._y), 2)
        z_sum = pow((z - self._z), 2)
        return np.sqrt(x_sum + y_sum + z_sum)


    @property
    def mass(self):
        """The atom's molar mass according to the Periodic Table, based on the
        atom's :py:meth:`element`. If the element doesn't match any symbol on
        the Periodic Table, a mass of 0 will be returned.

        The element lookup is case-insensitive.

        :rtype: ``float``"""

        return self.__data.PERIODIC_TABLE.get(self._element.upper(), 0)


    @property
    def is_metal(self):
        """Checks whether the atom's element matches a metal element.

        The element lookup is case-insensitive.

        :rtype: ``bool``"""

        return self._element.upper() in self.__data.METALS


    @property
    def structure(self):
        """Returns the :py:class:`.Residue` or :py:class:`.Ligand` the atom is
        part of, or ``None`` if it is not part of one.

        :rtype: ``Het```"""

        return self._structure


    @property
    def chain(self):
        """Returns the :py:class:`.Chain` the atom is part of, or ``None`` if
        it is not part of one.

        :rtype: ``Chain``"""

        if self._structure: return self._structure.chain


    @property
    def model(self):
        """Returns the :py:class:`.Model` the atom is part of, or ``None`` if
        it is not part of one.

        :rtype: ``Model``"""

        if self.chain: return self.chain.model


    def angle(self, atom1, atom2):
        """Gets the angle between two atom vectors with this atom as the origin.

        :param Atom atom1: The first atom.
        :param Atom atom2: Thne second atom."""

        vectors = [
         [v1 - v2 for v1, v2 in zip(atom.location, self.location)
        ] for atom in (atom1, atom2)]
        vectors = [v / np.linalg.norm(v) for v in vectors]
        return np.arccos(np.clip(np.dot(vectors[0], vectors[1]), -1.0, 1.0))


    def nearby_atoms(self, cutoff, *args, **kwargs):
        """Returns all atoms in the associated :py:class:`.Model` that are
        within a given distance (in the units of the atom coordinates) of this
        atom. If the atom is not part of a model, no atoms will be returned.

        :param float cutoff: The radius to search within.
        :rtype: ``set``"""

        if self.model:
            atoms =  self.model.atoms_in_sphere(
             self.location, cutoff, *args, **kwargs
            )
            try:
                atoms.remove(self)
            except: pass
            return atoms
        return set()


    def nearby_structures(self, *args, residues=True, ligands=True,
                          waters=False, **kwargs):
        """Returns all residues and ligandsin the associated :py:class:`.Model`
        that are within a given distance (in the units of the atom coordinates)
        of this atom. If the atom is not part of a model, no residues will be
        returned.

        :param float cutoff: the distance cutoff to use.
        :param bool residues: if ``False``, residues will not be returned.
        :param bool ligands: if ``False``, ligands will not be returned.
        :param bool waters: if ``True``, waters will be returned too.
        :rtype: ``set``"""

        atoms = self.nearby_atoms(*args, **kwargs)
        structures = set()
        for atom in atoms:
            structures.add(atom.structure)
        try:
            structures.remove(self.structure)
        except: pass
        if not residues:
            structures = {s for s in structures if not isinstance(s, Residue)}
        if not ligands:
            structures = {s for s in structures
             if not (isinstance(s, Ligand) and not s.water)}
        if not waters:
            structures = {s for s in structures
             if not (isinstance(s, Ligand) and s.water)}
        return structures


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
            _,_,_ = dx
            vector = dx
        except TypeError: vector = (dx, dy, dz)
        Atom.translate_atoms(vector, self)
        self.trim(trim)


    def transform(self, matrix, trim=12):
        """Transforms the atom using a 3x3 matrix supplied. This is useful if
        the :py:meth:`.rotate` method isn't powerful enough for your needs.

        :param array matrix: A NumPy matrix representing the transformation.\
        You can supply a list of lists if you like and it will be converted to\
        a NumPy matrix.
        :param int trim: The amount of rounding to do to the atom's coordinates\
        after transforming - the default is 12 decimal places but this can be\
        set to ``None`` if no rounding is to be done."""

        Atom.transform_atoms(matrix, self)
        self.trim(trim)


    def rotate(self, angle, axis, trim=12):
        """Rotates the atom by an angle in radians, around one of the the three
        axes.

        :param float angle: The angle to rotate by in radians.
        :param str axis: the axis to rotate around.
        :param int trim: The amount of rounding to do to the atom's coordinates\
        after rotating - the default is 12 decimal places but this can be\
        set to ``None`` if no rounding is to be done."""

        Atom.rotate_atoms(angle, axis, self)
        self.trim(trim)


    def move_to(self, x, y, z):
        """Moves the atom to the coordinates given.

        :param number x: The atom's new x coordinate.
        :param number y: The atom's new y coordinate.
        :param number z: The atom's new z coordinate."""

        self._x, self._y, self._z = x, y, z


    def trim(self, places):
        """Rounds the coordinate values to a given number of decimal places.
        Useful for removing floating point rounding errors after transformation.

        :param int places: The number of places to round the coordinates to. If\
        ``None``, no rounding will be done."""

        if places is not None:
            self._x = round(self._x, places)
            self._y = round(self._y, places)
            self._z = round(self._z, places)


    def equivalent_to(self, other):
        """Two atoms are equivalent if they occupy the same position and have
        the same properties.

        :param Atom other: the atom to compare with.
        :rtype: `bool``"""

        if not isinstance(other, Atom):
            raise TypeError("Can't compare Atom with object type '{}'".format(
             other.__class__.__name__
            ))
        for attr in ["_element", "location", "_id", "_name",
         "_charge", "_bvalue", "_anisotropy"]:
            if getattr(self, attr) != getattr(other, attr): return False
        return True


    def copy(self):
        """Returns a copy of the atom. The new atom will have the same element,
        location, name, charge, ID, bvalue etc. as the original, but will not
        be part of any model or other molecule.

        :rtype: ``Atom``"""

        args = [getattr(self, k) for k in self.__slots__[:-3]]
        return Atom(*args)
