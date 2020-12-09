"""Structure classes."""

import numpy as np
import rmsd
import math
import warnings
from collections import Counter, OrderedDict, defaultdict
from .base import StructureClass, query, StructureSet

class AtomStructure:
    """A structure made of atoms. This contains various useful methods that rely
    on a ``atoms()`` method, which the inheriting object must supply itself. All
    atomic structures also have IDs and names.

    Two atomic structures are equal if every pairwise atom in their pairing
    are equal.

    The class would never be instantiated directly."""

    def __init__(self, id=None, name=None):
        self._id, self._name = id, name


    def __eq__(self, other):
        try:
            mapping = self.pairing_with(other)
            for atom1, atom2 in mapping.items():
                if not atom1 == atom2: return False
            return True
        except: return False


    def __hash__(self):
        return id(self)


    @property
    def id(self):
        """The structure's unique ID.

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
        locations = np.array([a._location * a.mass for a in self.atoms()])
        return np.sum(locations, axis=0) / mass


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
             a._id, a._element, a._name, id(a)
            ))
        return {**pair, **{a1: a2 for a1, a2 in zip(atoms, other_atoms)}}


    def rmsd_with(self, structure):
        """Calculates the Root Mean Square Deviation between this structure and
        another.

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


    def create_grid(self, size=1, margin=0):
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


    def check_ids(self):
        """Looks through all the structure's sub-structures and raises a
        warning if they have duplicate ID."""

        for objects in ("chains", "ligands", "waters", "residues", "atoms"):
            try:
                ids = [obj.id for obj in getattr(self, objects)()]
                unique_ids = set(ids)
                if len(ids) != len(unique_ids):
                    warnings.warn(f"{objects} have duplicate IDs")
            except AttributeError: pass


    def save(self, path):
        """Saves the structure to file. The file extension given in the filename
        will be used to determine which file format to save in.

        If the structure you are saving has any duplicate IDs, a warning will be
        issued, as the file saved will likely be nonsensical.

        :param str path: the filename and location to save to."""

        from .utilities import save
        self.check_ids()
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


    def atoms_in_sphere(self, location, radius, *args, **kwargs):
        """Returns all the atoms in a given sphere within this structure. This
        will be a lot faster if the structure is a :py:class:`.Model` and if
        :py:meth:`.optimise_distances` has been called, as it won't have to
        search all atoms.

        :param tuple location: the centre of the sphere.
        :param float radius: the radius of the sphere.
        :rtype: ``set``"""

        if "_internal_grid" in self.__dict__ and self._internal_grid:
            r, atoms = math.ceil(radius / 10), set()
            x, y, z = [int(math.floor(n / 10)) * 10 for n in location]
            x_range, y_range, z_range = [
             [(val - (n * 10)) for n in range(1, r + 1)][::-1] + [val] + [
              (val + n * 10) for n in range(1, r + 1)
             ] for val in (x, y, z)
            ]
            for x in x_range:
                for y in y_range:
                    for z in z_range:
                        atoms = atoms.union(self._internal_grid[x][y][z])
            atoms = StructureSet(*atoms)
            atoms = query(lambda self: atoms)(self, *args, **kwargs)
        else:
            atoms = self.atoms(*args, **kwargs)
        return {a for a in atoms if a.distance_to(location) <= radius}


    def pairwise_atoms(self, *args, **kwargs):
        """A generator which yeilds all the pairwise atom combinations of the
        structure. There will be no duplicates in the returned generator, and
        the number of returned pairs will be a triangle number.

        :rtype: ``tuple``"""

        atoms = list(self.atoms(*args, **kwargs))
        for a_index in range(len(atoms) - 1):
            for o_index in range(a_index + 1, len(atoms)):
                yield {atoms[a_index], atoms[o_index]}


    def nearby_atoms(self, *args, **kwargs):
        """Returns all atoms within a given distance of this structure,
        excluding the structure's own atoms.

        This will be a lot faster if the model's
        :py:meth:`.optimise_distances` has been called, as it won't have to
        search all atoms.

        :param float cutoff: the distance cutoff to use.
        :rtype: ``set``"""

        atoms = set()
        for atom in self.atoms():
            atoms.update(atom.nearby_atoms(*args, **kwargs))
        return atoms - self.atoms()
    

    def nearby_hets(self, *args, **kwargs):
        """Returns all other het structures within a given distance of this
        structure, excluding itself.

        This will be a lot faster if the model's
        :py:meth:`.optimise_distances` has been called, as it won't have to
        search all atoms.

        :param float cutoff: the distance cutoff to use.
        :rtype: ``set``"""

        structures = set()
        hets = set()
        for atom in self.atoms():
            structures.update(atom.nearby_hets(*args, **kwargs))
            hets.add(atom.het)
        return structures - hets
    

    def nearby_chains(self, *args, **kwargs):
        """Returns all other chain structures within a given distance of this
        structure, excluding itself.

        :param float cutoff: the distance cutoff to use.
        :rtype: ``set``"""

        structures = set()
        chains = set()
        for atom in self.atoms():
            structures.update(atom.nearby_chains(*args, **kwargs))
            chains.add(atom.chain)
        return structures - chains


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



class Molecule(AtomStructure):
    """A molecule is a top-level constituent of a :py:class:`.Model` - a chain,
    a ligand, or a water molecule. They can have internal IDs, separate from the
    standard ID."""

    def __init__(self, id, name, internal_id):
        AtomStructure.__init__(self, id, name)
        self._internal_id = internal_id
        self._model = None


    @property
    def internal_id(self):
        """The molecule's internal ID - how it is refered to by atomium
        operations. This will be identical to regular IDs when the model comes
        from a .pdb file, but .cif and .mmtf files make this distinction.
        
        :rtype: ``str``"""

        return self._internal_id or self._id
       

    @property
    def model(self):
        """Returns the molecules :py:class:`.Model`.

        :rtype: ``Model``"""

        return self._model



class Het(AtomStructure):
    """A direct container of atoms, such as a residue or ligand. Though never
    instantiated directly, there is an initaliser method for setting up the
    atom dictionary."""

    from atomium import data as __data

    def __init__(self, id, name, full_name, *atoms):
        AtomStructure.__init__(self, id, name)
        self._full_name = full_name
        for atom in atoms: atom._het = self
        self._atoms = StructureSet(*atoms)


    def __contains__(self, atom):
        return atom in self._atoms.structures
    

    @property
    def full_name(self):
        """Returns the residue's full name, based on its three letter name - or
        just the three letter name if it doesn't match anything. Or you can just
        supply a full name when you instantiate the Het.

        :rtype: ``str``"""

        if self._full_name: return self._full_name
        return self.__data.FULL_NAMES.get(self._name, self._name)
    

    @full_name.setter
    def full_name(self, full_name):
        self._full_name = full_name
    

    @property
    def chain(self):
        """Returns the :py:class:`.Chain` the structure is part of (if a
        residue) or associated with (if a ligand).

        :rtype: ``Chain``"""

        return self._chain


    def atoms(self):
        """Returns the atoms in the ligand.

        :rtype: ``set``"""

        return self._atoms


    
class Model(AtomStructure, metaclass=StructureClass):
    """The universe in which all other molecules live, interact, and generally
    exist.

    It is a cotainer of its molecules, residues, and atoms.

    :param \*molecules: The chains, ligands, and waters that will inhabit the\
    model."""

    def __init__(self, *molecules, file=None):
        AtomStructure.__init__(self, None, None)
        self._chains = set()
        self._ligands = set()
        self._waters = set()
        for mol in molecules:
            mol._model = self
            d = (self._chains if isinstance(mol, Chain) else self._waters
             if mol._water else self._ligands)
            d.add(mol)
        self._chains = StructureSet(*self._chains)
        self._ligands = StructureSet(*self._ligands)
        self._waters = StructureSet(*self._waters)
        self._file = file
        self._internal_grid = None


    def __repr__(self):
        chains = "{} chains".format(len(self._chains))
        if len(self._chains) == 1: chains = chains[:-1]
        ligands = "{} ligands".format(len(self._ligands))
        if len(self._ligands) == 1: ligands = ligands[:-1]
        return "<Model ({}, {})>".format(chains, ligands)


    def __contains__(self, obj):
        return (obj in self.molecules() or obj in self.residues()
         or obj in self.atoms())


    @property
    def file(self):
        """The :py:class:`.File` the model comes from."""

        return self._file


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
    

    def optimise_distances(self):
        """Calling this method makes finding atoms within a sphere faster, and
        consequently makes all 'nearby' methods faster. It organises the atoms
        in the model into grids, so that only relevant atoms are checked for
        distances."""

        self._internal_grid = defaultdict(
         lambda: defaultdict(lambda: defaultdict(set))
        )
        for atom in self.atoms():
            x, y, z = [int(math.floor(n / 10)) * 10 for n in atom.location]
            self._internal_grid[x][y][z].add(atom)


    #TODO copy



class Chain(Molecule, metaclass=StructureClass):
    """A sequence of residues. Unlike other structures, they are iterable, and
    have a length.

    Residues can also be accessed using indexing.

    :param \*residues: The residues that will make up the chain.
    :param str id: the chain's unique ID.
    :param str internal_id: the internal ID used for transformations.
    :param str sequence: the actual sequence the chain should have.
    :param list helices: the alpha helices within the chain.
    :param list strands: the beta strands within the chain."""

    def __init__(self, *residues, sequence="", helices=None, strands=None, **kwargs):
        Molecule.__init__(
         self, kwargs.get("id"), kwargs.get("name"), kwargs.get("internal_id")
        )
        self._sequence = sequence
        for res in residues: res._chain = self
        self._residues = StructureSet(*residues)
        self._model = None
        self._helices = helices or []
        self._strands = strands or []


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
    def sequence(self):
        """Returns the sequence associated with the chain. Note that this is the
        sequence that the molecule actually has in real life - some may be
        missing from the actual sequence of residues in the structure.

        :rtype: ``str``"""

        return self._sequence


    @sequence.setter
    def sequence(self, sequence):
        self._sequence = sequence


    @property
    def helices(self):
        """The alpha helix residues in the chain

        :rtype: ``tuple``"""

        return tuple(self._helices)
    

    @property
    def strands(self):
        """The beta strand residues in the chain

        :rtype: ``tuple``"""

        return tuple(self._strands)


    @property
    def length(self):
        """Returns the number of residues in the chain.

        :rtype: ``int``"""

        return len(self)
    

    @property
    def present_sequence(self):
        """The sequence of residues actually present in the atoms present.

        :rtype: ``str``"""

        return "".join(r.code for r in self.residues())
    

    def copy(self, id=None, residue_ids=None, atom_ids=None):
        """Creates a copy of the chain, with new atoms and residues.

        :param str id: if given, the ID of the new chain.
        :param function residue_ids: a callable which, if given, will generate\
        new residue IDs.
        :param function atom_ids: a callable which, if given, will generate new\
        atom IDs.
        :rtype: ``Chain``"""

        residue_ids = residue_ids or (lambda i: i)
        residues = {r: r.copy(
         id=residue_ids(r.id), atom_ids=atom_ids
        ) for r in self.residues()}
        for r in self.residues():
            residues[r].next = residues[r.next] if r.next else None
        return Chain(
         *residues.values(), id=id or self._id, internal_id=self._internal_id,
         name=self._name, sequence=self._sequence,
         helices=[tuple(residues[r] for r in h) for h in self._helices],
         strands=[tuple(residues[r] for r in s) for s in self._strands]
        )
    

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



class Ligand(Molecule, Het, metaclass=StructureClass):
    """A small molecule, usually associated with a polymer chain.

    :param \*atoms: The atoms that will make up the ligand.
    :param str id: the ligand's unique ID.
    :param str name: the ligand's name.
    :param str internal_id: the internal ID used for transformations.
    :param Chain chain: the chain the ligand is associated with.
    :param bool water: if ``True``, the ligand will be treated as water."""

    def __init__(self, *atoms, chain=None, water=False, **kwargs):
        Het.__init__(
        self, kwargs.get("id"), kwargs.get("name"),
         kwargs.get("full_name"), *atoms)
        Molecule.__init__(self, kwargs.get("id"), kwargs.get("name"),
         kwargs.get("internal_id"))
        self._chain, self._water = chain, water


    def __repr__(self):
        return "<{} {} ({})>".format(
         "Water" if self._water else "Ligand", self._name, self._id
        )


    @property
    def is_water(self):
        """Returns ``True`` if the ligand is a water ligand.

        :rtype: ``bool``"""

        return self._water


    def copy(self, id=None, atom_ids=None):
        """Creates a copy of the ligand, with new atoms.

        :param str id: if given, the ID of the new ligand.
        :param function atom_ids: a callable which, if given, will generate new\
        atom IDs.
        :rtype: ``Ligand``"""

        atoms = list(self.atoms())
        if atom_ids:
            new_ids = [atom_ids(a.id) for a in atoms]
            atoms = [a.copy(id=id) for a, id in zip(atoms, new_ids)]
        else:
            atoms = [a.copy() for a in self.atoms()]
        return self.__class__(*atoms, id=id or self._id,
         name=self._name, internal_id=self._internal_id, water=self._water)



class Residue(Het, metaclass=StructureClass):
    """A small subunit within a chain.

    :param \*atoms: The atoms the residue is to be made of.
    :param str id: The residue's ID.
    :param str name: The residue's name."""

    from atomium import data as __data

    def __init__(self, *atoms, **kwargs):
        Het.__init__(self, kwargs.get("id"), kwargs.get("name"),
         kwargs.get("full_name"), *atoms)
        self._next, self._previous = None, None
        self._chain = None


    def __repr__(self):
        return "<Residue {} ({})>".format(self._name, self._id)


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
    def helix(self):
        """Returns ``True`` if the residue is part of an alpha helix.

        :rtype: ``bool``"""

        if self.chain:
            for helix in self.chain.helices:
                if self in helix: return True
        return False
    

    @property
    def strand(self):
        """Returns ``True`` if the residue is part of a beta strand.

        :rtype: ``bool``"""

        if self.chain:
            for strand in self.chain.strands:
                if self in strand: return True
        return False


    def copy(self, id=None, atom_ids=None):
        """Creates a copy of the residue, with new atoms.

        :param str id: if given, the ID of the new residue.
        :param function atom_ids: a callable which, if given, will\
        generate new atom IDs.
        :rtype: ``Residue``"""

        atoms = list(self.atoms())
        if atom_ids:
            new_ids = [atom_ids(a.id) for a in atoms]
            atoms = [a.copy(id=id) for a, id in zip(atoms, new_ids)]
        else:
            atoms = [a.copy() for a in self.atoms()]
        return self.__class__(*atoms, id=id or self._id, name=self._name)
    

    @property
    def model(self):
        """Returns the :py:class:`.Model` the residue is part of, via its
        chain.

        :rtype: ``Model``"""

        try:
            return self._chain._model
        except AttributeError: return None



class Atom:
    """An atom in space - a point particle with a location, element, charge etc.

    Atoms are the building blocks of all structures in atomium.

    Two atoms are equal if they have the same properties (not including ID).

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
     "_element", "_location", "_id", "_name", "_charge",
     "_bvalue", "_anisotropy", "_het", "_bonded_atoms", "_is_hetatm"
    ]

    def __init__(self, element, x, y, z, id, name, charge, bvalue, anisotropy, is_hetatm=False):
        self._location = np.array([x, y, z])
        self._element = element
        self._id, self._name, self._charge = id, name, charge
        self._bvalue, self._anisotropy = bvalue, anisotropy
        self._het, self._bonded_atoms, self._is_hetatm = None, set(), is_hetatm


    def __repr__(self):
        return "<Atom {} ({})>".format(self._id, self._name)


    def __iter__(self):
        return iter(self._location)


    def __eq__(self, other):
        if not isinstance(other, Atom): return False
        for attr in self.__slots__:
            if attr not in ("_id", "_het", "_bonded_atoms", "_location"):
                if getattr(self, attr) != getattr(other, attr): return False
            if list(self._location) != list(other._location): return False
        return True


    def __hash__(self):
        return id(self)


    @staticmethod
    def translate_atoms(vector, *atoms):
        """Translates multiple atoms using some vector.

        :param vector: the three values representing the delta position.
        :param \*atoms: the atoms to translate."""

        for atom in atoms:
            atom._location += np.array(vector)


    @staticmethod
    def transform_atoms(matrix, *atoms):
        """Transforms multiple atoms using some matrix.

        :param matrix: the transformation matrix.
        :param \*atoms: the atoms to transform."""

        locations = [list(a) for a in atoms]
        output = np.dot(np.array(matrix), np.array(locations).transpose())
        for atom, location in zip(atoms, output.transpose()):
            atom._location = location


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


    @property
    def location(self):
        """The atom's location.

        :rtype: ``tuple``"""

        return tuple(self._location)


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
    def bonded_atoms(self):
        """Returns the atoms this atom is bonded to.

        :rtype: ``set```"""

        return self._bonded_atoms


    @property
    def mass(self):
        """The atom's molar mass according to the Periodic Table, based on the
        atom's :py:meth:`element`. If the element doesn't match any symbol on
        the Periodic Table, a mass of 0 will be returned.

        The element lookup is case-insensitive.

        :rtype: ``float``"""

        return self.__data.PERIODIC_TABLE.get(self._element.upper(), 0)
    

    @property
    def covalent_radius(self):
        """The atom's covalent radius, based on the atom's :py:meth:`element`.
        If the element doesn't match any symbol on the Periodic Table, a radius
        of 0 will be returned.

        The element lookup is case-insensitive.

        :rtype: ``float``"""

        return self.__data.COVALENT_RADII.get(self._element.upper(), 0)


    @property
    def is_metal(self):
        """Checks whether the atom's element matches a metal element.

        The element lookup is case-insensitive.

        :rtype: ``bool``"""

        return self._element.upper() in self.__data.METALS


    @property
    def is_backbone(self):
        """Returns ``True`` if the atom has a backbone atom name.

        :rtype: ``bool``"""

        return isinstance(self._het, Residue) and \
         self._name in ["CA", "C", "N", "O"]


    @property
    def is_side_chain(self):
        """Returns ``True`` if the atom has a side chain atom name.

        :rtype: ``bool``"""

        return isinstance(self._het, Residue) and not self.is_backbone


    def distance_to(self, other):
        """Returns the distance (in whatever units the coordinates are defined
        in) between this atom and another. You can also give a (x, y, z) tuple
        instead of another atom if you so wish.

        :param Atom other: The other atom (or location tuple).
        :rtype: ``float``"""

        return np.linalg.norm(self._location - np.array(list(other)))


    def angle(self, atom1, atom2):
        """Gets the angle between two atom vectors with this atom as the origin.

        :param Atom atom1: The first atom.
        :param Atom atom2: Thne second atom."""

        vectors = [
         [v1 - v2 for v1, v2 in zip(atom.location, self.location)
        ] for atom in (atom1, atom2)]
        normalized = [np.linalg.norm(v) for v in vectors]
        if 0 in normalized: return 0
        vectors = [v / n for v, n in zip(vectors, normalized)]
        return np.arccos(np.clip(np.dot(vectors[0], vectors[1]), -1.0, 1.0))
    

    def copy(self, id=None):
        """Returns a copy of the atom. The new atom will have the same element,
        location, name, charge, ID, bvalue etc. as the original, but will not
        be part of any model or other molecule.

        :rtype: ``Atom``"""

        return Atom(
         self._element, *self._location, id or self._id, self._name,
         self._charge, self._bvalue, self._anisotropy
        )
    

    @property
    def het(self):
        """Returns the :py:class:`.Residue` or :py:class:`.Ligand` the atom is
        part of, or ``None`` if it is not part of one.

        :rtype: ``Het```"""

        return self._het
    

    @property
    def chain(self):
        """Returns the :py:class:`.Chain` the atom is part of, or ``None`` if
        it is not part of one.

        :rtype: ``Chain``"""

        if self._het: return self._het.chain


    @property
    def model(self):
        """Returns the :py:class:`.Model` the atom is part of, or ``None`` if
        it is not part of one.

        :rtype: ``Model``"""

        if self.chain: return self.chain.model


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


    def nearby_hets(self, *args, residues=True, ligands=True, **kwargs):
        """Returns all residues and ligands in the associated :py:class:`.Model`
        that are within a given distance (in the units of the atom coordinates)
        of this atom. If the atom is not part of a model, no residues will be
        returned.

        :param float cutoff: the distance cutoff to use.
        :param bool residues: if ``False``, residues will not be returned.
        :param bool ligands: if ``False``, ligands will not be returned.
        :rtype: ``set``"""

        atoms = self.nearby_atoms(*args, **kwargs)
        structures = set()
        for atom in atoms:
            if atom.het is not None: structures.add(atom.het)
        try:
            structures.remove(self.het)
        except: pass
        if not residues:
            structures = {s for s in structures if not isinstance(s, Residue)}
        if not ligands:
            structures = {s for s in structures if not (isinstance(s, Ligand))}
        return structures
    

    def nearby_chains(self, *args, **kwargs):
        """Returns all chain structures in the associated :py:class:`.Model`
        that are within a given distance (in the units of the atom coordinates)
        of this atom. If the atom is not part of a model, no chains will be
        returned.

        :param float cutoff: the distance cutoff to use.
        :rtype: ``set``"""

        atoms = self.nearby_atoms(*args, **kwargs)
        chains = set()
        for atom in atoms:
            if atom.chain is not None: chains.add(atom.chain)
        try:
            chains.remove(self.chain)
        except: pass
        return chains


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

        self._location[0], self._location[1], self._location[2] = x, y, z


    def trim(self, places):
        """Rounds the coordinate values to a given number of decimal places.
        Useful for removing floating point rounding errors after transformation.

        :param int places: The number of places to round the coordinates to. If\
        ``None``, no rounding will be done."""

        if places is not None:
            self._location = np.round(self._location, places)


    def bond(self, other):
        """Bonds the atom to some other atom.

        :param Atom other: the other atom to bond to."""
        
        self._bonded_atoms.add(other)
        other._bonded_atoms.add(self)
