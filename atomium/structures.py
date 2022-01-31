import math
import numpy as np
from scipy.spatial.distance import cdist
from collections import Counter, defaultdict
from .search import StructureSet, StructureClass, query
from .data import CODES, PERIODIC_TABLE, COVALENT_RADII, METALS, FULL_NAMES

class AtomStructure:
    
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
        locations = np.array([a.location * a.mass for a in self.atoms()])
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


    def create_grid(self, size=1, margin=0):
        """A generator which models a grid around the structure and returns the
        coordinates of all the points in that grid.

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
            while values[0] > min_: values.insert(0, round(values[0] - size, 12))
            while values[-1] < max_: values.append(round(values[-1] + size, 12))
            values = [
                v for v in values
                if (v > min_ or abs(v - min_) < size) and
                (v < max_ or abs(v - max_) < size)
            ]
            dimension_values.append(values)
        for x in dimension_values[0]:
            for y in dimension_values[1]:
                for z in dimension_values[2]:
                    yield (x, y, z)
    

    def pairwise_atoms(self, *args, **kwargs):
        """A generator which yields all the pairwise atom combinations of the
        structure. There will be no duplicates in the returned generator, and
        the number of returned pairs will be a triangle number.

        :rtype: ``tuple``"""

        atoms = list(self.atoms(*args, **kwargs))
        for a_index in range(len(atoms) - 1):
            for o_index in range(a_index + 1, len(atoms)):
                yield {atoms[a_index], atoms[o_index]}
    

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
            atoms = StructureSet(atoms)
            atoms = query(lambda self: atoms)(self, *args, **kwargs)
        else:
            atoms = self.atoms(*args, **kwargs)
        if not atoms: return atoms
        X = np.tile(location, [len(atoms), 1])
        Y = np.array([a.location for a in atoms])
        distances = cdist(X, Y)[0]
        return {a for index, a in enumerate(atoms) if distances[index] <= radius}
    

    def nearby_residues(self, *args, **kwargs):
        residues = set()
        for atom in self.atoms():
            residues.update(atom.nearby_residues(*args, **kwargs)) 
        return residues
    

    def nearby_molecules(self, *args, **kwargs):
        molecules = set()
        for atom in self.atoms():
            molecules.update(atom.nearby_molecules(*args, **kwargs)) 
        return molecules
    

    def pairing_with(self, structure, *args, **kwargs):
        self_atoms = list(self.atoms(*args, **kwargs))
        other_atoms = list(structure.atoms(*args, **kwargs))
        if len(self_atoms) != len(other_atoms):
            raise ValueError(f"{self} and {structure} have different numbers of atoms")
        for struct, atoms in ((self, self_atoms), (structure, other_atoms)):
            center = struct.center_of_mass
            atoms.sort(key=lambda a: (
                a.element, a.name, a.id, a.distance_to(center), id(a)
            ))
        return {atom1: atom2 for atom1, atom2 in zip(self_atoms, other_atoms)}
    

    def rmsd_with(self, structure, *args, **kwargs):
        pairing = self.pairing_with(structure, *args, **kwargs)
        coords1, coords2 = [[
            a.location for a in atoms
        ] for atoms in zip(*pairing.items())]
        coords1 = coords1 - np.array(self.center_of_mass)
        coords2 = coords2 - np.array(structure.center_of_mass)
        cov_matrix = np.dot(np.transpose(coords1), coords2)
        V, S, W = np.linalg.svd(cov_matrix)
        d = (np.linalg.det(V) * np.linalg.det(W)) < 0.0
        if d: S[-1] = -S[-1]
        if d: V[:, -1] = -V[:, -1]
        rot_matrix = np.dot(V, W)
        coords1_rotated = np.dot(coords1, rot_matrix)
        diff = coords1_rotated - coords2
        N = len(coords2)
        return np.sqrt((diff * diff).sum() / N)
    

    def trim(self, places):
        """Rounds the coordinate values to a given number of decimal places.
        Useful for removing floating point rounding errors after transformation.

        :param int places: The number of places to round the coordinates to."""

        for atom in self.atoms():
            atom.trim(places)



class Model(AtomStructure, metaclass=StructureClass):

    def __init__(self, molecules, file=None):
        for molecule in molecules: molecule.model = self
        self._molecules = StructureSet(molecules)
        self.file = file
    

    def __repr__(self):
        lists = [[
            m for m in self._molecules.structures if isinstance(m, Type)
        ] for Type in [Polymer, BranchedPolymer, NonPolymer]]
        strings = []
        if lists[0]:
            strings.append("{} polymer{}".format(len(lists[0]), "" if len(lists[0]) == 1 else "s"))
        if lists[1]:
            strings.append("{} branched-polymer{}".format(len(lists[1]), "" if len(lists[1]) == 1 else "s"))
        if lists[2]:
            strings.append("{} non-polymer{}".format(len(lists[2]), "" if len(lists[2]) == 1 else "s"))
        return "<Model ({})>".format(", ".join(strings))
    

    def __contains__(self, obj):
        return obj in self._molecules.structures or\
            obj in self.residues() or obj in self.atoms()
    

    def molecules(self):
        return self._molecules
    

    def polymers(self):
        return StructureSet([m for m in self._molecules.structures if isinstance(m, Polymer)])
    

    def branched_polymers(self):
        return StructureSet([m for m in self._molecules.structures if isinstance(m, BranchedPolymer)])
    

    def non_polymers(self):
        return StructureSet([m for m in self._molecules.structures if isinstance(m, NonPolymer)])
    
    
    def waters(self):
        return StructureSet([m for m in self._molecules.structures if isinstance(m, Water)])
    

    def residues(self):
        residues = []
        for chain in self.polymers():
            residues += chain.residues()
        for chain in self.branched_polymers():
            residues += chain.residues()
        return StructureSet(residues)
    

    def atoms(self):
        atoms = []
        for molecule in self.molecules():
            atoms += molecule.atoms()
        return StructureSet(atoms)
    

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
    

    def remove(self, structure):
        self._molecules.remove(structure)
        if isinstance(structure, (Residue, Atom)):
            for molecule in self._molecules.structures:
                molecule.remove(structure)




class Entity(AtomStructure, metaclass=StructureClass):
    
    def __init__(self, id, auth_id):
        self.id = id
        self.auth_id = auth_id
        self.model = None



class Polymer(Entity):

    def __init__(self, residues, helices, strands, *args, **kwargs):
        for residue in residues: residue.polymer = self
        self._residues = StructureSet(residues)
        self._helices, self._strands = helices, strands
        Entity.__init__(self, *args, **kwargs)

    
    def __repr__(self):
        return f"<Polymer {self.id} ({len(self._residues)} residue{'' if len(self._residues) == 1 else 's'})>"


    def __len__(self):
        return len(self._residues)
    

    def __iter__(self):
        return iter(self._residues.structures)


    def __getitem__(self, key):
        return self.residues()[key]


    def __contains__(self, obj):
        return obj in self._residues.structures or obj in self.atoms()


    def residues(self):
        return self._residues
    

    def atoms(self):
        atoms = []
        for residue in self.residues():
            atoms += residue.atoms()
        return StructureSet(atoms)
    

    @property
    def helices(self):
        residues = self.residues()
        ids = [r.id for r in residues]
        indices = [[
            ids.index(start), ids.index(end)
        ] for start, end in self._helices if start in ids and end in ids]
        return [residues[start:end + 1] for start, end in indices]
    

    @property
    def strands(self):
        residues = self.residues()
        ids = [r.id for r in residues]
        indices = [[
            ids.index(start), ids.index(end)
        ] for start, end in self._strands if start in ids and end in ids]
        return [residues[start:end + 1] for start, end in indices]
    

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
    

    def remove(self, structure):
        if isinstance(structure, Atom):
            for residue in self._residues.structures:
                residue.remove(structure)
        self._residues.remove(structure)
    


class BranchedPolymer(Entity):

    def __init__(self, residues, *args, **kwargs):
        for residue in residues: residue.branched_polymer = self
        self._residues = StructureSet(residues)
        Entity.__init__(self, *args, **kwargs)

    
    def __repr__(self):
        return f"<BranchedPolymer {self.id} ({len(self._residues)} residue{'' if len(self._residues) == 1 else 's'})>"


    def __contains__(self, obj):
        return obj in self._residues.structures or obj in self.atoms()

    
    def residues(self):
        return self._residues
    

    def atoms(self):
        atoms = []
        for residue in self.residues():
            atoms += residue.atoms()
        return StructureSet(atoms)
    

    def remove(self, structure):
        if isinstance(structure, Atom):
            for residue in self._residues.structures:
                residue.remove(structure)
        self._residues.remove(structure)



class NonPolymer(Entity):

    def __init__(self, name, atoms, *args, **kwargs):
        for atom in atoms: atom.non_polymer = self
        self._atoms = StructureSet(atoms)
        self.name = name
        Entity.__init__(self, *args, **kwargs)

    
    def __repr__(self):
        return f"<NonPolymer {self.id} ({self.name})>"
    

    def __contains__(self, atom):
        return atom in self._atoms.structures
    

    def atoms(self):
        return self._atoms
    

    def remove(self, atom):
        self._atoms.remove(atom)



class Water(Entity):

    def __init__(self, name, atoms, *args, **kwargs):
        for atom in atoms: atom.water = self
        self._atoms = StructureSet(atoms)
        self.name = name
        Entity.__init__(self, *args, **kwargs)

    
    def __repr__(self):
        return f"<Water {self.id} ({self.name})>"
    

    def __contains__(self, atom):
        return atom in self._atoms.structures
    

    def atoms(self):
        return self._atoms
    

    def remove(self, atom):
        self._atoms.remove(atom)



class Residue(AtomStructure, metaclass=StructureClass):

    def __init__(self, id, name, number, atoms):
        for atom in atoms: atom.residue = self
        self.id = id
        self.name = name
        self.number = number
        self._atoms = StructureSet(atoms)
        self.polymer, self.branched_polymer = None, None
    

    def __repr__(self):
        return f"<Residue {self.name} ({self.id})>"
    

    def __contains__(self, atom):
        return atom in self._atoms.structures
    

    def atoms(self):
        return self._atoms
    

    @property
    def model(self):
        if self.polymer: return self.polymer.model
        if self.branched_polymer: return self.branched_polymer.model
    

    @property
    def code(self):
        """Returns the single letter code, based on its three letter name - or
        just 'X' if it doesn't match anything.

        :rtype: ``str``"""

        return CODES.get(self.name, "X")
    

    @property
    def full_name(self):
        """Returns the residue's full name, based on its three letter name - or
        just the three letter name if it doesn't match anything. Or you can just
        supply a full name when you instantiate the Het.

        :rtype: ``str``"""

        return FULL_NAMES.get(self.name, self.name)
    

    @property
    def next(self):
        if not self.polymer: return None
        residues = self.polymer.residues()
        index = residues.index(self)
        if index == len(residues) - 1: return None
        return residues[index + 1]
    

    @property
    def previous(self):
        if not self.polymer: return None
        residues = self.polymer.residues()
        index = residues.index(self)
        if index == 0: return None
        return residues[index - 1]
    

    @property
    def in_helix(self):
        if not self.polymer: return None
        helices = self.polymer.helices
        for helix in helices:
            for residue in helix:
                if residue is self: return True
        return False
    

    @property
    def in_strand(self):
        if not self.polymer: return None
        strands = self.polymer.strands
        for strand in strands:
            for residue in strand:
                if residue is self: return True
        return False
    

    def remove(self, atom):
        self._atoms.remove(atom)



class Atom:
    """Represents a single atom - neutral or charged.
    
    Atoms have attributes for the residue, non-polymer or water they belong to,
    and dynamic properties for the polymer, branched-polymer and model they
    belong to.
    
    :param str element - the atom's element.
    :param float x - the atom's element.
    :param float y - the atom's element.
    :param float z - the atom's element.
    :param int id - the atom's ID.
    :param str name - the atom's name.
    :param float charge - the atom's charge.
    :param float bvalue - the atom's isotropic uncertainty.
    :param tuple anisotropy - the atom's anisotropic uncertainty.
    :param Residue residue - the atom's residue container.
    :param NonPolymer non_polymer - the atom's non-polymer container.
    :param Water water - the atom's water container."""

    __slots__ = [
        "element", "location", "id", "name", "charge",
        "bvalue", "anisotropy", "residue", "non_polymer", "water"
    ]

    def __init__(self, element, x, y, z, id, name, charge=0, bvalue=0, anisotropy=None):
        self.location = np.array([x, y, z])
        self.element, self.id, self.name = element, id, name
        self.charge, self.bvalue, self.anisotropy = charge, bvalue, anisotropy
        self.residue, self.non_polymer, self.water = None, None, None
    

    def __repr__(self):
        return f"<Atom {self.id} ({self.name})>"
    

    def __iter__(self):
        return iter(self.location)
    

    @property
    def polymer(self):
        """The atom's polymer."""

        if self.residue: return self.residue.polymer
    

    @property
    def branched_polymer(self):
        """The atom's branched-polymer."""

        if self.residue: return self.residue.branched_polymer
        

    @property
    def model(self):
        """The model the atom is in."""

        if self.residue: return self.residue.model
        if self.non_polymer: return self.non_polymer.model
        if self.water: return self.water.model
    

    @property
    def mass(self):
        """The atom's molar mass according to the Periodic Table, based on the
        atom's :py:meth:`element`. If the element doesn't match any symbol on
        the Periodic Table, a mass of 0 will be returned.

        The element lookup is case-insensitive.

        :rtype: ``float``"""

        return PERIODIC_TABLE.get(self.element.upper(), 0)
    

    @property
    def covalent_radius(self):
        """The atom's covalent radius, based on the atom's :py:meth:`element`.
        If the element doesn't match any symbol on the Periodic Table, a radius
        of 0 will be returned.

        The element lookup is case-insensitive.

        :rtype: ``float``"""

        return COVALENT_RADII.get(self.element.upper(), 0)


    @property
    def is_metal(self):
        """Checks whether the atom's element matches a metal element.

        The element lookup is case-insensitive.

        :rtype: ``bool``"""

        return self.element.upper() in METALS


    @property
    def is_backbone(self):
        """Returns ``True`` if the atom has a backbone atom name.

        :rtype: ``bool``"""

        if not self.residue: return False
        return self.name in ["CA", "C", "N", "O"]


    @property
    def is_side_chain(self):
        """Returns ``True`` if the atom has a side chain atom name.

        :rtype: ``bool``"""

        if not self.residue: return False
        return not self.is_backbone
    

    def distance_to(self, other):
        """Returns the distance (in whatever units the coordinates are defined
        in) between this atom and another. You can also give a (x, y, z) tuple
        instead of another atom if you so wish.

        :param Atom other: The other atom (or location tuple).
        :rtype: ``float``"""

        return np.linalg.norm(self.location - np.array(list(other)))
    

    def angle(self, atom1, atom2):
        """Gets the angle between two atom vectors with this atom as the origin.

        :param Atom atom1: The first atom.
        :param Atom atom2: Thne second atom."""

        vectors = [
            [v1 - v2 for v1, v2 in zip(list(atom), self.location)
        ] for atom in (atom1, atom2)]
        normalized = [np.linalg.norm(v) for v in vectors]
        if 0 in normalized: return 0
        vectors = [v / n for v, n in zip(vectors, normalized)]
        return np.arccos(np.clip(np.dot(vectors[0], vectors[1]), -1.0, 1.0))
    

    def nearby_atoms(self, cutoff, *args, **kwargs):
        """Returns all atoms in the associated :py:class:`.Model` that are
        within a given distance (in the units of the atom coordinates) of this
        atom. If the atom is not part of a model, no atoms will be returned.

        :param float cutoff: The radius to search within.
        :rtype: ``set``"""

        if self.model:
            atoms = self.model.atoms_in_sphere(
                self.location, cutoff, *args, **kwargs
            )
            try:
                atoms.remove(self)
            except: pass
            return atoms
        return set()
    

    def nearby_residues(self, cutoff, *args, **kwargs):
        atoms = self.nearby_atoms(cutoff, *args, **kwargs)
        residues = set()
        for atom in atoms:
            if atom.residue and atom.residue is not self.residue:
                residues.add(atom.residue)
        return residues
    

    def nearby_molecules(self, cutoff, *args, **kwargs):
        atoms = self.nearby_atoms(cutoff, *args, **kwargs)
        molecules = set()
        for atom in atoms:
            if atom.polymer and atom.polymer is not self.polymer:
                molecules.add(atom.polymer)
            if atom.branched_polymer and atom.branched_polymer is not self.branched_polymer:
                molecules.add(atom.branched_polymer)
            if atom.non_polymer and atom.non_polymer is not self.non_polymer:
                molecules.add(atom.non_polymer)
            if atom.water and atom.water is not self.water:
                molecules.add(atom.water)
        return molecules
    

    def trim(self, places):
        """Rounds the coordinate values to a given number of decimal places.
        Useful for removing floating point rounding errors after transformation.

        :param int places: The number of places to round the coordinates to."""

        self.location = np.round(self.location, places)