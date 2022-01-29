import numpy as np
from .search import StructureSet, StructureClass
from .data import CODES, PERIODIC_TABLE, COVALENT_RADII, METALS, FULL_NAMES

class Model(metaclass=StructureClass):

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




class Entity(metaclass=StructureClass):
    
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



class Residue(metaclass=StructureClass):

    def __init__(self, id, name, atoms):
        for atom in atoms: atom.residue = self
        self.id = id
        self.name = name
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