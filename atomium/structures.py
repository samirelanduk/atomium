import numpy as np
from .search import StructureSet, StructureClass

class Model(metaclass=StructureClass):

    def __init__(self, molecules):
        for molecule in molecules: molecule.model = self
        self._molecules = StructureSet(molecules)
    

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



class Atom:

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
        if self.residue: return self.residue.polymer
    

    @property
    def branched_polymer(self):
        if self.residue: return self.residue.branched_polymer
        

    @property
    def model(self):
        if self.residue: return self.residue.model
        if self.non_polymer: return self.non_polymer.model
        if self.water: return self.water.model