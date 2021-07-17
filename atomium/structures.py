import numpy as np
from .search import StructureSet, StructureClass
class Model(metaclass=StructureClass):

    def __init__(self, chains=None, ligands=None, waters=None):
        for mol in chains + ligands + waters:
            mol.model = self
        self._chains = StructureSet(chains)
        self._ligands = StructureSet(ligands)
        self._waters = StructureSet(waters)
    

    def __repr__(self):
        chains = "{} chains".format(len(self._chains))
        if len(self._chains) == 1: chains = chains[:-1]
        ligands = "{} ligands".format(len(self._ligands))
        if len(self._ligands) == 1: ligands = ligands[:-1]
        return "<Model ({}, {})>".format(chains, ligands)
    

    def chains(self):
        return self._chains

    
    def ligands(self):
        return self._ligands


    def waters(self):
        return self._waters
    

    def molecules(self):
        """Returns all of the model's molecules (chains, ligands, waters).

        :rtype: ``set``"""

        return self._chains + self._ligands + self._waters
    

    def residues(self):
        """Returns all of the model's residues in all its chains.

        :rtype: ``set``"""

        residues = []
        for chain in self._chains.structures:
            residues += chain.residues()
        return StructureSet(residues)
    

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
        return StructureSet(atoms)
    

    def dehydrate(self):
        """Removes all water ligands from the model."""

        self._waters = StructureSet([])



class Chain(metaclass=StructureClass):

    def __init__(self, *residues, id="", asym_id="", auth_asym_id="", sequence=None, helices=None, strands=None):
        for res in residues: res.chain = self
        self._residues = StructureSet(residues)
        self.id, self.asym_id, self.auth_asym_id = id, asym_id, auth_asym_id
        self.sequence = sequence
        self._helices = helices or []
        self._strands = strands or []
        self.model = None
    

    def __repr__(self):
        return f"<Chain {self.id} ({len(self._residues)} residues)>"
    

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
    

    def ligands(self):
        if not self.model: return None
        return StructureSet(self.model.ligands(auth_asym_id=self.auth_asym_id))
    

    def waters(self):
        if not self.model: return None
        return StructureSet(self.model.waters(auth_asym_id=self.auth_asym_id))
    

    def atoms(self):
        """Returns all the atoms in with the chain.

        :rtype: ``set``"""

        atoms = set()
        for res in self._residues.structures:
            atoms.update(res._atoms.structures)
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



class Residue(metaclass=StructureClass):

    def __init__(self, *atoms, id="", asym_id="", auth_asym_id="", name=""):
        for atom in atoms: atom.residue = self
        self._atoms = StructureSet(atoms)
        self.id, self.name = id, name
        self.asym_id, self.auth_asym_id = asym_id, auth_asym_id
        self.chain = None
    

    def __repr__(self):
        return f"<Residue {self.name} ({self.id})>"
    

    def __contains__(self, atom):
        return atom in self._atoms.structures

    
    def atoms(self):
        return self._atoms
    

    @property
    def model(self):
        """Returns the :py:class:`.Model` the residue is part of, via its
        chain.

        :rtype: ``Model``"""

        try:
            return self.chain.model
        except AttributeError: return None



class Ligand(metaclass=StructureClass):

    def __init__(self, *atoms, id="", asym_id="", auth_asym_id="", name="", water=False):
        for atom in atoms: atom.ligand = self
        self._atoms = StructureSet(atoms)
        self.id, self.name, self.is_water = id, name, water
        self.asym_id, self.auth_asym_id = asym_id, auth_asym_id
        self.model = None
    

    def __repr__(self):
        return f"<{'Water' if self.is_water else 'Ligand'} {self.name} ({self.id})>"
    

    def __contains__(self, atom):
        return atom in self._atoms.structures
    

    @property
    def chain(self):
        if not self.model: return None
        return self.model.chain(self.auth_asym_id)
    

    def atoms(self):
        return self._atoms



class Atom:

    def __init__(self, element, x, y, z, id, name, charge=0, bvalue=0, anisotropy=None):
        self.location = np.array([x, y, z])
        self.element, self.id, self.name = element, id, name
        self.charge, self.bvalue, self.anisotropy = charge, bvalue, anisotropy
        self.residue, self.ligand = None, None
    

    def __repr__(self):
        return f"<Atom {self.id} ({self.name})>"
    

    def __iter__(self):
        return iter(self.location)
    

    @property
    def chain(self):
        if self.residue: return self.residue.chain
        if self.ligand: return self.ligand.chain
    

    @property
    def model(self):
        if self.chain: return self.chain.model
        if self.ligand: return self.ligand.chain