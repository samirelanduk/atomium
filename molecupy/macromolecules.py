from . import molecules
from .exceptions import *

class Residue(molecules.Molecule):

    def __init__(self, residue_id, residue_name, *atoms):
        molecules.Molecule.__init__(self, *atoms, molecule_id=residue_id, molecule_name=residue_name)
        self.residue_id = self.molecule_id
        self.residue_name = self.molecule_name
        del self.__dict__["molecule_id"]
        del self.__dict__["molecule_name"]
        self.downstream_residue = None
        self.upstream_residue = None
        self.chain = None


    def __repr__(self):
        return "<Residue (%s)>" % self.residue_name


    def connect_to(self, other_residue, this_atom, their_atom):
        if this_atom not in self:
            raise InvalidAtomError("Atom %s not in residue %s" % (this_atom, self))
        if their_atom not in other_residue:
            raise InvalidAtomError("Atom %s not in residue %s" % (their_atom, other_residue))
        this_atom.covalent_bond_to(their_atom)
        self.downstream_residue = other_residue
        other_residue.upstream_residue = self


    def get_downstream_residues(self):
        downstream_residues = []
        downstream_residue = self.downstream_residue
        while downstream_residue and downstream_residue is not self:
            downstream_residues.append(downstream_residue)
            downstream_residue = downstream_residue.downstream_residue
        return tuple(downstream_residues)


    def get_upstream_residues(self):
        upstream_residues = []
        upstream_residue = self.upstream_residue
        while upstream_residue and upstream_residue is not self:
            upstream_residues.append(upstream_residue)
            upstream_residue = upstream_residue.upstream_residue
        return tuple(upstream_residues)


    def get_accessible_residues(self):
        return set(self.get_upstream_residues() + self.get_downstream_residues())




class ResiduicStructure(molecules.AtomicStructure):

    def __init__(self, *residues):
        if not all(isinstance(residue, Residue) for residue in residues):
            non_residues = [residue for residue in residues if not isinstance(residue, Residue)]
            raise TypeError("ResiduicStructure needs residues, not '%s'" % non_residues[0])
        if not residues:
            raise NoResiduesError("Cannot make ResiduicStructure with zero residues")
        residue_ids = [residue.residue_id for residue in residues if residue.residue_id is not None]
        if len(set(residue_ids)) < len(residue_ids):
            raise DuplicateResidueIdError("Cannot make ResiduicStructure with duplicate residue_ids")
        self.residues = set(residues)
        all_atoms = set()
        for residue in self.residues:
            all_atoms.update(residue.atoms)
        molecules.AtomicStructure.__init__(self, *all_atoms)


    def __repr__(self):
        return "<ResiduicStructure (%i residues)>" % len(self.residues)


    def __contains__(self, residue):
        return residue in self.residues



class Chain(ResiduicStructure, molecules.Molecule):

    def __init__(self, *residues, chain_id=None):
        ResiduicStructure.__init__(self, *residues)
        if len(residues) > 1 and set(residues[1:]) != residues[0].get_accessible_residues():
            raise BrokenChainError("Cannot make Chain with unconnected residues")
        self.residues = tuple(residues)
        all_atoms = set()
        for residue in self.residues:
            residue.chain = self
            all_atoms.update(residue.atoms)
        molecules.Molecule.__init__(self, *all_atoms)

        if not isinstance(chain_id, str) and chain_id is not None:
            raise TypeError("'%s' is not a valid chain_id" % str(chain_id))
        self.chain_id = chain_id


    def __repr__(self):
        return "<Chain (%i residues)>" % len(self.residues)


    def __getitem__(self, key):
        return self.residues[key]



class MacroModel(molecules.Model):

    def __init__(self):
        molecules.Model.__init__(self)
        self._chains = set()
        self._small_molecules = set()


    def __repr__(self):
        return "<MacroModel (%i atoms)>" % len(self.atoms)


    def add_chain(self, chain):
        if not isinstance(chain, Chain):
            raise TypeError(
             "Only chains can be added to a model with add_chain, not '%s'"
              % str(chain)
            )
        self._chains.add(chain)
        self.add_molecule(chain)


    def add_small_molecule(self, small_molecule):
        self._small_molecules.add(small_molecule)
        self.add_molecule(small_molecule)


    def get_chains(self):
        return self._chains


    def get_small_molecules(self):
        return self._small_molecules
