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


    def __contains__(self, obj):
        return obj in self.residues or obj in self.atoms


    def get_residue_by_name(self, residue_name):
        if not isinstance(residue_name, str):
            raise TypeError("Residue name search must be by str, not '%s'" % str(residue_name))
        for residue in self.residues:
            if residue.residue_name == residue_name:
                return residue


    def get_residues_by_name(self, residue_name):
        if not isinstance(residue_name, str):
            raise TypeError("Residue name search must be by str, not '%s'" % str(residue_name))
        return set([
         residue for residue in self.residues if residue.residue_name == residue_name
        ])



class Chain(ResiduicStructure, molecules.Molecule):

    def __init__(self, *residues, chain_id=None):
        ResiduicStructure.__init__(self, *residues)
        if len(residues) > 1 and not set(residues[1:]).issubset(residues[0].get_accessible_residues()):
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

        self.model = None
        self.complex = None


    def __repr__(self):
        return "<Chain (%i residues)>" % len(self.residues)


    def __getitem__(self, key):
        return self.residues[key]


    def get_residue_by_id(self, residue_id):
        if not isinstance(residue_id, int):
            raise TypeError("Residue ID search must be by int, not '%s'" % str(residue_id))
        for residue in self.residues:
            if residue.residue_id == residue_id:
                return residue



class Complex(molecules.AtomicStructure):

    def __init__(self, *chains, complex_id=None, complex_name=None):
        if not all(isinstance(chain, Chain) for chain in chains):
            non_chains = [chain for chain in chains if not isinstance(chain, Chain)]
            raise TypeError("Complex needs chains, not '%s'" % non_chains[0])
        if not chains:
            raise NoChainsError("Cannot make Complex with zero chains")
        self.chains = set(chains)
        self.model = None
        all_atoms = set()
        for chain in self.chains:
            chain.complex = self
            all_atoms.update(chain.atoms)
        molecules.AtomicStructure.__init__(self, *all_atoms)

        if not isinstance(complex_id, int) and complex_id is not None:
            raise TypeError("'%s' is not a valid complex_id" % str(complex_id))
        self.complex_id = complex_id

        if not isinstance(complex_name, str) and complex_name is not None:
            raise TypeError("'%s' is not a valid complex_name" % str(complex_name))
        self.complex_name = complex_name


    def __repr__(self):
        return "<Complex (%i chains)>" % len(self.chains)



class Site(ResiduicStructure):

    def __init__(self, *residues, site_id=None, site_name=None):
        ResiduicStructure.__init__(self, *residues)

        if not isinstance(site_id, int) and site_id is not None:
            raise TypeError("'%s' is not a valid site_id" % str(site_id))
        self.site_id = site_id

        if not isinstance(site_name, str) and site_name is not None:
            raise TypeError("'%s' is not a valid site_name" % str(site_name))
        self.site_name = site_name

        self.model = None
        self.ligand = None


    def __repr__(self):
        return "<Site (%i residues)>" % len(self.residues)



class MacroModel(molecules.Model):

    def __init__(self):
        molecules.Model.__init__(self)
        self._chains = set()
        self._small_molecules = set()
        self._complexes = set()
        self._sites = set()


    def __repr__(self):
        return "<MacroModel (%i atoms)>" % len(self.atoms)


    def __contains__(self, obj):
        return obj in self._molecules or obj in self._complexes or obj in self._sites


    def __getattr__(self, attribute):
        if attribute == "atoms":
            atoms = set()
            for molecule in self._molecules:
                atoms.update(molecule.atoms)
            for site in self._sites:
                atoms.update(site.atoms)
            return atoms
        else:
            return self.__getattribute__(attribute)


    def add_small_molecule(self, small_molecule):
        self._small_molecules.add(small_molecule)
        self.add_molecule(small_molecule)


    def add_chain(self, chain):
        if not isinstance(chain, Chain):
            raise TypeError(
             "Only chains can be added to a model with add_chain, not '%s'"
              % str(chain)
            )
        if chain.chain_id is not None:
            existing_chain_ids = [
             chain.chain_id for chain in self._chains
            ]
            if chain.chain_id in existing_chain_ids:
                raise DuplicateChainError(
                 "Cannot have two chains in a model with ID %s" % chain.chain_id
                )
        self._chains.add(chain)
        self.add_molecule(chain)


    def add_complex(self, complex_):
        if not isinstance(complex_, Complex):
            raise TypeError(
             "Only complexes can be added to a model with add_complex, not '%s'"
              % str(complex_)
            )
        if complex_.complex_id is not None:
            existing_complex_ids = [
             complex_.complex_id for complex_ in self._complexes
            ]
            if complex_.complex_id in existing_complex_ids:
                raise DuplicateComplexError(
                 "Cannot have two complexes in a model with ID %s" % complex_.complex_id
                )
        self._complexes.add(complex_)
        complex_.model = self
        for chain in complex_.chains:
            self.add_chain(chain)


    def add_site(self, site):
        if not isinstance(site, Site):
            raise TypeError(
             "Only sites can be added to a model with add_site, not '%s'"
              % str(site)
            )
        if site.site_id is not None:
            existing_site_ids = [
             site.site_id for site in self._sites
            ]
            if site.site_id in existing_site_ids:
                raise DuplicateSiteError(
                 "Cannot have two sites in a model with ID %s" % site.site_id
                )
        self._sites.add(site)
        site.model = self


    def get_small_molecules(self):
        return self._small_molecules


    def get_small_molecule_by_id(self, molecule_id):
        if not isinstance(molecule_id, int):
            raise TypeError(
             "Small Molecule ID search must be by int, not '%s'" % str(molecule_id)
            )
        for molecule in self._small_molecules:
            if molecule.molecule_id == molecule_id:
                return molecule


    def get_small_molecule_by_name(self, small_molecule_name):
        if not isinstance(small_molecule_name, str):
            raise TypeError("Small molecule name search must be by str, not '%s'" % str(small_molecule_name))
        for small_molecule in self._small_molecules:
            if small_molecule.molecule_name == small_molecule_name:
                return small_molecule


    def get_small_molecules_by_name(self, small_molecule_name):
        if not isinstance(small_molecule_name, str):
            raise TypeError("Small molecule name search must be by str, not '%s'" % str(small_molecule_name))
        return set([
         small_molecule for small_molecule in self._small_molecules
          if small_molecule.molecule_name == small_molecule_name
        ])


    def get_chains(self):
        return self._chains


    def get_chain_by_id(self, chain_id):
        if not isinstance(chain_id, str):
            raise TypeError(
             "Chain ID search must be by str, not '%s'" % str(chain_id)
            )
        for chain in self._chains:
            if chain.chain_id == chain_id:
                return chain


    def get_complexes(self):
        return self._complexes


    def get_sites(self):
        return self._sites


    def get_site_by_id(self, site_id):
        if not isinstance(site_id, int):
            raise TypeError(
             "Site ID search must be by int, not '%s'" % str(site_id)
            )
        for site in self._sites:
            if site.site_id == site_id:
                return site


    def get_site_by_name(self, site_name):
        if not isinstance(site_name, str):
            raise TypeError("Site name search must be by str, not '%s'" % str(site_name))
        for site in self._sites:
            if site.site_name == site_name:
                return site


    def get_sites_by_name(self, site_name):
        if not isinstance(site_name, str):
            raise TypeError("Site name search must be by str, not '%s'" % str(site_name))
        return set([
         site for site in self._sites
          if site.site_name == site_name
        ])
