"""Contains the Model class."""

import re
from .atoms import Atom, GhostAtom
from .molecules import AtomicStructure, SmallMolecule, Residue
from .chains import Chain, BindSite
from .complexes import Complex
from ..exceptions import DuplicateSmallMoleculesError, DuplicateChainsError
from ..exceptions import DuplicateBindSitesError, DuplicateComplexesError

class Model(AtomicStructure):
    """Base class: :py:class:`.AtomicStructure`

    Represents the structural environment in which the other structures exist."""

    def __init__(self):
        self._small_molecules = set()
        self._chains = set()
        self._bind_sites = set()
        self._complexes = set()
        self._source = None


    def __getattr__(self, attribute):
        if attribute == "_atoms":
            atoms = set()
            for molecule in self.small_molecules():
                atoms.update(molecule.atoms(atom_type="all"))
            for chain in self.chains():
                atoms.update(chain.atoms(atom_type="all"))
            for site in self.bind_sites():
                atoms.update(site.atoms(atom_type="all"))
            return atoms
        else:
            return self.__getattribute__(attribute)


    def source(self):
        """The object the Model was created from."""

        return self._source


    def small_molecules(self):
        """Returns all the :py:class:`.SmallMolecule` objects in this model.

        :rtype: ``set``"""

        return set(self._small_molecules)


    def add_small_molecule(self, small_molecule):
        """Adds a small molecule to the model.

        :param SmallMolecule small_molecule: The small molecule to add."""

        if not isinstance(small_molecule, SmallMolecule):
            raise TypeError(
             "Can only add SmallMolecule to Model, not '%s'" % str(small_molecule)
            )
        if small_molecule.molecule_id() in [mol.molecule_id() for mol in self.small_molecules()]:
             raise DuplicateSmallMoleculesError(
              "Cannot add small_molecule with ID %s to %s as there is already a small_molecule with that ID" % (
               small_molecule.molecule_id(), str(self)
              )
             )
        self._small_molecules.add(small_molecule)
        small_molecule._model = self


    def remove_small_molecule(self, small_molecule):
        """Removes a small molecule from the structure.

        :param SmallMolecule small_molecule: The small molecule to remove."""

        self._small_molecules.remove(small_molecule)
        small_molecule._model = None


    def get_small_molecule_by_id(self, molecule_id):
        """Returns the first small molecule that matches a given molecule ID.

        :param str molecule_id: The molecule ID to search by.
        :rtype: :py:class:`.SmallMolecule` or ``None``"""

        if not isinstance(molecule_id, str):
            raise TypeError(
             "Small molecule ID search must be by str, not '%s'" % str(molecule_id)
            )
        for molecule in self.small_molecules():
            if molecule.molecule_id() == molecule_id:
                return molecule


    def get_small_molecule_by_name(self, molecule_name):
        """Returns the first small molecules that matches a given name.

        :param str molecule_name: The name to search by.
        :rtype: :py:class:`.SmallMolecule` or ``None``"""

        if not isinstance(molecule_name, str):
            raise TypeError(
             "Small molecule name search must be by str, not '%s'" % str(molecule_name)
            )
        for molecule in self.small_molecules():
            if molecule.molecule_name() == molecule_name:
                return molecule


    def get_small_molecules_by_name(self, molecule_name):
        """Returns all the small molecules of a given name.

        :param str molecule_name: The name to search by.
        :rtype: ``set`` of :py:class:`.SmallMolecule` objects."""

        if not isinstance(molecule_name, str):
            raise TypeError(
             "Small molecule name search must be by str, not '%s'" % str(molecule_name)
            )
        return set([molecule for molecule in self.small_molecules()
         if molecule.molecule_name() == molecule_name])


    def duplicate_small_molecule(self, small_molecule, molecule_id=None):
        """Creates a copy of a small molecule in the Model. The coordinates will
        be identical but it will have a unique ID.

        :param SmalllMolecule small_molecule: The molecule to duplicate.
        :param str molecule_id: If given, this will determine the ID of the new\
        molecule."""

        if not isinstance(small_molecule, SmallMolecule):
            raise TypeError(
             "Can only duplicate SmallMolecule with this method, not '%s'" % str(
              small_molecule
             )
            )
        if small_molecule not in self.small_molecules():
            raise ValueError(
             "%s is not in this Model and so cannot be duplicated." % str(
              small_molecule
             )
            )
        current_molecule_ids = [
         mol.molecule_id() for mol in self.small_molecules()
        ]
        new_molecule = None
        new_atoms = set()
        next_id = sorted([atom.atom_id() for atom in self.atoms(atom_type="all")])[-1] + 1
        for atom in small_molecule.atoms(atom_type="all"):
            if isinstance(atom, Atom):
                new_atoms.add(Atom(
                 atom.x(), atom.y(), atom.z(), atom.element(), next_id, atom.atom_name()
                ))
            else:
                new_atoms.add(GhostAtom(atom.element(), next_id, atom.atom_name()))
            next_id += 1
        if molecule_id:
            if molecule_id in current_molecule_ids:
                raise ValueError(
                 "There is already a SmallMolecule with ID %s" % molecule_id
                )
            new_molecule = SmallMolecule(
             molecule_id, small_molecule.molecule_name(), *new_atoms
            )
        else:
            chain, residue = small_molecule.molecule_id()[0], int(small_molecule.molecule_id()[1:])
            id_ = residue
            while "%s%i" % (chain, id_) in current_molecule_ids:
                id_ += 1
            new_molecule = SmallMolecule(
             "%s%i" % (chain, id_),
             small_molecule.molecule_name(),
             *new_atoms
            )
        self.add_small_molecule(new_molecule)
        return new_molecule


    def chains(self):
        """Returns all the :py:class:`.Chain` objects in this model.

        :rtype: ``set``"""

        return set(self._chains)


    def add_chain(self, chain):
        """Adds a chain to the model.

        :param Chain chain: The chain to add."""

        if not isinstance(chain, Chain):
            raise TypeError(
             "Can only add Chain to Model, not '%s'" % str(chain)
            )
        if chain.chain_id() in [chain.chain_id() for chain in self.chains()]:
             raise DuplicateChainsError(
              "Cannot add chain with ID %s to %s as there is already a chain with that ID" % (
               chain.chain_id(), str(self)
              )
             )
        self._chains.add(chain)
        chain._model = self


    def remove_chain(self, chain):
        """Removes a chain from the structure.

        :param Chain chain: The chain to remove."""

        self._chains.remove(chain)
        chain._model = None


    def get_chain_by_id(self, chain_id):
        """Returns the first chain that matches a given chain ID.

        :param str chain_id: The chain ID to search by.
        :rtype: :py:class:`.Chain` or ``None``"""

        if not isinstance(chain_id, str):
            raise TypeError(
             "Chain ID search must be by str, not '%s'" % str(chain_id)
            )
        for chain in self.chains():
            if chain.chain_id() == chain_id:
                return chain


    def duplicate_chain(self, chain, chain_id=None):
        """Creates a copy of a chain in the Model. The coordinates will
        be identical but it will have a unique ID.

        :param Chain chain: The chain to duplicate.
        :param str chain_id: If given, this will determine the ID of the new\
        chain."""

        if not isinstance(chain, Chain):
            raise TypeError(
             "Can only duplicate Chain with this method, not '%s'" % str(
              chain
             )
            )
        if chain not in self.chains():
            raise ValueError(
             "%s is not in this Model and so cannot be duplicated." % str(
              chain
             )
            )
        current_chain_ids = [chain.chain_id() for chain in self.chains()]
        new_chain_id = ""
        if chain_id:
            if chain_id in current_chain_ids:
                raise ValueError(
                 "There is already a Chain with ID %s" % chain_id
                )
            new_chain_id = chain_id
        else:
            new_chain_id = chain.chain_id()
            while new_chain_id in current_chain_ids:
                new_chain_id = chr(ord(new_chain_id) + 1)
        new_residues = []
        next_atom_id = sorted(
         [atom.atom_id() for atom in self.atoms(atom_type="all")]
        )[-1] + 1
        for residue in chain.residues():
            new_atoms = set()
            for atom in residue.atoms(atom_type="all"):
                if isinstance(atom, Atom):
                    new_atoms.add(Atom(
                     atom.x(), atom.y(), atom.z(), atom.element(), next_atom_id, atom.atom_name()
                    ))
                else:
                    new_atoms.add(GhostAtom(atom.element(), next_atom_id, atom.atom_name()))
                next_atom_id += 1
            new_residues.append(Residue(
             new_chain_id + residue.residue_id()[1:],
             residue.residue_name(),
             *new_atoms
            ))

        new_chain = Chain(new_chain_id, *new_residues)
        self.add_chain(new_chain)
        return new_chain


    def bind_sites(self):
        """Returns all the :py:class:`.BindSite` objects in this model.

        :rtype: ``set``"""

        return set(self._bind_sites)


    def add_bind_site(self, site):
        """Adds a bind site to the model.

        :param BindSite site: The bind site to add."""

        if not isinstance(site, BindSite):
            raise TypeError(
             "Can only add BindSite to Model, not '%s'" % str(site)
            )
        if site.site_id() in [mol.site_id() for mol in self.bind_sites()]:
             raise DuplicateBindSitesError(
              "Cannot add site with ID %s to %s as there is already a site with that ID" % (
               site.site_id(), str(self)
              )
             )
        self._bind_sites.add(site)
        site._model = self


    def remove_bind_site(self, site):
        """Removes a bind site from the structure.

        :param BindSite site: The bind site to remove."""

        self._bind_sites.remove(site)
        site._model = None


    def get_bind_site_by_id(self, site_id):
        """Returns the first bind site that matches a given site ID.

        :param str site_id: The site ID to search by.
        :rtype: :py:class:`.BindSite` or ``None``"""

        if not isinstance(site_id, str):
            raise TypeError(
             "BindSite ID search must be by str, not '%s'" % str(site_id)
            )
        for site in self.bind_sites():
            if site.site_id() == site_id:
                return site


    def complexes(self):
        """Returns all the :py:class:`.Complex` objects in this model.

        :rtype: ``set``"""

        return set(self._complexes)


    def add_complex(self, complex_):
        """Adds a complex to the model.

        :param Complex complex: The complex to add."""

        if not isinstance(complex_, Complex):
            raise TypeError(
             "Can only add Complex to Model, not '%s'" % str(complex_)
            )
        if complex_.complex_id() in [mol.complex_id() for mol in self.complexes()]:
             raise DuplicateComplexesError(
              "Cannot add complex with ID %s to %s as there is already a complex with that ID" % (
               complex_.complex_id(), str(self)
              )
             )
        self._complexes.add(complex_)
        complex_._model = self


    def remove_complex(self, complex_):
        """Removes a complex from the model.

        :param Complex complex: The complex to remove."""

        self._complexes.remove(complex_)
        complex_._model = None


    def get_complex_by_id(self, complex_id):
        """Returns the first complex that matches a given complex ID.

        :param str complex_id: The complex ID to search by.
        :rtype: :py:class:`.Complex` or ``None``"""

        if not isinstance(complex_id, str):
            raise TypeError(
             "Complex ID search must be by str, not '%s'" % str(complex_id)
            )
        for complex_ in self.complexes():
            if complex_.complex_id() == complex_id:
                return complex_


    def get_complex_by_name(self, complex_name):
        """Returns the first complex that matches a given name.

        :param str complex_name: The name to search by.
        :rtype: :py:class:`.Complex` or ``None``"""

        if not isinstance(complex_name, str):
            raise TypeError(
             "Complex name search must be by str, not '%s'" % str(complex_name)
            )
        for complex_ in self.complexes():
            if complex_.complex_name() == complex_name:
                return complex_


    def get_complexes_by_name(self, complex_name):
        """Returns all the complexes of a given name.

        :param str complex_name: The name to search by.
        :rtype: ``set`` of :py:class:`.Complex` objects."""

        if not isinstance(complex_name, str):
            raise TypeError(
             "Complex name search must be by str, not '%s'" % str(complex_name)
            )
        return set([complex_ for complex_ in self.complexes()
         if complex_.complex_name() == complex_name])


    def duplicate_complex(self, complex_, complex_id=None, complex_name=None):
        """Creates a copy of a complex in the Model. The coordinates will
        be identical but it will have a unique ID.

        :param Complex complex: The complex to duplicate.
        :param str complex_id: If given, this will determine the ID of the new\
        complex.
        :param str complex_name: If given, this will determine the name of the\
        new complex."""

        if not isinstance(complex_, Complex):
            raise TypeError(
             "Can only duplicate Complex with this method, not '%s'" % str(
              complex_
             )
            )
        if complex_ not in self.complexes():
            raise ValueError(
             "%s is not in this Model and so cannot be duplicated." % str(
              complex_
             )
            )
        current_complex_ids = [complex_.complex_id() for complex_ in self.complexes()]
        new_complex_id = ""
        if complex_id:
            if complex_id in current_complex_ids:
                raise ValueError(
                 "There is already a Complex with ID %s" % complex_id
                )
            new_complex_id = complex_id
        elif re.match(r"(.+)_\d+$", complex_.complex_id()):
            current_iteration = int(complex_.complex_id().split("_")[-1])
            while "_".join(complex_.complex_id().split("_")[:-1]) + "_" +\
             str(current_iteration) in current_complex_ids:
                current_iteration += 1
            new_complex_id = "_".join(complex_.complex_id().split("_")[:-1]) +\
             "_" + str(current_iteration)
        elif complex_.complex_id() + "_1" in current_complex_ids:
            current_iteration = 1
            while "%s_%i" % (complex_.complex_id(), current_iteration) in current_complex_ids:
                current_iteration += 1
            new_complex_id = "%s_%i" % (complex_.complex_id(), current_iteration)
        else:
            new_complex_id = complex_.complex_id() + "_1"
        new_chains = []
        for chain in complex_.chains():
            new_chains.append(self.duplicate_chain(chain))
        new_complex = Complex(
         new_complex_id,
         complex_.complex_name(),
         *new_chains
        )
        self.add_complex(new_complex)
        return new_complex


    def to_pdb_data_file(self):
        """Converts the Model to a :py:class:`.PdbDataFile`."""

        from ..converters.model2pdbdatafile import pdb_data_file_from_model
        return pdb_data_file_from_model(self)


    def save_as_pdb(self, path):
        """Saves the Model to file as a PDB file.

        :param str path: The location and file name to save as."""

        data_file = self.to_pdb_data_file()
        pdb_file = data_file.to_pdb_file()
        with open(path, "w") as f:
            f.write(pdb_file.convert_to_string())
