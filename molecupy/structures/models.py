from .molecules import AtomicStructure, SmallMolecule
from .chains import Chain, BindSite
from ..exceptions import DuplicateSmallMoleculesError, DuplicateChainsError
from ..exceptions import DuplicateBindSitesError

class Model(AtomicStructure):
    """Base class: :py:class:`.AtomicStructure`

    Represents the structural environment in which the other structures exist."""

    def __init__(self):
        self._small_molecules = set()
        self._chains = set()
        self._bind_sites = set()


    def __getattr__(self, attribute):
        if attribute == "_atoms":
            atoms = set()
            for molecule in self.small_molecules():
                atoms.update(molecule.atoms())
            for chain in self.chains():
                atoms.update(chain.atoms())
            for site in self.bind_sites():
                atoms.update(site.atoms())
            return atoms
        else:
            return self.__getattribute__(attribute)


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
