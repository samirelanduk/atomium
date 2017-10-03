"""This module contains the Model class and its interfaces."""

from .molecules import AtomicStructure, Molecule, Residue
from .chains import ResidueStructure, Chain

class ChainStructure:
    """This is an interface class which confers upon an object the ability to
    retrieve the :py:class:`.Chain` objects it contains. It should not be
    instantiated directly.

    Only classes which inherit from :py:class:`.AtomicStructure` should inherit
    this class, because it requires the :py:meth:`~.AtomicStructure.atoms`
    method."""

    def chains(self, chain_id=None, name=None):
        """Returns the :py:class:`.Chain` objects in the structure. It can be
        given search criteria if you wish.

        :param str chain_id: Filter by chain ID.
        :param str name: Filter by name.
        :rtype: ``set``"""

        chains = set()
        for atom in self.atoms():
            chains.add(atom.chain())
        try:
            chains.remove(None)
        except KeyError: pass
        if chain_id:
            chains = set(filter(lambda c: c.chain_id() == chain_id, chains))
        if name:
            chains = set(filter(lambda c: c.name() == name, chains))
        return chains


    def chain(self, *args, **kwargs):
        """Returns the first :py:class:`.Chain` object in the structure which
        matches the given criteria.

        :param str chain_id: Filter by chain ID.
        :param str name: Filter by name.
        :rtype: ``Chain``"""

        chains = self.chains(*args, **kwargs)
        for chain in chains: return chain


    def add_chain(self, chain):
        """Adds a :py:class:`.Chain` to the structure.

        :param Chain chain: The Chain to add.
        :raises TypeError: if a non-Chain is given."""

        if not isinstance(chain, Chain):
            raise TypeError("{} is not a Chain".format(chain))
        for atom in chain.atoms():
            self.add_atom(atom)


    def remove_chain(self, chain):
        """Removes a :py:class:`.Chain` from the structure.

        :param Chain chain: The Chain to remove.
        :raises TypeError: if a non-Chain is given."""

        if not isinstance(chain, Chain):
            raise TypeError("{} is not a Chain".format(chain))
        for atom in chain.atoms():
            self.remove_atom(atom)



class MoleculeStructure:
    """This is an interface class which confers upon an object the ability to
    retrieve the :py:class:`.Molecule` objects it contains. It should not be
    instantiated directly.

    Only classes which inherit from :py:class:`.AtomicStructure` should inherit
    this class, because it requires the :py:meth:`~.AtomicStructure.atoms`
    method."""

    def molecules(self, molecule_id=None, name=None, generic=False, water=True):
        """Returns the :py:class:`.Molecule` objects in the structure. It can be
        given search criteria if you wish.

        :param str molecule_id: Filter by molecule ID.
        :param str name: Filter by name.
        :param bool generic: If ``True``, chains and residues will exlcuded, \
        and only isolated small molecules will be returned.
        :param bool water: If ``False``, water molecules will be excluded.
        :rtype: ``set``"""

        molecules = set()
        for atom in self.atoms():
            molecules.add(atom.molecule())
        try:
            molecules.remove(None)
        except KeyError: pass
        if generic:
            molecules = set(filter(
             lambda m: not isinstance(m, (Residue, Chain)), molecules
            ))
        if not water:
            molecules = set(filter(
             lambda m: m.name() not in ("HOH", "WAT"), molecules
            ))
        if molecule_id:
            molecules = set(filter(
             lambda m: m.molecule_id() == molecule_id, molecules
            ))
        if name:
            molecules = set(filter(lambda m: m.name() == name, molecules))
        return molecules


    def molecule(self, *args, **kwargs):
        """Returns the first :py:class:`.Molecule` object in the structure which
        matches the given criteria.

        :param str molecule_id: Filter by molecule ID.
        :param str name: Filter by name.
        :param bool generic: If ``True``, chains and residues will exlcuded, \
        and only isolated small molecules will be returned.
        :param bool water: If ``False``, water molecules will be excluded.
        :rtype: ``Molecule``"""

        molecules = self.molecules(*args, **kwargs)
        for molecule in molecules: return molecule


    def add_molecule(self, molecule):
        """Adds a :py:class:`.Molecule` to the structure.

        :param Molecule molecule: The Molecule to add.
        :raises TypeError: if a non-Molecule is given."""

        if not isinstance(molecule, Molecule):
            raise TypeError("{} is not a Molecule".format(molecule))
        for atom in molecule.atoms():
            self.add_atom(atom)


    def remove_molecule(self, molecule):
        """Removes a :py:class:`.Molecule` from the structure.

        :param Molecule molecule: The Molecule to remove.
        :raises TypeError: if a non-Molecule is given."""

        if not isinstance(molecule, Molecule):
            raise TypeError("{} is not a Molecule".format(molecule))
        for atom in molecule.atoms():
            self.remove_atom(atom)



class Model(AtomicStructure, ResidueStructure, ChainStructure,
            MoleculeStructure):
    """Base classes: :py:class:`.AtomicStructure` and
    :py:class:`.ResidueStructure` and :py:class:`.ChainStructure` and
    :py:class:`.MoleculeStructure`

    Represents molecular systems. These are essentially the isolated universes
    in which the other structures live.

    :param \*atoms: The atoms that make up the model. These can also be\
    :py:class:`.AtomicStructure` objects, in which case the atoms of that\
    structure will be used in its place."""

    def __init__(self, *atoms):
        AtomicStructure.__init__(self, *atoms)
        for atom in self._atoms:
            atom._model = self
        for atom in self._id_atoms:
            self._id_atoms[atom]._model = self


    def add_atom(self, atom, *args, **kwargs):
        """Adds an :py:class:`.Atom` to the model.

        :param Atom atom: The atom to add."""

        AtomicStructure.add_atom(self, atom, *args, **kwargs)
        atom._model = self


    def remove_atom(self, atom, *args, **kwargs):
        """Removes an :py:class:`.Atom` from the model.

        :param Atom atom: The atom to remove."""

        AtomicStructure.remove_atom(self, atom, *args, **kwargs)
        atom._model = None
