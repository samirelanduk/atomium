"""Contains the useful sub-classes of AtomStructure."""

from .structures import AtomStructure
from .data import CODES, FULL_NAMES
from .exceptions import SequenceConnectivityError

def lower(name):
    """A factory function which creates two functions for querying the lower
    structures a structure contains - one for returning all, one for returning
    the first match.

    :param name: The kind of structure to look for."""

    doc = """Returns the {} in this structure.

    :param str id: filter by ID.
    :param str id_regex: filter by ID regex.
    :param str name: filter by name.
    :param str name_regex: filter by name regex.

    :rtype: ``{}``"""

    def multi_func(self, *args, **kwargs):
        return self._get(name, *args, **kwargs)
    def single_func(self, *args, **kwargs):
        results = self._get(name, *args, **kwargs)
        for result in results: return result
    multi_func.__name__ = name + "s"
    single_func.__name__ = name
    multi_func.__doc__ = doc.format(name + "s", "set")
    single_func.__doc__ = doc.format("first matching " + name, name.title())
    return multi_func, single_func


def upper(name):
    """A factory function which creates a function for querying the upper
    structure a structure is a part of - or ``None`` if there is no single
    match.

    :param name: The kind of structure to look for."""

    doc = """Returns the {} the structure is part of.

    :rtype: ``{}``"""

    def func(self):
        objects = self._get(name)
        if len(objects) == 1: return list(objects)[0]
    func.__name__ = name
    func.__doc__ = doc.format(name, name.title())
    return func



class Model(AtomStructure):
    """The universe in which all other molecules live, interact, and generally
    exist.

    :param \*atoms: The atoms the structure is to be made of. The atoms will\
    be updated with awareness of the new structure they are part of if a\
    sub-class is used. You can also pass in other atom structures here, and\
    their atoms will be used.
    :param str id: The structure's ID.
    :param str name: The structure's name."""

    chains, chain = lower("chain")
    residues, residue = lower("residue")
    ligands, ligand = lower("ligand")

    def copy(self):
        """Creates a copy of the Model, as well as copies of its substructures
        such as chains etc.

        :rtype: ``Model``"""

        model = Model(*[chain.copy() for chain in self.chains()])
        atoms = [a for a in self._atoms if a._chain is None]
        atom_copies = [atom.copy() for atom in atoms]
        model._atoms.update(atom_copies)
        model._id, model._name = self._id, self._name
        return model



class Chain(AtomStructure):
    """A sequence of residues. Unlike other structures, they are iterables, and
    have a length.

    Residues can also be accessed using indexing.

    Residues *must* be connected in series when a chain is created from them.
    This is necessary because otherwise it won't know how to return a series of
    residues in order when you ask for them.

    :param \*atoms: The atoms the structure is to be made of. The atoms will\
    be updated with awareness of the new structure they are part of if a\
    sub-class is used. You can also pass in other atom structures here, and\
    their atoms will be used.
    :param str id: The structure's ID.
    :param str name: The structure's name."""

    model = property(upper("model"))
    residue = lower("residue")[1]
    ligands, ligand = lower("ligand")

    def __init__(self, *args, rep="", **kwargs):
        AtomStructure.__init__(self, *args, **kwargs)
        self._rep_sequence = rep
        self.verify()


    def __len__(self):
        return len(self.residues())


    def __getitem__(self, index):
        return self.residues()[index]


    def copy(self):
        """Creates a copy of the Chain, as well as copies of its substructures
        such as ligands etc.

        :rtype: ``Chain``"""

        ligands = [ligand.copy() for ligand in self.ligands()]
        residues = [residue.copy() for residue in self.residues()]
        for res1, res2 in zip(residues[:-1], residues[1:]): res1.next = res2
        substructures = ligands + residues
        chain = Chain(*substructures, rep=self._rep_sequence)
        atoms = [a for a in self._atoms
         if a._residue is None and a._ligand is None]
        atom_copies = [atom.copy() for atom in atoms]
        chain._atoms.update(atom_copies)
        chain._id, chain._name = self._id, self._name
        return chain


    def verify(self):
        """A static method for checking that the residues in a sequence are all
        connected together, and that there are no gaps in the sequence.

        :raises SequenceConnectivityError: if not properly connected.
        :returns: ``True`` if the test passes."""

        residues = set()
        for atom in self._atoms:
            if atom.residue: residues.add(atom.residue)
        if residues:
            seq = [list(residues)[0]]
            while seq[-1].next is not None and seq[-1].next in residues:
                seq.append(seq[-1].next)
            while seq[-1].previous is not None and seq[-1].previous in residues:
                seq.append(seq[-1].previous)
            if set(seq) != residues:
                raise SequenceConnectivityError(
                 "{} missing from sequence {} - check connections".format(
                  residues, residues - set(seq)
                 )
                )
        return True


    @property
    def length(self):
        """The number of :py:class:`.Residue` objects the chain has.

        :rtype: ``int``"""

        return len(self)


    def residues(self, *args, **kwargs):
        """Returns the :py:class:`.Residue` objects in the structure.

        .. WARNING::
            If the residues in the chain are not connected together with
            ``next`` and ``previous``, the method will not return all residues.
            Chains can only be made with connected residues but if you have
            modified the structure since then, this may produce unexpected
            results.

        :rtype: ``tuple``"""

        all_res = lower("residue")[0](self, *args, **kwargs)
        if len(all_res) == 0: return tuple()
        initial_res = list(all_res)[0]
        full_sequence = [initial_res]
        while full_sequence[-1].next:
            full_sequence.append(full_sequence[-1].next)
        while full_sequence[0].previous:
            full_sequence.insert(0, full_sequence[0].previous)
        return tuple(sorted(all_res, key=lambda k: full_sequence.index(k)))


    @property
    def sequence(self):
        """Returns the sequence of amino acids as represented in the residues
        contained.

        :rtype: ``str``"""

        return "".join(residue.code for residue in self.residues())


    @property
    def rep_sequence(self):
        """The sequence the chain represents, *not* the actual sequence of the
        residues present. A model can often have residues missing that the
        experiment it was created from didn't capture. This is the true sequence
        as it should be, not the sequence as it is in the model.

        :raises TypeError: if the sequence given is not str."""

        return self._rep_sequence


    @rep_sequence.setter
    def rep_sequence(self, rep):
        self._rep_sequence = rep



class Het(AtomStructure):
    """A component of a chain."""

    model = property(upper("model"))
    chain = property(upper("chain"))



class Ligand(Het):
    """A small molecule associated with a chain but not connected to it.

    :param \*atoms: The atoms the structure is to be made of. The atoms will\
    be updated with awareness of the new structure they are part of if a\
    sub-class is used. You can also pass in other atom structures here, and\
    their atoms will be used.
    :param str id: The structure's ID.
    :param str name: The structure's name."""



class Residue(Het):
    """A small subunit within a chain.

    :param \*atoms: The atoms the structure is to be made of. The atoms will\
    be updated with awareness of the new structure they are part of if a\
    sub-class is used. You can also pass in other atom structures here, and\
    their atoms will be used.
    :param str id: The structure's ID.
    :param str name: The structure's name."""

    def __init__(self, *atoms, **kwargs):
        Het.__init__(self, *atoms, **kwargs)
        self._next, self._previous = None, None


    @property
    def full_name(self):
        """Returns the full name of the reside if it is one of the 20 canonical
        amino acids - otherwise it just returns the name itself.

        :rtype: ``str``"""

        if self._name:
            return FULL_NAMES.get(self._name.upper(), self._name)


    @property
    def code(self):
        """Returns the amino acid code of the reside if it is one of the 20
        canonical amino acids - otherwise it just returns 'X'.

        :rtype: ``str``"""

        if self._name:
            return CODES.get(self._name.upper(), "X")


    @property
    def next(self):
        """Residues can be linked to each other in a linear chain. This property
        returns the :py:class:`.Residue` downstream of this one. Alternatively,
        if you supply a residue, that residue will be assigned as the 'next' one
        downstream to this, and this residue will be upstream to that.
        Note that is a separate concept from bonds. Creating a connection of
        this kind implies no, and requires no, explicit bonding.

        :param Residue residue: The residue to connect to. If ``None`` is\
        given, any existing connection downstream of this residue will be\
        broken.
        :raises ValueError: if you try to connect a residue to itself.
        :rtype: ``Residue``"""

        return self._next


    @next.setter
    def next(self, residue):
        if residue is None:
            if self._next: self._next._previous = None
            self._next = None
        elif residue is self:
            raise ValueError("Cannot link {} to itself".format(self))
        else:
            self._next = residue
            residue._previous = self


    @property
    def previous(self):
        """Residues can be linked to each other in a linear chain. This property
        returns the :py:class:`.Residue` upstream of this one. Alternatively,
        if you supply a residue, that residue will be assigned as the 'previous'
        one upstream to this, and this residue will be downstream from that.
        Note that is a separate concept from bonds. Creating a connection of
        this kind implies no, and requires no, explicit bonding.

        :param Residue residue: The residue to connect to. If ``None`` is\
        given, any existing connection upstream of this residue will be\
        broken.
        :raises ValueError: if you try to connect a residue to itself.
        :rtype: ``Residue``"""

        return self._previous


    @previous.setter
    def previous(self, residue):
        if residue is None:
            if self._previous: self._previous._next = None
            self._previous = None
        elif residue is self:
            raise ValueError("Cannot link {} to itself".format(self))
        else:
            self._previous = residue
            residue._next = self
