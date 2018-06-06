from .structures import AtomStructure
from .data import CODES, FULL_NAMES
from .exceptions import SequenceConnectivityError

def lower(name):
    """A factory function which creates two functions for querying the lower
    structures a structure contains - one for returning all, one for returning
    the first match.

    :param name: The kind of structure to look for."""

    def multi_func(self, *args, **kwargs):
        return self._get(name, *args, **kwargs)
    def single_func(self, *args, **kwargs):
        results = self._get(name, *args, **kwargs)
        for result in results: return result
    return multi_func, single_func


def upper(name):
    """A factory function which creates a function for querying the upper
    structure a structure is a part of - or ``None`` if there is no single
    match.

    :param name: The kind of structure to look for."""

    def func(self):
        objects = self._get(name)
        if len(objects) == 1: return list(objects)[0]
    return func



class Model(AtomStructure):
    """The universe in which all over molecules live."""

    complexes, complex = lower("complex")
    chains, chain = lower("chain")
    residues, residue = lower("residue")
    ligands, ligand = lower("ligand")



class Complex(AtomStructure):
    """An amalgamation of chains which form a functional unit."""

    model = property(upper("model"))
    chains, chain = lower("chain")
    residues, residue = lower("residue")
    ligands, ligand = lower("ligand")



class Chain(AtomStructure):
    """A sequence of residues."""

    def __init__(self, *args, rep="", **kwargs):
        AtomStructure.__init__(self, *args, **kwargs)
        self._rep_sequence = rep
        self.verify()


    def __len__(self):
        return len(self.residues())


    def __getitem__(self, index):
        return self.residues()[index]


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



    model = property(upper("model"))
    complex = property(upper("complex"))
    residue = lower("residue")[1]
    ligands, ligand = lower("ligand")



class Het(AtomStructure):

    model = property(upper("model"))
    complex = property(upper("complex"))
    chain = property(upper("chain"))



class Ligand(Het):
    """A small molecule associated with a chain but not connected to it."""



class Residue(Het):
    """A small subunit within a chain."""

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
