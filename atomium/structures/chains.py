"""Contains classes for Protein chains and other polymers."""

from .molecules import AtomicStructure, Residue

class ResidueSequence(AtomicStructure):
    """A ResidueSequence represents some linear sequence of :py:class:`.Residue`
    objects.

    :param \*residues: The residues that make up the sequence.
    :raises TypeError: if non-residues are given.
    :raises ValueError: if no residues are given."""

    def __init__(self, *residues):
        if not residues:
            raise ValueError("ResidueSequence needs at least one Residue")
        atoms = set()
        for residue in residues:
            if not isinstance(residue, Residue):
                raise TypeError("{} is not a Residue".format(residue))
            atoms.update(residue.atoms())
        AtomicStructure.__init__(self, *atoms)
        self._residues = list(residues)


    def __repr__(self):
        return "<ResidueSequence ({} residues)>".format(len(self._residues))


    def __len__(self):
        return len(self._residues)


    def __iter__(self):
        return iter(self._residues)


    def length(self):
        """Returns the number of residues in the sequence.

        :rtype: ``int``"""

        return len(self)


    def residues(self, *args, **kwargs):
        """Returns the :py:class:`.Residue` objects in the structure, in order.
        You can filter these by element if you wish.

        :param int residue_id: If given, only residues whose residue ID matches\
        this will be returned (this will only return one residue).
        :param str name: If given, only residues whose name matches this will\
        be returned.
        :rtype: ``tuple``"""

        if args or kwargs:
            residues = AtomicStructure.residues(self, *args, **kwargs)
            return tuple(sorted(residues, key=lambda r: self._residues.index(r)))
        return tuple(self._residues)


    def add_residue(self, residue):
        AtomicStructure.add_residue(self, residue)
        self._residues.append(residue)


    def remove_residue(self, residue):
        AtomicStructure.remove_residue(self, residue)
        self._residues.remove(residue)
