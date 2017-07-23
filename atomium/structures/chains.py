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
