"""This module contains the Model class and its interfaces."""

from .atoms import Atom
from .molecules import AtomicStructure, Molecule, Residue
from .chains import Chain

class Model(AtomicStructure):
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
        for cluster in self._atoms.values():
            for atom in cluster:
                atom._model = self
