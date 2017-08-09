"""Contains the Model class."""

from .molecules import AtomicStructure

class Model(AtomicStructure):
    """Base class: :py:class:`.AtomicStructure`

    Represents molecular systems.

    :param \*atoms: The atoms that make up the model."""

    def __init__(self, *atoms):
        AtomicStructure.__init__(self, *atoms)
        for atom in atoms:
            atom._model = self


    def add_atom(self, atom, *args, **kwargs):
        AtomicStructure.add_atom(self, atom, *args, **kwargs)
        atom._model = self


    def remove_atom(self, atom, *args, **kwargs):
        AtomicStructure.remove_atom(self, atom, *args, **kwargs)
        atom._model = None
