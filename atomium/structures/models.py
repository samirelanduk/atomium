"""This module contains the Model class and its interfaces."""

from .atoms import Atom
from .molecules import AtomicStructure, Molecule, Residue
from .chains import Chain

class Model(AtomicStructure):
    """Base class: :py:class:`.AtomicStructure`

    Represents molecular systems. These are essentially the isolated universes
    in which the other structures live.

    :param \*atoms: The atoms that make up the model. These can also be\
    :py:class:`.AtomicStructure` objects, in which case the atoms of that\
    structure will be used in its place."""

    def __init__(self, *atoms):
        AtomicStructure.__init__(self, *atoms)
        for atom in self._atoms:
            atom._model = self



class Complex(AtomicStructure):
    """Base class: :py:class:`.AtomicStructure`.

    Clusters of chains which form a single functional unit.

    :param \*atoms: The atoms that make up the model. These can also be\
    :py:class:`.AtomicStructure` objects, in which case the atoms of that\
    structure will be used in its place."""

    def __init__(self, *atoms, id=None, name=None):
        AtomicStructure.__init__(self, *atoms)
        if id is not None and not isinstance(id, str):
            raise TypeError("ID {} is not a string".format(id))
        if name is not None and not isinstance(name, str):
            raise TypeError("Complex name {} is not a string".format(name))
        self._id = id
        self._name = name
        for atom in self._atoms:
            atom._complex = self


    @property
    def id(self):
        """The complex's unique string ID.

        :rtype: ``str``"""

        return self._id


    @property
    def name(self):
        """The complex's name.

        :raises TypeError: if the name given is not str."""

        return self._name


    @name.setter
    def name(self, name):
        if not isinstance(name, str):
            raise TypeError("Complex name '{}' is not str".format(name))
        self._name = name
