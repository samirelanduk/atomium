"""Contains classes for structures made of atoms."""

class AtomicStructure:
    """Represents structures made of :py:class:`.Atom`s, which tends to be
    rather a lot of things in practice. This class would not generally be
    instantiated directly, and is here to be a parent class to other, more
    specific entities.

    :param \*atoms: The :py:class:`.Atom`s that make up the structure."""

    def __init__(self, *atoms):
        self._atoms = set(atoms)
