"""Contains the PdbDataFile class."""

class PdbDataFile:
    """A PdbDataFile is used to represent the parsed data from a
    :py:class:`.PdbFile`"""

    __slots__ = ["atoms", "heteroatoms", "connections"]

    def __repr__(self):
        return "<PdbDataFile>"
