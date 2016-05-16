"""molecuPy custom exceptions."""

class InvalidElementError(Exception):
    """The exception raised if an atom is created with an invalid element. The
    element doesn't need to be on the periodic table, but it does need to be a
    string either one or two characters long."""
    pass


class NoAtomsError(Exception):
    """The exception raised if an atomic structure is created without passing
    any atoms."""
    pass


class NoResiduesError(Exception):
    """The exception raised if a residuic structure is created without passing
    any residues."""
    pass


class DuplicateAtomIdError(Exception):
    """The exception raised if an atomic structure is created with two atoms of
    the same atom_id."""
    pass


class DuplicateResidueIdError(Exception):
    """The exception raised if a residuic structure is created with two residues
    of the same residue_id."""
    pass


class DuplicateSmallMoleculeError(Exception):
    """The exception raised if a PdbModel is given a small molecule when there
    is already a small molecule with that molecule_id."""
    pass


class DuplicateChainError(Exception):
    """The exception raised if a PdbModel is given a chain when there is already
     a small molecule with that chain_id."""
    pass


class InvalidPdbCodeError(Exception):
    """The exception raised when a PDB file is requested that does not seem to
    exist."""
    pass


class LongBondWarning(Warning):
    """The warning issued if a covalent bond is made between two atoms that
    is unrealistically long."""
    pass
