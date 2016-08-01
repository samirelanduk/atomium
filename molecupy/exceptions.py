"""molecuPy custom exceptions."""

class LongBondWarning(Warning):
    """The warning issued if a covalent bond is made between two atoms that
    is unrealistically long."""
    pass



class NoAtomsError(Exception):
    """The exception raised if an atomic structure is created without passing
    any atoms."""
    pass



class NoResiduesError(Exception):
    """The exception raised if a residuic structure is created without passing
    any residues."""
    pass



class MultipleResidueConnectionError(Exception):
    """The exception raised when a residue connection is made to a residue which
    is already connected to a residue in that fashion."""
    pass



class BrokenHelixError(Exception):
    """The exception raised when an alpha helix is created with residues on
    different chains."""
    pass



class BrokenStrandError(Exception):
    """The exception raised when a beta strand is created with residues on
    different chains."""
    pass



class DuplicateAtomsError(Exception):
    """The exception raised if an atomic structure is created with two atoms of
    the same atom_id."""
    pass



class DuplicateSmallMoleculesError(Exception):
    """The exception raised if a Model is given a small molecule when there
    is already a small molecule with that molecule_id."""
    pass



class DuplicateResiduesError(Exception):
    """The exception raised if a residuic structure is created with two residues
    of the same residue_id."""
    pass



class DuplicateChainsError(Exception):
    """The exception raised if a Model is given a chain when there is already
     a chain with that chain_id."""
    pass



class DuplicateBindSitesError(Exception):
    """The exception raised if a Model is given a bindsite when there is already
     a site with that site_id."""
    pass



class InvalidPdbCodeError(Exception):
    """The exception raised when a PDB file is requested that does not seem to
    exist."""
    pass
