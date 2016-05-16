class InvalidElementError(Exception):
    pass


class NoAtomsError(Exception):
    pass


class NoResiduesError(Exception):
    pass


class DuplicateAtomIdError(Exception):
    pass


class DuplicateResidueIdError(Exception):
    pass


class DuplicateSmallMoleculeError(Exception):
    pass


class InvalidPdbCodeError(Exception):
    pass


class LongBondWarning(Warning):
    pass
