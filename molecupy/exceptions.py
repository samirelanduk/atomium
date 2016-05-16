class InvalidElementError(Exception):
    pass


class NoAtomsError(Exception):
    pass


class DuplicateAtomIdError(Exception):
    pass


class DuplicateSmallMoleculeError(Exception):
    pass


class InvalidPdbCodeError(Exception):
    pass


class LongBondWarning(Warning):
    pass
