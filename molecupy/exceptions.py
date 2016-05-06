class InvalidElementError(Exception):
    pass


class NoAtomsError(Exception):
    pass


class DuplicateAtomIdError(Exception):
    pass


class BrokenMoleculeError(Exception):
    pass


class NoMonomersError(Exception):
    pass


class DuplicateMonomerIdError(Exception):
    pass


class InvalidAtomError(Exception):
    pass


class BrokenPolymerError(Exception):
    pass


class LongBondWarning(Warning):
    pass
