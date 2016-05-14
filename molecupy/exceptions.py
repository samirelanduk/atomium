class InvalidElementError(Exception):
    pass


class NoAtomsError(Exception):
    pass


class DuplicateAtomIdError(Exception):
    pass


class BrokenMoleculeError(Exception):
    pass


class NoResiduesError(Exception):
    pass


class DuplicateResidueIdError(Exception):
    pass


class InvalidAtomError(Exception):
    pass


class BrokenChainError(Exception):
    pass


class NoChainsError(Exception):
    pass


class DuplicateMoleculeError(Exception):
    pass


class DuplicateChainError(Exception):
    pass


class DuplicateComplexError(Exception):
    pass


class DuplicateSiteError(Exception):
    pass


class InvalidPdbCodeError(Exception):
    pass


class LongBondWarning(Warning):
    pass
