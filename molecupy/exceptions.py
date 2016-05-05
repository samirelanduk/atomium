class InvalidElementError(Exception):
    pass


class NoAtomsError(Exception):
    pass


class DuplicateAtomIdError(Exception):
    pass


class BrokenMoleculeError(Exception):
    pass

    
class LongBondWarning(Warning):
    pass
