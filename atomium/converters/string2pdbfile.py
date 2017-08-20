"""Contains the function for creating PdbFilees from strings."""

from .strings import string2lines
from ..files.pdbfile import PdbRecord, PdbFile

def string_to_pdb_file(s):
    """Converts a string taken from a .pdb file and turns it into a
    :py:class:`.PdbFile` object.

    :param str s: The string to convert.
    :rtype: ``PdbFile``"""

    pdb_file = PdbFile()
    lines = string2lines(s)
    pdb_file._records = [PdbRecord(line) for line in lines]
    return pdb_file
