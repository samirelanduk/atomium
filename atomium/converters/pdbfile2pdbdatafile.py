"""Contains the function for creating PdbDataFiles from PdbFiles."""

from atomium.parse.pdbdatafile import PdbDataFile

def pdb_file_to_pdb_data_file(pdb_file):
    """Converts a :py:class:`.PdbFile` to a :py:class:`.PdbDataFile`

    :param PdbFile pdb_file: The PdbFile.
    :rtype: ``PdbDataFile``"""

    data_file = PdbDataFile()
    return data_file
