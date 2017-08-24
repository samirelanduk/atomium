"""Contains the functions for creating PdbDataFiles from Pdbs."""

from .structure2pdbdatafile import structure_to_pdb_data_file

def pdb_to_pdb_data_file(pdb):
    """Converts a :py:class:`.Pdb` to a :py:class:`.PdbDataFile`

    :param Pdb pdb: The Pdb.
    :rtype: ``PdbDataFile``"""

    data_file = structure_to_pdb_data_file(pdb.model())
    return data_file
