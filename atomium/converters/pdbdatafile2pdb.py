"""Contains the function for creating Pdbs from PdbDataFiles."""

from .pdbdatafile2model import pdb_data_file_to_model
from ..files.pdb import Pdb

def pdb_data_file_to_pdb(data_file):
    """Converts a :py:class:`.PdbDataFile` to a :py:class:`.Pdb`

    :param PdbDataFile data_file: The PdbDataFile.
    :rtype: ``Pdb``"""

    pdb = Pdb()
    pdb._model = pdb_data_file_to_model(data_file)
    return pdb
