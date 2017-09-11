"""Contains the function for creating Pdbs from PdbDataFiles."""

from .pdbdatafile2models import pdb_data_file_to_models
from ..files.pdb import Pdb

def pdb_data_file_to_pdb(data_file):
    """Converts a :py:class:`.PdbDataFile` to a :py:class:`.Pdb`

    :param PdbDataFile data_file: The PdbDataFile.
    :rtype: ``Pdb``"""

    pdb = Pdb()
    pdb._models = pdb_data_file_to_models(data_file)
    return pdb
