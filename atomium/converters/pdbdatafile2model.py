"""Contains the function for creating Models from PdbDataFiles."""

from ..structures.models import Model

def pdb_data_file_to_model(data_file):
    """Converts a :py:class:`.PdbDataFile` to a :py:class:`.Model`

    :param PdbDataFile data_file: The PdbDataFile.
    :rtype: ``Model``"""

    model = Model()
    return model
