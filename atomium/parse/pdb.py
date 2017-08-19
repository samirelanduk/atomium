"""Contains the Pdb class and functions for opening them."""

from ..structures.models import Model

class Pdb:
    """A PdbDataFile is used to represent a fully processed PDB file."""

    def __init__(self):
        self._model = None


    def __repr__(self):
        return "<Pdb>"


    def model(self, model=None):
        """Returns the :py:class:`.Model` that the .pdb file contains. If a
        model is given, the model will be changed to the new model.

        :param Model model: If given, the model will be updated to this.
        :raises TypeError: if the model given is not a :py:class:`.Model`."""

        if model is None:
            return self._model
        else:
            if not isinstance(model, Model):
                raise TypeError("model must be Model, not '{}'".format(model))
            self._model = model



def pdb_from_file(path):
    """Opens a .pdb file at the specified path and creates a :py:class:`.Pdb`
    from it.

    :param str path: The path to open."""

    from ..converters.strings import string_from_file
    from ..converters.string2pdbfile import string_to_pdb_file
    from ..converters.pdbfile2pdbdatafile import pdb_file_to_pdb_data_file
    from ..converters.pdbdatafile2pdb import pdb_data_file_to_pdb
    s = string_from_file(path)
    pdb_file = string_to_pdb_file(s)
    data_file = pdb_file_to_pdb_data_file(pdb_file)
    pdb = pdb_data_file_to_pdb(data_file)
    return pdb
