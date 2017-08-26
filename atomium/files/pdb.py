"""Contains the Pdb class and functions for opening them."""

from ..structures.models import Model

class Pdb:
    """A Pdb is used to represent a fully processed PDB file."""

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


    def to_file_string(self):
        """Returns the file text that represents this Pdb.

        :rtype: ``str``"""

        from ..converters.pdb2pdbdatafile import pdb_to_pdb_data_file
        from ..converters.pdbdatafile2pdbfile import pdb_data_file_to_pdb_file
        from ..converters.pdbfile2pdbstring import pdb_file_to_pdb_string
        data_file = pdb_to_pdb_data_file(self)
        pdb_file = pdb_data_file_to_pdb_file(data_file)
        return pdb_file_to_pdb_string(pdb_file)


    def save(self, path):
        """Saves the Pdb as a .pdb file.

        :param str path: The path to save to."""

        from ..converters.strings import string_to_file
        string_to_file(self.to_file_string(), path)



def pdb_from_file(path):
    """Opens a .pdb file at the specified path and creates a :py:class:`.Pdb`
    from it.

    :param str path: The path to open.
    :rtype: ``Pdb``"""

    from ..converters.strings import string_from_file
    from ..converters.string2pdbfile import string_to_pdb_file
    from ..converters.pdbfile2pdbdatafile import pdb_file_to_pdb_data_file
    from ..converters.pdbdatafile2pdb import pdb_data_file_to_pdb
    s = string_from_file(path)
    pdb_file = string_to_pdb_file(s)
    data_file = pdb_file_to_pdb_data_file(pdb_file)
    pdb = pdb_data_file_to_pdb(data_file)
    return pdb


def fetch(code):
    """Gets a :py:class:`.Pdb` from the RCSB web services.

    :param str code: The PDB code to fetch.
    :raises TypeError: if the code is not a string.
    :raises ValueError: if the code is not four caracters long.
    :rtype: ``Pdb``"""

    if not isinstance(code, str):
        raise TypeError("PDB code {} is not string".format(code))
    if len(code) != 4:
        raise ValueError("PDB code {} is not of length 4".format(code))
    from requests import get
    from ..converters.string2pdbfile import string_to_pdb_file
    from ..converters.pdbfile2pdbdatafile import pdb_file_to_pdb_data_file
    from ..converters.pdbdatafile2pdb import pdb_data_file_to_pdb
    response = get("https://files.rcsb.org/view/{}.pdb".format(code))
    if response.status_code == 200:
        s = response.text
        pdb_file = string_to_pdb_file(s)
        data_file = pdb_file_to_pdb_data_file(pdb_file)
        pdb = pdb_data_file_to_pdb(data_file)
        return pdb
