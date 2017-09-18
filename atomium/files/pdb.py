"""Contains the Pdb class and functions for opening them."""

from ..structures.models import Model

class Pdb:
    """A Pdb is used to represent a fully processed PDB file."""

    def __init__(self):
        self._models = []


    def __repr__(self):
        num = len(self._models)
        return "<Pdb ({} model{})>".format(num, "" if num == 1 else "s")


    def models(self):
        """Returns the :py:class:`.Model` objects that the Pdb contains.

        :rtype: ``tuple``"""

        return tuple(self._models)


    def model(self, model=None):
        """Returns the first :py:class:`.Model` that the Pdb file contains."""

        return self._models[0] if self._models else None


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

    from ..files.pdbdatafile import pdb_data_file_from_file
    from ..converters.pdbdatafile2pdb import pdb_data_file_to_pdb
    data_file = pdb_data_file_from_file(path)
    return pdb_data_file_to_pdb(data_file)


def fetch(code, **kwargs):
    """Gets a :py:class:`.Pdb` from the RCSB web services.

    :param str code: The PDB code to fetch.
    :param bool pdbe: If ``True``, the PDB will instead be fetched from PDBe.
    :rtype: ``PdbFile``"""

    from ..files.pdbdatafile import fetch_data_file
    from ..converters.pdbdatafile2pdb import pdb_data_file_to_pdb
    data_file = fetch_data_file(code, **kwargs)
    return pdb_data_file_to_pdb(data_file)
