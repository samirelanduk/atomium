"""Contains the PdbDataFile class."""



class PdbDataFile:
    """A PdbDataFile is used to represent the parsed data from a
    :py:class:`.PdbFile`"""

    __slots__ = [
     "code", "title", "deposition_date", "atoms", "heteroatoms", "connections"
    ]

    def __repr__(self):
        if hasattr(self, "code") and self.code:
            return "<PdbDataFile ({})>".format(self.code)
        return "<PdbDataFile>"


def pdb_data_file_from_file(path):
    """Opens a .pdb file at the specified path and creates a
    :py:class:`.PdbdataFile` from it.

    :param str path: The path to open.
    :rtype: ``PdbDataFile``"""

    from ..files.pdbfile import pdb_file_from_file
    from ..converters.pdbfile2pdbdatafile import pdb_file_to_pdb_data_file
    pdb_file = pdb_file_from_file(path)
    return pdb_file_to_pdb_data_file(pdb_file)


def fetch_data_file(code, **kwargs):
    """Gets a :py:class:`.PdbdataFile` from the RCSB web services.

    :param str code: The PDB code to fetch.
    :param bool pdbe: If ``True``, the PDB will instead be fetched from PDBe.
    :rtype: ``PdbFile``"""

    from ..files.pdbfile import fetch_file
    from ..converters.pdbfile2pdbdatafile import pdb_file_to_pdb_data_file
    pdb_file = fetch_file(code, **kwargs)
    return pdb_file_to_pdb_data_file(pdb_file)
