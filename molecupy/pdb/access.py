"""This module contains the functions used to access PDB files themselves. These
are the only functions to be imported into the top level directory, and so are
all accesisble by importing molecupy itself."""

import requests
from .pdbfile import PdbFile
from .pdbdatafile import PdbDataFile
from .pdb import Pdb
from ..exceptions import InvalidPdbCodeError

def pdb_from_string(text):
    """Creates a :py:class:`.Pdb` object from the text of a PDB file.

    :param str string: The raw text of a PDB file.
    :rtype: :py:class:`.Pdb`"""

    return Pdb(PdbFile(text).to_pdb_data_file())


def pdb_data_file_from_string(text):
    return PdbFile(text).to_pdb_data_file()


def pdb_file_from_string(text):
    return PdbFile(text)


def get_pdb_from_file(path, processing="pdb"):
    if processing not in ("pdb", "datafile", "pdbfile"):
        raise ValueError(
         "Only valid processing arguments are 'pdb', 'datafile' and 'pdbfile'"
        )
    with open(path) as f:
        if processing == "pdbfile":
            return pdb_file_from_string(f.read())
        elif processing == "datafile":
            return pdb_data_file_from_string(f.read())
        else:
            return pdb_from_string(f.read())


def get_pdb_remotely(code, processing="pdb"):
    if processing not in ("pdb", "datafile", "pdbfile"):
        raise ValueError(
         "Only valid processing arguments are 'pdb', 'datafile' and 'pdbfile'"
        )
    response = requests.get(
     "http://www.ebi.ac.uk/pdbe/entry-files/pdb%s.ent" % code.lower()
    )
    if response.status_code == 200 and response.text[:6] == "HEADER":
        contents = response.text
        if processing == "pdbfile":
            return pdb_file_from_string(contents)
        elif processing == "datafile":
            return pdb_data_file_from_string(contents)
        else:
            return pdb_from_string(contents)
    else:
        raise InvalidPdbCodeError(
         "%s does not seem to be a valid PDB code." % code
        )
