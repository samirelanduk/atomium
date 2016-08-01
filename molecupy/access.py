"""This module contains the functions used to access PDB files themselves. These
are the only functions to be imported into the top level directory, and so are
all accesisble by importing molecupy itself."""

import requests
from .pdbfile import PdbFile
from .pdbdatafile import PdbDataFile
from .pdb import Pdb
from .exceptions import InvalidPdbCodeError

def pdb_from_string(text):
    """Creates a :py:class:`.Pdb` object from the text of a PDB file.

    :param str string: The raw text of a PDB file.
    :rtype: :py:class:`.Pdb`"""

    return Pdb(PdbDataFile(PdbFile(text)))


def get_pdb_from_file(path):
    """Gets a :py:class:`.Pdb` object from a PDB file stored on disk.

    :param str path: The path to the PDB file.
    :rtype: :py:class:`.Pdb`"""

    with open(path) as f:
        return pdb_from_string(f.read())


def get_pdb_remotely(code):
    """Gets a :py:class:`.Pdb` object from a PDB code.
    The file is requested from the RCSB servers via a HTTP request.
    
    :param str code: The PDB code required - e.g. '1NVQ'.
    :rtype: :py:class:`.Pdb`
    :raises: :class:`.InvalidPdbCodeError` if the PDB doesn't exist"""

    response = requests.get(
     "http://www.ebi.ac.uk/pdbe/entry-files/pdb%s.ent" % code.lower()
    )
    if response.status_code == 200 and response.text[:6] == "HEADER":
        contents = response.text
        return pdb_from_string(contents)
    else:
        raise InvalidPdbCodeError(
         "%s does not seem to be a valid PDB code." % code
        )
