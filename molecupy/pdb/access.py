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
    """Creates a :py:class:`.PdbDataFile` object from the text of a PDB file.

    :param str string: The raw text of a PDB file.
    :rtype: :py:class:`.PdbDataFile`"""

    return PdbFile(text).to_pdb_data_file()


def pdb_file_from_string(text):
    """Creates a :py:class:`.PdbFile` object from the text of a PDB file.

    :param str string: The raw text of a PDB file.
    :rtype: :py:class:`.PdbFile`"""

    return PdbFile(text)


def get_pdb_from_file(path, processing="pdb"):
    """Creates a :py:class:`.Pdb`, :py:class:`.PdbDataFile`, or
    :py:class:`.PdbFile` from a file path on disk - the default behaviour being
    to create a :py:class:`.Pdb`.

    :param str path: The location of the PDB file on disk.
    :param str processing: The level of processing you want the returned object\
    to have. Propviding ``"pdbfile"`` will just return a :py:class:`.PdbFile`,\
    ``"datafile"`` will return a :py:class:`.PdbDataFile`, and ``"pdb"`` (the\
    default) will return a fully processed :py:class:`.Pdb` object.
    :raises FileNotFoundError: if there is no file at the specified location."""

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
    """Creates a :py:class:`.Pdb`, :py:class:`.PdbDataFile`, or
    :py:class:`.PdbFile` from a 4-letter PDB code - the default behaviour being
    to create a :py:class:`.Pdb`.

    :param str code: The 4-letter PDB code.
    :param str processing: The level of processing you want the returned object\
    to have. Propviding ``"pdbfile"`` will just return a :py:class:`.PdbFile`,\
    ``"datafile"`` will return a :py:class:`.PdbDataFile`, and ``"pdb"`` (the\
    default) will return a fully processed :py:class:`.Pdb` object.
    :raises InvalidPdbCodeError: if there is no PDB with the given code."""

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
