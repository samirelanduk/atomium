import requests
from .pdbfile import PdbFile
from .pdbdatafile import PdbDataFile
from .pdb import Pdb
from .exceptions import InvalidPdbCodeError

def pdb_from_string(text):
    return Pdb(PdbDataFile(PdbFile(text)))


def get_pdb_from_file(path):
    with open(path) as f:
        return pdb_from_string(f.read())


def get_pdb_remotely(code):
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
