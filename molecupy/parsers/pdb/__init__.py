import requests
from .pdb_file import PdbFile
from .pdb_data_file import PdbDataFile
from .pdb import Pdb
from ...exceptions import InvalidPdbCodeError

def pdb_from_string(string):
    pdb_file = PdbFile(string)
    pdb_data_file = PdbDataFile(pdb_file)
    pdb = Pdb(pdb_data_file)
    return pdb


def get_pdb_from_file(path):
    with open(path) as f:
        return pdb_from_string(f.read())


def get_pdb_remotely(code):
    response = requests.get(
     "http://www.rcsb.org/pdb/files/%s.pdb" % code
    )
    if response.status_code == 200 and response.text[:6] == "HEADER":
        contents = response.text
        return pdb_from_string(contents)
    else:
        raise InvalidPdbCodeError(
         "%s does not seem to be a valid PDB code." % code
        )
