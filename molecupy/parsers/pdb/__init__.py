from .pdb_file import PdbFile
from .pdb_data_file import PdbDataFile
from .pdb import Pdb

def pdb_from_string(string):
    pdb_file = PdbFile(string)
    pdb_data_file = PdbDataFile(pdb_file)
    pdb = Pdb(pdb_data_file)
    return pdb


def get_pdb_from_file(path):
    with open(path) as f:
        return pdb_from_string(f.read())
