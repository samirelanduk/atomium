from .pdbfile import PdbFile
from .pdbdatafile import PdbDataFile
from .pdb import Pdb

def pdb_from_string(text):
    return Pdb(PdbDataFile(PdbFile(text)))


def get_pdb_from_file(path):
    with open(path) as f:
        return pdb_from_string(f.read())
