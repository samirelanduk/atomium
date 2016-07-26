from .pdbfile import PdbFile
from .pdbdatafile import PdbDataFile
from .pdb import Pdb

def pdb_from_string(text):
    return Pdb(PdbDataFile(PdbFile(text)))
