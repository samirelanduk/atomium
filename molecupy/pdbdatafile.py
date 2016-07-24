import datetime

class PdbDataFile:

    def __init__(self, pdb_file):
        self._pdb_file = pdb_file


    def __repr__(self):
        return "<PdbDataFile (????)>"


    def pdb_file(self):
        return self._pdb_file



def date_from_string(s):
    return datetime.datetime.strptime(
     s, "%d-%b-%y"
    ).date()
