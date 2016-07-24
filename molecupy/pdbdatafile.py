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


def merge_records(records, start, join=" ", dont_condense=""):
    string = join.join(
     str(record[start:]
    ) if record[start:] else "" for record in records)
    condense = [char for char in " ;:,-" if char not in dont_condense]
    for char in condense:
        string = string.replace(char + " ", char)
    return string
