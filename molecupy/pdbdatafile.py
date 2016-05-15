import datetime

class PdbDataFile:

    def __init__(self, pdb_file):
        self.pdb_file = pdb_file

        process_header(self)
        process_obslte(self)
        process_title(self)


    def __repr__(self):
        return "<PdbDataFile (????)>"



def date_from_string(s):
    return datetime.datetime.strptime(
     s, "%d-%b-%y"
    ).date()


def merge_records(records, start, join=" ", dont_condense=""):
    string = join.join(str(record[start:]) if record[start:] else "" for record in records)
    condense = [char for char in " ;:,-" if char not in dont_condense]
    for char in condense:
        string = string.replace(char + " ", char)
    return string


def process_header(data_file):
    header = data_file.pdb_file.get_record_by_name("HEADER")
    data_file.classification = header[10:50] if header else None
    data_file.deposition_date = date_from_string(header[50:59]) if header else None
    data_file.pdb_code = header[62:66] if header else None


def process_obslte(data_file):
    obslte = data_file.pdb_file.get_record_by_name("OBSLTE")
    data_file.is_obsolete = bool(obslte)
    data_file.obsolete_date = date_from_string(obslte[11:20]) if obslte else None
    data_file.replacement_code = obslte[31:35] if obslte else None


def process_title(data_file):
    titles = data_file.pdb_file.get_records_by_name("TITLE")
    title = merge_records(titles, 10, dont_condense=",;:-")
    data_file.title = title if title else None
