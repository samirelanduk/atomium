import datetime

class PdbDataFile:

    def __init__(self, pdb_file):
        self.pdb_file = pdb_file

        self.process_header()
        self.process_obslte()
        self.process_title()


    def __repr__(self):
        return "<PdbDataFile (????)>"


    def process_header(self):
        header = self.pdb_file.get_record_by_name("HEADER")
        self.classification = header[10:50] if header else None
        self.deposition_date = date_from_string(header[50:59]) if header else None
        self.pdb_code = header[62:66] if header else None


    def process_obslte(self):
        obslte = self.pdb_file.get_record_by_name("OBSLTE")
        self.is_obsolete = bool(obslte)
        self.obsolete_date = date_from_string(obslte[11:20]) if obslte else None
        self.replacement_code = obslte[31:35] if obslte else None


    def process_title(self):
        titles = self.pdb_file.get_records_by_name("TITLE")
        title = merge_records(titles, 10, dont_condense=",;:-")
        self.title = title if title else None



def date_from_string(s):
    return datetime.datetime.strptime(
     s, "%d-%b-%y"
    ).date()


def merge_records(records, start, join=" ", dont_condense=""):
        string = join.join(record[start:] for record in records)
        condense = [char for char in " ;:,-" if char not in dont_condense]
        for char in condense:
            string = string.replace(char + " ", char)
        return string
