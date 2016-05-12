import datetime

class PdbDataFile:

    def __init__(self, pdb_file):
        self.pdb_file = pdb_file

        self.process_header()


    def __repr__(self):
        return "<PdbDataFile (????)>"


    def process_header(self):
        header = self.pdb_file.get_record_by_name("HEADER")
        self.classification = header[10:50] if header else None
        self.deposition_date = date_from_string(header[50:59]) if header else None
        self.pdb_code = header[62:66] if header else None



def date_from_string(s):
    return datetime.datetime.strptime(
     s, "%d-%b-%y"
    ).date()
