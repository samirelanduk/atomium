class PdbFile:

    def __init__(self, file_string):
        self.records = []
        self.file_contents = file_string


    def __repr__(self):
        return "<PdbFile (%i records)>" % len(self.records)
