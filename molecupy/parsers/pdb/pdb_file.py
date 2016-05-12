class PdbFile:

    def __init__(self, file_string):
        self.records = []
        self.file_contents = "".join([
         char for char in file_string if 32 <= ord(char) <= 126 or char=="\n"
        ])


    def __repr__(self):
        return "<PdbFile (%i records)>" % len(self.records)
