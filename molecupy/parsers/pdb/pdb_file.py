class PdbFile:

    def __init__(self, file_string):
        self.records = []
        self.file_contents = "".join([
         char for char in file_string if 32 <= ord(char) <= 126 or char=="\n"
        ])


    def __repr__(self):
        return "<PdbFile (%i records)>" % len(self.records)



class PdbRecord:

    def __init__(self, line, number):
        self.number = number
        self.text = line.ljust(80)
        self.name = self.text[:6].strip()
        self.contents = self.text[6:]


    def __repr__(self):
        return "<PdbRecord (%s)>" % self.name


    def __getitem__(self, key):
        chunk = self.text[key].strip()
        try:
            return int(chunk)
        except ValueError:
            return chunk
