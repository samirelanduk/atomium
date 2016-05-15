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
        if not chunk: return None
        if chunk.count(".") == 1:
            try:
                return float(chunk)
            except ValueError:
                pass
        try:
            return int(chunk)
        except ValueError:
            return chunk


    def get_as_string(self, start, end):
        splice = self[start:end]
        return str(splice) if splice is not None else None



class PdbFile:

    def __init__(self, file_string):
        self.file_contents = "".join([
         char for char in file_string if 32 <= ord(char) <= 126 or char=="\n"
        ])
        self.records = [
         PdbRecord(line, i) for i, line in
          enumerate(self.file_contents.split("\n"), start=1) if line
        ]


    def __repr__(self):
        return "<PdbFile (%i records)>" % len(self.records)


    def get_record_by_name(self, record_name):
        for record in self.records:
            if record.name == record_name:
                return record


    def get_records_by_name(self, record_name):
        return [record for record in self.records if record.name == record_name]
