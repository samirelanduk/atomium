class PdbRecord:

    def __init__(self, text, number):
        if not isinstance(number, int):
            raise TypeError("PdbRecord number must be int, not '%s'" % str(number))
        self._number = number

        if not isinstance(text, str):
            raise TypeError("PdbRecord text must be str, not '%s'" % str(text))
        if len(text) == 0:
            raise ValueError("PdbRecord cannot be created wiyh empty string")
        expanded_text = text[:80] if len(text) >= 80 else text + (" " * (80 - len(text)))
        self._text = expanded_text
        self._name = expanded_text[:6].strip()
        self._content = expanded_text[6:]


    def __repr__(self):
        return "<PdbRecord %i (%s)>" % (self._number, self._name)


    def __getitem__(self, key):
        chunk = self._text[key].strip()
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


    def number(self):
        return self._number


    def name(self):
        return self._name


    def text(self):
        return self._text


    def content(self):
        return self._content



class PdbFile:

    def __init__(self, file_string):
        self._file_string = "".join([
         char for char in file_string if 32 <= ord(char) <= 126 or char=="\n"
        ])
        self._records = [
         PdbRecord(line, i) for i, line in
          enumerate(self._file_string.split("\n"), start=1) if line
        ]


    def __repr__(self):
        return "<PdbFile (%i Records)>" % len(self._records)


    def file_string(self):
        return self._file_string


    def records(self):
        return self._records


    def get_record_by_name(self, record_name):
        if not isinstance(record_name, str):
            raise TypeError(
             "Can only search for record by str, not '%s'" % str(record_name)
            )
        for record in self.records():
            if record.name() == record_name:
                return record


    def get_records_by_name(self, record_name):
        if not isinstance(record_name, str):
            raise TypeError(
             "Can only search for records by str, not '%s'" % str(record_name)
            )
        return [record for record in self.records() if record.name() == record_name]
