"""This module is used to provide a container to the PDB file itself and its
records - but not the data contained within them."""

class PdbRecord:
    """Represents the lines, or 'records' in a PDB file.

    Indexing a ``PdbRecord`` will get the equivalent slice of the record text,
    only stripped, and converted to ``int`` or ``float`` if possible. Empty
    strings will return ``None``.

    :param str text: The raw text of the record.
    :param int number: The line number in the file."""

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
        """Indexing a record will automatically convert the value to an integer
        or float if it can - using this method instead will force it to return a
        string.

        :param int start: The start of the subsection.
        :param int end: The end of the subsection.
        :rtype: ``str``"""

        splice = self[start:end]
        return str(splice) if splice is not None else None


    def number(self):
        """The record's line number.

        :rtype: ``int``"""

        return self._number


    def name(self):
        """The record's name (the first six characters).

        :rtype: ``str``"""

        return self._name


    def text(self):
        """The record's text, extended to 80 characters.

        :rtype: ``str``"""

        return self._text


    def content(self):
        """The record's text exlcuding the first six characters.

        :rtype: ``str``"""

        return self._content



class PdbFile:
    """A PDB File - a representation of the file itself, with no processing of
    the data it contains (other than reading record names from the start of each
    line).

    :param str file_string: The raw text of a PDB file."""


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
        """The file string from which the object was created.

        :rtype: ``str``"""

        return self._file_string


    def records(self):
        """A list of :py:class:`PdbRecord` objects.

        :returns: list of :py:class:`PdbRecord` objects."""

        return self._records


    def get_record_by_name(self, record_name):
        """Gets the first :py:class:`PdbRecord` of a given name.

        :param str record_name: record name to search by.
        :rtype: :py:class:`PdbRecord` or ``None`` if there is no match."""

        if not isinstance(record_name, str):
            raise TypeError(
             "Can only search for record by str, not '%s'" % str(record_name)
            )
        for record in self.records():
            if record.name() == record_name:
                return record


    def get_records_by_name(self, record_name):
        """Gets all :py:class:`PdbRecord` objects of a given name.

        :param str record_name: record name to search by.
        :returns: ``list`` of :py:class:`PdbRecord` objects."""

        if not isinstance(record_name, str):
            raise TypeError(
             "Can only search for records by str, not '%s'" % str(record_name)
            )
        return [record for record in self.records() if record.name() == record_name]
