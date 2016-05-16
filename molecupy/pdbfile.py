"""This module is used to provide a container to the PDB file itself and its
records - but not the data contained within them."""

class PdbRecord:
    """Represents the lines, or 'records' in a PDB file.

    Indexing a ``PdbRecord`` will get the equivalent slice of the record text,
    only stripped, and converted to ``int`` or ``float`` if possible. Empty
    strings will return ``None``.

    :param str line: The raw text of the record.

    :param int number: The line number in the file.

    .. py:attribute:: number:

        The record's line number.

    .. py:attribute:: text:

        The record's text, extended to 80 characters.

    .. py:attribute:: name:

        The record's name (the first six characters).

    .. py:attribute:: contents:

        The record's text exlcuding the first six characters."""

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
        """Indexing a record will automatically convert the value to an integer
        or float if it can - using this method instead will force it to return a
        string

        :param int start: The start of the subsection.
        :param int end: The end of the subsection.
        :rtype: ``str``"""

        splice = self[start:end]
        return str(splice) if splice is not None else None



class PdbFile:
    """A PDB File - a representation of the file itself, with no processing of
    the data it contains (other than reading record names from the start of each
    line).

    :param str file_string: The raw text of a PDB file.

    .. py:attribute:: file_contents:

        The file string from which the object was created.

    .. py:attribute:: records:

        A list of :py:class:`PdbRecord` objects."""

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
        """Gets the first :py:class:`PdbRecord` of a given name.

        :param str record_name: record name to search by.
        :rtype: :py:class:`PdbRecord` or ``None`` if there is no match."""

        for record in self.records:
            if record.name == record_name:
                return record


    def get_records_by_name(self, record_name):
        """Gets all :py:class:`PdbRecord` objects of a given name.

        :param str record_name: record name to search by.
        :returns: list of :py:class:`PdbRecord` objects."""

        return [record for record in self.records if record.name == record_name]
