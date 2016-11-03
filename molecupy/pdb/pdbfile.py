"""This module is used to provide a container to the PDB file itself and its
records - but not the data contained within them."""

class PdbRecord:
    """Represents the lines, or 'records' in a PDB file.

    Indexing a ``PdbRecord`` will get the equivalent slice of the record text,
    only stripped, and converted to ``int`` or ``float`` if possible. Empty
    sub-strings will return ``None``.

    :param str text: The raw text of the record.
    :param PdbFile pdb_file: Optional: a :py:class:`.PdbFile` that the record should be associated with."""

    def __init__(self, text, pdb_file=None):
        if not isinstance(text, str):
            raise TypeError("PdbRecord text must be str, not '%s'" % str(text))
        if len(text) == 0:
            raise ValueError("PdbRecord cannot be created wiyh empty string")
        expanded_text = text[:80] if len(text) >= 80 else text + (" " * (80 - len(text)))
        self._text = expanded_text
        self._name = expanded_text[:6].strip()
        self._content = expanded_text[6:]
        if pdb_file is not None and not isinstance(pdb_file, PdbFile):
            raise TypeError("pdb_file text must be PdbFile, not '%s'" % str(pdb_file))
        self._pdb_file = pdb_file


    def __repr__(self):
        if self.number():
            return "<PdbRecord %i (%s)>" % (self.number(), self._name)
        else:
            return "<PdbRecord (%s)>" % self._name


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
        """The record's line number in its associated :py:class:`.PdbFile`. If
        there is no file associated, this will return ``None``.

        :rtype: ``int``"""

        if self._pdb_file:
            return self._pdb_file.records().index(self) + 1
        else:
            return None


    def name(self, name=None):
        """The record's name (the first six characters). If a string value is
        supplied, the name will be set to the new value, and the text will also
        be updated.

        :param str name: (optional) A new name to change to.
        :rtype: ``str``"""

        if name:
            if not isinstance(name, str):
                raise TypeError("Record name must be str, not '%s'" % str(name))
            if len(name) > 6:
                raise ValueError(
                 "Record name must be <= 6 chars, which '%s' is not" % name
                )
            self._name = name
            self._text = "%-6s%s" % (self._name, self._content)
        else:
            return self._name


    def content(self, content=None):
        """The record's text exlcuding the first six characters. If a string
        value is supplied, the content will be set to the new value, and the
        text will also be updated.

        :param str content: (optional) A new content to change to.
        :rtype: ``str``"""

        if content:
            if not isinstance(content, str):
                raise TypeError(
                 "Record content must be str, not '%s'" % str(content)
                )
            if len(content) > 74:
                raise ValueError(
                 "Record content must be <= 74 chars, which '%s' is not" % content
                )
            self._content = content + (" " * (74 - len(content)))
            self._text = "%-6s%s" % (self._name, self._content)
        else:
            return self._content


    def text(self, text=None):
        """The record's text, extended to 80 characters. If a string
        value is supplied, the text will be set to the new value, and the
        name and content will also be updated.

        :param str text: (optional) A new text to change to.
        :rtype: ``str``"""

        if text:
            if not isinstance(text, str):
                raise TypeError(
                 "Record text must be str, not '%s'" % str(text)
                )
            if len(text) > 80:
                raise ValueError(
                 "Record text must be <= 80 chars, which '%s' is not" % text
                )
            self._text = text + (" " * (80 - len(text)))
            self._name = self._text[:6].strip()
            self._content = self._text[6:]
        else:
            return self._text


    def pdb_file(self, pdb_file=None):
        """The :py:class:`.PdbFile` that the record is associated with. This
        method can update the associated file by passing a :py:class:`.PdbFile`
        to it.

        :param PdbFile pdb_file: (optional) A new :py:class:`.PdbFile` to set.
        :rtype: ``PdbFile``"""

        if pdb_file:
            if not isinstance(pdb_file, PdbFile):
                raise TypeError(
                 "pdb_file text must be PdbFile, not '%s'" % str(pdb_file)
                )
            self._pdb_file = pdb_file
        else:
            return self._pdb_file



class PdbFile:
    """A PDB File - a representation of the file itself, with no processing of
    the data it contains (other than reading record names from the start of each
    line).

    :param str file_string: The raw text of a PDB file."""


    def __init__(self, file_string=""):
        self._file_string = "".join([
         char for char in file_string if 32 <= ord(char) <= 126 or char=="\n"
        ])
        self._records = [
         PdbRecord(line, self) for line in self._file_string.split("\n") if line
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

        return self._records[:]


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


    def add_record(self, record):
        if not isinstance(record, PdbRecord):
            raise TypeError(
             "Can only add PdbRecord objects to PdbFiles, not '%s'" % str(record)
            )
        self._records.append(record)
        record.pdb_file(self)


    def remove_record(self, record):
        if not isinstance(record, PdbRecord):
            raise TypeError(
             "Can only remove PdbRecord objects from PdbFiles, not '%s'" % str(record)
            )
        self._records.remove(record)
        record._pdb_file = None


    def convert_to_string(self):
        lines = [record.text() for record in self.records()]
        return "\n".join(lines)
