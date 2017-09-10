from numbers import Number

class PdbRecord:
    """A PdbRecord represents a single line in a PDB file. As such, it follows
    the same constraints that the PDB file formats impose on the lines in the
    actual files - namely that it cannot be more than 80 characters.

    Two PdbRecords are considered equal if they have the same text.

    PdbRecords are containers and iterables, in the same way that strings are.

    Also like strings you can index them using ``record[x:y]`` notation, but in
    order to aid parsing, they will process the substring before returning it
    It will strip whitespace from the substring, attempt to convert it to an
    int or a float, and return None if an empty string would be returned. It
    will also treat the record as an 80 character string, regardless of the
    actual length of the record. If this is more of an annoyance than a
    convenience, you can use :py:meth:`.get_as_string`, which works just like
    regualar indexing.

    :param str text: The string contents of the line in the PDB file.
    :raises ValueError: if a string longer than 80 characters is supplied."""

    def __init__(self, text):
        if not isinstance(text, str):
            raise TypeError("PdbRecord needs str, not '%s'" % str(text))
        if len(text) > 80:
            raise ValueError("'%s' is longer than 80 characters" % str(text))
        self._text = text.rstrip()
        self._number = None


    def __repr__(self):
        return "<{} record>".format(self.name())


    def __eq__(self, other):
        return isinstance(other, PdbRecord) and self._text == other._text


    def __len__(self):
        return len(self._text)


    def __contains__(self, member):
        return member in self._text


    def __iter__(self):
        return iter(self._text)


    def __getitem__(self, index):
        full_line = self._text + (" " * (80 - len(self._text)))
        chunk = full_line[index].strip()
        try:
            chunk = int(chunk)
        except ValueError:
            try:
                chunk = float(chunk)
            except ValueError:
                pass
        if chunk or isinstance(chunk, Number): return chunk
        return None


    def get_as_string(self, index, index2=None):
        """Indexing a PdbRecord will process the resultant substring in various
        ways, such as by converting it to an int or float. If you just want the
        string regardless, you can use this method to force a straight string
        return.

        :param int index: the index of the text you wish to get.
        :param int index2: if given, a slice will be taken using the two indeces."""

        full_line = self._text + (" " * (80 - len(self._text)))
        if index2:
            return full_line[index: index2]
        return full_line[index]


    def number(self):
        """Returns the line number of the record in its associated
        :py:class:`.PdbFile`.

        :rtype: ``int``"""

        return self._number


    def text(self, text=None):
        """The full string of the record. If a string is passed as an argument,
        the text will be updated, though any trailing whitespace will be
        removed.

        :param str text: If given, the record's text will be updated to this.
        :raises ValueError: if a string longer than 80 characters is supplied."""

        if text is None:
            return self._text
        else:
            if not isinstance(text, str):
                raise TypeError("PdbRecord needs str, not '%s'" % str(text))
            if len(text) > 80:
                raise ValueError("'%s' is longer than 80 characters" % str(text))
            self._text = text.rstrip()


    def name(self, name=None):
        """The name of the record, given by the first six characters of the line
        in the PDB file. The name will have any whitespace removed. If a string
        is passed as an argument, the name will be updated.

        :param str name: If given, the record's name will be updated to this.
        :raises ValueError: if a string longer than 6 characters is supplied."""

        if name is None:
            return self._text[:6].strip()
        else:
            if not isinstance(name, str):
                raise TypeError("PdbRecord needs str, not '%s'" % str(name))
            if len(name) > 6:
                raise ValueError("'%s' is longer than 6 characters" % str(name))
            self._text = name + (" " * (6 - len(name))) + self._text[6:]


    def body(self, body=None):
        """The body of the record - everything from column seven onwards (after
        the name). If a string is passed as an argument, the body will be
        updated, though any trailing whitespace will be removed.

        :param str body: If given, the record's body will be updated to this.
        :raises ValueError: if a string longer than 74 characters is supplied."""

        if body is None:
            return self._text[6:]
        else:
            if not isinstance(body, str):
                raise TypeError("PdbRecord needs str, not '%s'" % str(body))
            if len(body) > 74:
                raise ValueError("'%s' is longer than 74 characters" % str(body))
            self._text = self._text[:6] + body.rstrip()



class PdbFile:
    """Represents a Pdb file structure and the records it contains.

    Two PdbFiles are considered equal if they have the same number of records,
    and their records are equal. They do not need to share the same record
    objects in memory.

    PdbFiles are containers and iterables of the records they contain, and can
    be indexed and sliced to obtain them.

    :param \*records: The :py:class:`PdbRecord` objects that make up the file."""

    def __init__(self, *records):
        for record in records:
            if not isinstance(record, PdbRecord):
                raise TypeError("'%s' is not a PdbRecord" % str(record))
        self._records = list(records)
        for index, record in enumerate(self._records, start=1):
            record._number = index


    def __repr__(self):
        return "<PdbFile ({} records)>".format(self.length())


    def __len__(self):
        return self.length()


    def __eq__(self, other):
        if not isinstance(other, PdbFile) or self.length() != other.length():
            return False
        for record, other_record in zip(self._records, other._records):
            if record != other_record:
                return False
        return True


    def __contains__(self, member):
        return member in self._records


    def __iter__(self):
        return iter(self._records)


    def __getitem__(self, index):
        return self._records[index]


    def length(self):
        """Returns the number of records in this PdbFile.

        :rtype: ``int``"""

        return len(self._records)


    def records(self, name=None):
        """Returns the :py:class:`.PdbRecord` objects in this PdbFile. You can
        filter these by name if you wish.

        :param str name: The record name to filter by.
        :returns: ``list`` of ``PdbRecord``"""

        records = self._records
        if name:
            records = filter(lambda r: r.name().upper() == name.upper(), records)
        return tuple(records)


    def record(self, name):
        """Returns the first :py:class:`.PdbRecord` that matches the name given.

        :returns: ``PdbRecord``"""

        for record in self._records:
            if record.name().upper() == name.upper():
                return record


    def add_record(self, record):
        """Adds a :py:class:`.PdbRecord` to the end of the PdbFile.

        :param PdbRecord record: The record to add."""

        if not isinstance(record, PdbRecord):
            raise TypeError("'%s' is not a PdbRecord" % str(record))
        self._records.append(record)
        record._number = len(self._records)


    def remove_record(self, record):
        """Removes a :py:class:`.PdbRecord` from the PdbFile.

        :param PdbRecord record: The record to remove."""

        self._records.remove(record)
        record._number = None
        for index, record in enumerate(self._records, start=1):
            record._number = index
