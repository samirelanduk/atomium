class PdbRecord:

    def __init__(self, line):
        if not isinstance(line, str):
            raise TypeError("PdbRecord needs str, not '%s'" % str(line))
        if len(line) > 80:
            raise ValueError("'%s' is longer than 80 characters" % str(line))
        self._text = line
