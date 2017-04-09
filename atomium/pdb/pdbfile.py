class PdbRecord:

    def __init__(self, text):
        if not isinstance(text, str):
            raise TypeError("PdbRecord needs str, not '%s'" % str(text))
        if len(text) > 80:
            raise ValueError("'%s' is longer than 80 characters" % str(text))
        self._text = text


    def text(self, text=None):
        if text is None:
            return self._text
        else:
            if not isinstance(text, str):
                raise TypeError("PdbRecord needs str, not '%s'" % str(text))
            if len(text) > 80:
                raise ValueError("'%s' is longer than 80 characters" % str(text))
            self._text = text
