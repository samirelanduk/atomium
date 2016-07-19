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
