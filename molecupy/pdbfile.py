class PdbRecord:

    def __init__(self, text, number):
        self._number = number
        expanded_text = text[:80] if len(text) >= 80 else text + (" " * (80 - len(text)))
        self._text = expanded_text
        self._name = expanded_text[:6].strip()
        self._content = expanded_text[6:]
