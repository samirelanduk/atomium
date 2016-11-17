class Complex:

    def __init__(self, complex_id, complex_name, *chains):
        self._complex_id = complex_id
        self._complex_name = complex_name
        self._chains = set(chains)
