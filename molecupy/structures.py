class Atom:

    def __init__(self, element, atom_id, atom_name):
        if not isinstance(element, str):
            raise TypeError("element must be str, not '%s'" % str(element))
        if not 0 < len(element) <= 2:
            raise ValueError("element must be of length 1 or 2, not %s" % element)
        self._element = element
        if not isinstance(atom_id, int):
            raise TypeError("atom_id must be int, not '%s'" % str(atom_id))
        self._atom_id = atom_id
        if not isinstance(atom_name, str):
            raise TypeError("atom_name must be str, not '%s'" % str(atom_name))
        self._atom_name = atom_name


    def __repr__(self):
        return "<Atom %i (%s)>" % (self._atom_id, self._atom_name)
