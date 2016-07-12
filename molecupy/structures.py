class Atom:

    def __init__(self, element, atom_id, atom_name):
        self._element = element
        self._atom_id = atom_id
        self._atom_name = atom_name


    def __repr__(self):
        return "<Atom %i (%s)>" % (self._atom_id, self._atom_name)
