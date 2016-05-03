from .exceptions import InvalidElementError

class Atom:

    def __init__(self, x, y, z, element, atom_id=None, atom_name=None):
        for coordinate in (x, y, z):
            if not isinstance(coordinate, float):
                raise TypeError("'%s' is not a valid coordinate" % str(coordinate))
        self.x = x
        self.y = y
        self.z = z

        if not isinstance(element, str):
            raise TypeError("'%s' is not a valid element" % str(element))
        if len(element) == 0:
            raise InvalidElementError("Atom's element can't be an empty string")
        elif len(element) >= 2:
            raise InvalidElementError("'%s' is not a valid element" % element)
        self.element = element

        if not isinstance(atom_id, int) and atom_id is not None:
            raise TypeError("'%s' is not a valid atom_id" % str(atom_id))
        self.atom_id = atom_id
        
        if not isinstance(atom_name, str) and atom_name is not None:
            raise TypeError("'%s' is not a valid atom_name" % str(atom_name))
        self.atom_name = atom_name


    def __repr__(self):
        return "<Atom (%s)>" % self.element
