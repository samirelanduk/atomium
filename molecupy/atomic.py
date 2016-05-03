class Atom:

    def __init__(self, x, y, z, element):
        for coordinate in (x, y, z):
            if not isinstance(coordinate, float): raise TypeError
        self.x = x
        self.y = y
        self.z = z
        self.element = element


    def __repr__(self):
        return "<Atom (%s)>" % self.element
