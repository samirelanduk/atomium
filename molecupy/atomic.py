class Atom:

    def __init__(self, x, y, z, element):
        self.x = x
        self.y = y
        self.z = z
        self.element = element


    def __repr__(self):
        return "<Atom (%s)>" % self.element
