import math

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
        return "<%s %i (%s)>" % (
         self.__class__.__name__, self._atom_id, self._atom_name
        )


    def element(self, element=None):
        if element is None:
            return self._element
        else:
            if not isinstance(element, str):
                raise TypeError("element must be str, not '%s'" % str(element))
            if not 0 < len(element) <= 2:
                raise ValueError(
                 "element must be of length 1 or 2, not %s" % element
                )
            self._element = element


    def atom_id(self, atom_id=None):
        if atom_id is None:
            return self._atom_id
        else:
            if not isinstance(atom_id, int):
                raise TypeError("atom_id must be int, not '%s'" % str(atom_id))
            self._atom_id = atom_id


    def atom_name(self, atom_name=None):
        if atom_name is None:
            return self._atom_name
        else:
            if not isinstance(atom_name, str):
                raise TypeError("atom_name must be str, not '%s'" % str(atom_name))
            self._atom_name = atom_name


    def mass(self):
        return PERIODIC_TABLE.get(self.element().upper(), 0)



class PdbAtom(Atom):

    def __init__(self, x, y, z, *args):
        if not isinstance(x, float):
            raise TypeError("x coordinate must be float, not '%s'" % str(x))
        self._x = x
        if not isinstance(y, float):
            raise TypeError("y coordinate must be float, not '%s'" % str(y))
        self._y = y
        if not isinstance(z, float):
            raise TypeError("z coordinate must be float, not '%s'" % str(z))
        self._z = z
        Atom.__init__(self, *args)


    def x(self, x=None):
        if x is None:
            return self._x
        else:
            if not isinstance(x, float):
                raise TypeError("x coordinate must be float, not '%s'" % str(x))
            self._x = x


    def y(self, y=None):
        if y is None:
            return self._y
        else:
            if not isinstance(y, float):
                raise TypeError("y coordinate must be float, not '%s'" % str(y))
            self._y = y


    def z(self, z=None):
        if z is None:
            return self._z
        else:
            if not isinstance(z, float):
                raise TypeError("z coordinate must be float, not '%s'" % str(z))
            self._z = z


    def distance_to(self, other_atom):
        if not isinstance(other_atom, PdbAtom):
            raise TypeError(
             "Can only get distance between PdbAtoms, not '%s'" % str(other_atom)
            )
        x_sum = math.pow((other_atom.x() - self.x()), 2)
        y_sum = math.pow((other_atom.y() - self.y()), 2)
        z_sum = math.pow((other_atom.z() - self.z()), 2)
        distance = math.sqrt(x_sum + y_sum + z_sum)
        return distance



class CovalentBond:

    def __init__(self, atom1, atom2):
        if not isinstance(atom1, PdbAtom) or not isinstance(atom2, PdbAtom):
            raise TypeError(
             "Can only create covalent bond between PdbAtoms, not %s and %s" % (
              str(atom1), str(atom2)
             )
            )
        if atom1 is atom2:
            raise ValueError("Cannot covalent bond %s to itself." % str(atom1))
        self._atoms = set((atom1, atom2))



PERIODIC_TABLE = {
 "H": 1.0079, "HE": 4.0026, "LI": 6.941, "BE": 9.0122, "B": 10.811,
 "C": 12.0107, "N": 14.0067, "O": 15.9994, "F": 18.9984, "NE": 20.1797,
 "NA": 22.9897, "MG": 24.305, "AL": 26.9815, "SI": 28.0855, "P": 30.9738,
 "S": 32.065, "CL": 35.453, "K": 39.0983, "AR": 39.948, "CA": 40.078,
 "SC": 44.9559, "TI": 47.867, "V": 50.9415, "CR": 51.9961, "MN": 54.938,
 "FE": 55.845, "NI": 58.6934, "CO": 58.9332, "CU": 63.546, "ZN": 65.39,
 "GA": 69.723, "GE": 72.64, "AS": 74.9216, "SE": 78.96, "BR": 79.904,
 "KR": 83.8, "RB": 85.4678, "SR": 87.62, "Y": 88.9059, "ZR": 91.224,
 "NB": 92.9064, "MO": 95.94, "TC": 98, "RU": 101.07, "RH": 102.9055,
 "PD": 106.42, "AG": 107.8682, "CD": 112.411, "IN": 114.818, "SN": 118.71,
 "SB": 121.76, "I": 126.9045, "TE": 127.6, "XE": 131.293, "CS": 132.9055,
 "BA": 137.327, "LA": 138.9055, "CE": 140.116, "PR": 140.9077, "ND": 144.24,
 "PM": 145, "SM": 150.36, "EU": 151.964, "GD": 157.25, "TB": 158.9253,
 "DY": 162.5, "HO": 164.9303, "ER": 167.259, "TM": 168.9342, "YB": 173.04,
 "LU": 174.967, "HF": 178.49, "TA": 180.9479, "W": 183.84, "RE": 186.207,
 "OS": 190.23, "IR": 192.217, "PT": 195.078, "AU": 196.9665, "HG": 200.59,
 "TL": 204.3833, "PB": 207.2, "BI": 208.9804, "PO": 209, "AT": 210, "RN": 222,
 "FR": 223, "RA": 226, "AC": 227, "PA": 231.0359, "TH": 232.0381, "NP": 237,
 "U": 238.0289, "AM": 243, "PU": 244, "CM": 247, "BK": 247, "CF": 251,
 "ES": 252, "FM": 257, "MD": 258, "NO": 259, "RF": 261, "LR": 262, "DB": 262,
 "BH": 264, "SG": 266, "MT": 268, "RG": 272, "HS": 277
}
