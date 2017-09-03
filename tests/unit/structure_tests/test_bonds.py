from unittest import TestCase
from unittest.mock import Mock, patch
from atomium.structures.atoms import Bond, Atom

class BondTest(TestCase):

    def setUp(self):
        self.atom1 = Mock(Atom)
        self.atom2 = Mock(Atom)
        self.atom1.element.return_value = "C"
        self.atom2.element.return_value = "N"
        self.atom1._bonds, self.atom2._bonds = set(), set()
        self.atoms = set([self.atom1, self.atom2])



class BondCreationTests(BondTest):

    def test_can_create_bond(self):
        bond = Bond(self.atom1, self.atom2)
        self.assertEqual(bond._atoms, self.atoms)


    def test_bond_atoms_must_be_atoms(self):
        with self.assertRaises(TypeError):
            Bond(self.atom1, "self.atom2")
        with self.assertRaises(TypeError):
            Bond("self.atom1", self.atom2)


    def test_bond_atoms_must_be_different(self):
        with self.assertRaises(ValueError):
            Bond(self.atom1, self.atom1)


    def test_creating_bonds_updates_atoms(self):
        bond = Bond(self.atom1, self.atom2)
        self.assertEqual(self.atom1._bonds, set([bond]))
        self.assertEqual(self.atom2._bonds, set([bond]))



class BondReprTests(BondTest):

    def test_bond_repr(self):
        bond = Bond(self.atom1, self.atom2)
        self.assertIn(str(bond), ("<C-N Bond>", "<N-C Bond>"))



class BondAtomsTests(BondTest):

    def test_bond_atoms(self):
        bond = Bond(self.atom1, self.atom2)
        self.assertEqual(bond.atoms(), bond._atoms)
        self.assertIsNot(bond.atoms(), bond._atoms)



class BondLengthTests(BondTest):

    def test_bond_length(self):
        self.atom1.distance_to.return_value = 10
        self.atom2.distance_to.return_value = 10
        bond = Bond(self.atom1, self.atom2)
        self.assertEqual(bond.length(), 10)



class BondVectorTests(BondTest):

    def test_bond_to_vector(self):
        bond = Bond(self.atom1, self.atom2)
        self.atom1.location.return_value = (1, 2, 3)
        self.atom2.location.return_value = (-5, 32, 9)
        vector = bond.vector(self.atom1)
        self.assertEqual(vector.values(), (6, -30, -6))
        vector = bond.vector(self.atom2)
        self.assertEqual(vector.values(), (-6, 30, 6))


    def test_vector_needs_atom(self):
        bond = Bond(self.atom1, self.atom2)
        with self.assertRaises(TypeError):
            bond.vector("atom")


    def test_vector_needs_present_atom(self):
        bond = Bond(self.atom1, self.atom2)
        with self.assertRaises(ValueError):
            bond.vector(Mock(Atom))



class BondDestructionTests(BondTest):

    def test_destroying_bond_removes_from_atoms(self):
        self.atom1._bonds.add("some bond")
        bond = Bond(self.atom1, self.atom2)
        self.assertEqual(self.atom1._bonds, set([bond, "some bond"]))
        self.assertEqual(self.atom2._bonds, set([bond]))
        bond.destroy()
        self.assertEqual(self.atom1._bonds, set(["some bond"]))
        self.assertEqual(self.atom2._bonds, set())
