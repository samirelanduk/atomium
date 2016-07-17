from unittest import TestCase
from molecupy.structures import Atom, PdbAtom

class PdbAtomCreationTests(TestCase):

    def test_can_create_pdb_atom(self):
        atom = PdbAtom(10.0, 20.0, 15.0, "C", 100, "CA")
        self.assertIsInstance(atom, Atom)
        self.assertEqual(atom._x, 10.0)
        self.assertEqual(atom._y, 20.0)
        self.assertEqual(atom._z, 15.0)
        self.assertEqual(atom._element, "C")
        self.assertEqual(atom._atom_id, 100)
        self.assertEqual(atom._atom_name, "CA")
        self.assertEqual(atom._bonds, set())


    def test_repr(self):
        atom = PdbAtom(10.0, 20.0, 15.0, "C", 100, "CA")
        self.assertEqual(str(atom), "<PdbAtom 100 (CA)>")


    def test_coordinates_must_be_float(self):
        with self.assertRaises(TypeError):
            atom = PdbAtom("10.0", 20.0, 15.0, "C", 100, "CA")
        with self.assertRaises(TypeError):
            atom = PdbAtom(10.0, "20.0", 15.0, "C", 100, "CA")
        with self.assertRaises(TypeError):
            atom = PdbAtom(10.0, 20.0, "15.0", "C", 100, "CA")



class PdbAtomPropertyTests(TestCase):

    def setUp(self):
        self.atom = PdbAtom(10.0, 20.0, 15.0, "C", 100, "CA")


    def test_basic_properties(self):
        self.assertEqual(self.atom.x(), 10.0)
        self.assertEqual(self.atom.y(), 20.0)
        self.assertEqual(self.atom.z(), 15.0)


    def test_can_set_coordinates(self):
        self.atom.x(1000.0)
        self.assertEqual(self.atom.x(), 1000.0)
        self.atom.y(1000.0)
        self.assertEqual(self.atom.y(), 1000.0)
        self.atom.z(1000.0)
        self.assertEqual(self.atom.z(), 1000.0)


    def test_coordinates_must_be_float(self):
        with self.assertRaises(TypeError):
            self.atom.x("10.0")
        with self.assertRaises(TypeError):
            self.atom.y("10.0")
        with self.assertRaises(TypeError):
            self.atom.z("10.0")



class PdbAtomDistanceTests(TestCase):

    def test_can_get_inter_atomic_distance(self):
        atom1 = PdbAtom(-0.791, 64.789, 30.59, "O", 2621, "OD1") # Atom 2621 in 1LOL
        atom2 = PdbAtom(5.132, 63.307, 56.785, "C", 1011, "CD") # Atom 1011 in 1LOL
        pymol_calculated_distance = 26.9
        self.assertAlmostEqual(
         atom1.distance_to(atom2),
         pymol_calculated_distance,
         delta=0.01
        )


    def test_can_only_measure_distance_between_pdb_atoms(self):
        atom1 = PdbAtom(-0.791, 64.789, 30.59, "O", 2621, "OD1")
        atom2 = Atom("C", 1011, "CD")
        with self.assertRaises(TypeError):
            atom1.distance_to(atom2)



class PdbAtomBondtests(TestCase):

    def setUp(self):
        self.atom1 = PdbAtom(10.0, 20.0, 15.0, "C", 100, "CA")
        self.atom2 = PdbAtom(10.0, 20.0, 17.0, "C", 100, "CA")
        self.atom3 = PdbAtom(10.0, 20.0, 13.0, "C", 100, "CA")
        self.atom4 = PdbAtom(10.0, 22.0, 17.0, "C", 100, "CA")
        self.atom5 = PdbAtom(10.0, 18.0, 17.0, "C", 100, "CA")
        self.atom6 = PdbAtom(10.0, 20.0, 19.0, "C", 100, "CA")


    def test_can_bond_atoms(self):
        self.atom1.bond_to(self.atom2)
        self.assertEqual(self.atom1.bonds(), self.atom2.bonds())
        self.assertEqual(len(self.atom1.bonds()), 1)
        self.assertEqual(self.atom1.bonded_atoms(), set((self.atom2,)))
        self.assertEqual(self.atom2.bonded_atoms(), set((self.atom1,)))


    def test_can_bond_many_atoms(self):
        self.atom1.bond_to(self.atom2)
        self.atom1.bond_to(self.atom3)
        self.atom1.bond_to(self.atom4)
        self.atom1.bond_to(self.atom5)
        self.atom2.bond_to(self.atom6)
        self.assertEqual(len(self.atom1.bonds()), 4)
        self.assertIn(self.atom2, self.atom1.bonded_atoms())
        self.assertIn(self.atom3, self.atom1.bonded_atoms())
        self.assertIn(self.atom4, self.atom1.bonded_atoms())
        self.assertIn(self.atom5, self.atom1.bonded_atoms())
        self.assertEqual(len(self.atom2.bonds()), 2)
        self.assertIn(self.atom1, self.atom2.bonded_atoms())
        self.assertIn(self.atom6, self.atom2.bonded_atoms())
        self.assertEqual(len(self.atom3.bonds()), 1)
        self.assertIn(self.atom1, self.atom3.bonded_atoms())
        self.assertEqual(len(self.atom4.bonds()), 1)
        self.assertIn(self.atom1, self.atom4.bonded_atoms())
        self.assertEqual(len(self.atom5.bonds()), 1)
        self.assertIn(self.atom1, self.atom5.bonded_atoms())
        self.assertEqual(len(self.atom6.bonds()), 1)
        self.assertIn(self.atom2, self.atom6.bonded_atoms())


    def test_can_only_bond_pdb_atoms(self):
        with self.assertRaises(TypeError):
            self.atom1.bond_to("atom")


    def test_can_get_bond_between_atoms(self):
        self.assertEqual(self.atom1.get_bond_with(self.atom2), None)
        self.assertEqual(self.atom2.get_bond_with(self.atom1), None)
        self.atom1.bond_to(self.atom2)
        self.assertIn(self.atom1, self.atom2.bonded_atoms())
        self.assertIs(
         self.atom1.get_bond_with(self.atom2),
         self.atom2.get_bond_with(self.atom1)
        )
        self.assertNotEqual(self.atom1.get_bond_with(self.atom2), None)


    def test_can_delete_bonds(self):
        self.atom1.bond_to(self.atom2)
        self.atom2.bond_to(self.atom6)
        self.assertIn(self.atom1, self.atom2.bonded_atoms())
        self.assertIn(self.atom6, self.atom2.bonded_atoms())
        self.atom2.break_bond_with(self.atom1)
        self.assertIn(self.atom6, self.atom2.bonded_atoms())
        self.assertNotIn(self.atom1, self.atom2.bonded_atoms())
        self.assertNotIn(self.atom2, self.atom1.bonded_atoms())
