from unittest import TestCase
from unittest.mock import Mock, patch
from molecupy.structures import GhostAtom, Atom, AtomicStructure

class AtomCreationTests(TestCase):

    def test_can_create_atom(self):
        atom = Atom(10.0, 20.0, 15.0, "C", 100, "CA")
        self.assertIsInstance(atom, GhostAtom)
        self.assertEqual(atom._x, 10.0)
        self.assertEqual(atom._y, 20.0)
        self.assertEqual(atom._z, 15.0)
        self.assertEqual(atom._element, "C")
        self.assertEqual(atom._atom_id, 100)
        self.assertEqual(atom._atom_name, "CA")
        self.assertEqual(atom._bonds, set())


    @patch("molecupy.structures.atoms.GhostAtom.__init__")
    def test_atom_uses_ghost_atom_initialisation(self, mock):
        Atom(10.0, 20.0, 15.0, "C", 100, "CA")
        self.assertTrue(mock.called)


    def test_coordinates_must_be_float(self):
        with self.assertRaises(TypeError):
            Atom("10.0", 20.0, 15.0, "C", 100, "CA")
        with self.assertRaises(TypeError):
            Atom(10.0, "20.0", 15.0, "C", 100, "CA")
        with self.assertRaises(TypeError):
            Atom(10.0, 20.0, "15.0", "C", 100, "CA")


    def test_repr(self):
        atom = Atom(10.0, 20.0, 15.0, "C", 100, "CA")
        self.assertEqual(str(atom), "<Atom 100 (CA)>")



class AtomLocationTests(TestCase):

    def setUp(self):
        self.atom = Atom(10.0, 20.0, 15.0, "C", 100, "CA")


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


    def test_can_get_all_coordinates(self):
        self.assertEqual(
         self.atom.location(),
         (self.atom.x(), self.atom.y(), self.atom.z())
        )



class AlternateLocationTests(TestCase):

    def test_by_default_atom_has_one_location(self):
        atom = Atom(10.0, 20.0, 15.0, "C", 100, "CA")
        self.assertEqual(
         atom.locations(),
         {(10.0, 20.0, 15.0): 1}
        )



class AtomDistanceTests(TestCase):

    def test_can_get_inter_atomic_distance(self):
        atom1 = Atom(-0.791, 64.789, 30.59, "O", 2621, "OD1") # Atom 2621 in 1LOL
        atom2 = Atom(5.132, 63.307, 56.785, "C", 1011, "CD") # Atom 1011 in 1LOL
        pymol_calculated_distance = 26.9
        self.assertAlmostEqual(
         atom1.distance_to(atom2),
         pymol_calculated_distance,
         delta=0.01
        )
        self.assertEqual(atom1.distance_to(atom2), atom2.distance_to(atom1))


    def test_can_only_measure_distance_between_localised_atoms(self):
        atom1 = Atom(-0.791, 64.789, 30.59, "O", 2621, "OD1")
        atom2 = Mock(GhostAtom)
        with self.assertRaises(TypeError):
            atom1.distance_to(atom2)


    def test_can_distance_to_atomic_structure_uses_center_of_mass(self):
        atom1 = Atom(-0.791, 64.789, 30.59, "O", 2621, "OD1")
        atomic_structure = Mock(AtomicStructure)
        atomic_structure.center_of_mass.return_value = (5.132, 63.307, 56.785)
        pymol_calculated_distance = 26.9
        self.assertAlmostEqual(
         atom1.distance_to(atomic_structure),
         pymol_calculated_distance,
         delta=0.01
        )



class AtomBondtests(TestCase):

    def setUp(self):
        self.atom1 = Atom(10.0, 20.0, 15.0, "C", 100, "CA")
        self.atom2 = Atom(10.0, 20.0, 17.0, "C", 100, "CA")
        self.atom3 = Atom(10.0, 20.0, 13.0, "C", 100, "CA")
        self.atom4 = Atom(10.0, 22.0, 17.0, "C", 100, "CA")
        self.atom5 = Atom(10.0, 18.0, 17.0, "C", 100, "CA")
        self.atom6 = Atom(10.0, 20.0, 19.0, "C", 100, "CA")


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


    def test_can_only_bond_atoms(self):
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


    def test_atom_can_get_all_atoms_covalently_accessible(self):
        atoms = [Atom(x / 10, x / 10, x / 10, "C", 1, "C") for x in range(11)]
        atoms[0].bond_to(atoms[1])
        atoms[1].bond_to(atoms[2])
        atoms[2].bond_to(atoms[3])
        atoms[3].bond_to(atoms[4])
        atoms[4].bond_to(atoms[5])
        atoms[5].bond_to(atoms[0])
        atoms[6].bond_to(atoms[3])
        atoms[7].bond_to(atoms[6])
        atoms[7].bond_to(atoms[8])
        atoms[7].bond_to(atoms[9])
        atoms[10].bond_to(atoms[8])
        atoms[10].bond_to(atoms[9])
        for atom in atoms:
            self.assertEqual(
             atom.accessible_atoms(),
             set([a for a in atoms if a is not atom])
            )



class NearbyAtomTests(TestCase):

    def setUp(self):
        self.atoms = [
         Atom(10.0, 20.0, float(i), "C", int(i), "CA") for i in range(10)
        ]
        molecule = Mock()
        self.model = Mock()
        molecule.model.return_value = self.model
        self.model.atoms.return_value = set(self.atoms)
        for atom in self.atoms:
            atom._molecule = molecule
            atom._element = "H" if atom.atom_id() == 3 else "C"


    def test_can_get_nearby_atoms(self):
        self.assertEqual(
         self.atoms[4].local_atoms(1),
         set([self.atoms[3], self.atoms[5]])
        )
        self.assertEqual(
         self.atoms[4].local_atoms(1.5),
         set([self.atoms[3], self.atoms[5]])
        )
        self.assertEqual(
         self.atoms[4].local_atoms(2),
         set([self.atoms[3], self.atoms[5], self.atoms[2], self.atoms[6]])
        )


    def test_can_exclude_hydrogens_from_nearby_atoms(self):
        self.assertEqual(
         self.atoms[4].local_atoms(1, include_hydrogens=False),
         set([self.atoms[5]])
        )
