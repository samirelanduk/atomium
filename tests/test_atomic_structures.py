from unittest import TestCase
from molecupy import exceptions
from molecupy.atomic import AtomicStructure, Atom

class AtomicStructureTest(TestCase):

    def setUp(self):
        self.atom1 = Atom(1.0, 1.0, 1.0, "H", atom_id=1, atom_name="H1")
        self.atom2 = Atom(1.0, 1.0, 2.0, "C", atom_id=2, atom_name="CA")
        self.atom3 = Atom(1.0, 1.0, 3.0, "O", atom_id=3, atom_name="OX1")


    def check_valid_atomic_structure(self, atomic_structure):
        self.assertIsInstance(atomic_structure, AtomicStructure)
        self.assertIsInstance(atomic_structure.atoms, set)
        for atom in atomic_structure.atoms:
            self.assertIsInstance(atom, Atom)
        self.assertRegex(str(atomic_structure), r"<AtomicStructure \((\d+) atoms\)>")



class AtomicStructureCreationTests(AtomicStructureTest):

    def test_can_create_atomic_structure(self):
        atomic_structure = AtomicStructure(self.atom1, self.atom2, self.atom3)
        self.check_valid_atomic_structure(atomic_structure)


    def test_atomic_structure_needs_atoms(self):
        with self.assertRaises(TypeError):
            atomic_structure = AtomicStructure(self.atom1, self.atom2, "atom3")


    def test_atomic_structure_needs_at_least_one_atom(self):
        with self.assertRaises(exceptions.NoAtomsError):
            atomic_structure = AtomicStructure()


    def test_atom_structure_needs_unique_atom_ids(self):
        atom4 = Atom(1.0, 1.0, 4.0, "H", atom_id=1, atom_name="H1")
        atom5 = Atom(1.0, 1.0, 5.0, "H", atom_name="H1")
        atom6 = Atom(1.0, 1.0, 6.0, "H", atom_name="H1")
        with self.assertRaises(exceptions.DuplicateAtomIdError):
            atomic_structure = AtomicStructure(
             self.atom1, self.atom2, self.atom3, atom4, atom5, atom6
            )
        atom4.atom_id = 4
        atomic_structure = AtomicStructure(
         self.atom1, self.atom2, self.atom3, atom4, atom5, atom6
        )


    def test_atoms_in_atomic_structure(self):
        atomic_structure = AtomicStructure(self.atom1, self.atom2, self.atom3)
        self.assertEqual(len(atomic_structure.atoms), 3)
        self.assertIn(self.atom1, atomic_structure)
        self.assertIn(self.atom2, atomic_structure)
        self.assertIn(self.atom3, atomic_structure)



class AtomicStructureBehaviorTests(AtomicStructureTest):

    def test_can_get_atomic_structure_mass(self):
        atomic_structure = AtomicStructure(self.atom1, self.atom2, self.atom3)
        self.assertAlmostEqual(atomic_structure.get_mass(), 1 + 12 + 16, delta=0.5)



class AtomicStructureAtomRetrievalTests(AtomicStructureTest):

    def test_can_get_atom_by_id(self):
        atomic_structure = AtomicStructure(self.atom1, self.atom2, self.atom3)
        self.assertIs(atomic_structure.get_atom_by_id(1), self.atom1)
        self.assertIs(atomic_structure.get_atom_by_id(2), self.atom2)
        self.assertIs(atomic_structure.get_atom_by_id(3), self.atom3)
        self.assertIs(atomic_structure.get_atom_by_id(4), None)


    def test_can_only_search_by_numeric_id(self):
        atomic_structure = AtomicStructure(self.atom1, self.atom2, self.atom3)
        with self.assertRaises(TypeError):
            atomic_structure.get_atom_by_id(None)


    def test_can_get_atom_by_name(self):
        atom4 = Atom(1.0, 1.0, 4.0, "H", atom_id=4, atom_name="H1")
        atomic_structure = AtomicStructure(self.atom1, self.atom2, self.atom3, atom4)
        self.assertIn(atomic_structure.get_atom_by_name("H1"), (self.atom1, atom4))
        self.assertIs(atomic_structure.get_atom_by_name("xxx"), None)


    def test_can_get_multiple_atoms_by_name(self):
        atom4 = Atom(1.0, 1.0, 4.0, "H", atom_id=4, atom_name="H1")
        atomic_structure = AtomicStructure(self.atom1, self.atom2, self.atom3, atom4)
        self.assertEqual(
         atomic_structure.get_atoms_by_name("H1"),
         set([self.atom1, atom4])
        )
        self.assertEqual(atomic_structure.get_atoms_by_name("XXX"), set())


    def test_can_only_search_by_string_name(self):
        atomic_structure = AtomicStructure(self.atom1, self.atom2, self.atom3)
        with self.assertRaises(TypeError):
            atomic_structure.get_atom_by_name(None)
        with self.assertRaises(TypeError):
            atomic_structure.get_atoms_by_name(None)


    def test_can_get_multiple_atoms_by_element(self):
        atom4 = Atom(1.0, 1.0, 4.0, "H", atom_id=4, atom_name="H2")
        atomic_structure = AtomicStructure(self.atom1, self.atom2, self.atom3, atom4)
        self.assertEqual(
         atomic_structure.get_atoms_by_element("H"),
         set([self.atom1, atom4])
        )
        self.assertEqual(
         atomic_structure.get_atoms_by_element("C"),
         set([self.atom2])
        )
        self.assertEqual(atomic_structure.get_atoms_by_element("XXX"), set())


    def test_can_only_search_by_string_element(self):
        atomic_structure = AtomicStructure(self.atom1, self.atom2, self.atom3)
        with self.assertRaises(TypeError):
            atomic_structure.get_atoms_by_element(None)
