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
