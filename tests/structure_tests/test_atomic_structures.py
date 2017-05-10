from unittest import TestCase
from unittest.mock import Mock, patch
from atomium.structures.atoms import Atom
from atomium.structures.molecules import AtomicStructure

class AtomicStructureTest(TestCase):

    def setUp(self):
        self.atom1 = Mock(Atom)
        self.atom2 = Mock(Atom)
        self.atom3 = Mock(Atom)
        self.atoms = [self.atom1, self.atom2, self.atom3]



class AtomicStructureCreationTests(AtomicStructureTest):

    def test_can_create_atomic_structure(self):
        structure = AtomicStructure(self.atom1, self.atom2, self.atom3)
        self.assertEqual(structure._atoms, set(self.atoms))


    def test_atomic_structure_needs_atoms(self):
        with self.assertRaises(TypeError):
            AtomicStructure("self.atom1", self.atom2, self.atom3)



class AtomicStructureReprTests(AtomicStructureTest):

    def test_atomic_structure_repr(self):
        structure = AtomicStructure(self.atom1, self.atom2, self.atom3)
        self.assertEqual(str(structure), "<AtomicStructure (3 atoms)>")



class AtomicStructureContainerTests(AtomicStructureTest):

    def test_atomic_structure_is_container_of_its_atoms(self):
        structure = AtomicStructure(self.atom1, self.atom2)
        self.assertIn(self.atom1, structure)
        self.assertIn(self.atom2, structure)
        self.assertNotIn(self.atom3, structure)



class AtomicStructureAtomsTests(AtomicStructureTest):

    def test_can_get_all_atoms(self):
        structure = AtomicStructure(self.atom1, self.atom2, self.atom3)
        self.assertEqual(structure.atoms(), structure._atoms)
        self.assertIsNot(structure.atoms(), structure._atoms)


    def test_can_get_atoms_by_element(self):
        structure = AtomicStructure(self.atom1, self.atom2, self.atom3)
        self.atom1.element.return_value = "A"
        self.atom2.element.return_value = "B"
        self.atom3.element.return_value = "B"
        self.assertEqual(structure.atoms(element="A"), set(self.atoms[:1]))
        self.assertEqual(structure.atoms(element="B"), set(self.atoms[1:]))
        self.assertEqual(structure.atoms(element="C"), set())



class AtomicStructureAtomTest(AtomicStructureTest):

    @patch("atomium.structures.molecules.AtomicStructure.atoms")
    def test_atom_calls_atoms(self, mock_atoms):
        mock_atoms.return_value = set(self.atoms[:1])
        structure = AtomicStructure(self.atom1, self.atom2, self.atom3)
        atom = structure.atom(element="A")
        mock_atoms.assert_called_with(element="A")
        self.assertIs(atom, self.atom1)


    @patch("atomium.structures.molecules.AtomicStructure.atoms")
    def test_atom_can_return_none(self, mock_atoms):
        structure = AtomicStructure(self.atom1, self.atom2, self.atom3)
        mock_atoms.return_value = set()
        self.assertIs(structure.atom(element="C"), None)



class AtomicStructureAtomAdditionTests(AtomicStructureTest):

    def test_can_add_atoms_to_structure(self):
        structure = AtomicStructure(self.atom1, self.atom2)
        structure.add_atom(self.atom3)
        self.assertEqual(structure._atoms, set(self.atoms))


    def test_can_only_add_atoms(self):
        structure = AtomicStructure(self.atom1, self.atom2)
        with self.assertRaises(TypeError):
            structure.add_atom("atom")



class AtomicStructureAtomRemovalTests(AtomicStructureTest):

    def test_can_remove_atoms_from_structure(self):
        structure = AtomicStructure(self.atom1, self.atom2)
        structure.remove_atom(self.atom1)
        self.assertEqual(structure._atoms, set(self.atoms[1:2]))
