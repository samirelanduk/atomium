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
