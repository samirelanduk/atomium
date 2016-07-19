from unittest import TestCase
import unittest.mock
from molecupy.structures import AtomicStructure, PdbAtom

class AtomicStructureTest(TestCase):

    def setUp(self):
        self.pdb_atoms = [unittest.mock.Mock(spec=PdbAtom) for _ in range(10)]



class AtomicStructureCreationTests(AtomicStructureTest):

    def test_can_create_atomic_structure_with_pdb_atoms(self):
        atomic_structure = AtomicStructure(*self.pdb_atoms)
        self.assertEqual(atomic_structure._atoms, set(self.pdb_atoms))
