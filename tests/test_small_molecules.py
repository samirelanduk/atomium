
from unittest import TestCase
import unittest.mock
from molecupy.structures import SmallMolecule, AtomicStructure, PdbAtom

class SmallMoleculeTest(TestCase):

    def setUp(self):
        self.atoms = [unittest.mock.Mock(spec=PdbAtom) for _ in range(10)]



class SmallMoleculeCreationTests(SmallMoleculeTest):

    def test_can_create_small_molecule(self):
        small_molecule = SmallMolecule("A500", "MOL", *self.atoms)
        self.assertIsInstance(small_molecule, AtomicStructure)
        self.assertEqual(small_molecule._molecule_id, "A500")
        self.assertEqual(small_molecule._molecule_name, "MOL")
        self.assertEqual(small_molecule._atoms, set(self.atoms))
