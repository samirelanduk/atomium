from unittest import TestCase
import unittest.mock
from molecupy.structures import ResiduicStructure, Residue, AtomicStructure

class ResiduicStructureTest(TestCase):

    def setUp(self):
        self.residues = [unittest.mock.Mock(spec=Residue) for _ in range(10)]



class ResiduicStructureCreationTests(ResiduicStructureTest):

    def test_can_create_residuic_structure(self):
        residuic_structure = ResiduicStructure(*self.residues)
        self.assertIsInstance(residuic_structure, AtomicStructure)
        self.assertEqual(residuic_structure._residues, set(self.residues))


    def test_residuic_structure_repr(self):
        residuic_structure = ResiduicStructure(*self.residues)
        self.assertEqual(str(residuic_structure), "<ResiduicStructure (10 residues)>")
