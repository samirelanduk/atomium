from unittest import TestCase
import unittest.mock
from molecupy.structures import ResiduicSequence, ResiduicStructure, Residue, Atom

class ResiduicSequenceTest(TestCase):

    def setUp(self):
        self.residues = [unittest.mock.Mock(spec=Residue) for _ in range(10)]



class ResiduicSequenceCreationTests(ResiduicSequenceTest):

    def test_can_create_residuic_structure(self):
        residuic_sequence = ResiduicSequence(*self.residues)
        self.assertIsInstance(residuic_sequence, ResiduicStructure)
        self.assertEqual(residuic_sequence._residues, self.residues)


    def test_can_only_create_residuic_sequence_with_residues(self):
        with self.assertRaises(TypeError):
            ResiduicSequence("Atom1", "Atom2")


    def test_residuic_structure_repr(self):
        residuic_sequence = ResiduicSequence(*self.residues)
        self.assertEqual(str(residuic_sequence), "<ResiduicSequence (10 residues)>")
