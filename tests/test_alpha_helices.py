from unittest import TestCase
import unittest.mock
from molecupy.structures import AlphaHelix, Chain, ResiduicSequence, Residue
from molecupy import BrokenHelixError

class HelixTest(TestCase):

    def setUp(self):
        self.residues = [unittest.mock.Mock(spec=Residue) for _ in range(10)]
        self.chain = Chain("A", *self.residues[:5])
        self.chain = Chain("B", *self.residues[5:])
        for residue in self.residues:
            residue.chain.return_value = residue._chain



class HelixCreationTests(HelixTest):

    def test_can_create_helix(self):
        helix = AlphaHelix("AA", *self.residues[1:4])
        self.assertIsInstance(helix, ResiduicSequence)
        self.assertEqual(helix._helix_id, "AA")
        self.assertEqual(helix._residues, self.residues[1:4])
        self.assertEqual(helix._helix_class, None)
        self.assertEqual(helix._comment, None)


    def test_can_create_with_class_and_comment(self):
        helix = AlphaHelix("AA", *self.residues[1:4], helix_class=".", comment="..")
        self.assertEqual(helix._helix_class, ".")
        self.assertEqual(helix._comment, "..")


    def test_helix_id_must_be_str(self):
        with self.assertRaises(TypeError):
            AlphaHelix(100, *self.residues[1:4])


    def test_helix_residues_must_be_on_same_chain(self):
        with self.assertRaises(BrokenHelixError):
            AlphaHelix("AA", *self.residues[1:9])


    def test_helix_class_must_be_str(self):
        with self.assertRaises(TypeError):
            AlphaHelix("100", *self.residues[1:4], helix_class=1)


    def test_helix_comment_must_be_str(self):
        with self.assertRaises(TypeError):
            AlphaHelix("100", *self.residues[1:4], comment=1)



class HelixPropertyTests(HelixTest):

    def test_can_get_helix_properties(self):
        helix = AlphaHelix("AA", *self.residues[1:4], helix_class=".", comment="..")
        self.assertEqual(helix.helix_id(), "AA")
        self.assertEqual(helix.helix_class(), ".")
        self.assertEqual(helix.comment(), "..")
