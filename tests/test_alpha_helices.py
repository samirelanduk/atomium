from unittest import TestCase
import unittest.mock
from molecupy.structures import AlphaHelix, Chain, ResiduicSequence, Residue

class HelixTest(TestCase):

    def setUp(self):
        self.residues = [unittest.mock.Mock(spec=Residue) for _ in range(10)]



class HelixCreationTests(HelixTest):

    def test_can_create_helix(self):
        helix = AlphaHelix("AA", *self.residues)
        self.assertIsInstance(helix, ResiduicSequence)
        self.assertEqual(helix._helix_id, "AA")
        self.assertEqual(helix._residues, self.residues)
        self.assertEqual(helix._helix_class, None)
        self.assertEqual(helix._comment, None)


    def test_can_create_with_class_and_comment(self):
        helix = AlphaHelix("AA", *self.residues, helix_class=".", comment="..")
        self.assertEqual(helix._helix_class, ".")
        self.assertEqual(helix._comment, "..")
