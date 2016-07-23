from unittest import TestCase
import unittest.mock
from molecupy.structures import Chain, ResiduicSequence, Residue

class ChainTest(TestCase):

    def setUp(self):
        self.residues = [unittest.mock.Mock(spec=Residue) for _ in range(10)]



class ChainCreationTests(ChainTest):

    def test_can_create_small_molecule(self):
        chain = Chain("A", *self.residues)
        self.assertIsInstance(chain, ResiduicSequence)
        self.assertEqual(chain._chain_id, "A")
        self.assertEqual(chain._residues, self.residues)
