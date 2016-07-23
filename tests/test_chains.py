from unittest import TestCase
import unittest.mock
from molecupy.structures import Chain, ResiduicSequence, Residue

class ChainTest(TestCase):

    def setUp(self):
        self.residues = [unittest.mock.Mock(spec=Residue) for _ in range(10)]



class ChainCreationTests(ChainTest):

    def test_can_create_chain(self):
        chain = Chain("A", *self.residues)
        self.assertIsInstance(chain, ResiduicSequence)
        self.assertEqual(chain._chain_id, "A")
        self.assertEqual(chain._residues, self.residues)


    def test_chain_id_must_be_str(self):
        with self.assertRaises(TypeError):
            Chain(200, *self.residues)


    def test_chain_repr(self):
        chain = Chain("A", *self.residues)
        self.assertEqual(str(chain), "<Chain A (10 residues)>")
