import math
from omnicanvas.canvas import Canvas
import omnicanvas.graphics
from unittest import TestCase
import unittest.mock
from molecupy.structures import Chain, ResiduicSequence, Residue, Atom, GhostAtom
from molecupy.structures import AlphaHelix, BetaStrand

class ChainTest(TestCase):

    def setUp(self):
        self.residues = [unittest.mock.Mock(spec=Residue) for _ in range(10)]



class ChainCreationTests(ChainTest):

    def test_can_create_chain(self):
        chain = Chain("A", *self.residues)
        self.assertIsInstance(chain, ResiduicSequence)
        self.assertEqual(chain._chain_id, "A")
        self.assertEqual(chain._residues, self.residues)
        self.assertEqual(chain._alpha_helices, set())
        self.assertEqual(chain._beta_strands, set())
        self.assertEqual(chain._model, None)


    def test_chain_updates_residues(self):
        for residue in self.residues:
            residue._chain = None
        chain = Chain("A", *self.residues)
        for residue in self.residues:
            self.assertIs(residue._chain, chain)


    def test_chain_id_must_be_str(self):
        with self.assertRaises(TypeError):
            Chain(200, *self.residues)


    def test_chain_repr(self):
        chain = Chain("A", *self.residues)
        self.assertEqual(str(chain), "<Chain A (10 residues)>")



class ChainPropertyTests(ChainTest):

    def test_chain_properties(self):
        chain = Chain("A", *self.residues)
        self.assertEqual(chain.chain_id(), "A")
        self.assertEqual(chain.alpha_helices(), set())
        self.assertEqual(chain.beta_strands(), set())
        self.assertEqual(chain.model(), None)


    def test_can_add_residues_and_update_them(self):
        chain = Chain("A", *self.residues)
        residue = unittest.mock.Mock(spec=Residue)
        residue._chain = None
        chain.add_residue(residue)
        self.assertIs(chain.residues()[-1], residue)
        self.assertIs(residue._chain, chain)


    def test_can_remove_residues_and_update_them(self):
        chain = Chain("A", *self.residues)
        chain.remove_residue(self.residues[5])
        self.assertNotIn(self.residues[5], chain.residues())
        self.assertIs(self.residues[5]._chain, None)



class SecondaryStructureRetrievalTests(ChainTest):

    def setUp(self):
        ChainTest.setUp(self)
        self.chain = Chain("A", *self.residues)
        self.alpha_helices = [unittest.mock.Mock(spec=AlphaHelix) for _ in range(5)]
        for index, helix in enumerate(self.alpha_helices):
            helix.helix_id.return_value = str(index + 1)
        self.chain._alpha_helices = set(self.alpha_helices)
        self.beta_strands = [unittest.mock.Mock(spec=BetaStrand) for _ in range(5)]
        for index, strand in enumerate(self.beta_strands):
            strand.strand_id.return_value = str(index + 1)
        self.chain._beta_strands = set(self.beta_strands)


    def test_can_get_helix_by_id(self):
        self.assertIs(self.chain.get_helix_by_id("1"), self.alpha_helices[0])
        self.assertIs(self.chain.get_helix_by_id("2"), self.alpha_helices[1])
        self.assertIs(self.chain.get_helix_by_id("3"), self.alpha_helices[2])
        self.assertIs(self.chain.get_helix_by_id("4"), self.alpha_helices[3])
        self.assertIs(self.chain.get_helix_by_id("5"), self.alpha_helices[4])
        self.assertIs(self.chain.get_helix_by_id("6"), None)


    def test_can_only_search_for_string_helix_id(self):
        with self.assertRaises(TypeError):
            self.chain.get_helix_by_id(1)


    def test_can_get_strand_by_id(self):
        self.assertIs(self.chain.get_strand_by_id("1"), self.beta_strands[0])
        self.assertIs(self.chain.get_strand_by_id("2"), self.beta_strands[1])
        self.assertIs(self.chain.get_strand_by_id("3"), self.beta_strands[2])
        self.assertIs(self.chain.get_strand_by_id("4"), self.beta_strands[3])
        self.assertIs(self.chain.get_strand_by_id("5"), self.beta_strands[4])
        self.assertIs(self.chain.get_strand_by_id("6"), None)


    def test_can_only_search_for_string_strand_id(self):
        with self.assertRaises(TypeError):
            self.chain.get_strand_by_id(1)
