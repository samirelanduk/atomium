from unittest import TestCase
import unittest.mock
from molecupy.structures import BetaStrand, Chain, ResiduicSequence, Residue
from molecupy import BrokenStrandError

class StrandTest(TestCase):

    def setUp(self):
        self.residues = [unittest.mock.Mock(spec=Residue) for _ in range(10)]
        self.chainA = Chain("A", *self.residues[:5])
        self.chainB = Chain("B", *self.residues[5:])
        for residue in self.residues:
            residue.chain.return_value = residue._chain



class StrandCreationTests(StrandTest):

    def test_can_create_strand(self):
        strand = BetaStrand("1", -1, *self.residues[1:4])
        self.assertIsInstance(strand, ResiduicSequence)
        self.assertEqual(strand._strand_id, "1")
        self.assertEqual(strand._sense, -1)
        self.assertEqual(strand._residues, self.residues[1:4])


    def test_strand_id_must_be_str(self):
        with self.assertRaises(TypeError):
            BetaStrand(100, -1, *self.residues[1:4])


    def test_strand_residues_must_be_on_same_chain(self):
        with self.assertRaises(BrokenStrandError):
            BetaStrand("100", -1, *self.residues[1:9])


    def test_sense_must_be_int(self):
        with self.assertRaises(TypeError):
            BetaStrand("100", -1.5, *self.residues[1:4])


    def test_sense_must_be_within_range(self):
        with self.assertRaises(ValueError):
            BetaStrand("100", -2, *self.residues[1:4])
        with self.assertRaises(ValueError):
            BetaStrand("100", 2, *self.residues[1:4])


    def test_strand_repr(self):
        strand = BetaStrand("1", -1, *self.residues[1:4])
        self.assertEqual(str(strand), "<BetaStrand 1 (3 residues)>")



class StrandPropertyTests(StrandTest):

    def test_can_get_strand_properties(self):
        strand = BetaStrand("1", -1, *self.residues[1:4])
        self.assertEqual(strand.strand_id(), "1")
        self.assertEqual(strand.sense(), -1)


    def test_cannot_add_residue_from_other_chain(self):
        strand = BetaStrand("1", -1, *self.residues[1:4])
        strand.add_residue(self.residues[4])
        with self.assertRaises(BrokenStrandError):
            strand.add_residue(self.residues[5])


    def test_can_update_strand_sense(self):
        strand = BetaStrand("1", -1, *self.residues[1:4])
        self.assertEqual(strand.sense(), -1)
        strand.sense(1)
        self.assertEqual(strand.sense(), 1)


    def test_strand_sense_can_only_be_set_to_int(self):
        strand = BetaStrand("1", -1, *self.residues[1:4])
        with self.assertRaises(TypeError):
            strand.sense(0.5)


    def test_strand_sense_can_only_be_set_to_valid_int(self):
        strand = BetaStrand("1", -1, *self.residues[1:4])
        with self.assertRaises(ValueError):
            strand.sense(-2)
        with self.assertRaises(ValueError):
            strand.sense(2)


    def test_can_get_chain(self):
        strand = BetaStrand("1", -1, *self.residues[1:4])
        self.assertIs(strand.chain(), self.chainA)
        strand = BetaStrand("1", -1, *self.residues[6:9])
        self.assertIs(strand.chain(), self.chainB)
