from unittest import TestCase
from omnicanvas.canvas import Canvas
from molecupy import exceptions
from molecupy.structures import PdbChain, ResiduicSequence, PdbAtom, PdbResidue

class ChainTest(TestCase):

    def setUp(self):
        self.atom1 = PdbAtom(1.0, 1.0, 1.0, "H", 1, "H1")
        self.atom2 = PdbAtom(1.0, 1.0, 2.0, "C", 2, "CA")
        self.atom3 = PdbAtom(1.0, 1.0, 3.0, "O", 3, "OX1")
        self.residue1 = PdbResidue("A1", "ARG", self.atom1, self.atom2, self.atom3)
        self.atom4 = PdbAtom(1.0, 1.0, 4.0, "H", 4, "H1")
        self.atom5 = PdbAtom(1.0, 1.0, 5.0, "C", 5, "CA")
        self.atom6 = PdbAtom(1.0, 1.0, 6.0, "O", 6, "OX1")
        self.residue2 = PdbResidue("A2", "HST", self.atom4, self.atom5, self.atom6)
        self.atom7 = PdbAtom(1.0, 1.0, 7.0, "H", 7, "H1")
        self.atom8 = PdbAtom(1.0, 1.0, 8.0, "C", 8, "CA")
        self.atom9 = PdbAtom(1.0, 1.0, 9.0, "O", 9, "OX1")
        self.residue3 = PdbResidue("A3", "TRP", self.atom7, self.atom8, self.atom9)


    def check_valid_chain(self, chain):
        self.assertIsInstance(chain, PdbChain)
        self.assertIsInstance(chain, ResiduicSequence)
        self.assertIsInstance(chain.chain_id, str)
        self.assertIsInstance(chain.residues, list)
        for residue in chain.residues:
            self.assertIs(residue.chain, chain)
        chain.model
        self.assertIsInstance(chain.missing_residues, list)
        self.assertRegex(str(chain), r"<Chain [A-Z] \((\d+) residues\)>")



class ChainCreationTests(ChainTest):

    def test_can_create_chain(self):
        chain = PdbChain("A", self.residue1, self.residue2, self.residue3)
        self.check_valid_chain(chain)


    def test_chain_id_must_be_str(self):
        with self.assertRaises(TypeError):
            chain = PdbChain(1, self.residue1, self.residue2, self.residue3)
        with self.assertRaises(TypeError):
            chain = PdbChain(None, self.residue1, self.residue2, self.residue3)



class ChainSequenceTests(ChainTest):

    def test_can_get_sequence_with_missing_residues(self):
        chain = PdbChain("A", self.residue1, self.residue2, self.residue3)
        chain.missing_residues = ["A2A", "A2B", "A3A", "A4", "A5", "A5A"]
        self.assertEqual(
         chain.get_residue_ids_including_missing(),
         ["A1", "A2", "A2A", "A2B", "A3", "A3A", "A4", "A5", "A5A"]
        )



class ChainMatrixTests(ChainTest):

    def setUp(self):
        self.atom1 = PdbAtom(3.696, 33.898, 63.219, "N", 1, "N")
        self.atom2 = PdbAtom(3.198, 33.218, 61.983, "C", 2, "CA")
        self.atom3 = PdbAtom(3.914, 31.863, 61.818, "O", 3, "C")
        self.residue1 = PdbResidue("A11", "VAL", self.atom1, self.atom2, self.atom3)
        self.atom4 = PdbAtom(3.155, 30.797, 61.557, "N", 4, "N")
        self.atom5 = PdbAtom(3.728, 29.464, 61.400, "C", 5, "CA")
        self.atom6 = PdbAtom(4.757, 29.459, 60.275, "O", 6, "C")
        self.residue2 = PdbResidue("A12", "MET", self.atom4, self.atom5, self.atom6)
        self.atom7 = PdbAtom(5.980, 29.039, 60.600, "N", 7, "N")
        self.atom8 = PdbAtom(7.092, 28.983, 59.649, "C", 8, "CA")
        self.atom9 = PdbAtom(7.092, 30.310, 58.987, "O", 9, "C")
        self.residue3 = PdbResidue("A13", "ASN", self.atom7, self.atom8, self.atom9)
        self.chain = PdbChain("A", self.residue1, self.residue2, self.residue3)
        self.chain.missing_residues = ["A12A", "A14"]


    def check_valid_matrix(self, matrix):
        self.assertIsInstance(matrix, Canvas)


    def test_can_generate_basic_matrix(self):
        matrix = self.chain.generate_residue_distance_matrix()
        self.check_valid_matrix(matrix)
