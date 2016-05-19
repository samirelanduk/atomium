from unittest import TestCase
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
