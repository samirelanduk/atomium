from unittest import TestCase
from molecupy import exceptions
from molecupy.structures import PdbAlphaHelix, ResiduicSequence, PdbAtom, PdbResidue, PdbChain

class HelixTest(TestCase):

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
        self.chain = PdbChain("A", self.residue1, self.residue2, self.residue3)


    def check_valid_helix(self, helix):
        self.assertIsInstance(helix, PdbAlphaHelix)
        self.assertIsInstance(helix, ResiduicSequence)
        self.assertIsInstance(helix.helix_id, str)
        self.assertIsInstance(helix.residues, list)
        if helix.helix_class is not None:
            self.assertIsInstance(helix.helix_class, str)
        if helix.comment is not None:
            self.assertIsInstance(helix.comment, str)
        self.assertIsInstance(helix.get_chain(), PdbChain)
        self.assertIn(helix, helix.get_chain().alpha_helices)
        self.assertRegex(str(helix), r"<AlphaHelix (.+) \((\d+) residues\)>")



class HelixCreationTests(HelixTest):

    def test_can_create_helix(self):
        helix = PdbAlphaHelix("AA", self.residue1, self.residue2)
        self.check_valid_helix(helix)
        self.assertIs(helix.get_chain(), self.chain)


    def test_helix_id_must_be_str(self):
        with self.assertRaises(TypeError):
            helix = PdbAlphaHelix(1, self.residue1, self.residue2)
        with self.assertRaises(TypeError):
            helix = PdbAlphaHelix(None, self.residue1, self.residue2)


    def test_all_helix_residues_must_be_on_same_chain(self):
        new_residue = PdbResidue("B1", "VAL", PdbAtom(1.0, 1.0, 1.0, "H", 10, "H1"))
        new_chain = PdbChain("B", new_residue)
        with self.assertRaises(exceptions.BrokenHelixError):
            helix = PdbAlphaHelix("AH", self.residue1, new_residue)


    def test_can_make_helix_with_classification(self):
        helix = PdbAlphaHelix("AA", self.residue1, self.residue2, helix_class=".")
        self.assertEqual(helix.helix_class, ".")


    def test_helix_class_must_be_str(self):
        with self.assertRaises(TypeError):
            helix = PdbAlphaHelix("AA", self.residue1, self.residue2, helix_class=1)


    def test_can_make_helix_with_comment(self):
        helix = PdbAlphaHelix("AA", self.residue1, self.residue2, comment=".")
        self.assertEqual(helix.comment, ".")


    def test_helix_comment_must_be_str(self):
        with self.assertRaises(TypeError):
            helix = PdbAlphaHelix("AA", self.residue1, self.residue2, comment=1)
