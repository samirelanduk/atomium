from unittest import TestCase
from molecupy import exceptions
from molecupy.structures import ResiduicStructure, PdbSite, PdbAtom, PdbResidue

class SiteTest(TestCase):

    def setUp(self):
        self.atom1 = PdbAtom(1.0, 1.0, 1.0, "H", 1, "H1")
        self.atom2 = PdbAtom(1.0, 1.0, 2.0, "C", 2, "CA")
        self.atom3 = PdbAtom(1.0, 1.0, 3.0, "O", 3, "OX1")
        self.residue1 = PdbResidue(1, "ARG", self.atom1, self.atom2, self.atom3)
        self.atom4 = PdbAtom(1.0, 1.0, 4.0, "H", 4, "H1")
        self.atom5 = PdbAtom(1.0, 1.0, 5.0, "C", 5, "CA")
        self.atom6 = PdbAtom(1.0, 1.0, 6.0, "O", 6, "OX1")
        self.residue2 = PdbResidue(2, "HIS", self.atom4, self.atom5, self.atom6)
        self.atom7 = PdbAtom(1.0, 1.0, 7.0, "H", 7, "H1")
        self.atom8 = PdbAtom(1.0, 1.0, 8.0, "C", 8, "CA")
        self.atom9 = PdbAtom(1.0, 1.0, 9.0, "O", 9, "OX1")
        self.residue3 = PdbResidue(3, "TRP", self.atom7, self.atom8, self.atom9)


    def check_valid_site(self, site):
        self.assertIsInstance(site, PdbSite)
        self.assertIsInstance(site, ResiduicStructure)
        self.assertIsInstance(site.site_id, str)
        site.ligand
        self.assertRegex(str(site), r"<Site (.+) \((\d+) residues\)>")



class SiteTests(SiteTest):

    def test_can_create_site(self):
        site = PdbSite("AB1", self.residue1, self.residue2, self.residue3)
        self.check_valid_site(site)


    def test_site_id_must_be_str(self):
        with self.assertRaises(TypeError):
            site = PdbSite(1, self.residue1, self.residue2, self.residue3)
