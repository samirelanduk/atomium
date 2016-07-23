from unittest import TestCase
import unittest.mock
from molecupy.structures import BindSite, ResiduicStructure, Residue

class BindSiteTest(TestCase):

    def setUp(self):
        self.residues = [unittest.mock.Mock(spec=Residue) for _ in range(10)]



class SiteCreationTests(BindSiteTest):

    def test_can_create_bind_site(self):
        site = BindSite("A1", *self.residues)
        self.assertIsInstance(site, ResiduicStructure)
        self.assertEqual(site._site_id, "A1")
        self.assertEqual(site._ligand, None)
        self.assertEqual(site._residues, set(self.residues))


    def test_site_id_must_be_str(self):
        with self.assertRaises(TypeError):
            BindSite(200, *self.residues)


    def test_site_repr(self):
        site = BindSite("A1", *self.residues)
        self.assertEqual(str(site), "<BindSite A1 (10 residues)>")



class SitePropertyTests(BindSiteTest):

    def test_site_properties(self):
        site = BindSite("A1", *self.residues)
        self.assertEqual(site.site_id(), "A1")
        self.assertEqual(site.ligand(), None)
