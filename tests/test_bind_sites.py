from unittest import TestCase
import unittest.mock
from molecupy.structures import BindSite, ResiduicStructure, ResiduicSequence
from molecupy.structures import Residue, SmallMolecule, Chain

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


    def test_can_assign_ligand(self):
        site = BindSite("A1", *self.residues)
        ligand = unittest.mock.Mock(spec=SmallMolecule)
        ligand._bind_site = None
        self.assertEqual(site.ligand(), None)
        site.ligand(ligand)
        self.assertEqual(site.ligand(), ligand)
        self.assertEqual(ligand._bind_site, site)


    def test_ligand_assigning_must_small_molecule(self):
        site = BindSite("A1", *self.residues)
        with self.assertRaises(TypeError):
            site.ligand("ligand")



class ContinuousSequenceTests(BindSiteTest):

    def test_can_make_continuous_sequence(self):
        chain = Chain("A", *self.residues)
        for residue in self.residues:
            residue.chain.return_value = residue._chain
        site = BindSite("A1", self.residues[2], self.residues[8], self.residues[4])
        bind_sequence = site.continuous_sequence()
        self.assertIsInstance(bind_sequence, ResiduicSequence)
        self.assertEqual(
         bind_sequence.residues(),
         self.residues[2:8]
        )
