from unittest import TestCase
from unittest.mock import Mock, patch
from atomium.structures.chains import Site
from atomium.structures.molecules import Residue, AtomicStructure, Molecule

class SiteTest(TestCase):

    def setUp(self):
        self.residues = [Mock(Residue), Mock(Residue), Mock(Residue), Mock(Residue)]
        def mock_init(obj, *args, **kwargs):
            obj._atoms = {i: {arg} for i, arg in enumerate(args)}
        self.patch1 = patch("atomium.structures.molecules.AtomicStructure.__init__")
        self.mock_init = self.patch1.start()
        self.mock_init.side_effect = mock_init
        self.ligand = Mock(Molecule)


    def tearDown(self):
        self.patch1.stop()


class SiteCreationTests(SiteTest):

    def test_can_create_site(self):
        site = Site(*self.residues)
        self.assertIsInstance(site, AtomicStructure)
        self.mock_init.assert_called_with(site, *self.residues)
        self.assertIsNone(site._ligand)


    def test_can_create_site_with_ligand(self):
        site = Site(*self.residues, ligand=self.ligand)
        self.mock_init.assert_called_with(site, *self.residues)
        self.assertIs(site._ligand, self.ligand)


    def test_ligand_must_be_molecule(self):
        with self.assertRaises(TypeError):
            Site(*self.residues, ligand="ligand")



class SiteReprTests(SiteTest):

    @patch("atomium.structures.chains.Site.residues")
    def test_site_repr_no_ligand(self, mock_res):
        mock_res.return_value = [1, 2, 3]
        site = Site(*self.residues)
        self.assertEqual(str(site), "<Site (3 residues)>")


    @patch("atomium.structures.chains.Site.residues")
    def test_site_repr_with_ligand(self, mock_res):
        mock_res.return_value = [1, 2, 3]
        self.ligand._id = "LIGID"
        site = Site(*self.residues, ligand=self.ligand)
        self.assertEqual(str(site), "<'LIGID' Site (3 residues)>")



class SiteLigandTests(SiteTest):

    def test_can_get_site_ligand(self):
        site = Site(*self.residues, ligand=self.ligand)
        self.assertIs(site._ligand, site.ligand)


    def test_can_update_ligand(self):
        site = Site(*self.residues, ligand=self.ligand)
        ligand = Mock(Molecule)
        site.ligand = ligand
        self.assertIs(site._ligand, ligand)


    def test_ligand_must_be_molecule(self):
        site = Site(*self.residues, ligand=self.ligand)
        with self.assertRaises(TypeError):
            site.ligand = "ligand"



class SiteResiduesTests(SiteTest):

    @patch("atomium.structures.chains.Site.residues")
    def test_can_get_normal_residues(self, mock_res):
        mock_res.return_value = {1, 2, 3}
        site = Site(*self.residues)
        self.assertEqual(site.residues(), {1, 2, 3})


    @patch("atomium.structures.chains.AtomicStructure.residues")
    @patch("atomium.structures.chains.Site.atoms")
    def test_can_get_water_residues(self, mock_atoms, mock_res):
        mock_res.return_value = {1, 2, 3}
        site = Site(*self.residues)
        atoms = [Mock(), Mock(), Mock()]
        site.atoms.return_value = atoms
        water = Mock()
        water.residue_name.return_value = "HOH"
        atoms[0].molecule = water
        atoms[1].molecule = 1
        atoms[2].molecule = 2
        self.assertEqual(site.residues(), {1, 2, 3, water})
