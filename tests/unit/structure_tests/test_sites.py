from unittest import TestCase
from unittest.mock import Mock, patch
from atomium.structures.chains import ResidueStructure, Site
from atomium.structures.molecules import Residue, AtomicStructure, Molecule

class SiteTest(TestCase):

    def setUp(self):
        self.residues = [Mock(Residue), Mock(Residue), Mock(Residue), Mock(Residue)]
        def mock_init(obj, *args, **kwargs):
            obj._atoms = set()
            obj._id_atoms = {}
        self.mock_init = mock_init



class SiteCreationTests(SiteTest):

    @patch("atomium.structures.molecules.AtomicStructure.__init__")
    def test_can_create_site(self, mock_init):
        mock_init.side_effect = self.mock_init
        site = Site(*self.residues)
        self.assertIsInstance(site, AtomicStructure)
        self.assertIsInstance(site, ResidueStructure)
        mock_init.assert_called_with(site, *self.residues)
        self.assertIsNone(site._ligand)


    @patch("atomium.structures.molecules.AtomicStructure.__init__")
    def test_can_create_site_with_ligand(self, mock_init):
        mock_init.side_effect = self.mock_init
        ligand = Mock(Molecule)
        site = Site(*self.residues, ligand=ligand)
        mock_init.assert_called_with(site, *self.residues)
        self.assertIs(site._ligand, ligand)


    def test_ligand_must_be_molecule(self):
        with self.assertRaises(TypeError):
            Site(*self.residues, ligand="ligand")
