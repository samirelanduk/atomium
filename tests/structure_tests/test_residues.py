from unittest import TestCase
from unittest.mock import Mock, patch
from atomium.structures.molecules import Molecule, Residue
from atomium.structures.atoms import Atom

class ResidueTest(TestCase):

    def setUp(self):
        self.atom1, self.atom2, self.atom3 = Mock(Atom), Mock(Atom), Mock(Atom)
        self.atoms = [self.atom1, self.atom2, self.atom3]



class ResidueCreationTests(ResidueTest):

    @patch("atomium.structures.molecules.Molecule.__init__")
    def test_can_create_residue(self, mock_init):
        mock_init.return_value = None
        res = Residue(self.atom1, self.atom2, self.atom3)
        self.assertIsInstance(res, Molecule)
        mock_init.assert_called_with(res, self.atom1, self.atom2, self.atom3)


    @patch("atomium.structures.molecules.Molecule.__init__")
    def test_can_create_residue_with_name(self, mock_init):
        mock_init.return_value = None
        res = Residue(self.atom1, self.atom2, self.atom3, name="GLY")
        mock_init.assert_called_with(res, self.atom1, self.atom2, self.atom3, name="GLY")


    @patch("atomium.structures.molecules.Molecule.__init__")
    def test_can_residue_init_doesnt_send_residue_id(self, mock_init):
        mock_init.return_value = None
        res = Residue(self.atom1, self.atom2, self.atom3, residue_id="A100")
        mock_init.assert_called_with(res, self.atom1, self.atom2, self.atom3)


    def test_can_create_residue_with_id(self):
        res = Residue(self.atom1, self.atom2, self.atom3, residue_id="A101")
        self.assertEqual(res._id, "A101")


    def test_residue_id_must_be_str(self):
        with self.assertRaises(TypeError):
            Residue(self.atom1, self.atom2, self.atom3, residue_id=100)


    def test_residue_id_must_be_unique(self):
        res = Residue(self.atom1, self.atom2, self.atom3, residue_id="A200")
        with self.assertRaises(ValueError):
            Residue(self.atom1, self.atom2, self.atom3, residue_id="A200")
