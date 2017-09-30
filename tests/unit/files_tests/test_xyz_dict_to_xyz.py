from unittest import TestCase
from unittest.mock import patch, Mock, MagicMock
from atomium.files.xyzdict2xyz import *

class XyzDictToXyzTests(TestCase):

    @patch("atomium.files.xyzdict2xyz.Xyz")
    @patch("atomium.files.xyzdict2xyz.xyz_dict_to_model")
    def test_can_convert_pdb_dict_to_pdb(self, mock_model, mock_xyz):
        xyz_dict = {"title": "The xyz title"}
        mock_model.return_value = "MODEL"
        xyz = Mock()
        mock_xyz.return_value = xyz
        returned_xyz = xyz_dict_to_xyz(xyz_dict)
        mock_model.assert_called_with(xyz_dict)
        self.assertIs(returned_xyz, xyz)
        self.assertEqual(xyz._model, "MODEL")
        self.assertEqual(xyz._title, "The xyz title")



class XyzDictToXyzTests(TestCase):

    @patch("atomium.files.xyzdict2xyz.Model")
    @patch("atomium.files.xyzdict2xyz.atom_dict_to_atom")
    def test_can_convert_xyz_dict_to_model(self, mock_atom, mock_model):
        model = Mock()
        mock_model.return_value = model
        atoms = [Mock(), Mock(), Mock()]
        mock_atom.side_effect = atoms
        xyz_dict = {"atoms": ["atom1", "atom2", "atom3"]}
        returned_model = xyz_dict_to_model(xyz_dict)
        mock_atom.assert_any_call("atom1")
        mock_atom.assert_any_call("atom2")
        mock_atom.assert_any_call("atom3")
        mock_model.assert_called_with(*atoms)
        self.assertIs(returned_model, model)



class AtomDictToAtomTests(TestCase):

    @patch("atomium.files.xyzdict2xyz.Atom")
    def test_can_convert_atom_dict_to_atom(self, mock_atom):
        mock_atom.return_value = "ATOM"
        atom_dict = {"element": "C", "x": 1.2, "y": -4.5, "z": 3.3}
        atom = atom_dict_to_atom(atom_dict)
        mock_atom.assert_called_with("C", 1.2, -4.5, 3.3)
        self.assertEqual(atom, "ATOM")
