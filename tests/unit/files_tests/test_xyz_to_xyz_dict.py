from unittest import TestCase
from unittest.mock import patch, Mock, MagicMock
from atomium.files.xyz2xyzdict import *

class XyzToXyzDictTests(TestCase):

    @patch("atomium.files.xyz2xyzdict.structure_to_xyz_dict")
    def test_can_convert_xyz_to_xyz_dict(self, mock_dict):
        xyz = Mock()
        xyz._title = "T"
        xyz._model = "model1"
        mock_dict.return_value = {"atoms": ["a1", "a2"]}
        xyz_dict = xyz_to_xyz_dict(xyz)
        mock_dict.assert_called_with("model1")
        self.assertEqual(xyz_dict, {"atoms": ["a1", "a2"], "title": "T"})



class StructureToXyzDictTests(TestCase):

    @patch("atomium.files.xyz2xyzdict.atom_to_atom_dict")
    def test_can_convert_structure_to_xyz_dict(self, mock_atom):
        structure = Mock()
        atoms = [Mock(), Mock(), Mock()]
        mock_atom.side_effect = ["a1", "a2", "a3"]
        structure.atoms.return_value = atoms
        xyz_dict = structure_to_xyz_dict(structure)
        for atom in atoms:
            mock_atom.assert_any_call(atom)
        self.assertEqual(xyz_dict, {"title": None, "atoms": ["a1", "a2", "a3"]})



class AtomToAtomDictTests(TestCase):

    def test_can_convert_atom_to_atom_dict(self):
        atom = Mock()
        atom.element = "E"
        atom.x = "X"
        atom.y = "Y"
        atom.z = "Z"
        self.assertEqual(
         atom_to_atom_dict(atom), {"element": "E", "x": "X", "y": "Y", "z": "Z"}
        )
