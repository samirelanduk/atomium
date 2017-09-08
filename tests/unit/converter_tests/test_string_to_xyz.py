from unittest import TestCase
from unittest.mock import patch, Mock
from atomium.converters.string2xyz import *
from atomium.files.xyz import Xyz
from atomium.structures.atoms import Atom
from atomium.structures.models import Model

class StringToXyzTests(TestCase):

    def setUp(self):
        self.patcher1 = patch("atomium.converters.string2xyz.string2lines")
        self.patcher2 = patch("atomium.converters.string2xyz.remove_atom_num")
        self.patcher3 = patch("atomium.converters.string2xyz.parse_atom")
        self.patcher4 = patch("atomium.converters.string2xyz.extract_comment")
        self.mock_lines = self.patcher1.start()
        self.mock_remove = self.patcher2.start()
        self.mock_parse = self.patcher3.start()
        self.mock_extract = self.patcher4.start()
        self.mock_lines.return_value = ["12", "comment", "...", "line1", "line2"]
        self.mock_extract.return_value = "comment"
        self.atoms = [Mock(Atom) for _ in range(5)]
        for atom in self.atoms:
            atom.atom_id.return_value = None
        self.mock_parse.side_effect = self.atoms


    def tearDown(self):
        self.patcher1.stop()
        self.patcher2.stop()
        self.patcher3.stop()
        self.patcher4.stop()


    def test_converter_produces_xyz(self):
        xyz = string_to_xyz("teststring")
        self.assertIsInstance(xyz, Xyz)


    def test_converter_creates_model(self):
        xyz = string_to_xyz("teststring")
        self.assertIsInstance(xyz._model, Model)


    def test_xyz_has_correct_comment(self):
        xyz = string_to_xyz("teststring")
        self.assertEqual(xyz._comment, self.mock_extract.return_value)


    def test_model_has_correct_atoms(self):
        xyz = string_to_xyz("teststring")
        self.assertEqual(xyz._model._atoms, set(self.atoms))


    def test_converter_function_calls(self):
        xyz = string_to_xyz("teststring")
        self.mock_lines.assert_called_with("teststring")
        self.mock_remove.assert_called_with(self.mock_lines.return_value)
        self.mock_extract.assert_called_with(self.mock_lines.return_value)
        for line in self.mock_lines.return_value:
            self.mock_parse.assert_any_call(line)



class AtomNumRemovalTests(TestCase):

    def test_can_remove_atom_num(self):
        lines = ["34", "line1", "line2"]
        remove_atom_num(lines)
        self.assertEqual(lines, ["line1", "line2"])


    def test_can_ignore_non_atom_num(self):
        lines = ["34.7", "line1", "line2"]
        remove_atom_num(lines)
        self.assertEqual(lines, ["34.7", "line1", "line2"])


    def test_can_handle_empty_lines(self):
        lines = []
        remove_atom_num(lines)
        self.assertEqual(lines, [])



class AtomParserTests(TestCase):

    def test_can_parse_atom(self):
        atom = parse_atom("C  35.884  30.895  49.120")
        self.assertIsInstance(atom, Atom)
        self.assertEqual(atom._element, "C")
        self.assertEqual(atom._x, 35.884)
        self.assertEqual(atom._y, 30.895)
        self.assertEqual(atom._z, 49.12)


    def test_unparseable_line_returns_none(self):
        atom = parse_atom("A comment line")
        self.assertIs(atom, None)



class CommentExtractorTests(TestCase):

    @patch("atomium.converters.string2xyz.parse_atom")
    def test_can_extract_comment(self, mock_parse):
        mock_parse.return_value = None
        lines = ["comment", "atom1", "atom2"]
        self.assertEqual(extract_comment(lines), "comment")
        self.assertEqual(lines, ["atom1", "atom2"])


    @patch("atomium.converters.string2xyz.parse_atom")
    def test_stops_at_parseable_line(self, mock_parse):
        mock_parse.return_value = Mock(Atom)
        lines = ["atom1", "atom2"]
        self.assertEqual(extract_comment(lines), "")
        self.assertEqual(lines, ["atom1", "atom2"])


    @patch("atomium.converters.string2xyz.parse_atom")
    def test_can_handle_empty_lines(self, mock_parse):
        mock_parse.return_value = None
        lines = []
        self.assertEqual(extract_comment(lines), "")
        self.assertEqual(lines, [])
