from unittest import TestCase
from unittest.mock import Mock, MagicMock, patch
from atomium.converters.strings import string2lines, string_from_file
from atomium.converters.strings import string_to_file, fetch_string

class String2LinesTests(TestCase):

    def test_line_becomes_line(self):
        self.assertEqual(string2lines("A line."), ["A line."])


    def test_can_break_string_on_newline(self):
        self.assertEqual(string2lines("A line.\nB line."), ["A line.", "B line."])


    def test_can_handle_windows_linebreaks(self):
        self.assertEqual(string2lines("A line.\r\nB line."), ["A line.", "B line."])


    def test_trailing_empty_lines_ignored(self):
        self.assertEqual(
         string2lines("A line.\nB line.\n\n\n\n"),
         ["A line.", "B line."]
        )



class StringFromFileTests(TestCase):

    @patch("builtins.open")
    def test_gets_string_from_file(self, mock_open):
        open_return = MagicMock()
        mock_file = Mock()
        open_return.__enter__.return_value = mock_file
        mock_file.read.return_value = "returnstring"
        mock_open.return_value = open_return
        string = string_from_file("path/to/file")
        mock_open.assert_called_with("path/to/file")
        self.assertEqual(string, "returnstring")



class StringToFileTests(TestCase):

    @patch("builtins.open")
    def test_saves_string_to_file(self, mock_open):
        open_return = MagicMock()
        mock_file = Mock()
        mock_write = MagicMock()
        mock_file.write = mock_write
        open_return.__enter__.return_value = mock_file
        mock_open.return_value = open_return
        string_to_file("filestring", "filename")
        mock_open.assert_called_once_with("filename", "w")
        mock_write.assert_called_once_with("filestring")



class StringFromWebServicesTests(TestCase):

    @patch("requests.get")
    def test_can_fetch_string(self, mock_get):
        response = Mock()
        response.status_code = 200
        response.text = "filestring"
        mock_get.return_value = response
        returned_string = fetch_string("1XXX")
        mock_get.assert_called_with("https://files.rcsb.org/view/1xxx.pdb")
        self.assertEqual(returned_string, "filestring")


    @patch("requests.get")
    def test_can_fetch_pdbe_string(self, mock_get):
        response = Mock()
        response.status_code = 200
        response.text = "filestring"
        mock_get.return_value = response
        returned_string = fetch_string("1XXX", pdbe=True)
        mock_get.assert_called_with(
         "https://www.ebi.ac.uk/pdbe/entry-files/pdb1xxx.ent"
        )
        self.assertEqual(returned_string, "filestring")


    @patch("requests.get")
    def test_can_get_none_if_no_file_found(self, mock_get):
        response = Mock()
        response.status_code = 404
        self.assertIsNone(fetch_string("1XXX"))


    def test_pdb_code_must_be_string(self):
        with self.assertRaises(TypeError):
            fetch_string(4000)


    def test_pdb_code_must_be_len_4(self):
        with self.assertRaises(ValueError):
            fetch_string("1xx")
