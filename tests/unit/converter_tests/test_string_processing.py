from unittest import TestCase
from unittest.mock import Mock, MagicMock, patch
from atomium.converters.strings import string2lines, string_from_file, string_to_file

class String2LinesTests(TestCase):

    def test_line_becomes_line(self):
        self.assertEqual(string2lines("A line."), ["A line."])


    def test_can_break_string_on_newline(self):
        self.assertEqual(string2lines("A line.\nB line."), ["A line.", "B line."])


    def test_can_handle_windows_linebreaks(self):
        self.assertEqual(string2lines("A line.\r\nB line."), ["A line.", "B line."])



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
