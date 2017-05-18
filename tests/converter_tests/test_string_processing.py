from unittest import TestCase
from atomium.converters.strings import string2lines

class String2LinesTests(TestCase):

    def test_line_becomes_line(self):
        self.assertEqual(string2lines("A line."), ["A line."])


    def test_can_break_string_on_newline(self):
        self.assertEqual(string2lines("A line.\nB line."), ["A line.", "B line."])


    def test_can_handle_windows_linebreaks(self):
        self.assertEqual(string2lines("A line.\r\nB line."), ["A line.", "B line."])
