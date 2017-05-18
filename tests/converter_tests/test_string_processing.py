from unittest import TestCase
from atomium.converters.strings import string2lines

class String2LinesTests(TestCase):

    def test_line_becomes_line(self):
        self.assertEqual(string2lines("A line."), ["A line."])
