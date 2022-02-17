from datetime import date
from unittest import TestCase
import atomium

class DictParsingTests(TestCase):

    def test_1lol_bcif(self):
        d = atomium.open("tests/integration/files/1lol.bcif", dictionary=True)
        self.assertEqual(len(d.keys()), 66)
        self.assertEqual(d["entry"], [{"id": "1LOL"}])