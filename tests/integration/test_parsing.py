from unittest import TestCase
import atomium

class DictParsingTests(TestCase):

    def test_1lol_mmcif(self):
        d = atomium.open("tests/integration/files/1lol.cif", dictionary=True)
        self.assertEqual(len(d.keys()), 66)
        print(len(d.keys()))