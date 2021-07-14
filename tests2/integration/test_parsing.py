from unittest import TestCase
import atomium

class Test(TestCase):

    def test(self):
        f = atomium.open("tests/integration/files/1lol.cif")
        print(f)