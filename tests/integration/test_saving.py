import os
import shutil
import atomium
from unittest import TestCase

class DictSavingTests(TestCase):

    def setUp(self):
        try:
            os.mkdir("tests/integration/files/output")
        except: pass


    def tearDown(self):
        '''if os.path.exists("tests/integration/files/output/"):
            shutil.rmtree("tests/integration/files/output/")'''
    

    def save(self, code):
        original = atomium.open(f"tests/integration/files/{code}.cif", dictionary=True)
        atomium.save_mmcif_dict(original, f"tests/integration/files/output/{code}.cif")
        saved = atomium.open(f"tests/integration/files/output/{code}.cif", dictionary=True)
        self.assertEqual(original, saved)


    def test_1lol_saving(self):
        self.save("1lol")