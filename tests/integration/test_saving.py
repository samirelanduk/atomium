import os
import shutil
import atomium
from unittest import TestCase

class SaveTest(TestCase):

    def setUp(self):
        os.mkdir("tests/integration/files/output")


    def tearDown(self):
        if os.path.exists("tests/integration/files/output/"):
            shutil.rmtree("tests/integration/files/output/")



class MmcifDictSavingTests(SaveTest):

    def save(self, code):
        original = atomium.open(f"tests/integration/files/{code}.cif", dictionary=True)
        atomium.save_mmcif_dict(original, f"tests/integration/files/output/{code}.cif")
        saved = atomium.open(f"tests/integration/files/output/{code}.cif", dictionary=True)
        self.assertEqual(original, saved)


    def test_1lol_saving(self):
        self.save("1lol")
    

    def test_1xda_saving(self):
        self.save("1xda")
    

    def test_1m4x_saving(self):
        self.save("1m4x")
    

    def test_12ca_saving(self):
        self.save("12ca")
    

    def test_1cbn_saving(self):
        self.save("1cbn")
    

    def test_2bfb_saving(self):
        self.save("2bfb")
    

    def test_2igd_saving(self):
        self.save("2igd")
    

    def test_3jbp_saving(self):
        self.save("3jbp")
    

    def test_6xlu_saving(self):
        self.save("6xlu")


    def test_6uaj_saving(self):
        self.save("6uaj")
    
    
    def test_4jtd_saving(self):
        self.save("4jtd")
    

    def test_3cvx_saving(self):
        self.save("3cvx")
    

    def test_2fur_saving(self):
        self.save("2fur")



class BcifDictSavingTests(SaveTest):

    def save(self, code):
        original = atomium.open(f"tests/integration/files/{code}.bcif", dictionary=True)
        atomium.save_mmcif_dict(original, f"tests/integration/files/output/{code}.bcif")
        saved = atomium.open(f"tests/integration/files/output/{code}.bcif", dictionary=True)
        self.assertEqual(original, saved)


    def test_1lol_saving(self):
        self.save("1lol")