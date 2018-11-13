import atomium
import os
from unittest import TestCase

class SavingTest(TestCase):

    def setUp(self):
        self.files_at_start = os.listdir("tests/integration/files")


    def tearDown(self):
        files_at_end = os.listdir("tests/integration/files")
        to_remove = [f for f in files_at_end if f not in self.files_at_start]
        for f in to_remove:
            os.remove("tests/integration/files/" + f)


    def check_file_saving(self, filename):
        f = atomium.open("tests/integration/files/" + filename)
        f.model.save("tests/integration/files/saved_" + filename)
        f2 = atomium.open("tests/integration/files/saved_" + filename)
        self.assertTrue(f.model.equivalent_to(f2.model))



class MmcifFileSavingTests(SavingTest):

    def test_can_save_1lol(self):
        self.check_file_saving("1lol.cif")


    def test_can_save_1cbn(self):
        self.check_file_saving("1cbn.cif")


    def test_can_save_1m4x(self):
        self.check_file_saving("1m4x.cif")


    def test_can_save_1xda(self):
        self.check_file_saving("1xda.cif")


    def test_can_save_5xme(self):
        self.check_file_saving("5xme.cif")


    def test_can_save_4y60(self):
        self.check_file_saving("4y60.cif")



class MmtfFileSavingTests(SavingTest):

    def test_can_save_1lol(self):
        self.check_file_saving("1lol.mmtf")


    def test_can_save_1cbn(self):
        self.check_file_saving("1cbn.mmtf")


    def test_can_save_1m4x(self):
        self.check_file_saving("1m4x.mmtf")


    def test_can_save_1xda(self):
        self.check_file_saving("1xda.mmtf")


    def test_can_save_5xme(self):
        self.check_file_saving("5xme.mmtf")


    def test_can_save_4y60(self):
        self.check_file_saving("4y60.mmtf")
