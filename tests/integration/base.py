import os
from unittest import TestCase

class IntegratedTest(TestCase):

    def setUp(self):
        self.files_at_start = os.listdir("tests/integration/files")


    def tearDown(self):
        files_at_end = os.listdir("tests/integration/files")
        to_remove = [f for f in files_at_end if f not in self.files_at_start]
        for f in to_remove:
            os.remove("tests/integration/files/%s" % f)
