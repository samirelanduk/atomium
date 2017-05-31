import os
from unittest import TestCase

class IntegratedTest(TestCase):

    def setUp(self):
        self.files_at_start = os.listdir("itests/files")


    def tearDown(self):
        files_at_end = os.listdir("itests/files")
        to_remove = [f for f in files_at_end if f not in self.files_at_start]
        for f in to_remove:
            os.remove("itests/files/%s" % f)
