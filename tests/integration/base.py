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


    def check_files_the_same(self, new, ref):
        with open("tests/integration/files/" + new) as f:
            new = [l.strip() for l in f.readlines() if l.strip()]
        with open("tests/integration/files/" + ref) as f:
            ref = [l.strip() for l in f.readlines() if l.strip()]
        for new_line, ref_line in zip(new, ref):
            self.assertEqual(new_line, ref_line)
