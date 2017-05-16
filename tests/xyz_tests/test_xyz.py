from unittest import TestCase
from atomium.xyz.xyz import Xyz

class XyzCreationTests(TestCase):

    def test_can_create_xyz(self):
        xyz = Xyz()
        self.assertEqual(xyz._model, None)
        self.assertEqual(xyz._comment, "")


    def test_can_create_xyz_with_comment(self):
        xyz = Xyz("Glucose molecule")
        self.assertEqual(xyz._model, None)
        self.assertEqual(xyz._comment, "Glucose molecule")


    def test_xyz_comment_must_be_str(self):
        with self.assertRaises(TypeError):
            Xyz(100)
