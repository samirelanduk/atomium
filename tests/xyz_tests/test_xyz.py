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



class XyzReprTests(TestCase):

    def test_xyz_repr(self):
        xyz = Xyz("Glucose molecule")
        self.assertEqual(str(xyz), "<Xyz (Glucose molecule)>")



class XyzCommentTests(TestCase):

    def test_comment_property(self):
        xyz = Xyz("Glucose molecule")
        self.assertIs(xyz._comment, xyz.comment())


    def test_can_change_comment(self):
        xyz = Xyz("Glucose molecule")
        xyz.comment("Fructose molecule")
        self.assertEqual(xyz._comment, "Fructose molecule")


    def test_xyz_comment_must_be_str(self):
        xyz = Xyz("Glucose molecule")
        with self.assertRaises(TypeError):
            xyz.comment(100)
