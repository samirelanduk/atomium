from unittest import TestCase
from unittest.mock import patch, Mock
from atomium.converters.string2xyz import string_to_xyz
from atomium.xyz.xyz import Xyz

class StringToXyzTests(TestCase):

    def test_empty_string_is_blank_xyz(self):
        xyz = string_to_xyz("")
        self.assertIsInstance(xyz, Xyz)
        self.assertEqual(xyz._comment, "")
        self.assertEqual(xyz._model._atoms, set())
