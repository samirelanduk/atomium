from datetime import date
from unittest import TestCase
from unittest.mock import Mock, patch, PropertyMock, MagicMock
from atomium.data import *

class DataDictToFileTests(TestCase):

    @patch("atomium.data.File")
    def test_can_convert_data_dict_to_file(self, mock_file):
        d = {"K1": {"A": 1, "B": 2}, "K2": {"C": 3}, "models": []}
        f = data_dict_to_file(d, "abc")
        mock_file.assert_called_with("abc")
        self.assertEqual(f, mock_file.return_value)
        self.assertEqual(f._A, 1)
        self.assertEqual(f._B, 2)
        self.assertEqual(f._C, 3)
