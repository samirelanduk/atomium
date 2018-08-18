from copy import deepcopy
from unittest import TestCase
from unittest.mock import Mock, patch, PropertyMock, MagicMock
from atomium.files.xyz import *

class XyzStringToXyzDictTests(TestCase):

    def setUp(self):
        self.lines = [
         "",
         "C       -0.180226841      0.360945118     -1.120304970",
         "CU      -0.180226841      1.559292118     -0.407860970",
         "C       -0.180226841      1.503191118      0.986935030"
        ]

    def test_can_turn_xyz_string_to_xyz_dict_just_atoms(self):
        self.assertEqual(xyz_string_to_xyz_dict("\n".join(self.lines)), {
         "header_lines": [], "atom_lines": self.lines[1:]
        })


    def test_can_turn_xyz_string_to_xyz_dict_with_count(self):
        self.lines.insert(0, "11")
        self.assertEqual(xyz_string_to_xyz_dict("\n".join(self.lines)), {
         "header_lines": ["11"], "atom_lines": self.lines[2:]
        })


    def test_can_turn_xyz_string_to_xyz_dict_with_header(self):
        self.lines.insert(0, "11")
        self.lines.insert(1, "name of molecule")
        self.assertEqual(xyz_string_to_xyz_dict("\n".join(self.lines)), {
         "header_lines": ["11", "name of molecule"], "atom_lines": self.lines[3:]
        })



class XyzDictToDataDictTests(TestCase):

    def test_can_convert_xyz_dict_to_data_dict(self):
        d = deepcopy(DATA_DICT)
        d["description"]["title"] = "name"
        d["models"].append(deepcopy(MODEL_DICT))
        d["models"][0]["atoms"] = [{
         "id": 0, "element": "R", "name": None,
         "x": 1.1, "y": 2.2, "z": -4.5,
         "bfactor": None, "charge": 0,
         "residue_id": None, "residue_name": None, "residue_insert": "",
         "chain_id": None, "occupancy": 1, "alt_loc": None,
         "anisotropy": [], "polymer": False, "full_res_id": None
        }, {
         "id": 0, "element": "CA", "name": None,
         "x": 0.435, "y": 2.0, "z": -19.234,
         "bfactor": None, "charge": 0,
         "residue_id": None, "residue_name": None, "residue_insert": "",
         "chain_id": None, "occupancy": 1, "alt_loc": None,
         "anisotropy": [], "polymer": False, "full_res_id": None
        }]
        self.assertEqual(xyz_dict_to_data_dict({
         "header_lines": ["11", "name"], "atom_lines": [
          "R      1.1        2.2       -4.5",
          "CA     0.435      2.0      -19.234"
         ]
        }), d)
