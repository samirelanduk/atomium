from unittest import TestCase
from unittest.mock import patch, Mock
from atomium.converters.pdbdatafile2model import *
from atomium.parse.pdbdatafile import PdbDataFile
from atomium.structures.models import Model
from atomium.structures.atoms import Atom

class PdbDataFileToModelTests(TestCase):

    def test_can_get_model_from_data_file(self):
        data_file = Mock(PdbDataFile)
        model = pdb_data_file_to_model(data_file)
        self.assertIsInstance(model, Model)



class AtomDictToAtomTest(TestCase):

    def setUp(self):
        self.atom_dict = {
         "atom_id": 107, "atom_name": "N1", "alt_loc": None,
         "residue_name": "GLY",
         "chain_id": "A", "residue_id": 13, "insert_code": "A",
         "x": 12.681, "y": 37.302, "z": -25.211,
         "occupancy": 1.0, "temperature_factor": 15.56,
         "element": "N", "charge": -2
        }


    def test_can_make_atom_from_dict(self):
        atom = atom_dict_to_atom(self.atom_dict)
        self.assertIsInstance(atom, Atom)
        self.assertEqual(atom._x, 12.681)
        self.assertEqual(atom._y, 37.302)
        self.assertEqual(atom._z, -25.211)
        self.assertEqual(atom._id, 107)
        self.assertEqual(atom._name, "N1")
        self.assertEqual(atom._element, "N")
        self.assertEqual(atom._charge, -2)
        self.assertEqual(atom.temp_chain, "A")
        self.assertEqual(atom.temp_residue_id, "A13A")
        self.assertEqual(atom.temp_residue_name, "GLY")


    def test_can_make_atom_from_minimal_dict(self):
        self.atom_dict["insert_code"] = None
        self.atom_dict["occupancy"] = None
        self.atom_dict["temperature_factor"] = None
        self.atom_dict["charge"] = None
        self.atom_dict["element"] = None
        atom = atom_dict_to_atom(self.atom_dict)
        self.assertIsInstance(atom, Atom)
        self.assertEqual(atom._element, "X")
        self.assertEqual(atom._charge, 0)
        self.assertEqual(atom.temp_residue_id, "A13")
