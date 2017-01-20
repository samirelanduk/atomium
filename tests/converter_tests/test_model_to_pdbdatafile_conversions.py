import os
from unittest import TestCase
import unittest.mock
from molecupy.converters.model2pdbdatafile import pdb_data_file_from_model
from molecupy.structures.models import Model
from molecupy.structures.complexes import Complex
from molecupy.structures.chains import Chain
from molecupy.structures.molecules import SmallMolecule, Residue
from molecupy.structures.atoms import Atom
from molecupy.pdb.pdbdatafile import PdbDataFile

class Model2PdbDataFile(TestCase):

    def setUp(self):
        self.model = Model()


class BasicDataFileCreationTests(Model2PdbDataFile):

    def test_model_can_create_pdbdatafile(self):
        data_file = pdb_data_file_from_model(self.model)
        self.assertIsInstance(data_file, PdbDataFile)


    def test_can_only_convert_models(self):
        with self.assertRaises(TypeError):
            pdb_data_file_from_model("PDB file")


    def test_data_file_knows_source(self):
        data_file = pdb_data_file_from_model(self.model)
        self.assertIs(data_file.source(), self.model)


    '''def test_can_save_to_file(self):
        model = Model()
        model.add_small_molecule(
         SmallMolecule("A1", "AA", Atom(1.0, 1.0, 1.0, "A", 1, "A"))
        )
        model.save_as_pdb("temp.pdb")
        try:
            with open("temp.pdb") as f:
                self.assertEqual(
                 f.read(),
                 model.pdb_data_file().generate_pdb_file().convert_to_string()
                )
        finally:
            os.remove("temp.pdb")'''



'''class ComplexConversionTests(TestCase):

    def test_can_add_complexes_to_pdb_data_file(self):
        chain1 = unittest.mock.Mock(Chain)
        chain2 = unittest.mock.Mock(Chain)
        chain3 = unittest.mock.Mock(Chain)
        chain1.chain_id.return_value = "A"
        chain2.chain_id.return_value = "B"
        chain3.chain_id.return_value = "C"
        complex1 = unittest.mock.Mock(Complex)
        complex2 = unittest.mock.Mock(Complex)
        complex1.chains.return_value = set([chain1, chain2])
        complex2.chains.return_value = set([chain3])
        complex1.complex_id.return_value = "1"
        complex2.complex_id.return_value = "2"
        complex1.complex_name.return_value = "FIRST COMPLEX"
        complex2.complex_name.return_value = "SECOND COMPLEX"
        model = Model()
        model.add_complex(complex1)
        model.add_complex(complex2)
        data_file = model.pdb_data_file()
        self.assertEqual(
         data_file.compounds(),
         [{
          "MOL_ID": 1,
          "MOLECULE": "FIRST COMPLEX",
          "CHAIN": ["A", "B"]
         }, {
          "MOL_ID": 2,
          "MOLECULE": "SECOND COMPLEX",
          "CHAIN": ["C"]
         }]
        )



class AtomConversionTests(TestCase):

    def test_can_add_atoms_to_pdb_data_file(self):
        model = Model()
        model.add_chain(Chain(
         "A",
         Residue(
          "A1",
          "PRO",
          Atom(1.2, 2.3, 3.4, "G", 23, "GX")
         ), Residue(
          "A1A",
          "VAL",
          Atom(11.2, 11.3, 34.4, "Y", 38, "YT")
         )
        ))
        data_file = model.pdb_data_file()
        self.assertEqual(
         data_file.atoms(),
         [
          {
           "atom_id": 23,
           "atom_name": "GX",
           "alt_loc": None,
           "residue_name": "PRO",
           "chain_id": "A",
           "residue_id": 1,
           "insert_code": None,
           "x": 1.2,
           "y": 2.3,
           "z": 3.4,
           "occupancy": 1.0,
           "temperature_factor": 0.0,
           "element": "G",
           "charge": None,
           "model_id": 1
          }, {
           "atom_id": 38,
           "atom_name": "YT",
           "alt_loc": None,
           "residue_name": "VAL",
           "chain_id": "A",
           "residue_id": 1,
           "insert_code": "A",
           "x": 11.2,
           "y": 11.3,
           "z": 34.4,
           "occupancy": 1.0,
           "temperature_factor": 0.0,
           "element": "Y",
           "charge": None,
           "model_id": 1
          }
         ]
        )
        self.assertEqual(data_file.heteroatoms(), [])


    def test_can_add_heteroatoms_to_pdb_data_file(self):
        model = Model()
        model.add_small_molecule(SmallMolecule(
         "A1001A",
         "MOL",
         Atom(1.2, 2.3, 3.4, "G", 23, "GX"),
         Atom(11.2, 11.3, 34.4, "Y", 38, "YT")
        ))
        data_file = model.pdb_data_file()
        self.assertEqual(
         data_file.heteroatoms(),
         [
          {
           "atom_id": 23,
           "atom_name": "GX",
           "alt_loc": None,
           "residue_name": "MOL",
           "chain_id": "A",
           "residue_id": 1001,
           "insert_code": "A",
           "x": 1.2,
           "y": 2.3,
           "z": 3.4,
           "occupancy": 1.0,
           "temperature_factor": 0.0,
           "element": "G",
           "charge": None,
           "model_id": 1
          }, {
           "atom_id": 38,
           "atom_name": "YT",
           "alt_loc": None,
           "residue_name": "MOL",
           "chain_id": "A",
           "residue_id": 1001,
           "insert_code": "A",
           "x": 11.2,
           "y": 11.3,
           "z": 34.4,
           "occupancy": 1.0,
           "temperature_factor": 0.0,
           "element": "Y",
           "charge": None,
           "model_id": 1
          }
         ]
        )
        self.assertEqual(data_file.atoms(), [])


    def test_can_add_heteroatom_bonds_to_pdb_data_file(self):
        model = Model()
        atom1 = Atom(1.0, 1.0, 1.0, "A", 1, "1")
        atom2 = Atom(1.0, 1.0, 1.0, "A", 2, "1")
        atom3 = Atom(1.0, 1.0, 1.0, "A", 3, "1")
        atom4 = Atom(1.0, 1.0, 1.0, "A", 4, "1")
        atom5 = Atom(1.0, 1.0, 1.0, "A", 5, "1")
        atom1.bond_to(atom2)
        atom1.bond_to(atom2)
        atom1.bond_to(atom3)
        atom1.bond_to(atom4)
        atom5.bond_to(atom4)
        molecule = SmallMolecule("A1", "MOL", atom1, atom2, atom3, atom4, atom5)
        model.add_small_molecule(molecule)
        data_file = model.pdb_data_file()
        self.assertEqual(
         data_file.connections(),
         [{
          "atom_id": 1,
          "bonded_atoms": [2, 3, 4]
         }, {
          "atom_id": 2,
          "bonded_atoms": [1]
         }, {
          "atom_id": 3,
          "bonded_atoms": [1]
         }, {
          "atom_id": 4,
          "bonded_atoms": [1, 5]
         }, {
          "atom_id": 5,
          "bonded_atoms": [4]
         }]
        )'''
