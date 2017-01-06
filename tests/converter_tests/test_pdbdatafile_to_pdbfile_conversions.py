from molecupy.pdb.pdbfile import PdbFile
from molecupy.pdb.pdbdatafile import PdbDataFile

from unittest import TestCase

class PdbDataFileTest(TestCase):

    def setUp(self):
        self.empty = PdbDataFile(PdbFile(""))
        self.blank = PdbDataFile()



class CompoundsConversionTests(PdbDataFileTest):

    def test_can_produce_compnd_records(self):
        compounds = [
         {
          "MOL_ID": 1,
          "MOLECULE": "OROTIDINE 5'-MONOPHOSPHATE DECARBOXYLASE",
          "CHAIN": ["A", "B"],
          "SYNONYM": [
           "OMP DECARBOXYLASE",
           "OMPDCASE",
           "OMPDECASE"
          ],
          "EC": "4.1.1.23",
          "ENGINEERED": True
         }, {
          "MOL_ID": 2,
          "MOLECULE": "OROTIDINE 5'-MONOPHOSPHATE VERYPHOSPHATE INDEED DECARBOXYLASE PLUS"
         }
        ]
        for compound in compounds:
            self.blank.compounds().append(compound)
        pdb_file = self.blank.generate_pdb_file()
        self.assertEqual(len(pdb_file.records()), 9)
        self.assertEqual(
         "\n".join([line.rstrip() for line in pdb_file.convert_to_string().split("\n")]),
         "COMPND    MOL_ID: 1;\n"
         "COMPND   2 MOLECULE: OROTIDINE 5'-MONOPHOSPHATE DECARBOXYLASE;\n"
         "COMPND   3 CHAIN: A, B;\n"
         "COMPND   4 SYNONYM: OMP DECARBOXYLASE, OMPDCASE, OMPDECASE;\n"
         "COMPND   5 EC: 4.1.1.23;\n"
         "COMPND   6 ENGINEERED: YES;\n"
         "COMPND   7 MOL_ID: 2;\n"
         "COMPND   8 MOLECULE: OROTIDINE 5'-MONOPHOSPHATE VERYPHOSPHATE INDEED\n"
         "COMPND   9 DECARBOXYLASE PLUS;"
        )



class AtomsConversionTests(PdbDataFileTest):

    def test_can_produce_atom_records(self):
        atoms = [
         {
          "atom_id": 107,
          "atom_name": "N",
          "alt_loc": None,
          "residue_name": "GLY",
          "chain_id": "A",
          "residue_id": 13,
          "insert_code": "",
          "x": 12.681,
          "y": 37.302,
          "z": -25.211,
          "occupancy": 1.0,
          "temperature_factor": 15.56,
          "element": "N",
          "charge": None,
          "model_id": 1
         }, {
          "atom_id": 108,
          "atom_name": "CA",
          "alt_loc": None,
          "residue_name": "GLY",
          "chain_id": "A",
          "residue_id": 13,
          "insert_code": "A",
          "x": 11.982,
          "y": 37.996,
          "z": -26.241,
          "occupancy": 1.2,
          "temperature_factor": 16.92,
          "element": "C",
          "charge": -1,
          "model_id": 1
         }
        ]
        for atom in atoms:
            self.blank.atoms().append(atom)
        pdb_file = self.blank.generate_pdb_file()
        self.assertEqual(len(pdb_file.records()), 2)
        self.assertEqual(
         pdb_file.records()[0].text(),
         "ATOM    107 N    GLY A  13    12.681  37.302  -25.211 1.0   15.56           N   "
        )
        self.assertEqual(
         pdb_file.records()[1].text(),
         "ATOM    108 CA   GLY A  13A   11.982  37.996  -26.241 1.2   16.92           C -1"
        )


    def test_can_handle_missing_atom_properties(self):
        atoms = [
         {
          "atom_id": 107,
          "atom_name": "N",
          "alt_loc": None,
          "residue_name": None,
          "chain_id": None,
          "residue_id": None,
          "insert_code": None,
          "x": 12.681,
          "y": 37.302,
          "z": -25.211,
          "occupancy": None,
          "temperature_factor": None,
          "element": "N",
          "charge": None,
          "model_id": 1
         }
        ]
        for atom in atoms:
            self.blank.atoms().append(atom)
        pdb_file = self.blank.generate_pdb_file()
        self.assertEqual(len(pdb_file.records()), 1)
        self.assertEqual(
         pdb_file.records()[0].text(),
         "ATOM    107 N                 12.681  37.302  -25.211                       N   "
        )


    def test_can_handle_awkward_coordinates(self):
        atoms = [
         {
          "atom_id": 107,
          "atom_name": "N",
          "alt_loc": None,
          "residue_name": None,
          "chain_id": None,
          "residue_id": None,
          "insert_code": None,
          "x": 3.9999999999998,
          "y": -7.9999999999999,
          "z": -8.881784197001252e-16,
          "occupancy": None,
          "temperature_factor": None,
          "element": "N",
          "charge": None,
          "model_id": 1
         }, {
          "atom_id": 107,
          "atom_name": "N",
          "alt_loc": None,
          "residue_name": None,
          "chain_id": None,
          "residue_id": None,
          "insert_code": None,
          "x": -5.000022760448198,
          "y": -0.7071067811865478,
          "z": 8.881784197001252e-16,
          "occupancy": None,
          "temperature_factor": None,
          "element": "N",
          "charge": None,
          "model_id": 1
         }
        ]
        for atom in atoms:
            self.blank.atoms().append(atom)
        pdb_file = self.blank.generate_pdb_file()
        self.assertEqual(len(pdb_file.records()), 2)
        self.assertEqual(
         pdb_file.records()[0].text(),
         "ATOM    107 N                 4.0     -8.0    0.0                           N   "
        )
        self.assertEqual(
         pdb_file.records()[1].text(),
         "ATOM    107 N                 -5.0    -0.707110.0                           N   "
        )



class HeteroatomsConversionTests(PdbDataFileTest):

    def test_can_produce_hetatm_records(self):
        atoms = [
         {
          "atom_id": 107,
          "atom_name": "N",
          "alt_loc": None,
          "residue_name": "GLY",
          "chain_id": "A",
          "residue_id": 13,
          "insert_code": "",
          "x": 12.681,
          "y": 37.302,
          "z": -25.211,
          "occupancy": 1.0,
          "temperature_factor": 15.56,
          "element": "N",
          "charge": None,
          "model_id": 1
         }, {
          "atom_id": 108,
          "atom_name": "CA",
          "alt_loc": None,
          "residue_name": "GLY",
          "chain_id": "A",
          "residue_id": 13,
          "insert_code": "A",
          "x": 11.982,
          "y": 37.996,
          "z": -26.241,
          "occupancy": 1.2,
          "temperature_factor": 16.92,
          "element": "C",
          "charge": -1,
          "model_id": 1
         }
        ]
        for atom in atoms:
            self.blank.heteroatoms().append(atom)
        pdb_file = self.blank.generate_pdb_file()
        self.assertEqual(len(pdb_file.records()), 2)
        self.assertEqual(
         pdb_file.records()[0].text(),
         "HETATM  107 N    GLY A  13    12.681  37.302  -25.211 1.0   15.56           N   "
        )
        self.assertEqual(
         pdb_file.records()[1].text(),
         "HETATM  108 CA   GLY A  13A   11.982  37.996  -26.241 1.2   16.92           C -1"
        )


    def test_can_handle_missing_hetatom_properties(self):
        atoms = [
         {
          "atom_id": 107,
          "atom_name": "N",
          "alt_loc": None,
          "residue_name": None,
          "chain_id": None,
          "residue_id": None,
          "insert_code": None,
          "x": 12.681,
          "y": 37.302,
          "z": -25.211,
          "occupancy": None,
          "temperature_factor": None,
          "element": "N",
          "charge": None,
          "model_id": 1
         }
        ]
        for atom in atoms:
            self.blank.heteroatoms().append(atom)
        pdb_file = self.blank.generate_pdb_file()
        self.assertEqual(len(pdb_file.records()), 1)
        self.assertEqual(
         pdb_file.records()[0].text(),
         "HETATM  107 N                 12.681  37.302  -25.211                       N   "
        )


    def test_can_handle_awkward_coordinates(self):
        atoms = [
         {
          "atom_id": 107,
          "atom_name": "N",
          "alt_loc": None,
          "residue_name": None,
          "chain_id": None,
          "residue_id": None,
          "insert_code": None,
          "x": 3.9999999999998,
          "y": -7.99999999999998,
          "z": 8.881784197001252e-16,
          "occupancy": None,
          "temperature_factor": None,
          "element": "N",
          "charge": None,
          "model_id": 1
         }, {
          "atom_id": 107,
          "atom_name": "N",
          "alt_loc": None,
          "residue_name": None,
          "chain_id": None,
          "residue_id": None,
          "insert_code": None,
          "x": -5.000022760448198,
          "y": -0.7071067811865478,
          "z": 8.881784197001252e-16,
          "occupancy": None,
          "temperature_factor": None,
          "element": "N",
          "charge": None,
          "model_id": 1
         }
        ]
        for atom in atoms:
            self.blank.heteroatoms().append(atom)
        pdb_file = self.blank.generate_pdb_file()
        self.assertEqual(len(pdb_file.records()), 2)
        self.assertEqual(
         pdb_file.records()[0].text(),
         "HETATM  107 N                 4.0     -8.0    0.0                           N   "
        )
        self.assertEqual(
         pdb_file.records()[1].text(),
         "HETATM  107 N                 -5.0    -0.707110.0                           N   "
        )



class ConnectionsConversionTests(PdbDataFileTest):

    def test_can_produce_conect_records(self):
        connections = [{
         "atom_id": 1179,
         "bonded_atoms": [746, 1184, 1195, 1203, 1211, 1222]
        }, {
         "atom_id": 11,
         "bonded_atoms": [746, 1184]
        }]
        for connection in connections:
            self.blank.connections().append(connection)
        pdb_file = self.blank.generate_pdb_file()
        self.assertEqual(len(pdb_file.records()), 3)
        self.assertEqual(
         pdb_file.records()[0].text(),
         "CONECT 1179  746 1184 1195 1203" + (" " * 49)
        )
        self.assertEqual(
         pdb_file.records()[1].text(),
         "CONECT 1179 1211 1222" + (" " * 59)
        )
        self.assertEqual(
         pdb_file.records()[2].text(),
         "CONECT   11  746 1184" + (" " * 59)
        )
