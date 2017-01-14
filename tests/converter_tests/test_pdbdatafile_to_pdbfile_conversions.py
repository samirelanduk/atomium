from unittest import TestCase
from molecupy.pdb.pdbfile import PdbFile
from molecupy.pdb.pdbdatafile import PdbDataFile
from molecupy.converters.pdbdatafile2pdbfile import pdb_file_from_pdb_data_file

class PdbDataFile2PdbFileTest(TestCase):
    pass



class BasicPdbFileCreationTests(PdbDataFile2PdbFileTest):

    def test_can_create_pdb_file(self):
        pdb_file = pdb_file_from_pdb_data_file(PdbDataFile())
        self.assertIsInstance(pdb_file, PdbFile)


    def test_can_only_convert_pdb_data_files(self):
        with self.assertRaises(TypeError):
            pdb_file_from_pdb_data_file("PDB file")


    def test_pdb_file_knows_source(self):
        data_file = PdbDataFile()
        pdb_file = pdb_file_from_pdb_data_file(data_file)
        self.assertIs(pdb_file.source(), data_file)

'''class PdbDataFileTest(TestCase):

    def setUp(self):
        self.empty = PdbDataFile(PdbFile(""))
        self.blank = PdbDataFile()


    def add_compounds_to_blank(self, compounds):
        for compound in compounds:
            self.blank.compounds().append(compound)


    def add_atoms_to_blank(self, atoms):
        for atom in atoms:
            self.blank.atoms().append(atom)


    def add_heteroatoms_to_blank(self, heteroatoms):
        for atom in heteroatoms:
            self.blank.heteroatoms().append(atom)


    def add_connections_to_blank(self, connections):
        for connection in connections:
            self.blank.connections().append(connection)



class PdbFileCreationTests(PdbDataFileTest):

    def test_can_make_basic_pdb_file(self):
        pdb_file = self.blank.generate_pdb_file()
        self.assertIsInstance(pdb_file, PdbFile)
        self.assertEqual(pdb_file.records(), [])



class CompoundsConversionTests(PdbDataFileTest):

    def test_can_convert_simple_compound(self):
        self.add_compounds_to_blank([{
         "MOL_ID": 1,
         "MOLECULE": "OROTIDINE 5'-MONOPHOSPHATE DECARBOXYLASE",
         "CHAIN": ["A", "B"]
        }])
        pdb_file = self.blank.generate_pdb_file()
        self.assertEqual(
         "\n".join([line.rstrip() for line in pdb_file.convert_to_string().split("\n")]),
         "COMPND    MOL_ID: 1;\n"
         "COMPND   2 MOLECULE: OROTIDINE 5'-MONOPHOSPHATE DECARBOXYLASE;\n"
         "COMPND   3 CHAIN: A, B;"
        )


    def test_can_convert_longer_compound(self):
        self.add_compounds_to_blank([{
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
        }])
        pdb_file = self.blank.generate_pdb_file()
        self.assertEqual(
         "\n".join([line.rstrip() for line in pdb_file.convert_to_string().split("\n")]),
         "COMPND    MOL_ID: 1;\n"
         "COMPND   2 MOLECULE: OROTIDINE 5'-MONOPHOSPHATE DECARBOXYLASE;\n"
         "COMPND   3 CHAIN: A, B;\n"
         "COMPND   4 SYNONYM: OMP DECARBOXYLASE, OMPDCASE, OMPDECASE;\n"
         "COMPND   5 EC: 4.1.1.23;\n"
         "COMPND   6 ENGINEERED: YES;"
        )


    def test_can_convert_multiple_compounds(self):
        self.add_compounds_to_blank([{
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
        }])
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

    def setUp(self):
        PdbDataFileTest.setUp(self)
        self.atom = {
         "atom_id": 107,
         "atom_name": "N",
         "alt_loc": None,
         "residue_name": "GLY",
         "chain_id": "A",
         "residue_id": 13,
         "insert_code": "A",
         "x": 12.681,
         "y": 37.302,
         "z": -25.211,
         "occupancy": 1.0,
         "temperature_factor": 15.56,
         "element": "N",
         "charge": 1,
         "model_id": 1
        }


    def test_can_convert_atom(self):
        self.add_atoms_to_blank([self.atom])
        pdb_file = self.blank.generate_pdb_file()
        self.assertEqual(len(pdb_file.records()), 1)
        self.assertEqual(
         pdb_file.records()[0].text(),
         "ATOM    107 N    GLY A  13A   12.681  37.302  -25.211 1.0   15.56           N 1 "
        )


    def test_can_convert_multiple_atoms(self):
        self.add_atoms_to_blank([self.atom, {
         "atom_id": 108,
         "atom_name": "CA",
         "alt_loc": None,
         "residue_name": "GLY",
         "chain_id": "A",
         "residue_id": 13,
         "insert_code": "",
         "x": 11.982,
         "y": 37.996,
         "z": -26.241,
         "occupancy": 1.2,
         "temperature_factor": 16.92,
         "element": "C",
         "charge": -1,
         "model_id": 1
        }])
        pdb_file = self.blank.generate_pdb_file()
        self.assertEqual(len(pdb_file.records()), 2)
        self.assertEqual(
         pdb_file.records()[0].text(),
         "ATOM    107 N    GLY A  13A   12.681  37.302  -25.211 1.0   15.56           N 1 "
        )
        self.assertEqual(
         pdb_file.records()[1].text(),
         "ATOM    108 CA   GLY A  13    11.982  37.996  -26.241 1.2   16.92           C -1"
        )


    def test_can_handle_missing_atom_properties(self):
        self.atom["residue_name"] = None
        self.atom["chain_id"] = None
        self.atom["residue_id"] = None
        self.atom["insert_code"] = None
        self.atom["occupancy"] = None
        self.atom["temperature_factor"] = None
        self.atom["charge"] = None
        self.add_atoms_to_blank([self.atom])
        pdb_file = self.blank.generate_pdb_file()
        self.assertEqual(
         pdb_file.records()[0].text(),
         "ATOM    107 N                 12.681  37.302  -25.211                       N   "
        )


    def test_atom_conversion_can_round_coordinates(self):
        self.atom["x"] = 3.9999999999998
        self.atom["y"] = -7.9999999999999
        self.atom["z"] = -5.000022760448198
        self.add_atoms_to_blank([self.atom])
        pdb_file = self.blank.generate_pdb_file()
        self.assertEqual(
         pdb_file.records()[0].text(),
         "ATOM    107 N    GLY A  13A   4.0     -8.0    -5.0    1.0   15.56           N 1 "
        )


    def test_can_handle_scientific_notation_in_coordinates(self):
        self.atom["x"] = 8.881784197001252e-16
        self.atom["y"] = -0.7071067811865478
        self.atom["z"] = -8.881784197001252e-16
        self.add_atoms_to_blank([self.atom])
        pdb_file = self.blank.generate_pdb_file()
        self.assertEqual(
         pdb_file.records()[0].text(),
         "ATOM    107 N    GLY A  13A   0.0     -0.707110.0     1.0   15.56           N 1 "
        )


    def test_atom_conversion_will_produce_hetatm_records(self):
        self.add_heteroatoms_to_blank([self.atom])
        pdb_file = self.blank.generate_pdb_file()
        self.assertEqual(len(pdb_file.records()), 1)
        self.assertEqual(
         pdb_file.records()[0].text(),
         "HETATM  107 N    GLY A  13A   12.681  37.302  -25.211 1.0   15.56           N 1 "
        )



class ConnectionsConversionTests(PdbDataFileTest):

    def test_can_convert_connection(self):
        self.add_connections_to_blank([{
         "atom_id": 11,
         "bonded_atoms": [746, 1184]
        }])
        pdb_file = self.blank.generate_pdb_file()
        self.assertEqual(len(pdb_file.records()), 1)
        self.assertEqual(
         pdb_file.records()[0].text(),
         "CONECT   11  746 1184" + (" " * 59)
        )


    def test_can_convert_connection_with_more_than_four_bonds(self):
        self.add_connections_to_blank([{
         "atom_id": 1179,
         "bonded_atoms": [746, 1184, 1195, 1203, 1211, 1222]
        }])
        pdb_file = self.blank.generate_pdb_file()
        self.assertEqual(len(pdb_file.records()), 2)
        self.assertEqual(
         pdb_file.records()[0].text(),
         "CONECT 1179  746 1184 1195 1203" + (" " * 49)
        )
        self.assertEqual(
         pdb_file.records()[1].text(),
         "CONECT 1179 1211 1222" + (" " * 59)
        )


    def test_can_convert_multiple_connection(self):
        self.add_connections_to_blank([{
         "atom_id": 1179,
         "bonded_atoms": [746, 1184, 1195, 1203, 1211, 1222]
        }, {
         "atom_id": 11,
         "bonded_atoms": [746, 1184]
        }])
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
        )'''
