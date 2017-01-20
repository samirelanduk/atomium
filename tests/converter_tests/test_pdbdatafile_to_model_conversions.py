from unittest import TestCase
from unittest.mock import Mock
from molecupy.converters.pdbdatafile2model import model_from_pdb_data_file
from molecupy.converters.pdbdatafile2model import _residue_id_is_greater_than_residue_id
from molecupy.converters.pdbdatafile2model import _residue_id_to_int
from molecupy.converters.pdbdatafile2model import _create_missing_residue_from_id
from molecupy.converters.pdbdatafile2model import _get_missing_residue_info
from molecupy.converters.pdbdatafile2model import _mol_id_from_atom
from molecupy.converters.pdbdatafile2model import _add_residues_to_chain
from molecupy.converters.pdbdatafile2model import _add_missing_residues_to_chain
from molecupy.converters.pdbdatafile2model import _get_top_atom_id
from molecupy.pdb.pdbdatafile import PdbDataFile
from molecupy.structures import Model, SmallMolecule, Chain, Residue, BindSite
from molecupy.structures import AlphaHelix, BetaStrand, Complex

class PdbDataFile2ModelTest(TestCase):

    def setUp(self):
        self.data_file = Mock(PdbDataFile)
        self.data_file.models.return_value = [{
         "model_id": 1,
         "start_record": -1,
         "end_record": -1
        }]
        self.data_file.compounds.return_value = []
        self.data_file.helices.return_value = []
        self.data_file.sheets.return_value = []
        self.data_file.ss_bonds.return_value = []
        self.data_file.links.return_value = []
        self.data_file.sites.return_value = []
        self.data_file.atoms.return_value = []
        self.data_file.heteroatoms.return_value = []
        self.data_file.connections.return_value = []
        self.data_file.get_remark_by_number.return_value = {}
        self.atom1 = {
         "atom_id": 8237,
         "atom_name": "CA",
         "alt_loc": None,
         "residue_name": "123",
         "chain_id": "A",
         "residue_id": 1001,
         "insert_code": "A",
         "x": 13.872,
         "y": -2.555,
         "z": -29.045,
         "occupancy": 1.0,
         "temperature_factor": 27.36,
         "element": "C",
         "charge": None,
         "model_id": 1
        }
        self.atom2 = {
         "atom_id": 8238,
         "atom_name": "MG",
         "alt_loc": None,
         "residue_name": "123",
         "chain_id": "A",
         "residue_id": 1001,
         "insert_code": "A",
         "x": 13.872,
         "y": -2.555,
         "z": -29.045,
         "occupancy": 1.0,
         "temperature_factor": 27.36,
         "element": "MG",
         "charge": None,
         "model_id": 1
        }
        self.atom3 = {
         "atom_id": 8239,
         "atom_name": "CA",
         "alt_loc": None,
         "residue_name": "MOL",
         "chain_id": "A",
         "residue_id": 1002,
         "insert_code": "",
         "x": 13.872,
         "y": -2.555,
         "z": -29.045,
         "occupancy": 1.0,
         "temperature_factor": 27.36,
         "element": "C",
         "charge": None,
         "model_id": 1
        }
        self.atom4 = {
         "atom_id": 8240,
         "atom_name": "MG",
         "alt_loc": None,
         "residue_name": "MOL",
         "chain_id": "A",
         "residue_id": 1002,
         "insert_code": "",
         "x": 13.872,
         "y": -2.555,
         "z": -29.045,
         "occupancy": 1.0,
         "temperature_factor": 27.36,
         "element": "MG",
         "charge": None,
         "model_id": 1
        }



class BasicModelCreationTests(PdbDataFile2ModelTest):

    def test_can_create_model(self):
        model = model_from_pdb_data_file(self.data_file)
        self.assertIsInstance(model, Model)


    def test_can_only_convert_pdb_data_files(self):
        with self.assertRaises(TypeError):
            model_from_pdb_data_file("PDB file")


    def test_data_file_knows_source(self):
        model = model_from_pdb_data_file(self.data_file)
        self.assertIs(model.source(), self.data_file)


    def test_can_get_specific_model(self):
        self.data_file.models.return_value = [{
         "model_id": 2,
         "start_record": 3,
         "end_record": 5
        }]
        model2 = model_from_pdb_data_file(self.data_file, model_id=2)
        self.assertIsInstance(model2, Model)
        with self.assertRaises(ValueError):
            model_from_pdb_data_file(self.data_file)
        with self.assertRaises(ValueError):
            model_from_pdb_data_file(self.data_file, model_id=1)



class SmallMoleculeCreationTests(PdbDataFile2ModelTest):

    def test_single_small_molecule(self):
        self.data_file.heteroatoms.return_value = [self.atom1, self.atom2]
        model = model_from_pdb_data_file(self.data_file, 1)
        self.assertEqual(len(model.small_molecules()), 1)
        self.assertIsInstance(list(model.small_molecules())[0], SmallMolecule)
        self.assertEqual(len(list(model.small_molecules())[0].atoms()), 2)


    def test_multiple_small_molecules(self):
        self.data_file.heteroatoms.return_value = [
         self.atom1, self.atom2, self.atom3, self.atom4
        ]
        model = model_from_pdb_data_file(self.data_file, 1)
        self.assertEqual(len(model.small_molecules()), 2)
        self.assertEqual(
         set([mol.molecule_name() for mol in model.small_molecules()]),
         set(["123", "MOL"])
        )
        self.assertEqual(
         set([mol.molecule_id() for mol in model.small_molecules()]),
         set(["A1002", "A1001A"])
        )


    def test_single_small_molecules_in_multiple_models(self):
        self.data_file.models.return_value = [
         {"model_id": 1, "start_record": 0, "end_record": 1},
         {"model_id": 2, "start_record": 2, "end_record": 3}
        ]
        self.atom3["model_id"] = 2
        self.atom4["model_id"] = 2
        self.data_file.heteroatoms.return_value = [
         self.atom1, self.atom2, self.atom3, self.atom4
        ]
        model1 = model_from_pdb_data_file(self.data_file, 1)
        model2 = model_from_pdb_data_file(self.data_file, 2)
        self.assertEqual(len(model1.small_molecules()), 1)
        self.assertEqual(len(model2.small_molecules()), 1)
        self.assertIsNot(
         [model1.small_molecules()][0],
         [model2.small_molecules()][0]
        )



class ChainTests(PdbDataFile2ModelTest):

    def test_single_chain(self):
        self.data_file.atoms.return_value = [
         self.atom1, self.atom2, self.atom3, self.atom4
        ]
        model = model_from_pdb_data_file(self.data_file, 1)
        self.assertEqual(len(model.chains()), 1)
        chain = list(model.chains())[0]
        self.assertIsInstance(chain, Chain)
        self.assertEqual(chain.chain_id(), "A")
        self.assertEqual(len(chain.residues()), 2)
        self.assertIsInstance(chain.residues()[0], Residue)
        self.assertEqual(chain.residues()[0].residue_id(), "A1001A")
        self.assertEqual(chain.residues()[0].residue_name(), "123")
        self.assertEqual(chain.residues()[1].residue_id(), "A1002")
        self.assertEqual(chain.residues()[1].residue_name(), "MOL")


    def test_multiple_chains(self):
        self.data_file.atoms.return_value = [
         self.atom1, self.atom2, self.atom3, self.atom4
        ]
        self.atom3["chain_id"] = self.atom4["chain_id"] = "B"
        model = model_from_pdb_data_file(self.data_file, 1)
        self.assertEqual(len(model.chains()), 2)
        self.assertEqual(
         set([chain.chain_id() for chain in model.chains()]),
         set(["A", "B"])
        )
        self.assertEqual(
         set([len(chain.residues()) for chain in model.chains()]),
         set([1, 1])
        )


    def test_multiple_models(self):
        self.data_file.models.return_value = [
         {"model_id": 1, "start_record": 0, "end_record": 1},
         {"model_id": 2, "start_record": 2, "end_record": 3}
        ]
        self.atom3["model_id"] = 2
        self.atom4["model_id"] = 2
        self.data_file.atoms.return_value = [
         self.atom1, self.atom2, self.atom3, self.atom4
        ]
        model1 = model_from_pdb_data_file(self.data_file, 1)
        model2 = model_from_pdb_data_file(self.data_file, 2)
        self.assertEqual(len(model1.chains()), 1)
        self.assertEqual(len(model2.chains()), 1)
        self.assertIsNot(
         [model1.chains()][0],
         [model2.chains()][0]
        )


    def test_missing_residues(self):
        self.data_file.get_remark_by_number.return_value = {
         'content': 'MISSING RESIDUES\n'
         'THE FOLLOWING RESIDUES WERE NOT LOCATED IN THE\n'
         'EXPERIMENT. (M=MODEL NUMBER; RES=RESIDUE NAME; C=CHAIN\n'
         'IDENTIFIER; SSSEQ=SEQUENCE NUMBER; I=INSERTION CODE.)\n\n'
         'M RES C SSSEQI\n'
         'LEU A     12\n'
         'ARG A     14\n',
         'number': 465
        }

        self.data_file.atoms.return_value = [
         self.atom1, self.atom2, self.atom3, self.atom4
        ]
        for atom in self.data_file.atoms.return_value:
            atom["residue_id"] = 13
        self.atom1["insert_code"] = self.atom2["insert_code"] = ""
        self.atom3["insert_code"] = self.atom4["insert_code"] = "A"
        model = model_from_pdb_data_file(self.data_file)
        chain = list(model.chains())[0]
        self.assertEqual(len(chain.residues()), 4)
        self.assertEqual(
         [res.residue_id() for res in chain.residues()],
         ["A12", "A13", "A13A", "A14"]
        )
        self.assertTrue(chain.residues()[0].is_missing())
        self.assertFalse(chain.residues()[1].is_missing())
        self.assertFalse(chain.residues()[2].is_missing())
        self.assertTrue(chain.residues()[3].is_missing())
        self.assertEqual(len(chain.residues()[0].atoms(atom_type="all")), 22)
        self.assertEqual(len(chain.residues()[1].atoms(atom_type="all")), 2)
        self.assertEqual(len(chain.residues()[2].atoms(atom_type="all")), 2)
        self.assertEqual(len(chain.residues()[3].atoms(atom_type="all")), 27)


    def test_missing_residues_on_multiple_chains(self):
        self.data_file.get_remark_by_number.return_value = {
         'content': 'MISSING RESIDUES\n'
         'THE FOLLOWING RESIDUES WERE NOT LOCATED IN THE\n'
         'EXPERIMENT. (M=MODEL NUMBER; RES=RESIDUE NAME; C=CHAIN\n'
         'IDENTIFIER; SSSEQ=SEQUENCE NUMBER; I=INSERTION CODE.)\n\n'
         'M RES C SSSEQI\n'
         'ALA A     12\n'
         'XXX A     14\n'
         'ALA B     12\n'
         'XXX B     14\n',
         'number': 465
        }
        self.data_file.atoms.return_value = [
         self.atom1, self.atom2, self.atom3, self.atom4
        ]
        self.atom3["chain_id"] = self.atom4["chain_id"] = "B"
        self.atom2["residue_id"] = 1002
        self.atom3["residue_id"] = 1001
        self.atom1["insert_code"] = self.atom2["insert_code"] = ""
        model = model_from_pdb_data_file(self.data_file)
        chainA, chainB = model.get_chain_by_id("A"), model.get_chain_by_id("B")
        self.assertEqual(
         set([a.atom_id() for a in chainA.get_residue_by_id("A1001").atoms(atom_type="all")]),
         set([8237])
        )
        self.assertEqual(
         set([a.atom_id() for a in chainA.get_residue_by_id("A1002").atoms(atom_type="all")]),
         set([8238])
        )
        self.assertEqual(
         set([a.atom_id() for a in chainB.get_residue_by_id("B1001").atoms(atom_type="all")]),
         set([8239])
        )
        self.assertEqual(
         set([a.atom_id() for a in chainB.get_residue_by_id("B1002").atoms(atom_type="all")]),
         set([8240])
        )
        self.assertEqual(
         set([a.atom_id() for a in chainA.get_residue_by_id("A12").atoms(atom_type="all")]),
         set([8241 + n for n in range(13)])
        )
        self.assertEqual(
         set([a.atom_id() for a in chainA.get_residue_by_id("A14").atoms(atom_type="all")]),
         set([8254 + n for n in range(3)])
        )
        self.assertEqual(
         set([a.atom_id() for a in chainB.get_residue_by_id("B12").atoms(atom_type="all")]),
         set([8257 + n for n in range(13)])
        )
        self.assertEqual(
         set([a.atom_id() for a in chainB.get_residue_by_id("B14").atoms(atom_type="all")]),
         set([8270 + n for n in range(3)])
        )


    def test_handling_of_duplicate_missing_residues(self):
        self.data_file.get_remark_by_number.return_value = {
         'content': 'MISSING RESIDUES\n'
         'THE FOLLOWING RESIDUES WERE NOT LOCATED IN THE\n'
         'EXPERIMENT. (M=MODEL NUMBER; RES=RESIDUE NAME; C=CHAIN\n'
         'IDENTIFIER; SSSEQ=SEQUENCE NUMBER; I=INSERTION CODE.)\n\n'
         'M RES C SSSEQI\n'
         'LEU A     12\n'
         'ARG A     14\n'
         'PRO A     14\n'
         'CYS A     15\n',
         'number': 465
        }
        self.data_file.atoms.return_value = [
         self.atom1, self.atom2, self.atom3, self.atom4
        ]
        for atom in self.data_file.atoms.return_value:
            atom["residue_id"] = 13
        self.atom1["insert_code"] = self.atom2["insert_code"] = ""
        self.atom3["insert_code"] = self.atom4["insert_code"] = "A"
        model = model_from_pdb_data_file(self.data_file)
        chain = list(model.chains())[0]
        self.assertEqual(len(chain.residues()), 5)
        self.assertEqual(
         [res.residue_id() for res in chain.residues()],
         ["A12", "A13", "A13A", "A14", "A15"]
        )
        self.assertNotIn("PRO", [res.residue_name() for res in chain.residues()])



class BondTests(PdbDataFile2ModelTest):

    def setUp(self):
        PdbDataFile2ModelTest.setUp(self)
        self.data_file.atoms.return_value = [{
         "atom_id": 131,
         "atom_name": "N",
         "alt_loc": None,
         "residue_name": "ALA",
         "chain_id": "A",
         "residue_id": 27,
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
         "atom_id": 132,
         "atom_name": "CA",
         "alt_loc": None,
         "residue_name": "ALA",
         "chain_id": "A",
         "residue_id": 27,
         "insert_code": "",
         "x": 11.982,
         "y": 37.996,
         "z": -26.241,
         "occupancy": 1.0,
         "temperature_factor": 16.92,
         "element": "C",
         "charge": None,
         "model_id": 1
        }, {
         "atom_id": 133,
         "atom_name": "C",
         "alt_loc": None,
         "residue_name": "ALA",
         "chain_id": "A",
         "residue_id": 27,
         "insert_code": "",
         "x": 12.681,
         "y": 37.302,
         "z": -25.211,
         "occupancy": 1.0,
         "temperature_factor": 15.56,
         "element": "C",
         "charge": None,
         "model_id": 1
        }, {
         "atom_id": 134,
         "atom_name": "O",
         "alt_loc": None,
         "residue_name": "ALA",
         "chain_id": "A",
         "residue_id": 27,
         "insert_code": "",
         "x": 11.982,
         "y": 37.996,
         "z": -26.241,
         "occupancy": 1.0,
         "temperature_factor": 16.92,
         "element": "O",
         "charge": None,
         "model_id": 1
        }, {
         "atom_id": 135,
         "atom_name": "CB",
         "alt_loc": None,
         "residue_name": "ALA",
         "chain_id": "A",
         "residue_id": 27,
         "insert_code": "",
         "x": 11.982,
         "y": 37.996,
         "z": -26.241,
         "occupancy": 1.0,
         "temperature_factor": 16.92,
         "element": "C",
         "charge": None,
         "model_id": 1
        }, {
         "atom_id": 136,
         "atom_name": "N",
         "alt_loc": None,
         "residue_name": "ALA",
         "chain_id": "A",
         "residue_id": 28,
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
         "atom_id": 137,
         "atom_name": "CA",
         "alt_loc": None,
         "residue_name": "ALA",
         "chain_id": "A",
         "residue_id": 28,
         "insert_code": "",
         "x": 11.982,
         "y": 37.996,
         "z": -26.241,
         "occupancy": 1.0,
         "temperature_factor": 16.92,
         "element": "C",
         "charge": None,
         "model_id": 1
        }, {
         "atom_id": 138,
         "atom_name": "C",
         "alt_loc": None,
         "residue_name": "ALA",
         "chain_id": "A",
         "residue_id": 28,
         "insert_code": "",
         "x": 12.681,
         "y": 37.302,
         "z": -25.211,
         "occupancy": 1.0,
         "temperature_factor": 15.56,
         "element": "C",
         "charge": None,
         "model_id": 1
        }, {
         "atom_id": 139,
         "atom_name": "O",
         "alt_loc": None,
         "residue_name": "ALA",
         "chain_id": "A",
         "residue_id": 28,
         "insert_code": "",
         "x": 11.982,
         "y": 37.996,
         "z": -26.241,
         "occupancy": 1.0,
         "temperature_factor": 16.92,
         "element": "O",
         "charge": None,
         "model_id": 1
        }, {
         "atom_id": 140,
         "atom_name": "CB",
         "alt_loc": None,
         "residue_name": "ALA",
         "chain_id": "A",
         "residue_id": 28,
         "insert_code": "",
         "x": 11.982,
         "y": 37.996,
         "z": -26.241,
         "occupancy": 1.0,
         "temperature_factor": 16.92,
         "element": "C",
         "charge": None,
         "model_id": 1
        }]
        self.data_file.heteroatoms.return_value = [{
         "atom_id": 8237,
         "atom_name": "CA",
         "alt_loc": None,
         "residue_name": "123",
         "chain_id": "A",
         "residue_id": 1001,
         "insert_code": "A",
         "x": 13.872,
         "y": -2.555,
         "z": -29.045,
         "occupancy": 1.0,
         "temperature_factor": 27.36,
         "element": "C",
         "charge": None,
         "model_id": 1
        }, {
         "atom_id": 8238,
         "atom_name": "MG",
         "alt_loc": None,
         "residue_name": "123",
         "chain_id": "A",
         "residue_id": 1001,
         "insert_code": "A",
         "x": 13.872,
         "y": -2.555,
         "z": -29.045,
         "occupancy": 1.0,
         "temperature_factor": 27.36,
         "element": "MG",
         "charge": None,
         "model_id": 1
        }, {
         "atom_id": 8239,
         "atom_name": "CA",
         "alt_loc": None,
         "residue_name": "MOL",
         "chain_id": "A",
         "residue_id": 1002,
         "insert_code": "",
         "x": 13.872,
         "y": -2.555,
         "z": -29.045,
         "occupancy": 1.0,
         "temperature_factor": 27.36,
         "element": "C",
         "charge": None,
         "model_id": 1
        }, {
         "atom_id": 8240,
         "atom_name": "MG",
         "alt_loc": None,
         "residue_name": "MOL",
         "chain_id": "A",
         "residue_id": 1002,
         "insert_code": "",
         "x": 13.872,
         "y": -2.555,
         "z": -29.045,
         "occupancy": 1.0,
         "temperature_factor": 27.36,
         "element": "MG",
         "charge": None,
         "model_id": 1
        }]


    def test_small_molecules_bonded_together(self):
        self.data_file.connections.return_value = [{
          "atom_id": 8237,
          "bonded_atoms": [8238]
         }, {
          "atom_id": 8238,
          "bonded_atoms": [8237]
         }, {
          "atom_id": 8239,
          "bonded_atoms": [8240]
         }, {
          "atom_id": 8240,
          "bonded_atoms": [8239]
         }]
        model = model_from_pdb_data_file(self.data_file)
        atom1 = model.get_atom_by_id(8237)
        atom2 = model.get_atom_by_id(8238)
        atom3 = model.get_atom_by_id(8239)
        atom4 = model.get_atom_by_id(8240)
        self.assertEqual(atom1.bonded_atoms(), set([atom2]))
        self.assertEqual(atom2.bonded_atoms(), set([atom1]))
        self.assertEqual(atom3.bonded_atoms(), set([atom4]))
        self.assertEqual(atom4.bonded_atoms(), set([atom3]))


    def test_can_handle_self_conect_records(self):
        self.data_file.connections.return_value = [{
          "atom_id": 8237,
          "bonded_atoms": [8238]
         }, {
          "atom_id": 8238,
          "bonded_atoms": [8237]
         }, {
          "atom_id": 8239,
          "bonded_atoms": [8239]
         }]
        model_from_pdb_data_file(self.data_file)


    def test_residues_are_connected_internally(self):
        model = model_from_pdb_data_file(self.data_file)
        atom1 = model.get_atom_by_id(131)
        atom2 = model.get_atom_by_id(132)
        atom3 = model.get_atom_by_id(133)
        atom4 = model.get_atom_by_id(134)
        atom5 = model.get_atom_by_id(135)
        atom6 = model.get_atom_by_id(136)
        atom7 = model.get_atom_by_id(137)
        atom8 = model.get_atom_by_id(138)
        atom9 = model.get_atom_by_id(139)
        atom10 = model.get_atom_by_id(140)
        self.assertEqual(atom1.bonded_atoms(), set([atom2]))
        self.assertEqual(atom2.bonded_atoms(), set([atom1, atom3, atom5]))
        self.assertEqual(atom3.bonded_atoms(), set([atom2, atom4, atom6]))
        self.assertEqual(atom4.bonded_atoms(), set([atom3]))
        self.assertEqual(atom5.bonded_atoms(), set([atom2]))
        self.assertEqual(atom6.bonded_atoms(), set([atom7, atom3]))
        self.assertEqual(atom7.bonded_atoms(), set([atom6, atom8, atom10]))
        self.assertEqual(atom8.bonded_atoms(), set([atom7, atom9]))
        self.assertEqual(atom9.bonded_atoms(), set([atom8]))
        self.assertEqual(atom10.bonded_atoms(), set([atom7]))


    def test_residues_are_connected_to_each_other(self):
        model = model_from_pdb_data_file(self.data_file)
        chain = list(model.chains())[0]
        self.assertIs(chain.residues()[0].upstream_residue(), None)
        self.assertIs(chain.residues()[0].downstream_residue(), chain.residues()[1])
        self.assertIs(chain.residues()[1].upstream_residue(), chain.residues()[0])
        self.assertIs(chain.residues()[1].downstream_residue(), None)
        atom1 = model.get_atom_by_id(133)
        atom2 = model.get_atom_by_id(136)
        self.assertIn(atom1, atom2.bonded_atoms())
        self.assertIn(atom2, atom1.bonded_atoms())


    def test_disulphide_bonds_are_present(self):
        self.data_file.ss_bonds.return_value = [{
         "serial_num": 1,
         "residue_name_1": "ALA",
         "chain_id_1": "A",
         "residue_id_1": 27,
         "insert_code_1": "",
         "residue_name_2": "ALA",
         "chain_id_2": "A",
         "residue_id_2": 28,
         "insert_code_2": "",
         "symmetry_1": "1555",
         "symmetry_2": "1555",
         "length": 2.04
        }]
        atoms = self.data_file.atoms()
        atoms[4]["element"] = atoms[9]["element"] = "S"
        self.data_file.atoms.return_value = atoms
        model = model_from_pdb_data_file(self.data_file)
        atom1 = model.get_atom_by_id(135)
        atom2 = model.get_atom_by_id(140)
        self.assertIn(atom1, atom2.bonded_atoms())
        self.assertIn(atom2, atom1.bonded_atoms())


    def test_can_account_for_disulphide_bonds_incorrectly_assigned_to_same_atom(self):
        self.data_file.ss_bonds.return_value = [{
         "serial_num": 1,
         "residue_name_1": "ALA",
         "chain_id_1": "A",
         "residue_id_1": 27,
         "insert_code_1": "",
         "residue_name_2": "ALA",
         "chain_id_2": "A",
         "residue_id_2": 27,
         "insert_code_2": "",
         "symmetry_1": "1555",
         "symmetry_2": "1555",
         "length": 2.04
        }]
        atoms = self.data_file.atoms()
        atoms[4]["element"] = atoms[9]["element"] = "S"
        self.data_file.atoms.return_value = atoms
        model = model_from_pdb_data_file(self.data_file)


    def test_link_bonds_are_present(self):
        self.data_file.links.return_value = [{
         "atom_1": "CB",
         "alt_loc_1": None,
         "residue_name_1": "ALA",
         "chain_id_1": "A",
         "residue_id_1": 27,
         "insert_code_1": "",
         "atom_2": "CB",
         "alt_loc_2": None,
         "residue_name_2": "ALA",
         "chain_id_2": "A",
         "residue_id_2": 28,
         "insert_code_2": "",
         "symmetry_1": "1555",
         "symmetry_2": "1555",
         "length": 2.75
        }]
        model = model_from_pdb_data_file(self.data_file)
        atom1 = model.get_atom_by_id(135)
        atom2 = model.get_atom_by_id(140)
        self.assertIn(atom1, atom2.bonded_atoms())
        self.assertIn(atom2, atom1.bonded_atoms())



class BindSiteTests(PdbDataFile2ModelTest):

    def setUp(self):
        PdbDataFile2ModelTest.setUp(self)
        self.data_file.sites.return_value = [{
         "site_id": "AC1",
         "residue_count": 2,
         "residues": [
          {"residue_name": "ALA", "chain_id": "A", "residue_id": 1001, "insert_code": "A"},
          {"residue_name": "ALA", "chain_id": "A", "residue_id": 1002, "insert_code": ""},
         ]
        }]
        self.data_file.atoms.return_value = [
         self.atom1, self.atom2, self.atom3, self.atom4
        ]
        self.data_file.heteroatoms.return_value = [self.atom3, self.atom4]


    def test_bind_sites_present(self):
        model = model_from_pdb_data_file(self.data_file)
        self.assertEqual(len(model.bind_sites()), 1)
        site = list(model.bind_sites())[0]
        self.assertIsInstance(site, BindSite)
        self.assertEqual(site.site_id(), "AC1")
        self.assertEqual(
         site.residues(),
         set(list(model.chains())[0].residues())
        )


    def test_bind_site_and_ligand_know_about_each_other(self):
        self.data_file.get_remark_by_number.return_value = {
         'content': 'SITE\n'
         'SITE_IDENTIFIER: AC1\n'
         'EVIDENCE_CODE: SOFTWARE\n'
         'SITE_DESCRIPTION: BINDING SITE FOR RESIDUE MOL A1002\n\n',
         'number': 465
        }
        model = model_from_pdb_data_file(self.data_file)
        site = list(model.bind_sites())[0]
        molecule = model.get_small_molecule_by_id("A1002")
        self.assertIs(molecule.bind_site(), site)
        self.assertIs(site.ligand(), molecule)


    def test_site_can_discard_residues_on_absent_chain(self):
        self.data_file.sites.return_value = [{
         "site_id": "AC1",
         "residue_count": 3,
         "residues": [
          {"residue_name": "ALA", "chain_id": "A", "residue_id": 1001, "insert_code": "A"},
          {"residue_name": "ALA", "chain_id": "A", "residue_id": 1002, "insert_code": ""},
          {"residue_name": "ALA", "chain_id": "Q", "residue_id": 1, "insert_code": ""},
         ]
        }]
        model = model_from_pdb_data_file(self.data_file)
        site = list(model.bind_sites())[0]
        self.assertEqual(
         site.residues(),
         set(list(model.chains())[0].residues())
        )



class MolIdFromAtomTests(TestCase):

    def test_can_get_basic_id(self):
        self.assertEqual(
         _mol_id_from_atom({"chain_id": "A", "residue_id": 1, "insert_code": ""}),
         "A1"
        )


    def test_can_get_id_with_insert_code(self):
        self.assertEqual(
         _mol_id_from_atom({"chain_id": "A", "residue_id": 1, "insert_code": "A"}),
         "A1A"
        )



class SafeAtomIdTests(PdbDataFile2ModelTest):

    def test_can_get_safe_id_from_atoms(self):
        self.data_file.atoms.return_value = [self.atom1, self.atom3]
        self.assertEqual(_get_top_atom_id(self.data_file.atoms()), 8239)


    def test_can_get_safe_id_from_atoms(self):
        self.data_file.heteroatoms.return_value = [self.atom2, self.atom4]
        self.assertEqual(_get_top_atom_id(
         self.data_file.atoms(), self.data_file.heteroatoms()
        ), 8240)


    def test_can_get_safe_id_from_all_atoms(self):
        self.data_file.atoms.return_value = [self.atom1, self.atom3]
        self.data_file.heteroatoms.return_value = [self.atom2, self.atom4]
        self.assertEqual(_get_top_atom_id(
         self.data_file.atoms(), self.data_file.heteroatoms()
        ), 8240)



class ResidueAddingToChainTests(PdbDataFile2ModelTest):

    def setUp(self):
        PdbDataFile2ModelTest.setUp(self)
        self.atoms = [self.atom1, self.atom2, self.atom3, self.atom4]


    def test_can_add_residues(self):
        residues = []
        _add_residues_to_chain(residues, "A", self.atoms)
        self.assertEqual(len(residues), 2)
        self.assertIsInstance(residues[0], Residue)
        self.assertEqual(residues[0].residue_id(), "A1001A")
        self.assertEqual(residues[0].residue_name(), "123")
        self.assertEqual(residues[1].residue_id(), "A1002")
        self.assertEqual(residues[1].residue_name(), "MOL")


    def test_can_ignore_atoms_from_other_chains(self):
        self.atom3["chain_id"] = self.atom4["chain_id"] = "B"
        residues = []
        _add_residues_to_chain(residues, "A", self.atoms)
        self.assertEqual(len(residues), 1)



class MissingResidueAddingToChainTests(PdbDataFile2ModelTest):

    def setUp(self):
        PdbDataFile2ModelTest.setUp(self)
        residue1 = Mock(Residue)
        residue1.residue_id.return_value = "A98"
        residue2 = Mock(Residue)
        residue2.residue_id.return_value = "A99"
        self.residues = [residue1, residue2]
        self.missing_info = [("MET", "A100"), ("VAL", "A101")]


    def test_can_add_missing_residues(self):
        _add_missing_residues_to_chain(self.residues, self.missing_info, 100)
        self.assertEqual(len(self.residues), 4)


    def test_missing_residues_will_order_correctly(self):
        _add_missing_residues_to_chain(self.residues, self.missing_info, 100)
        self.assertEqual(len(self.residues), 4)
        self.assertEqual(self.residues[2].residue_id(), "A100")
        self.assertEqual(self.residues[3].residue_id(), "A101")


    def test_missing_residues_will_order_correctly_when_slotting_needed(self):
        self.residues[-1].residue_id.return_value = "A105"
        _add_missing_residues_to_chain(self.residues, self.missing_info, 100)
        self.assertEqual(len(self.residues), 4)
        self.assertEqual(self.residues[1].residue_id(), "A100")
        self.assertEqual(self.residues[2].residue_id(), "A101")


    def test_missing_residues_will_order_correctly_with_insert_codes(self):
        self.residues[-1].residue_id.return_value = "A101A"
        _add_missing_residues_to_chain(self.residues, self.missing_info, 100)
        self.assertEqual(len(self.residues), 4)
        self.assertEqual(self.residues[1].residue_id(), "A100")
        self.assertEqual(self.residues[2].residue_id(), "A101")


    def test_missing_residues_atoms_ordered_properly(self):
        _add_missing_residues_to_chain(self.residues, self.missing_info, 100)
        self.assertEqual(
         min([a.atom_id() for a in self.residues[2].atoms(atom_type="all")]),
         101
        )
        self.assertEqual(
         max([a.atom_id() for a in self.residues[2].atoms(atom_type="all")]),
         min([a.atom_id() for a in self.residues[3].atoms(atom_type="all")]) - 1
        )



class MissingResidueInfoRetrievalTests(PdbDataFile2ModelTest):

    def test_data_file_with_no_relevant_remark_returns_none(self):
        self.assertIs(_get_missing_residue_info(self.data_file, "A"), None)


    def test_can_parse_remark(self):
        self.data_file.get_remark_by_number.return_value = {
         'content': 'MISSING RESIDUES\n'
         'THE FOLLOWING RESIDUES WERE NOT LOCATED IN THE\n'
         'EXPERIMENT. (M=MODEL NUMBER; RES=RESIDUE NAME; C=CHAIN\n'
         'IDENTIFIER; SSSEQ=SEQUENCE NUMBER; I=INSERTION CODE.)\n\n'
         'M RES C SSSEQI\n'
         'LEU A     12\n'
         'ARG A     14A\n',
         'number': 465
        }
        self.assertEqual(
         _get_missing_residue_info(self.data_file, "A"),
         [("LEU", "A12"), ("ARG", "A14A")]
        )


    def test_can_handle_duplicate_residue_ids(self):
        self.data_file.get_remark_by_number.return_value = {
         'content': 'MISSING RESIDUES\n'
         'THE FOLLOWING RESIDUES WERE NOT LOCATED IN THE\n'
         'EXPERIMENT. (M=MODEL NUMBER; RES=RESIDUE NAME; C=CHAIN\n'
         'IDENTIFIER; SSSEQ=SEQUENCE NUMBER; I=INSERTION CODE.)\n\n'
         'M RES C SSSEQI\n'
         'LEU A     12\n'
         'ARG A     14\n'
         'PRO A     14\n'
         'CYS A     15\n',
         'number': 465
        }
        self.assertEqual(
         _get_missing_residue_info(self.data_file, "A"),
         [("LEU", "A12"), ("ARG", "A14"), ("CYS", "A15")]
        )


    def test_can_filter_out_other_chains_residues(self):
        self.data_file.get_remark_by_number.return_value = {
         'content': 'MISSING RESIDUES\n'
         'THE FOLLOWING RESIDUES WERE NOT LOCATED IN THE\n'
         'EXPERIMENT. (M=MODEL NUMBER; RES=RESIDUE NAME; C=CHAIN\n'
         'IDENTIFIER; SSSEQ=SEQUENCE NUMBER; I=INSERTION CODE.)\n\n'
         'M RES C SSSEQI\n'
         'LEU A     12\n'
         'ARG A     14A\n'
         'ARG B     14\n',
         'number': 465
        }
        self.assertEqual(
         _get_missing_residue_info(self.data_file, "A"),
         [("LEU", "A12"), ("ARG", "A14A")]
        )



class ResidueFromMissingIdTests(TestCase):

    def test_can_create_missing_residue(self):
        residue = _create_missing_residue_from_id(("VAL", "A1"), 100)
        self.assertIsInstance(residue, Residue)
        self.assertTrue(residue.is_missing())


    def test_known_residues_are_correct(self):
        residue = _create_missing_residue_from_id(("GLY", "A1"), 100)
        self.assertEqual(len(residue.atoms(atom_type="all")), 10)
        self.assertEqual(residue.residue_name(), "GLY")


    def test_unknown_residues_have_skeleton(self):
        residue = _create_missing_residue_from_id(("XXX", "A1"), 100)
        self.assertEqual(len(residue.atoms(atom_type="all")), 3)
        self.assertEqual(residue.residue_name(), "XXX")
        self.assertEqual(
         set([a.atom_name() for a in residue.atoms(atom_type="all")]),
         set(["C", "CA", "N"])
        )


    def test_atom_ids_count_on_from_starting_id(self):
        residue = _create_missing_residue_from_id(("GLY", "A1"), 100)
        self.assertEqual(
         set([a.atom_id() for a in residue.atoms(atom_type="all")]),
         set([100, 101, 102, 103, 104, 105, 106, 107, 108, 109])
        )



class ResidueIdToIntTests(TestCase):

    def test_can_convert_residue_ids_to_integers(self):
        self.assertEqual(_residue_id_to_int("A10"), 1000)
        self.assertEqual(_residue_id_to_int("A1"), 100)
        self.assertEqual(_residue_id_to_int("A-30"), -3000)


    def test_can_convert_insert_codes_to_integers(self):
        self.assertEqual(_residue_id_to_int("A10A"), 1001)
        self.assertEqual(_residue_id_to_int("A104Z"), 10426)



class ResidueIdOrderingTests(TestCase):

    def test_can_order_basic_ids_correctly(self):
        self.assertTrue(_residue_id_is_greater_than_residue_id("A10", "A9"))
        self.assertFalse(_residue_id_is_greater_than_residue_id("A9", "A10"))


    def test_can_order_negative_ids(self):
        self.assertTrue(_residue_id_is_greater_than_residue_id("A5", "A-5"))
        self.assertTrue(_residue_id_is_greater_than_residue_id("A-2", "A-5"))
        self.assertFalse(_residue_id_is_greater_than_residue_id("A-5", "A5"))
        self.assertFalse(_residue_id_is_greater_than_residue_id("A-5", "A-2"))


    def test_can_order_insert_codes(self):
        self.assertTrue(_residue_id_is_greater_than_residue_id("A9A", "A9"))
        self.assertTrue(_residue_id_is_greater_than_residue_id("A9B", "A9A"))
        self.assertTrue(_residue_id_is_greater_than_residue_id("A10", "A9C"))
        self.assertFalse(_residue_id_is_greater_than_residue_id("A9", "A9A"))
        self.assertFalse(_residue_id_is_greater_than_residue_id("A9A", "A9B"))
        self.assertFalse(_residue_id_is_greater_than_residue_id("A9C", "A10"))
