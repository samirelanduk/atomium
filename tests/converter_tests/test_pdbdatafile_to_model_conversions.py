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
