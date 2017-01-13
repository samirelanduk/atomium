from unittest import TestCase
from unittest.mock import Mock
from molecupy.converters.pdbdatafile2model import model_from_pdb_data_file
from molecupy.pdb.pdbdatafile import PdbDataFile
from molecupy.structures import Model, SmallMolecule, Chain, Residue, BindSite
from molecupy.structures import AlphaHelix, BetaStrand, Complex

class ModelCreationTest(TestCase):

    def setUp(self):
        self.data_file = Mock(spec=PdbDataFile)
        self.data_file.models.return_value = [
         {"model_id": 1, "start_record": 0, "end_record": 0}
        ]
        self.data_file.heteroatoms.return_value = []
        self.data_file.atoms.return_value = []
        self.data_file.get_remark_by_number.return_value = {}
        self.data_file.connections.return_value = []
        self.data_file.ss_bonds.return_value = []
        self.data_file.links.return_value = []
        self.data_file.sites.return_value = []
        self.data_file.helices.return_value = []
        self.data_file.sheets.return_value = []
        self.data_file.compounds.return_value = []

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



class SmallMoleculeTests(ModelCreationTest):

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



class ChainTests(ModelCreationTest):

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
