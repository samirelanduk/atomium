from unittest import TestCase
import unittest.mock
from molecupy.pdb.pdb import Pdb
from molecupy.pdb.pdbdatafile import PdbDataFile
from molecupy.structures import Model, SmallMolecule, Chain, Residue, BindSite
from molecupy.structures import AlphaHelix, BetaStrand, Complex

class PdbTest(TestCase):

    def setUp(self):
        self.data_file = unittest.mock.Mock(spec=PdbDataFile)
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


class PdbCreationTests(PdbTest):

    def test_can_create_pdb(self):
        pdb = Pdb(self.data_file)
        self.assertIs(pdb.data_file(), self.data_file)


    def test_can_get_data_attributes(self):
        self.data_file.classification.return_value = "val1",
        self.data_file.deposition_date.return_value = "val2",
        self.data_file.pdb_code.return_value = "val3",
        self.data_file.is_obsolete.return_value = "val4",
        self.data_file.obsolete_date.return_value = "val5",
        self.data_file.replacement_code.return_value = "val6",
        self.data_file.title.return_value = "val7",
        self.data_file.split_codes.return_value = "val8",
        self.data_file.caveat.return_value = "val9",
        self.data_file.keywords.return_value = "val10",
        self.data_file.experimental_techniques.return_value = "val11",
        self.data_file.model_count.return_value = "val12",
        self.data_file.model_annotations.return_value = "val13",
        self.data_file.authors.return_value = "val14",
        self.data_file.revisions.return_value = "val15",
        self.data_file.supercedes.return_value = "val16",
        self.data_file.supercede_date.return_value = "val17",
        self.data_file.journal.return_value = "val18"
        pdb = Pdb(self.data_file)
        self.assertIs(
         pdb.classification(),
         self.data_file.classification()
        )
        self.assertIs(
         pdb.deposition_date(),
         self.data_file.deposition_date()
        )
        self.assertIs(
         pdb.pdb_code(),
         self.data_file.pdb_code(),
        )
        self.assertIs(
         pdb.is_obsolete(),
         self.data_file.is_obsolete()
        )
        self.assertIs(
         pdb.obsolete_date(),
         self.data_file.obsolete_date()
        )
        self.assertIs(
         pdb.replacement_code(),
         self.data_file.replacement_code()
        )
        self.assertIs(
         pdb.title(),
         self.data_file.title()
        )
        self.assertIs(
         pdb.split_codes(),
         self.data_file.split_codes()
        )
        self.assertIs(
         pdb.caveat(),
         self.data_file.caveat()
        )
        self.assertIs(
         pdb.keywords(),
         self.data_file.keywords()
        )
        self.assertIs(
         pdb.experimental_techniques(),
         self.data_file.experimental_techniques()
        )
        self.assertIs(
         pdb.model_count(),
         self.data_file.model_count()
        )
        self.assertIs(
         pdb.model_annotations(),
         self.data_file.model_annotations()
        )
        self.assertIs(
         pdb.revisions(),
         self.data_file.revisions()
        )
        self.assertIs(
         pdb.supercedes(),
         self.data_file.supercedes()
        )
        self.assertIs(
         pdb.supercede_date(),
         self.data_file.supercede_date()
        )
        self.assertIs(
         pdb.journal(),
         self.data_file.journal()
        )


    def test_pdb_repr(self):
        pdb = Pdb(self.data_file)
        self.data_file.pdb_code.return_value = None
        self.assertEqual(str(pdb), "<Pdb (????)>")
        self.data_file.pdb_code.return_value = "1SAM"
        self.assertEqual(str(pdb), "<Pdb (1SAM)>")



class PdbModelsTests(PdbTest):

    def test_single_model(self):
        pdb = Pdb(self.data_file)
        self.assertEqual(len(pdb.models()), 1)
        self.assertIsInstance(pdb.models()[0], Model)


    def test_multiple_models(self):
        self.data_file.models.return_value = [
         {"model_id": 1, "start_record": 1, "end_record": 3},
         {"model_id": 2, "start_record": 4, "end_record": 6}
        ]
        pdb = Pdb(self.data_file)
        self.assertEqual(len(pdb.models()), 2)
        self.assertIsInstance(pdb.models()[0], Model)
        self.assertIsInstance(pdb.models()[1], Model)


    def test_one_model_access(self):
        self.data_file.models.return_value = [
         {"model_id": 1, "start_record": 1, "end_record": 3},
         {"model_id": 2, "start_record": 4, "end_record": 6}
        ]
        pdb = Pdb(self.data_file)
        self.assertIs(pdb.models()[0], pdb.model())



class PdbSmallMoleculeTests(PdbTest):

    def test_single_small_molecule(self):
        self.data_file.heteroatoms.return_value = [
         {
          "atom_id": 8237,
          "atom_name": "CA",
          "alt_loc": None,
          "residue_name": "123",
          "chain_id": "A",
          "residue_id": 1001,
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
          "atom_id": 8238,
          "atom_name": "MG",
          "alt_loc": None,
          "residue_name": "123",
          "chain_id": "A",
          "residue_id": 1001,
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
        ]
        pdb = Pdb(self.data_file)
        self.assertEqual(len(pdb.model().small_molecules()), 1)
        self.assertIsInstance(list(pdb.model().small_molecules())[0], SmallMolecule)
        self.assertEqual(len(list(pdb.model().small_molecules())[0].atoms()), 2)


    def test_multiple_small_molecules(self):
        self.data_file.heteroatoms.return_value = [
         {
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
         }
        ]
        pdb = Pdb(self.data_file)
        self.assertEqual(len(pdb.model().small_molecules()), 2)
        self.assertEqual(
         set([mol.molecule_name() for mol in pdb.model().small_molecules()]),
         set(["123", "MOL"])
        )
        self.assertEqual(
         set([mol.molecule_id() for mol in pdb.model().small_molecules()]),
         set(["A1002", "A1001A"])
        )


    def test_single_small_molecules_in_multiple_models(self):
        self.data_file.models.return_value = [
         {"model_id": 1, "start_record": 0, "end_record": 1},
         {"model_id": 2, "start_record": 2, "end_record": 3}
        ]
        self.data_file.heteroatoms.return_value = [
         {
          "atom_id": 8237,
          "atom_name": "CA",
          "alt_loc": None,
          "residue_name": "123",
          "chain_id": "A",
          "residue_id": 1001,
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
          "atom_id": 8238,
          "atom_name": "MG",
          "alt_loc": None,
          "residue_name": "123",
          "chain_id": "A",
          "residue_id": 1001,
          "insert_code": "",
          "x": 13.872,
          "y": -2.555,
          "z": -29.045,
          "occupancy": 1.0,
          "temperature_factor": 27.36,
          "element": "MG",
          "charge": None,
          "model_id": 1
         }, {
          "atom_id": 8237,
          "atom_name": "CA",
          "alt_loc": None,
          "residue_name": "123",
          "chain_id": "A",
          "residue_id": 1001,
          "insert_code": "",
          "x": 13.872,
          "y": -2.555,
          "z": -29.045,
          "occupancy": 1.0,
          "temperature_factor": 27.36,
          "element": "C",
          "charge": None,
          "model_id": 2
         }, {
          "atom_id": 8238,
          "atom_name": "MG",
          "alt_loc": None,
          "residue_name": "123",
          "chain_id": "A",
          "residue_id": 1001,
          "insert_code": "",
          "x": 13.872,
          "y": -2.555,
          "z": -29.045,
          "occupancy": 1.0,
          "temperature_factor": 27.36,
          "element": "MG",
          "charge": None,
          "model_id": 2
         }
        ]
        pdb = Pdb(self.data_file)
        self.assertEqual(len(pdb.models()[0].small_molecules()), 1)
        self.assertEqual(len(pdb.models()[1].small_molecules()), 1)
        self.assertIsNot(
         [pdb.models()[0].small_molecules()][0],
         [pdb.models()[1].small_molecules()][0]
        )



class PdbChainTests(PdbTest):

    def setUp(self):
        PdbTest.setUp(self)
        self.data_file.atoms.return_value = [
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
          "atom_id": 109,
          "atom_name": "N",
          "alt_loc": None,
          "residue_name": "MET",
          "chain_id": "A",
          "residue_id": 13,
          "insert_code": "A",
          "x": 12.681,
          "y": 37.302,
          "z": -25.211,
          "occupancy": 1.0,
          "temperature_factor": 15.56,
          "element": "N",
          "charge": None,
          "model_id": 1
         }, {
          "atom_id": 110,
          "atom_name": "CA",
          "alt_loc": None,
          "residue_name": "MET",
          "chain_id": "A",
          "residue_id": 13,
          "insert_code": "A",
          "x": 11.982,
          "y": 37.996,
          "z": -26.241,
          "occupancy": 1.0,
          "temperature_factor": 16.92,
          "element": "C",
          "charge": None,
          "model_id": 1
         }
        ]


    def test_single_chain(self):
        pdb = Pdb(self.data_file)
        self.assertEqual(len(pdb.model().chains()), 1)
        chain = list(pdb.model().chains())[0]
        self.assertIsInstance(chain, Chain)
        self.assertEqual(chain.chain_id(), "A")
        self.assertEqual(len(chain.residues()), 2)
        self.assertIsInstance(chain.residues()[0], Residue)
        self.assertEqual(chain.residues()[0].residue_id(), "A13")
        self.assertEqual(chain.residues()[0].residue_name(), "GLY")
        self.assertEqual(chain.residues()[1].residue_id(), "A13A")
        self.assertEqual(chain.residues()[1].residue_name(), "MET")


    def test_multiple_chains(self):
        atoms = self.data_file.atoms()
        atoms[2]["chain_id"] = atoms[3]["chain_id"] = "B"
        self.data_file.atoms.return_value = atoms
        pdb = Pdb(self.data_file)
        self.assertEqual(len(pdb.model().chains()), 2)
        self.assertEqual(
         set([chain.chain_id() for chain in pdb.model().chains()]),
         set(["A", "B"])
        )
        self.assertEqual(
         set([len(chain.residues()) for chain in pdb.model().chains()]),
         set([1, 1])
        )


    def test_multiple_models(self):
        self.data_file.models.return_value = [
         {"model_id": 1, "start_record": 0, "end_record": 1},
         {"model_id": 2, "start_record": 2, "end_record": 3}
        ]
        atoms = self.data_file.atoms()
        atoms[2]["model_id"] = atoms[3]["model_id"] = 2
        self.data_file.atoms.return_value = atoms
        pdb = Pdb(self.data_file)
        self.assertEqual(len(pdb.models()[0].chains()), 1)
        self.assertEqual(len(pdb.models()[1].chains()), 1)
        self.assertIsNot(
         [pdb.models()[0].chains()][0],
         [pdb.models()[1].chains()][0]
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
        pdb = Pdb(self.data_file)
        chain = list(pdb.model().chains())[0]
        self.assertEqual(len(chain.residues()), 5)
        self.assertEqual(
         [res.residue_id() for res in chain.residues()],
         ["A12", "A13", "A13A", "A14", "A15"]
        )
        self.assertNotIn("PRO", [res.residue_name() for res in chain.residues()])



class PdbBondTests(PdbTest):

    def setUp(self):
        PdbTest.setUp(self)
        self.data_file.atoms.return_value = [
         {
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
         }
        ]
        self.data_file.heteroatoms.return_value = [
         {
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
         }
        ]


    def test_small_molecules_bonded_together(self):
        self.data_file.connections.return_value = [
         {
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
         }
        ]
        pdb = Pdb(self.data_file)
        atom1 = pdb.model().get_atom_by_id(8237)
        atom2 = pdb.model().get_atom_by_id(8238)
        atom3 = pdb.model().get_atom_by_id(8239)
        atom4 = pdb.model().get_atom_by_id(8240)
        self.assertEqual(atom1.bonded_atoms(), set([atom2]))
        self.assertEqual(atom2.bonded_atoms(), set([atom1]))
        self.assertEqual(atom3.bonded_atoms(), set([atom4]))
        self.assertEqual(atom4.bonded_atoms(), set([atom3]))


    def test_can_handle_self_conect_records(self):
        self.data_file.connections.return_value = [
         {
          "atom_id": 8237,
          "bonded_atoms": [8238]
         }, {
          "atom_id": 8238,
          "bonded_atoms": [8237]
         }, {
          "atom_id": 8239,
          "bonded_atoms": [8239]
         }
        ]
        pdb = Pdb(self.data_file)


    def test_residues_are_connected_internally(self):
        pdb = Pdb(self.data_file)
        atom1 = pdb.model().get_atom_by_id(131)
        atom2 = pdb.model().get_atom_by_id(132)
        atom3 = pdb.model().get_atom_by_id(133)
        atom4 = pdb.model().get_atom_by_id(134)
        atom5 = pdb.model().get_atom_by_id(135)
        atom6 = pdb.model().get_atom_by_id(136)
        atom7 = pdb.model().get_atom_by_id(137)
        atom8 = pdb.model().get_atom_by_id(138)
        atom9 = pdb.model().get_atom_by_id(139)
        atom10 = pdb.model().get_atom_by_id(140)
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
        pdb = Pdb(self.data_file)
        chain = list(pdb.model().chains())[0]
        self.assertIs(chain.residues()[0].upstream_residue(), None)
        self.assertIs(chain.residues()[0].downstream_residue(), chain.residues()[1])
        self.assertIs(chain.residues()[1].upstream_residue(), chain.residues()[0])
        self.assertIs(chain.residues()[1].downstream_residue(), None)
        atom1 = pdb.model().get_atom_by_id(133)
        atom2 = pdb.model().get_atom_by_id(136)
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
        pdb = Pdb(self.data_file)
        atom1 = pdb.model().get_atom_by_id(135)
        atom2 = pdb.model().get_atom_by_id(140)
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
        pdb = Pdb(self.data_file)


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
        pdb = Pdb(self.data_file)
        atom1 = pdb.model().get_atom_by_id(135)
        atom2 = pdb.model().get_atom_by_id(140)
        self.assertIn(atom1, atom2.bonded_atoms())
        self.assertIn(atom2, atom1.bonded_atoms())



class PdbBindSiteTests(PdbBondTests):

    def setUp(self):
        PdbBondTests.setUp(self)
        self.data_file.sites.return_value = [{
         "site_id": "AC1",
         "residue_count": 2,
         "residues": [
          {"residue_name": "ALA", "chain_id": "A", "residue_id": 27, "insert_code": ""},
          {"residue_name": "ALA", "chain_id": "A", "residue_id": 28, "insert_code": ""},
         ]
        }]


    def test_bind_sites_present(self):
        pdb = Pdb(self.data_file)
        self.assertEqual(len(pdb.model().bind_sites()), 1)
        site = list(pdb.model().bind_sites())[0]
        self.assertIsInstance(site, BindSite)
        self.assertEqual(site.site_id(), "AC1")
        self.assertEqual(
         site.residues(),
         set(list(pdb.model().chains())[0].residues())
        )


    def test_bind_site_and_ligand_know_about_each_other(self):
        self.data_file.get_remark_by_number.return_value = {
         'content': 'SITE\n'
         'SITE_IDENTIFIER: AC1\n'
         'EVIDENCE_CODE: SOFTWARE\n'
         'SITE_DESCRIPTION: BINDING SITE FOR RESIDUE MOL A1002\n\n',
         'number': 465
        }
        pdb = Pdb(self.data_file)
        site = list(pdb.model().bind_sites())[0]
        molecule = pdb.model().get_small_molecule_by_id("A1002")
        self.assertIs(molecule.bind_site(), site)
        self.assertIs(site.ligand(), molecule)


    def test_site_can_discard_residues_on_absent_chain(self):
        self.data_file.sites.return_value = [{
         "site_id": "AC1",
         "residue_count": 3,
         "residues": [
          {"residue_name": "ALA", "chain_id": "A", "residue_id": 27, "insert_code": ""},
          {"residue_name": "ALA", "chain_id": "A", "residue_id": 28, "insert_code": ""},
          {"residue_name": "ALA", "chain_id": "Q", "residue_id": 1, "insert_code": ""},
         ]
        }]
        pdb = Pdb(self.data_file)
        site = list(pdb.model().bind_sites())[0]
        self.assertEqual(
         site.residues(),
         set(list(pdb.model().chains())[0].residues())
        )



class PdbSecondaryStructureTests(PdbBondTests):

    def setUp(self):
        PdbBondTests.setUp(self)


    def test_pdb_has_alpha_helices(self):
        self.data_file.helices.return_value = [{
         "helix_id": 1,
         "helix_name": "1",
         "start_residue_name": "ALA",
         "start_residue_chain_id": "A",
         "start_residue_id": 27,
         "start_residue_insert": "",
         "end_residue_name": "ALA",
         "end_residue_chain_id": "A",
         "end_residue_id": 28,
         "end_residue_insert": "",
         "helix_class": 5,
         "comment": None,
         "length": 3
        }]
        pdb = Pdb(self.data_file)
        self.assertEqual(len(pdb.model().get_chain_by_id("A").alpha_helices()), 1)
        helix = list(pdb.model().get_chain_by_id("A").alpha_helices())[0]
        self.assertIsInstance(helix, AlphaHelix)
        self.assertEqual(helix.helix_id(), "1")
        self.assertEqual(
         helix.residues(),
         list(pdb.model().chains())[0].residues()
        )
        self.assertEqual(helix.helix_class(), "Right-handed 3 - 10")


    def test_pdb_has_beta_strands(self):
        self.data_file.sheets.return_value = [{
         "sheet_id": "A",
         "strand_count": 1,
         "strands": [{
          "strand_id": 1,
          "start_residue_name": "ALA",
          "start_residue_chain_id": "A",
          "start_residue_id": 27,
          "start_residue_insert": "",
          "end_residue_name": "ALA",
          "end_residue_chain_id": "A",
          "end_residue_id": 28,
          "end_residue_insert": "",
          "sense": 0,
          "current_atom": None,
          "current_residue_name": None,
          "current_chain_id": None,
          "current_residue_id": None,
          "current_insert": "",
          "previous_atom": None,
          "previous_residue_name": None,
          "previous_chain_id": None,
          "previous_residue_id": None,
          "previous_insert": ""
         }]
        }]
        pdb = Pdb(self.data_file)
        self.assertEqual(len(pdb.model().get_chain_by_id("A").beta_strands()), 1)
        strand = list(pdb.model().get_chain_by_id("A").beta_strands())[0]
        self.assertIsInstance(strand, BetaStrand)
        self.assertEqual(strand.strand_id(), "1")
        self.assertEqual(
         strand.residues(),
         list(pdb.model().chains())[0].residues()
        )
        self.assertEqual(strand.sense(), 0)



class PdbComplexTests(PdbTest):

    def setUp(self):
        PdbTest.setUp(self)
        self.data_file.atoms.return_value = [
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
          "atom_id": 109,
          "atom_name": "N",
          "alt_loc": None,
          "residue_name": "MET",
          "chain_id": "B",
          "residue_id": 13,
          "insert_code": "A",
          "x": 12.681,
          "y": 37.302,
          "z": -25.211,
          "occupancy": 1.0,
          "temperature_factor": 15.56,
          "element": "N",
          "charge": None,
          "model_id": 1
         }, {
          "atom_id": 110,
          "atom_name": "CA",
          "alt_loc": None,
          "residue_name": "MET",
          "chain_id": "B",
          "residue_id": 13,
          "insert_code": "A",
          "x": 11.982,
          "y": 37.996,
          "z": -26.241,
          "occupancy": 1.0,
          "temperature_factor": 16.92,
          "element": "C",
          "charge": None,
          "model_id": 1
         }
        ]
        self.data_file.compounds.return_value = [{
         "MOL_ID": 1,
         "MOLECULE": "COMPLEX1",
         "CHAIN": ["A", "B"]
        }]


    def test_can_add_single_complex(self):
        pdb = Pdb(self.data_file)
        self.assertEqual(len(pdb.model().complexes()), 1)
        complex_ = list(pdb.model().complexes())[0]
        self.assertIsInstance(complex_, Complex)
        self.assertEqual(complex_.complex_id(), "1")
        self.assertEqual(complex_.complex_name(), "COMPLEX1")
        self.assertEqual(len(complex_.chains()), 2)
