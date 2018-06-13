from unittest import TestCase
from unittest.mock import patch, Mock, MagicMock
from atomium.files.pdb2pdbdict import *

class PdbToPdbDictTests(TestCase):

    @patch("atomium.files.pdb2pdbdict.structure_to_pdb_dict")
    @patch("atomium.files.pdb2pdbdict.sequences_from_model")
    def test_can_convert_pdb_to_pdb_dict_one_model(self, mock_seq, mock_dict):
        mock_seq.return_value = "SEQ"
        pdb = Mock()
        pdb._deposition_date = "D"
        pdb._code = "C"
        pdb._title = "T"
        pdb._resolution = 1.5
        pdb._rfactor = 1.8
        pdb._rfree = 2.4
        pdb._rcount = 18
        pdb._organism = "O"
        pdb._expression_system = "E"
        pdb._technique = "T"
        pdb._classification = "CLASS"
        pdb._keywords = ["a", "b"]
        pdb._models = ["model1"]
        mock_dict.return_value = {"models": ["m1"], "connections": ["c1", "c2"]}
        pdb_dict = pdb_to_pdb_dict(pdb)
        mock_dict.assert_called_with("model1")
        mock_seq.assert_called_with("model1")
        self.assertEqual(pdb_dict, {
         "deposition_date": "D", "code": "C", "title": "T", "resolution": 1.5,
         "organism": "O", "expression_system": "E", "technique": "T",
         "classification": "CLASS", "rfactor": 1.8, "keywords": ["a", "b"],
         "rfree": 2.4, "rcount": 18,
         "models": ["m1"], "connections": ["c1", "c2"], "sequences": "SEQ"
        })


    @patch("atomium.files.pdb2pdbdict.structure_to_pdb_dict")
    @patch("atomium.files.pdb2pdbdict.sequences_from_model")
    def test_can_convert_pdb_to_pdb_dict_two_models(self, mock_seq, mock_dict):
        mock_seq.return_value = "SEQ"
        pdb = Mock()
        pdb._deposition_date = "D"
        pdb._code = "C"
        pdb._title = "T"
        pdb._resolution = 1.5
        pdb._rfactor = 1.8
        pdb._rfree = 2.4
        pdb._rcount = 18
        pdb._organism = "O"
        pdb._expression_system = "E"
        pdb._technique = "T"
        pdb._classification = "CLASS"
        pdb._keywords = ["a", "b"]
        pdb._models = ["model1", "model2"]
        mock_dict.side_effect = [
         {"models": ["m1"], "connections": ["c1", "c2"]},
         {"models": ["m2"], "connections": ["c1", "c2"]}
        ]
        pdb_dict = pdb_to_pdb_dict(pdb)
        mock_dict.assert_any_call("model1")
        mock_dict.assert_any_call("model2")
        mock_seq.assert_called_with("model1")
        self.assertEqual(pdb_dict, {
         "deposition_date": "D", "code": "C", "title": "T", "resolution": 1.5,
         "organism": "O", "expression_system": "E", "technique": "T",
         "classification": "CLASS", "rfactor": 1.8, "keywords": ["a", "b"],
         "rfree": 2.4, "rcount": 18,
         "models": ["m1", "m2"], "connections": ["c1", "c2"], "sequences": "SEQ"
        })



class StructureToPdbDictTests(TestCase):

    @patch("atomium.files.pdb2pdbdict.atom_to_atom_dict")
    @patch("atomium.files.pdb2pdbdict.atoms_to_chains")
    @patch("atomium.files.pdb2pdbdict.structure_to_connections")
    def test_can_convert_model_to_pdb_dict(self, mock_con, mock_chain, mock_atom):
        structure = Mock()
        atoms = [Mock(), Mock(), Mock(), Mock(), Mock(), Mock()]
        mock_atom.side_effect = lambda a: "a" + str(a.id)
        chains = [Mock(), Mock()]
        for index, atom in enumerate(atoms):
            atom.id = index + 1
            atom.residue = "res" if index < 4 else None
        mock_chain.return_value = ["chain1", "chain2"]
        structure.atoms.return_value = set(atoms)
        mock_con.return_value = ["c1", "c2"]
        pdb_dict = structure_to_pdb_dict(structure)
        for atom in atoms:
            mock_atom.assert_any_call(atom)
        mock_chain.assert_called_with(["a1", "a2", "a3", "a4"], ["a5", "a6"])
        mock_con.assert_called_with(structure)
        self.assertEqual(pdb_dict, {
         "deposition_date": None,
         "code": None,
         "title": None,
         "resolution": None,
         "rfactor": None,
         "rfree": None,
         "rcount": None,
         "organism": None,
         "expression_system": None,
         "technique": None,
         "classification": None,
         "keywords": [],
         "sequences": {},
         "models": [{
          "chains": ["chain1", "chain2"]
         }],
         "connections": ["c1", "c2"]
        })



class AtomToAtomDictTests(TestCase):

    def setUp(self):
        self.atom = Mock()
        self.atom.id = 107
        self.atom.name = "N1"
        self.atom.x = 12.681
        self.atom.y = 37.302
        self.atom.z = -25.211
        self.atom.element = "N"
        self.atom.charge = -2
        self.atom.bfactor = 12.5
        self.atom.anisotropy = [1, 2, 0, 3, 3, 4]
        self.residue = Mock()
        self.residue.name = "GLY"
        self.residue.id = "A13B"
        self.chain = Mock()
        self.chain.id = "A"
        self.atom.residue = self.residue
        self.atom.chain = self.chain


    def test_can_convert_full_atom_to_dict(self):
        d = atom_to_atom_dict(self.atom)
        self.assertEqual(d, {
         "atom_id": 107, "atom_name": "N1", "alt_loc": None,
         "residue_name": "GLY",
         "chain_id": "A", "residue_id": 13, "insert_code": "B", "full_id": "A13B",
         "x": 12.681, "y": 37.302, "z": -25.211,
         "occupancy": 1.0, "temp_factor": 12.5, "anisotropy": [1, 2, 0, 3, 3, 4],
         "element": "N", "charge": -2,
        })


    def test_can_convert_full_heteroatom_to_dict(self):
        self.atom.residue = None
        mol = Mock()
        mol.id = "A200"
        mol.name = "SUC"
        self.atom.ligand = mol
        d = atom_to_atom_dict(self.atom)
        self.assertEqual(d, {
         "atom_id": 107, "atom_name": "N1", "alt_loc": None,
         "residue_name": "SUC",
         "chain_id": "A", "residue_id": 200, "insert_code": "", "full_id": "A200",
         "x": 12.681, "y": 37.302, "z": -25.211,
         "occupancy": 1.0, "temp_factor": 12.5, "anisotropy": [1, 2, 0, 3, 3, 4],
         "element": "N", "charge": -2,
        })


    def test_can_convert_atom_with_no_insert_code(self):
        self.residue.id = "A13"
        d = atom_to_atom_dict(self.atom)
        self.assertEqual(d["residue_id"], 13)
        self.assertEqual(d["insert_code"], "")
        self.assertEqual(d["full_id"], "A13")


    def test_can_convert_atom_with_no_chain_in_residue_id(self):
        self.residue.id = "13"
        d = atom_to_atom_dict(self.atom)
        self.assertEqual(d["residue_id"], 13)
        self.assertEqual(d["chain_id"], "A")
        self.assertEqual(d["insert_code"], "")
        self.assertEqual(d["full_id"], "13")



class StructureToConnectionsTests(TestCase):

    def test_can_convert_structure_to_connections(self):
        atoms = [Mock(), Mock(), Mock(), Mock(), Mock(), Mock(), Mock()]
        for i, atom in enumerate(atoms):
            atom.residue = "residue" if i < 4 else None
            atom.id = i + 1
        atoms[4].bonded_atoms = set(atoms[5:])
        atoms[5].bonded_atoms = set([atoms[4], atoms[6]])
        atoms[6].bonded_atoms = set(atoms[4:6])
        model = Mock()
        model.atoms.return_value = set(atoms)
        connections = structure_to_connections(model)
        self.assertEqual(connections, [{
         "atom": 5, "bond_to": [6, 7]
        }, {
         "atom": 6, "bond_to": [5, 7]
        }, {
         "atom": 7, "bond_to": [5, 6]
        }])



class SequencesFromModelTests(TestCase):

    def test_can_get_sequences_from_model(self):
        chaina, chainb = Mock(id="A", rep_sequence="TXV"), Mock(id="B", rep_sequence="UM")
        model = Mock()
        model.chains.return_value = [chaina, Mock(rep_sequence=""), chainb]
        sequences = sequences_from_model(model)
        self.assertEqual(sequences, {
         "A": ["THR", "???", "VAL"], "B": ["???", "MET"]
        })
