from unittest import TestCase
from unittest.mock import patch, Mock, MagicMock
from atomium.files.pdb2pdbdict import *

class PdbToPdbDictTests(TestCase):

    @patch("atomium.files.pdb2pdbdict.structure_to_pdb_dict")
    def test_can_convert_pdb_to_pdb_dict_one_model(self, mock_dict):
        pdb = Mock()
        pdb._deposition_date = "D"
        pdb._code = "C"
        pdb._title = "T"
        pdb._resolution = 1.5
        pdb._rfactor = 1.8
        pdb._organism = "O"
        pdb._expression_system = "E"
        pdb._technique = "T"
        pdb._classification = "CLASS"
        pdb._models = ["model1"]
        mock_dict.return_value = {"models": ["m1"], "connections": ["c1", "c2"]}
        pdb_dict = pdb_to_pdb_dict(pdb)
        mock_dict.assert_called_with("model1")
        self.assertEqual(pdb_dict, {
         "deposition_date": "D", "code": "C", "title": "T", "resolution": 1.5,
         "organism": "O", "expression_system": "E", "technique": "T",
         "classification": "CLASS", "rfactor": 1.8,
         "models": ["m1"], "connections": ["c1", "c2"]
        })


    @patch("atomium.files.pdb2pdbdict.structure_to_pdb_dict")
    def test_can_convert_pdb_to_pdb_dict_two_models(self, mock_dict):
        pdb = Mock()
        pdb._deposition_date = "D"
        pdb._code = "C"
        pdb._title = "T"
        pdb._resolution = 1.5
        pdb._rfactor = 1.8
        pdb._organism = "O"
        pdb._expression_system = "E"
        pdb._technique = "T"
        pdb._classification = "CLASS"
        pdb._models = ["model1", "model2"]
        mock_dict.side_effect = [
         {"models": ["m1"], "connections": ["c1", "c2"]},
         {"models": ["m2"], "connections": ["c1", "c2"]}
        ]
        pdb_dict = pdb_to_pdb_dict(pdb)
        mock_dict.assert_any_call("model1")
        mock_dict.assert_any_call("model2")
        self.assertEqual(pdb_dict, {
         "deposition_date": "D", "code": "C", "title": "T", "resolution": 1.5,
         "organism": "O", "expression_system": "E", "technique": "T",
         "classification": "CLASS", "rfactor": 1.8,
         "models": ["m1", "m2"], "connections": ["c1", "c2"]
        })



class StructureToPdbDictTests(TestCase):

    @patch("atomium.files.pdb2pdbdict.atom_to_atom_dict")
    @patch("atomium.files.pdb2pdbdict.atoms_to_chains")
    @patch("atomium.files.pdb2pdbdict.atoms_to_residues")
    @patch("atomium.files.pdb2pdbdict.structure_to_connections")
    def test_can_convert_model_to_pdb_dict(self, mock_con, mock_res, mock_chain, mock_atom):
        structure = Mock()
        atoms = [Mock(), Mock(), Mock(), Mock(), Mock(), Mock()]
        mock_atom.side_effect = lambda a: "a" + str(a.atom_id())
        chains = [Mock(), Mock()]
        for index, atom in enumerate(atoms):
            atom.atom_id.return_value = index + 1
            atom.chain.return_value = chains[index // 2] if index < 4 else None
        mock_chain.return_value = ["chain1", "chain2"]
        mock_res.return_value = ["mol1", "mol2"]
        structure.atoms.return_value = set(atoms)
        mock_con.return_value = ["c1", "c2"]
        pdb_dict = structure_to_pdb_dict(structure)
        for atom in atoms:
            mock_atom.assert_any_call(atom)
        mock_chain.assert_called_with(["a1", "a2", "a3", "a4"])
        mock_res.assert_called_with(["a5", "a6"])
        mock_con.assert_called_with(structure)
        self.assertEqual(pdb_dict, {
         "deposition_date": None,
         "code": None,
         "title": None,
         "resolution": None,
         "rfactor": None,
         "organism": None,
         "expression_system": None,
         "technique": None,
         "classification": None,
         "models": [{
          "chains": ["chain1", "chain2"],
          "molecules": ["mol1", "mol2"]
         }],
         "connections": ["c1", "c2"]
        })


    @patch("atomium.files.pdb2pdbdict.atom_to_atom_dict")
    @patch("atomium.files.pdb2pdbdict.atoms_to_chains")
    @patch("atomium.files.pdb2pdbdict.atoms_to_residues")
    @patch("atomium.files.pdb2pdbdict.structure_to_connections")
    def test_can_convert_model_to_pdb_dict_no_atom_id(self, mcok_con,  mock_res, mock_chain, mock_atom):
        structure = Mock()
        atoms = [Mock(), Mock(), Mock(), Mock(), Mock(), Mock()]
        mock_atom.side_effect = lambda a: "a" + str(a.atom_id()) if a.atom_id() else "N"
        chains = [Mock(), Mock()]
        for index, atom in enumerate(atoms):
            atom.atom_id.return_value = index + 1 if index else None
            atom.chain.return_value = chains[index // 2] if index < 4 else None
        mock_chain.return_value = ["chain1", "chain2"]
        mock_res.return_value = ["mol1", "mol2"]
        structure.atoms.return_value = set(atoms)
        pdb_dict = structure_to_pdb_dict(structure)
        for atom in atoms:
            mock_atom.assert_any_call(atom)
        mock_chain.assert_called_with(["a2", "a3", "a4", "N"])



class AtomToAtomDictTests(TestCase):

    def setUp(self):
        self.atom = Mock()
        self.atom.atom_id.return_value = 107
        self.atom.name.return_value = "N1"
        self.atom.x.return_value = 12.681
        self.atom.y.return_value = 37.302
        self.atom.z.return_value = -25.211
        self.atom.element.return_value = "N"
        self.atom.charge.return_value = -2
        self.atom.bfactor.return_value = 12.5
        self.residue = Mock()
        self.residue.name.return_value = "GLY"
        self.residue.residue_id.return_value = "A13B"
        self.chain = Mock()
        self.chain.chain_id.return_value = "A"
        self.atom.residue.return_value = self.residue
        self.atom.chain.return_value = self.chain


    def test_can_convert_full_atom_to_dict(self):
        d = atom_to_atom_dict(self.atom)
        self.assertEqual(d, {
         "atom_id": 107, "atom_name": "N1", "alt_loc": None,
         "residue_name": "GLY",
         "chain_id": "A", "residue_id": 13, "insert_code": "B", "full_id": "A13B",
         "x": 12.681, "y": 37.302, "z": -25.211,
         "occupancy": 1.0, "temp_factor": 12.5,
         "element": "N", "charge": -2,
        })


    def test_can_convert_full_heteroatom_to_dict(self):
        self.atom.residue.return_value = None
        self.atom.chain.return_value = None
        mol = Mock()
        mol.molecule_id.return_value = "A200"
        mol.name.return_value = "SUC"
        self.atom.molecule.return_value = mol
        d = atom_to_atom_dict(self.atom)
        self.assertEqual(d, {
         "atom_id": 107, "atom_name": "N1", "alt_loc": None,
         "residue_name": "SUC",
         "chain_id": "A", "residue_id": 200, "insert_code": "", "full_id": "A200",
         "x": 12.681, "y": 37.302, "z": -25.211,
         "occupancy": 1.0, "temp_factor": 12.5,
         "element": "N", "charge": -2,
        })


    def test_can_convert_atom_with_no_insert_code(self):
        self.residue.residue_id.return_value = "A13"
        d = atom_to_atom_dict(self.atom)
        self.assertEqual(d["residue_id"], 13)
        self.assertEqual(d["insert_code"], "")
        self.assertEqual(d["full_id"], "A13")


    def test_can_convert_atom_with_no_chain_in_residue_id(self):
        self.residue.residue_id.return_value = "13"
        d = atom_to_atom_dict(self.atom)
        self.assertEqual(d["residue_id"], 13)
        self.assertEqual(d["chain_id"], "A")
        self.assertEqual(d["insert_code"], "")
        self.assertEqual(d["full_id"], "13")



class StructureToConnectionsTests(TestCase):

    def test_can_convert_structure_to_connections(self):
        atoms = [Mock(), Mock(), Mock(), Mock(), Mock(), Mock(), Mock()]
        for i, atom in enumerate(atoms):
            atom.residue.return_value = "residue" if i < 4 else None
            atom.atom_id.return_value = i + 1
        atoms[4].bonded_atoms.return_value = set(atoms[5:])
        atoms[5].bonded_atoms.return_value = set([atoms[4], atoms[6]])
        atoms[6].bonded_atoms.return_value = set(atoms[4:6])
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
