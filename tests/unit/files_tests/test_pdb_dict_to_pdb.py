from unittest import TestCase
from unittest.mock import patch, Mock, MagicMock
from atomium.files.pdbdict2pdb import *
from atomium.structures.reference import bonds

class PdbDictToPdbTests(TestCase):

    @patch("atomium.files.pdbdict2pdb.Pdb")
    @patch("atomium.files.pdbdict2pdb.model_dict_to_model")
    def test_can_convert_pdb_dict_to_pdb(self, mock_model, mock_pdb):
        pdb = Mock()
        mock_pdb.return_value = pdb
        mock_model.side_effect = ["model1", "model2", "model3"]
        pdb_dict = {
         "deposition_date": "D", "code": "C", "title": "T", "resolution": 1.4,
         "organism": "H. sap", "expression_system": "M. mus",
         "technique": "TECHNIQUE", "classification": "CLASS", "rfactor": 4.5,
         "keywords": ["A", "B"],
         "models": ["1", "2", "3"],
         "connections": ["c1", "c2"]
        }
        returned_pdb = pdb_dict_to_pdb(pdb_dict)
        mock_model.assert_any_call("1", ["c1", "c2"])
        mock_model.assert_any_call("2", ["c1", "c2"])
        mock_model.assert_any_call("3", ["c1", "c2"])
        self.assertIs(returned_pdb, pdb)
        self.assertEqual(returned_pdb._deposition_date, "D")
        self.assertEqual(returned_pdb._code, "C")
        self.assertEqual(returned_pdb._title, "T")
        self.assertEqual(returned_pdb._resolution, 1.4)
        self.assertEqual(returned_pdb._organism, "H. sap")
        self.assertEqual(returned_pdb._expression_system, "M. mus")
        self.assertEqual(returned_pdb._technique, "TECHNIQUE")
        self.assertEqual(returned_pdb._classification, "CLASS")
        self.assertEqual(returned_pdb._rfactor, 4.5)
        self.assertEqual(returned_pdb._keywords, ["A", "B"])
        self.assertEqual(returned_pdb._models, ["model1", "model2", "model3"])



class ModelDictToModelTests(TestCase):

    @patch("atomium.files.pdbdict2pdb.Model")
    @patch("atomium.files.pdbdict2pdb.chain_dict_to_chain")
    @patch("atomium.files.pdbdict2pdb.residue_dict_to_residue")
    @patch("atomium.files.pdbdict2pdb.bond_atoms")
    def test_can_convert_model_dict_to_model(self, mock_bond, mock_res, mock_chain, mock_model):
        model = Mock()
        mock_model.return_value = model
        mock_chain.side_effect = ["chain1", "chain2"]
        mock_res.side_effect = ["mol1", "mol2", "mol3"]
        model_dict = {
         "molecules": ["m1", "m2", "m3"], "chains": ["c1", "c2"]
        }
        returned_model = model_dict_to_model(model_dict, ["c1", "c2"])
        mock_chain.assert_any_call("c1")
        mock_chain.assert_any_call("c2")
        mock_res.assert_any_call("m1", molecule=True)
        mock_res.assert_any_call("m2", molecule=True)
        mock_res.assert_any_call("m3", molecule=True)
        self.assertIs(returned_model, model)
        model.add.assert_any_call("chain1")
        model.add.assert_any_call("chain2")
        model.add.assert_any_call("mol1")
        model.add.assert_any_call("mol2")
        model.add.assert_any_call("mol3")
        mock_bond.assert_called_with(model, ["c1", "c2"])



class ChainDictToChainTests(TestCase):

    @patch("atomium.files.pdbdict2pdb.Chain")
    @patch("atomium.files.pdbdict2pdb.residue_dict_to_residue")
    def test_can_convert_chain_dict_to_chain(self, mock_res, mock_chain):
        chain = Mock()
        mock_chain.return_value = chain
        residues = [Mock(), Mock(), Mock()]
        mock_res.side_effect = residues
        chain_dict = {"chain_id": "A", "residues": ["r1", "r2", "r3"]}
        returned_chain = chain_dict_to_chain(chain_dict)
        self.assertIs(returned_chain, chain)
        mock_res.assert_any_call("r1")
        mock_res.assert_any_call("r2")
        mock_res.assert_any_call("r3")
        self.assertIs(residues[0].next, residues[1])
        self.assertIs(residues[1].next, residues[2])
        mock_chain.assert_called_with(*residues, id="A")



class MoleculeDictToMoleculeTests(TestCase):

    def setUp(self):
        self.atom_objects = [Mock(Atom), Mock(Atom), Mock(Atom), Mock(Atom)]
        for atom in self.atom_objects:
            atom._residue = None
        self.atom_objects[0].name = "C"
        self.atom_objects[1].name = "N"
        self.atom_objects[2].name = "CB"
        self.atom_objects[3].name = "CC"
        self.atom_objects[0].id = 1
        self.atom_objects[1].id = 2
        self.atom_objects[2].id = 3
        self.atom_objects[3].id = 4
        self.atom_dicts = [{
         "alt_loc": None, "atom_id": n, "occupancy": 1
        } for n in range(1, 5)]


    @patch("atomium.files.pdbdict2pdb.Residue")
    @patch("atomium.files.pdbdict2pdb.atom_dict_to_atom")
    def test_can_convert_residue_dict_to_residue(self, mock_atom, mock_res):
        mock_atom.side_effect = self.atom_objects
        residue = Mock()
        mock_res.return_value = residue
        res_dict = {"id": "A12", "name": "VAL", "atoms": self.atom_dicts}
        returned_residue = residue_dict_to_residue(res_dict)
        mock_atom.assert_any_call({"alt_loc": None, "atom_id": 1, "occupancy": 1})
        mock_atom.assert_any_call({"alt_loc": None, "atom_id": 2, "occupancy": 1})
        mock_atom.assert_any_call({"alt_loc": None, "atom_id": 3, "occupancy": 1})
        mock_atom.assert_any_call({"alt_loc": None, "atom_id": 4, "occupancy": 1})
        self.assertEqual(mock_atom.call_count, 4)
        self.assertIs(returned_residue, residue)
        self.assertEqual(residue._id, "A12")
        mock_res.assert_called_with(*self.atom_objects, name="VAL")


    @patch("atomium.files.pdbdict2pdb.Residue")
    @patch("atomium.files.pdbdict2pdb.atom_dict_to_atom")
    def test_can_convert_residue_dict_to_residue_alt_loc(self, mock_atom, mock_res):
        self.atom_dicts[2]["alt_loc"] = "A"
        self.atom_dicts[3]["alt_loc"] = "B"
        self.atom_dicts[2]["occupancy"] = 0.8
        self.atom_dicts[3]["occupancy"] = 0.2
        mock_atom.side_effect = self.atom_objects
        residue = Mock()
        mock_res.return_value = residue
        res_dict = {"id": "A12", "name": "VAL", "atoms": self.atom_dicts}
        returned_residue = residue_dict_to_residue(res_dict)
        mock_atom.assert_any_call({"alt_loc": None, "atom_id": 1, "occupancy": 1})
        mock_atom.assert_any_call({"alt_loc": None, "atom_id": 2, "occupancy": 1})
        mock_atom.assert_any_call({"alt_loc": "A", "atom_id": 3, "occupancy": 0.8})
        self.assertEqual(mock_atom.call_count, 3)
        self.assertIs(returned_residue, residue)
        self.assertEqual(residue._id, "A12")
        mock_res.assert_called_with(*self.atom_objects[:-1], name="VAL")


    @patch("atomium.files.pdbdict2pdb.Molecule")
    @patch("atomium.files.pdbdict2pdb.atom_dict_to_atom")
    def test_can_convert_molecule_dict_to_molecule(self, mock_atom, mock_mol):
        mock_atom.side_effect = self.atom_objects
        molecule = Mock()
        mock_mol.return_value = molecule
        mol_dict = {"id": "A500", "name": "XMP", "atoms": self.atom_dicts}
        returned_molecule = residue_dict_to_residue(mol_dict, molecule=True)
        mock_atom.assert_any_call({"alt_loc": None, "atom_id": 1, "occupancy": 1})
        mock_atom.assert_any_call({"alt_loc": None, "atom_id": 2, "occupancy": 1})
        mock_atom.assert_any_call({"alt_loc": None, "atom_id": 3, "occupancy": 1})
        mock_atom.assert_any_call({"alt_loc": None, "atom_id": 4, "occupancy": 1})
        self.assertEqual(mock_atom.call_count, 4)
        self.assertIs(returned_molecule, molecule)
        self.assertEqual(molecule._id, "A500")
        mock_mol.assert_called_with(*self.atom_objects, name="XMP")



class AtomDictToAtomTests(TestCase):

    def setUp(self):
        self.atom_dict = {
         "atom_id": 107, "atom_name": "N1", "alt_loc": "A",
         "residue_name": "GLY",
         "chain_id": "B", "residue_id": 13, "insert_code": "C",
         "x": 12.681, "y": 37.302, "z": -25.211,
         "occupancy": 0.5, "temp_factor": 15.56,
         "element": "N", "charge": -2
        }

    @patch("atomium.files.pdbdict2pdb.Atom")
    def test_can_convert_atom_dict_to_atom(self, mock_atom):
        atom = Mock()
        mock_atom.return_value = atom
        returned_atom = atom_dict_to_atom(self.atom_dict)
        self.assertIs(returned_atom, atom)
        mock_atom.assert_called_with(
         "N", 12.681, 37.302, -25.211,
         id=107, name="N1", charge=-2, bfactor=15.56
        )


    @patch("atomium.files.pdbdict2pdb.Atom")
    def test_can_convert_atom_dict_to_atom_with_no_temp_factor(self, mock_atom):
        atom = Mock()
        mock_atom.return_value = atom
        self.atom_dict["temp_factor"] = None
        returned_atom = atom_dict_to_atom(self.atom_dict)
        mock_atom.assert_called_with(
         "N", 12.681, 37.302, -25.211,
         id=107, name="N1", charge=-2, bfactor=0
        )



class AtomBondingTests(TestCase):

    @patch("atomium.files.pdbdict2pdb.make_intra_residue_bonds")
    @patch("atomium.files.pdbdict2pdb.make_inter_residue_bonds")
    @patch("atomium.files.pdbdict2pdb.make_connections_bonds")
    def test_can_bond_atoms(self, mock_con, mock_inter, mock_intra):
        residues = [Mock(), Mock()]
        model = Mock()
        model.residues.return_value = set(residues)
        connections = "ccc"
        bond_atoms(model, connections)
        mock_intra.assert_called_with(set(residues), bonds)
        mock_inter.assert_called_with(set(residues))
        mock_con.assert_called_with(model, "ccc")



class IntraResidueConnectionTests(TestCase):

    def test_can_connect_residue(self):
        residues = [Mock(), Mock(), Mock(), Mock()]
        atoms = [Mock(), Mock(), Mock(), Mock(), Mock(), Mock()]
        atoms += [Mock(), Mock(), Mock(), Mock(), Mock(), Mock()]
        for i, residue in enumerate(residues):
            residue.atoms.return_value = set(atoms[i * 3: i * 3 + 3])
            residue.name = ["CYS", "TYR", "MET", "VAL"][i]
            atoms[i * 3].name = "A"
            atoms[i * 3 + 1].name = "B"
            atoms[i * 3 + 2].name = ["C", "P", "U", "D"][i]
        d = {
         "CYS": {"A": ["B"], "B": ["A"]},
         "TYR": {"A": ["B", "P"], "B": ["A"], "P": ["A"]},
         "MET": {"A": ["B", "X"], "X": ["A"], "B": ["A"]}
        }
        make_intra_residue_bonds(residues, d)
        atoms[0].bond_to.assert_called_with(atoms[1])
        atoms[1].bond_to.assert_called_with(atoms[0])
        atoms[3].bond_to.assert_any_call(atoms[4])
        atoms[3].bond_to.assert_any_call(atoms[5])
        atoms[4].bond_to.assert_called_with(atoms[3])
        atoms[5].bond_to.assert_called_with(atoms[3])
        atoms[6].bond_to.assert_called_with(atoms[7])
        atoms[7].bond_to.assert_called_with(atoms[6])
        self.assertFalse(atoms[2].bond_to.called)
        self.assertFalse(atoms[8].bond_to.called)
        self.assertFalse(atoms[9].bond_to.called)
        self.assertFalse(atoms[10].bond_to.called)
        self.assertFalse(atoms[11].bond_to.called)



class InterResidueConnectionTests(TestCase):

    def test_can_connect_residues(self):
        residues = [Mock(), Mock(), Mock(), Mock()]
        atoms = [Mock(), Mock(), Mock(), Mock(), Mock(), Mock()]
        atoms += [Mock(), Mock(), Mock(), Mock(), Mock(), Mock()]
        for atom in atoms:
            atom.distance_to.return_value = 4.9
        def get_atom1(name=None):
            return atoms[0] if name == "N" else atoms[2]
        def get_atom2(name=None):
            return atoms[3] if name == "N" else atoms[5]
        def get_atom3(name=None):
            return atoms[6] if name == "N" else atoms[8]
        def get_atom4(name=None):
            return atoms[9] if name == "N" else atoms[11]
        residues[0].atom.side_effect = get_atom1
        residues[1].atom.side_effect = get_atom2
        residues[2].atom.side_effect = get_atom3
        residues[3].atom.side_effect = get_atom4
        for i, residue in enumerate(residues):
            residue.next = residues[i + 1] if i != 3 else None
        make_inter_residue_bonds(residues)
        atoms[2].bond_to.assert_called_with(atoms[3])
        atoms[5].bond_to.assert_called_with(atoms[6])
        atoms[8].bond_to.assert_called_with(atoms[9])
        self.assertFalse(atoms[0].bond_to.called)
        self.assertFalse(atoms[1].bond_to.called)
        self.assertFalse(atoms[3].bond_to.called)
        self.assertFalse(atoms[4].bond_to.called)
        self.assertFalse(atoms[6].bond_to.called)
        self.assertFalse(atoms[7].bond_to.called)
        self.assertFalse(atoms[9].bond_to.called)
        self.assertFalse(atoms[10].bond_to.called)
        self.assertFalse(atoms[11].bond_to.called)


    def test_can_skip_bond_where_atom_not_present(self):
        residue1, residue2 = Mock(), Mock()
        residue1.next = residue2
        residue2.next = None
        atom = Mock()
        atom.distance_to.return_value = 4.9
        residue1.atom.return_value = atom
        residue2.atom.return_value = None
        make_inter_residue_bonds([residue1, residue2])
        self.assertFalse(atom.bond_to.called)
        residue2.atom.return_value = atom
        residue1.atom.return_value = None
        make_inter_residue_bonds([residue1, residue2])
        self.assertFalse(atom.bond_to.called)


    def test_can_skip_bond_where_distance_too_great(self):
        residue1, residue2 = Mock(), Mock()
        residue1.next.return_value = residue2
        residue2.next.return_value = None
        atom = Mock()
        atom.bond = MagicMock()
        atom.distance_to.return_value = 5.1
        residue1.atom.return_value = atom
        residue2.atom.return_value = atom
        make_inter_residue_bonds([residue1, residue2])
        self.assertFalse(atom.bond.called)



class ConnectionBondsTests(TestCase):

    def test_can_bond_from_connections(self):
        model = Mock()
        atoms = [Mock(), Mock(), Mock(), Mock()]
        model.atom.side_effect = [
         atoms[1], atoms[2], atoms[3], None, atoms[2], atoms[1], atoms[3], atoms[1], None
        ]
        for index, atom in enumerate(atoms):
            atom.id = index
        connections = [{
         "atom": 1, "bond_to": [2, 3, 4]
        }, {
         "atom": 2, "bond_to": [1]
        }, {
         "atom": 3, "bond_to": [1]
        }, {
         "atom": 4, "bond_to": [1]
        }]
        make_connections_bonds(model, connections)
        atoms[1].bond_to.assert_any_call(atoms[2])
        atoms[1].bond_to.assert_any_call(atoms[3])
        atoms[2].bond_to.assert_called_with(atoms[1])
        atoms[3].bond_to.assert_called_with(atoms[1])
        self.assertFalse(atoms[0].called)


    def test_can_ignore_self_bonding(self):
        model = Mock()
        atoms = [Mock(), Mock(), Mock(), Mock()]
        model.atom.side_effect = [
         atoms[0], atoms[1], atoms[2], atoms[3]
        ]
        for index, atom in enumerate(atoms):
            atom.id = index + 1
        connections = [{
         "atom": 1, "bond_to": [1, 2, 3, 4]
        }]
        make_connections_bonds(model, connections)
        atoms[0].bond_to.assert_any_call(atoms[1])
        atoms[0].bond_to.assert_any_call(atoms[2])
        atoms[0].bond_to.assert_called_with(atoms[3])
