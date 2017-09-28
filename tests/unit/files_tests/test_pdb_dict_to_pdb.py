from unittest import TestCase
from unittest.mock import patch, Mock, MagicMock
from atomium.files.pdbdict2pdb import *

class PdbDictToPdbTests(TestCase):

    @patch("atomium.files.pdbdict2pdb.Pdb")
    @patch("atomium.files.pdbdict2pdb.model_dict_to_model")
    def test_can_convert_pdb_dict_to_pdb(self, mock_model, mock_pdb):
        pdb = Mock()
        mock_pdb.return_value = pdb
        mock_model.side_effect = ["model1", "model2", "model3"]
        pdb_dict = {
         "deposition_date": "D", "code": "C", "title": "T",
         "models": ["1", "2", "3"]
        }
        returned_pdb = pdb_dict_to_pdb(pdb_dict)
        self.assertIs(returned_pdb, pdb)
        self.assertEqual(returned_pdb._deposition_date, "D")
        self.assertEqual(returned_pdb._code, "C")
        self.assertEqual(returned_pdb._title, "T")
        self.assertEqual(returned_pdb._models, ["model1", "model2", "model3"])



class ModelDictToModelTests(TestCase):

    @patch("atomium.files.pdbdict2pdb.Model")
    @patch("atomium.files.pdbdict2pdb.chain_dict_to_chain")
    @patch("atomium.files.pdbdict2pdb.residue_dict_to_residue")
    def test_can_convert_model_dict_to_model(self, mock_res, mock_chain, mock_model):
        model = Mock()
        mock_model.return_value = model
        mock_chain.side_effect = ["chain1", "chain2"]
        mock_res.side_effect = ["mol1", "mol2", "mol3"]
        model_dict = {
         "molecules": ["m1", "m2", "m3"], "chains": ["c1", "c2"]
        }
        returned_model = model_dict_to_model(model_dict)
        mock_chain.assert_any_call("c1")
        mock_chain.assert_any_call("c2")
        mock_res.assert_any_call("m1", molecule=True)
        mock_res.assert_any_call("m2", molecule=True)
        mock_res.assert_any_call("m3", molecule=True)
        self.assertIs(returned_model, model)
        model.add_chain.assert_any_call("chain1")
        model.add_chain.assert_any_call("chain2")
        model.add_molecule.assert_any_call("mol1")
        model.add_molecule.assert_any_call("mol2")
        model.add_molecule.assert_any_call("mol3")



class ChainDictToChainTests(TestCase):

    @patch("atomium.files.pdbdict2pdb.Chain")
    @patch("atomium.files.pdbdict2pdb.residue_dict_to_residue")
    def test_can_convert_chain_dict_to_chain(self, mock_res, mock_chain):
        chain = Mock()
        mock_chain.return_value = chain
        mock_res.side_effect = ["residue1", "residue2"]
        chain_dict = {"chain_id": "A", "residues": ["r1", "r2"]}
        returned_chain = chain_dict_to_chain(chain_dict)
        self.assertIs(returned_chain, chain)
        mock_res.assert_any_call("r1")
        mock_res.assert_any_call("r2")
        mock_chain.assert_called_with("residue1", "residue2", chain_id="A")



class MoleculeDictToMoleculeTests(TestCase):

    @patch("atomium.files.pdbdict2pdb.Residue")
    @patch("atomium.files.pdbdict2pdb.atom_dict_to_atom")
    def test_can_convert_residue_dict_to_residue(self, mock_atom, mock_res):
        mock_atom.side_effect = ["atom1", "atom2"]
        residue = Mock()
        mock_res.return_value = residue
        res_dict = {"id": "A12", "name": "VAL", "atoms": ["a1", "a2"]}
        returned_residue = residue_dict_to_residue(res_dict)
        mock_atom.assert_any_call("a1")
        mock_atom.assert_any_call("a2")
        self.assertIs(returned_residue, residue)
        self.assertEqual(residue._id, "A12")
        mock_res.assert_called_with("atom1", "atom2", name="VAL")


    @patch("atomium.files.pdbdict2pdb.Molecule")
    @patch("atomium.files.pdbdict2pdb.atom_dict_to_atom")
    def test_can_convert_molecule_dict_to_molecule(self, mock_atom, mock_mol):
        mock_atom.side_effect = ["atom1", "atom2"]
        molecule = Mock()
        mock_mol.return_value = molecule
        mol_dict = {"id": "A500", "name": "XMP", "atoms": ["a1", "a2"]}
        returned_molecule = residue_dict_to_residue(mol_dict, molecule=True)
        mock_atom.assert_any_call("a1")
        mock_atom.assert_any_call("a2")
        self.assertIs(returned_molecule, molecule)
        self.assertEqual(molecule._id, "A500")
        mock_mol.assert_called_with("atom1", "atom2", name="XMP")



class AtomDictToAtomTests(TestCase):

    @patch("atomium.files.pdbdict2pdb.Atom")
    def test_can_convert_atom_dict_to_atom(self, mock_atom):
        atom = Mock()
        mock_atom.return_value = atom
        atom_dict = {
         "atom_id": 107, "atom_name": "N1", "alt_loc": "A",
         "residue_name": "GLY",
         "chain_id": "B", "residue_id": 13, "insert_code": "C",
         "x": 12.681, "y": 37.302, "z": -25.211,
         "occupancy": 0.5, "temp_factor": 15.56,
         "element": "N", "charge": -2
        }
        returned_atom = atom_dict_to_atom(atom_dict)
        self.assertIs(returned_atom, atom)
        mock_atom.assert_called_with(
         "N", 12.681, 37.302, -25.211,
         atom_id=107, name="N1", charge=-2, bfactor=15.56
        )
