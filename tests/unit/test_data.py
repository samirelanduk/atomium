from datetime import date
from unittest import TestCase
from unittest.mock import Mock, patch, PropertyMock, MagicMock
from atomium.data import *

class DataDictToFileTests(TestCase):

    @patch("atomium.data.File")
    @patch("atomium.data.model_dict_to_model")
    def test_can_convert_data_dict_to_file(self, mock_mod, mock_file):
        d = {"K1": {"A": 1, "B": 2}, "K2": {"C": 3}, "models": [1, 2]}
        mock_mod.side_effect = "XY"
        f = data_dict_to_file(d, "abc")
        mock_file.assert_called_with("abc")
        mock_mod.assert_any_call(1)
        mock_mod.assert_any_call(2)
        self.assertEqual(f, mock_file.return_value)
        self.assertEqual(f._A, 1)
        self.assertEqual(f._B, 2)
        self.assertEqual(f._C, 3)
        self.assertEqual(f._models, ["X", "Y"])



class ModelDictToModelTests(TestCase):

    @patch("atomium.data.create_chains")
    @patch("atomium.data.create_ligands")
    @patch("atomium.data.Model")
    def test_can_convert_model_dict_to_model(self, mock_mod, mock_lig, mock_ch):
        m = {1: 2}
        mock_ch.return_value = [1, 2]
        mock_lig.side_effect = [[3, 4], [5, 6]]
        model = model_dict_to_model(m)
        mock_ch.assert_called_with(m)
        mock_lig.assert_any_call(m, [1, 2])
        mock_lig.assert_any_call(m, [1, 2], water=True)
        mock_mod.assert_called_with(1, 2, 3, 4, 5, 6)
        self.assertEqual(model, mock_mod.return_value)



class ChainCreationTests(TestCase):

    @patch("atomium.data.create_het")
    @patch("atomium.data.Chain")
    def test_can_create_chains(self, mock_chain, mock_het):
        residues = [Mock() for _ in range(6)]
        mock_het.side_effect = residues
        chains = create_chains({
         "polymer": {
          "A": {"residues": {1: {"number": 4}, 3: {"number": 2}}, "internal_id": "1", "sequence": "GH"},
          "B": {"residues": {5: {"number": 6}, 7: {"number": 8}}, "internal_id": "2", "sequence": "AL"},
         }
        })
        mock_het.assert_any_call({"number": 4}, 1)
        mock_het.assert_any_call({"number": 2}, 3)
        mock_het.assert_any_call({"number": 6}, 5)
        mock_het.assert_any_call({"number": 8}, 7)
        self.assertIs(residues[0]._next, residues[1])
        self.assertIs(residues[1]._previous, residues[0])
        self.assertIs(residues[2]._next, residues[3])
        self.assertIs(residues[3]._previous, residues[2])
        mock_chain.assert_any_call(*residues[:2], id="A", internal_id="1", sequence="GH")
        mock_chain.assert_any_call(*residues[2:4], id="B", internal_id="2", sequence="AL")
        self.assertEqual(chains, [mock_chain.return_value, mock_chain.return_value])



class LigandCreationTests(TestCase):

    @patch("atomium.data.create_het")
    def test_can_create_ligands(self, mock_het):
        chains = [Mock(_id="1"), Mock(_id="2"), Mock(_id="3"), Mock(_id="4")]
        ligands = create_ligands({
         "non-polymer": {"A": {"polymer": "1"}, "B": {"polymer": "2"}},
         "water": {"C": {"polymer": "3"}, "D": {"polymer": "4"}}
        }, chains)
        mock_het.assert_any_call({"polymer": "1"}, "A", ligand=True, chain=chains[0], water=False)
        mock_het.assert_any_call({"polymer": "2"}, "B", ligand=True, chain=chains[1], water=False)
        self.assertEqual(ligands, [mock_het.return_value, mock_het.return_value])


    @patch("atomium.data.create_het")
    def test_can_create_water_ligand(self, mock_het):
        chains = [Mock(_id="1"), Mock(_id="2"), Mock(_id="3"), Mock(_id="4")]
        ligands = create_ligands({
         "non-polymer": {"A": {"polymer": "1"}, "B": {"polymer": "2"}},
         "water": {"C": {"polymer": "3"}, "D": {"polymer": "4"}}
        }, chains, water=True)
        mock_het.assert_any_call({"polymer": "3"}, "C", ligand=True, chain=chains[2], water=True)
        mock_het.assert_any_call({"polymer": "4"}, "D", ligand=True, chain=chains[3], water=True)
        self.assertEqual(ligands, [mock_het.return_value, mock_het.return_value])



class HetCreationTests(TestCase):

    @patch("atomium.data.atom_dict_to_atom")
    @patch("atomium.data.Residue")
    def test_can_create_residue(self, mock_res, mock_atom):
        mock_atom.side_effect = "IJ"
        create_het({"atoms": {
         1: {"occupancy": 1, "name": "C"}, 3: {"occupancy": 1, "name": "D"}
        }, "name": "VAL"}, "A")
        mock_atom.assert_any_call({"occupancy": 1, "name": "C"}, 1)
        mock_atom.assert_any_call({"occupancy": 1, "name": "D"}, 3)
        mock_res.assert_called_with("I", "J", id="A", name="VAL")


    @patch("atomium.data.atom_dict_to_atom")
    @patch("atomium.data.Ligand")
    def test_can_create_ligand(self, mock_lig, mock_atom):
        mock_atom.side_effect = "IJ"
        create_het({"atoms": {
         1: {"occupancy": 1, "name": "C"}, 3: {"occupancy": 1, "name": "D"}
        }, "name": "XMP", "internal_id": "Q"}, "A", ligand=True, chain=2, water=True)
        mock_atom.assert_any_call({"occupancy": 1, "name": "C"}, 1)
        mock_atom.assert_any_call({"occupancy": 1, "name": "D"}, 3)
        mock_lig.assert_called_with("I", "J", id="A", name="XMP", chain=2, water=True, internal_id="Q")


    @patch("atomium.data.atom_dict_to_atom")
    @patch("atomium.data.Residue")
    def test_can_create_residue_with_some_atoms(self, mock_res, mock_atom):
        mock_atom.side_effect = "IJK"
        create_het({"atoms": {
         1: {"occupancy": 1, "name": "C", "alt_loc": None},
         3: {"occupancy": 1, "name": "D", "alt_loc": None},
         10: {"occupancy": 0.8, "name": "C", "alt_loc": "A"},
         13: {"occupancy": 0.2, "name": "D", "alt_loc": "B"}
        }, "name": "VAL"}, "A")
        mock_atom.assert_any_call({"occupancy": 1, "name": "C", "alt_loc": None}, 1)
        mock_atom.assert_any_call({"occupancy": 1, "name": "D", "alt_loc": None}, 3)
        mock_atom.assert_any_call({"occupancy": 0.8, "name": "C", "alt_loc": "A"}, 10)
        mock_res.assert_called_with("I", "J", "K", id="A", name="VAL")



class AtomCreationTests(TestCase):

    @patch("atomium.data.Atom")
    def test_can_make_atom(self, mock_atom):
        a = atom_dict_to_atom({
         "element": 1, "x": 2, "y": 3, "z": 4, "name": 5,
         "charge": 6, "bvalue": 7, "anisotropy": 8
        }, 100)
        mock_atom.assert_called_with(1, 2, 3, 4, 100, 5, 6, 7, 8)
        self.assertIs(a, mock_atom.return_value)
