from datetime import date
from unittest import TestCase
from unittest.mock import Mock, patch, PropertyMock, MagicMock
from atomium.files.data import *

class HigherStructureGenerationTests(TestCase):

    def test_can_generate_higher_structures(self):
        model = {"atoms": [
         {"chain_id": "A", "full_res_id": "A1", "polymer": True, "residue_name": "TYR"},
         {"chain_id": "A", "full_res_id": "A1", "polymer": True, "residue_name": "TYR"},
         {"chain_id": "A", "full_res_id": "A2", "polymer": True, "residue_name": "VAL"},
         {"chain_id": "A", "full_res_id": "A2", "polymer": True, "residue_name": "VAL"},
         {"chain_id": "A", "full_res_id": "A100", "polymer": False, "residue_name": "HOH"},
         {"chain_id": "A", "full_res_id": "A100", "polymer": False, "residue_name": "HOH"},
         {"chain_id": "A", "full_res_id": "A101", "polymer": False, "residue_name": "XYZ"},
         {"chain_id": "A", "full_res_id": "A101", "polymer": False, "residue_name": "XYZ"},
         {"chain_id": "B", "full_res_id": "B1", "polymer": True, "residue_name": "TYR"},
         {"chain_id": "B", "full_res_id": "B1", "polymer": True, "residue_name": "TYR"},
         {"chain_id": "B", "full_res_id": "B2", "polymer": True, "residue_name": "VAL"},
         {"chain_id": "B", "full_res_id": "B2", "polymer": True, "residue_name": "VAL"},
         {"chain_id": "B", "full_res_id": "B100", "polymer": False, "residue_name": "HOH"},
         {"chain_id": "B", "full_res_id": "B100", "polymer": False, "residue_name": "HOH"},
         {"chain_id": "B", "full_res_id": "B101", "polymer": False, "residue_name": "XYZ"},
         {"chain_id": "B", "full_res_id": "B101", "polymer": False, "residue_name": "XYZ"},
        ], "chains": [], "residues": [], "ligands": []}
        models = [deepcopy(model), deepcopy(model)]
        generate_higher_structures(models)
        for model in models:
            self.assertEqual(model["chains"], [
             {"id": "A", "full_sequence": []}, {"id": "B", "full_sequence": []}
            ])
            self.assertEqual(model["residues"], [
             {"id": "A1", "name": "TYR", "chain_id": "A"},
             {"id": "A2", "name": "VAL", "chain_id": "A"},
             {"id": "B1", "name": "TYR", "chain_id": "B"},
             {"id": "B2", "name": "VAL", "chain_id": "B"},
            ])
            self.assertEqual(model["ligands"], [
             {"id": "A100", "name": "HOH", "chain_id": "A"},
             {"id": "A101", "name": "XYZ", "chain_id": "A"},
             {"id": "B100", "name": "HOH", "chain_id": "B"},
             {"id": "B101", "name": "XYZ", "chain_id": "B"},
            ])



class DataDictToFileConversionTests(TestCase):

    @patch("atomium.files.data.File")
    @patch("atomium.files.data.model_dict_to_model")
    def test_can_convert_data_dict_to_file_object(self, mock_model, mock_file):
        d = {
         "description": {
          "code": 1, "title": 2, "deposition_date": 3,
          "classification": 4, "keywords": 5, "authors": 6
         }, "experiment": {
          "technique": 10, "source_organism": 11, "expression_system": 12
         }, "quality": {
          "resolution": 7, "rvalue": 8, "rfree": 9
         }, "geometry": {
          "assemblies": 13
         }, "models": ["1", "2", "3"]
        }
        mock_model.side_effect = [10, 20, 30]
        f = data_dict_to_file(d)
        self.assertEqual(f._code, 1)
        self.assertEqual(f._title, 2)
        self.assertEqual(f._deposition_date, 3)
        self.assertEqual(f._classification, 4)
        self.assertEqual(f._keywords, 5)
        self.assertEqual(f._authors, 6)
        self.assertEqual(f._technique, 10)
        self.assertEqual(f._source_organism, 11)
        self.assertEqual(f._expression_system, 12)
        self.assertEqual(f._resolution, 7)
        self.assertEqual(f._rvalue, 8)
        self.assertEqual(f._rfree, 9)
        self.assertEqual(f._assemblies, 13)
        self.assertEqual(f._models, [10, 20, 30])
        mock_model.assert_any_call("1")
        mock_model.assert_any_call("2")
        mock_model.assert_any_call("3")



class ModelDictToModelTests(TestCase):

    @patch("atomium.files.data.atom_dict_to_atom")
    @patch("atomium.files.data.Model")
    @patch("atomium.files.data.create_het")
    @patch("atomium.files.data.Chain")
    @patch("atomium.files.data.bond_atoms")
    def test_can_convert_model_dict_to_model(self, mock_bond, mock_chain, mock_het, mock_mod, mock_atom):
        m = {
         "atoms": [1, 2, 3],  "connections": [],
         "residues": [
          {"id": "A1", "name": "MET", "chain_id": "A"},
          {"id": "A2", "name": "TYR", "chain_id": "A"},
          {"id": "A3", "name": "PRO", "chain_id": "A"},
          {"id": "B1", "name": "ASP", "chain_id": "B"},
          {"id": "B2", "name": "VAL", "chain_id": "B"},
          {"id": "B3", "name": "GLU", "chain_id": "B"},
         ],
         "ligands": [
          {"id": "A100", "name": "XMP", "chain_id": "A"},
          {"id": "A200", "name": "LOS", "chain_id": "A"},
          {"id": "B100", "name": "XYZ", "chain_id": "B"},
         ],
         "chains": [
          {"id": "A", "full_sequence": ["MET", "VAL"]},
          {"id": "B", "full_sequence": ["PLO", "HIS"]}
         ]
        }
        mock_atom.side_effect = lambda a: str(a)
        mod = Mock(atoms=lambda: [])
        mock_mod.return_value = mod
        residues = [Mock() for _ in range(6)]
        ligands = [Mock() for _ in range(3)]
        mock_het.side_effect = residues[:3] + ligands[:2] + residues[3:] + ligands[-1:]
        model = model_dict_to_model(m)
        for n in (1, 2, 3): mock_atom.assert_any_call(n)
        mock_mod.assert_called_with("1", "2", "3")
        for res in m["residues"]: mock_het.assert_any_call(res, [1, 2, 3], mod)
        for res in m["ligands"]: mock_het.assert_any_call(res, [1, 2, 3], mod)
        self.assertIs(residues[0].next, residues[1])
        self.assertIs(residues[1].next, residues[2])
        self.assertIs(residues[3].next, residues[4])
        self.assertIs(residues[4].next, residues[5])
        mock_chain.assert_any_call(*(residues[:3] + ligands[:2]), id="A", rep="MV")
        mock_chain.assert_any_call(*(residues[3:] + ligands[-1:]), id="B", rep="XH")



class AtomDictToAtomTests(TestCase):

    def setUp(self):
        self.atom_dict = {
         "id": 107, "name": "N1", "alt_loc": "A",
         "residue_name": "GLY",
         "chain_id": "B", "residue_id": 13, "residue_insert": "C",
         "x": 12.681, "y": 37.302, "z": -25.211,
         "occupancy": 0.5, "bfactor": 15.56,
         "element": "N", "charge": -2, "anisotropy": [1, 2, 3, 4, 5, 6]
        }

    @patch("atomium.files.data.Atom")
    def test_can_convert_atom_dict_to_atom(self, mock_atom):
        atom = Mock()
        mock_atom.return_value = atom
        returned_atom = atom_dict_to_atom(self.atom_dict)
        self.assertIs(returned_atom, atom)
        mock_atom.assert_called_with(
         "N", 12.681, 37.302, -25.211, anisotropy=[1, 2, 3, 4, 5, 6],
         id=107, name="N1", charge=-2, bfactor=15.56
        )


    @patch("atomium.files.data.Atom")
    def test_can_convert_atom_dict_to_atom_with_no_temp_factor(self, mock_atom):
        atom = Mock()
        mock_atom.return_value = atom
        self.atom_dict["bfactor"] = None
        returned_atom = atom_dict_to_atom(self.atom_dict)
        mock_atom.assert_called_with(
         "N", 12.681, 37.302, -25.211, anisotropy=[1, 2, 3, 4, 5, 6],
         id=107, name="N1", charge=-2, bfactor=0
        )



class HetCreationTests(TestCase):

    @patch("atomium.files.data.Residue")
    def test_can_create_residue(self, mock_res):
        res = {"id": "A1", "name": "TYR"}
        model = Mock()
        model.atom.side_effect = [1, 2, 3]
        atoms = [
         {"full_res_id": "A0", "polymer": True, "occupancy": 1, "alt_loc": None, "id": 1},
         {"full_res_id": "A1", "polymer": True, "occupancy": 1, "alt_loc": None, "id": 2},
         {"full_res_id": "A1", "polymer": True, "occupancy": 1, "alt_loc": None, "id": 3},
         {"full_res_id": "A1", "polymer": True, "occupancy": 1, "alt_loc": None, "id": 4},
         {"full_res_id": "A2", "polymer": True, "occupancy": 1, "alt_loc": None, "id": 5},
        ]
        het = create_het(res, atoms, model)
        for n in (2, 3, 4): model.atom.assert_any_call(n)
        mock_res.assert_called_with(1, 2, 3, id="A1", name="TYR")
        self.assertIs(het, mock_res.return_value)


    @patch("atomium.files.data.Ligand")
    def test_can_create_ligand(self, mock_lig):
        res = {"id": "A1", "name": "TYR"}
        model = Mock()
        model.atom.side_effect = [1, 2, 3]
        atoms = [
         {"full_res_id": "A0", "polymer": True, "occupancy": 1, "alt_loc": None, "id": 1},
         {"full_res_id": "A1", "polymer": False, "occupancy": 1, "alt_loc": None, "id": 2},
         {"full_res_id": "A1", "polymer": False, "occupancy": 1, "alt_loc": None, "id": 3},
         {"full_res_id": "A1", "polymer": False, "occupancy": 1, "alt_loc": None, "id": 4},
         {"full_res_id": "A2", "polymer": True, "occupancy": 1, "alt_loc": None, "id": 5},
        ]
        het = create_het(res, atoms, model)
        for n in (2, 3, 4): model.atom.assert_any_call(n)
        mock_lig.assert_called_with(1, 2, 3, id="A1", name="TYR")
        self.assertIs(het, mock_lig.return_value)


    @patch("atomium.files.data.Residue")
    def test_can_create_residue_with_low_occupancy_atoms(self, mock_res):
        res = {"id": "A1", "name": "TYR"}
        model = Mock()
        model.atom.side_effect = [1, 2, 3, 4]
        atoms = [
         {"full_res_id": "A0", "polymer": True, "occupancy": 1, "alt_loc": None, "id": 1},
         {"full_res_id": "A1", "polymer": True, "occupancy": 1, "alt_loc": None, "id": 2},
         {"full_res_id": "A1", "polymer": True, "occupancy": 1, "alt_loc": "A", "id": 3},
         {"full_res_id": "A1", "polymer": True, "occupancy": 0.7, "alt_loc": "A", "id": 4},
         {"full_res_id": "A1", "polymer": True, "occupancy": 0.3, "alt_loc": "B", "id": 5},
         {"full_res_id": "A1", "polymer": True, "occupancy": 0.4, "alt_loc": "A", "id": 6},
         {"full_res_id": "A1", "polymer": True, "occupancy": 0.6, "alt_loc": "B", "id": 7},
         {"full_res_id": "A2", "polymer": True, "occupancy": 1, "alt_loc": None, "id": 8},
        ]
        het = create_het(res, atoms, model)
        for n in (2, 3, 4, 6): model.atom.assert_any_call(n)
        mock_res.assert_called_with(1, 2, 3, 4, id="A1", name="TYR")
        self.assertIs(het, mock_res.return_value)



class AtomBondingTests(TestCase):

    @patch("atomium.files.data.make_intra_residue_bonds")
    @patch("atomium.files.data.make_inter_residue_bonds")
    @patch("atomium.files.data.make_connections_bonds")
    def test_can_bond_atoms(self, mock_con, mock_inter, mock_intra):
        residues = [Mock(), Mock()]
        model = Mock()
        model.residues.return_value = set(residues)
        connections = "ccc"
        bond_atoms(model, connections)
        mock_intra.assert_called_with(set(residues), BONDS)
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
