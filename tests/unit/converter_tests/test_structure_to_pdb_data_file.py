from unittest import TestCase
from unittest.mock import patch, Mock, MagicMock
from atomium.converters.structure2pdbdatafile import *
from atomium.files.pdbdatafile import PdbDataFile

class StructureToPdbDataFileTests(TestCase):

    @patch("atomium.converters.structure2pdbdatafile.atom_to_atom_dict")
    def test_can_convert_structure_to_data_file(self, mock_atom):
        mock_atom.side_effect = ["a1", "a2", "a3", "a4", "h1", "h2", "h3"]
        structure = Mock()
        atoms = [Mock(), Mock(), Mock(), Mock(), Mock(), Mock(), Mock()]
        structure.atoms.return_value = atoms
        for atom in atoms[:4]:
            atom.residue.return_value = "residue"
        for atom in atoms[4:]:
            atom.residue.return_value = None
        data_file = structure_to_pdb_data_file(structure)
        self.assertIsInstance(data_file, PdbDataFile)
        self.assertEqual(data_file.atoms, ["a1", "a2", "a3", "a4"])
        self.assertEqual(data_file.heteroatoms, ["h1", "h2", "h3"])
        for atom in atoms:
            mock_atom.assert_any_call(atom)



class AtomToAtomDictTest(TestCase):

    def setUp(self):
        self.atom = Mock()
        self.atom.atom_id.return_value = 107
        self.atom.atom_name.return_value = "N1"
        self.atom.x.return_value = 12.681
        self.atom.y.return_value = 37.302
        self.atom.z.return_value = -25.211
        self.atom.element.return_value = "N"
        self.atom.charge.return_value = -2
        self.residue = Mock()
        self.residue.name.return_value = "GLY"
        self.residue.residue_id.return_value = "A13A"
        self.chain = Mock()
        self.chain.chain_id.return_value = "A"
        self.atom.residue.return_value = self.residue
        self.atom.chain.return_value = self.chain


    def test_can_convert_full_atom(self):
        d = atom_to_atom_dict(self.atom)
        self.assertEqual(d, {
         "atom_id": 107, "atom_name": "N1", "alt_loc": None,
         "residue_name": "GLY",
         "chain_id": "A", "residue_id": 13, "insert_code": "A",
         "x": 12.681, "y": 37.302, "z": -25.211,
         "occupancy": 1.0, "temperature_factor": None,
         "element": "N", "charge": -2
        })


    def test_can_convert_atom_with_no_insert_code(self):
        self.residue.residue_id.return_value = "A13"
        d = atom_to_atom_dict(self.atom)
        self.assertEqual(d["residue_id"], 13)
        self.assertEqual(d["insert_code"], None)


    def test_can_convert_atom_with_no_chain_in_residue_id(self):
        self.residue.residue_id.return_value = "13"
        d = atom_to_atom_dict(self.atom)
        self.assertEqual(d["residue_id"], 13)
        self.assertEqual(d["chain_id"], "A")
        self.assertEqual(d["insert_code"], None)


    def test_can_convert_heteroatom(self):
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
         "chain_id": "A", "residue_id": 200, "insert_code": None,
         "x": 12.681, "y": 37.302, "z": -25.211,
         "occupancy": 1.0, "temperature_factor": None,
         "element": "N", "charge": -2
        })
