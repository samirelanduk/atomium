from unittest import TestCase
from unittest.mock import patch, Mock, MagicMock
from atomium.converters.structure2pdbdatafile import *
from atomium.files.pdbdatafile import PdbDataFile

class StructureToPdbDataFileTests(TestCase):

    @patch("atomium.converters.structure2pdbdatafile.atom_to_atom_dict")
    @patch("atomium.converters.structure2pdbdatafile.atoms_to_connections")
    def test_can_convert_structure_to_data_file(self, mock_con, mock_atom):
        mock_atom.side_effect = ["a1", "a2", "a4", "a3", "h2", "h1", "h3"]
        mock_atom.side_effect = [{"atom_id": a} for a in mock_atom.side_effect]
        mock_con.return_value = ["c1", "c2", "c3"]
        structure = Mock()
        atoms = [Mock(), Mock(), Mock(), Mock(), Mock(), Mock(), Mock()]
        structure.atoms.return_value = atoms
        for atom in atoms[:4]:
            atom.residue.return_value = "residue"
        for atom in atoms[4:]:
            atom.residue.return_value = None
        data_file = structure_to_pdb_data_file(structure, model=3)
        self.assertIsInstance(data_file, PdbDataFile)
        self.assertEqual(data_file.atoms, [
         {"atom_id": "a1"}, {"atom_id": "a2"}, {"atom_id": "a3"}, {"atom_id": "a4"}
        ])
        self.assertEqual(data_file.heteroatoms, [
         {"atom_id": "h1"}, {"atom_id": "h2"}, {"atom_id": "h3"}
        ])
        self.assertEqual(data_file.connections, ["c1", "c2", "c3"])
        for atom in atoms:
            mock_atom.assert_any_call(atom, model=3)
        mock_con.assert_called_with(atoms)



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
        self.residue.residue_id.return_value = "A13A"
        self.chain = Mock()
        self.chain.chain_id.return_value = "A"
        self.atom.residue.return_value = self.residue
        self.atom.chain.return_value = self.chain


    def test_can_convert_full_atom(self):
        d = atom_to_atom_dict(self.atom, model=2)
        self.assertEqual(d, {
         "atom_id": 107, "atom_name": "N1", "alt_loc": None,
         "residue_name": "GLY",
         "chain_id": "A", "residue_id": 13, "insert_code": "A",
         "x": 12.681, "y": 37.302, "z": -25.211,
         "occupancy": 1.0, "temperature_factor": 12.5,
         "element": "N", "charge": -2, "model": 2
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
         "occupancy": 1.0, "temperature_factor": 12.5,
         "element": "N", "charge": -2, "model": 1
        })



class AtomsToConnectionsTests(TestCase):

    def test_can_convert_atoms_to_connections(self):
        atoms = [Mock(), Mock(), Mock(), Mock(), Mock(), Mock(), Mock()]
        for i, atom in enumerate(atoms):
            atom.residue.return_value = "residue" if i < 4 else None
            atom.atom_id.return_value = i + 1
        atoms[4].bonded_atoms.return_value = set(atoms[5:])
        atoms[5].bonded_atoms.return_value = set([atoms[4], atoms[6]])
        atoms[6].bonded_atoms.return_value = set(atoms[4:6])
        connections = atoms_to_connections(set(atoms))
        self.assertEqual(connections, [{
         "atom": 5, "bond_to": [6, 7]
        }, {
         "atom": 6, "bond_to": [5, 7]
        }, {
         "atom": 7, "bond_to": [5, 6]
        }])
