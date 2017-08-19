from unittest import TestCase
from unittest.mock import patch, Mock, MagicMock
from atomium.converters.pdbdatafile2model import *
from atomium.parse.pdbdatafile import PdbDataFile
from atomium.structures.models import Model
from atomium.structures.molecules import Residue
from atomium.structures.atoms import Atom

class PdbDataFileToModelTests(TestCase):

    @patch("atomium.converters.pdbdatafile2model.load_chains")
    def test_can_get_model_from_data_file(self, mock_chains):
        data_file = Mock(PdbDataFile)
        data_file.atoms = "aaa"
        model = pdb_data_file_to_model(data_file)
        self.assertIsInstance(model, Model)
        mock_chains.assert_called_with("aaa", model)



class ChainLoadingTests(TestCase):

    @patch("atomium.converters.pdbdatafile2model.residues_to_chains")
    @patch("atomium.converters.pdbdatafile2model.atoms_to_residues")
    @patch("atomium.converters.pdbdatafile2model.atom_dict_to_atom")
    def test_can_load_chains(self, mock_atom, mock_residues, mock_chains):
        atom1, atom2, atom3 = Mock(), Mock(), Mock()
        mock_atom.side_effect = [atom1, atom2, atom3]
        residue1, residue2 = Mock(), Mock()
        mock_residues.return_value = [residue1, residue2]
        chain1, chain2 = Mock(), Mock()
        mock_chains.return_value = [chain1, chain2]
        atomdict1, atomdict2, atomdict3 = Mock(), Mock(), Mock()
        model = Mock()
        model.add_chain = MagicMock()
        load_chains([atomdict1, atomdict2, atomdict3], model)
        mock_atom.assert_any_call(atomdict1)
        mock_atom.assert_any_call(atomdict2)
        mock_atom.assert_any_call(atomdict3)
        mock_residues.assert_called_with([atom1, atom2, atom3])
        mock_chains.assert_called_with([residue1, residue2])
        model.add_chain.assert_any_call(chain1)
        model.add_chain.assert_any_call(chain2)



class AtomDictToAtomTest(TestCase):

    def setUp(self):
        self.atom_dict = {
         "atom_id": 107, "atom_name": "N1", "alt_loc": None,
         "residue_name": "GLY",
         "chain_id": "A", "residue_id": 13, "insert_code": "A",
         "x": 12.681, "y": 37.302, "z": -25.211,
         "occupancy": 1.0, "temperature_factor": 15.56,
         "element": "N", "charge": -2
        }


    def test_can_make_atom_from_dict(self):
        atom = atom_dict_to_atom(self.atom_dict)
        self.assertIsInstance(atom, Atom)
        self.assertEqual(atom._x, 12.681)
        self.assertEqual(atom._y, 37.302)
        self.assertEqual(atom._z, -25.211)
        self.assertEqual(atom._id, 107)
        self.assertEqual(atom._name, "N1")
        self.assertEqual(atom._element, "N")
        self.assertEqual(atom._charge, -2)
        self.assertEqual(atom.temp_chain_id, "A")
        self.assertEqual(atom.temp_residue_id, "A13A")
        self.assertEqual(atom.temp_residue_name, "GLY")


    def test_can_make_atom_from_minimal_dict(self):
        self.atom_dict["insert_code"] = None
        self.atom_dict["occupancy"] = None
        self.atom_dict["temperature_factor"] = None
        self.atom_dict["charge"] = None
        self.atom_dict["element"] = None
        atom = atom_dict_to_atom(self.atom_dict)
        self.assertIsInstance(atom, Atom)
        self.assertEqual(atom._element, "X")
        self.assertEqual(atom._charge, 0)
        self.assertEqual(atom.temp_residue_id, "A13")



class AtomsToResiduesTests(TestCase):

    def test_can_get_residues_from_atoms(self):
        atoms = []
        for n in range(12):
            atoms.append(Mock(Atom))
            atoms[-1].temp_residue_id = str(n // 2)
            atoms[-1].temp_residue_name = str(n // 2) + "nm"
            atoms[-1].temp_chain_id = "A" if n < 6 else "B"
        residues = atoms_to_residues(atoms)
        self.assertEqual(len(residues), 6)
        for index, residue in enumerate(residues):
            self.assertIsInstance(residue, Residue)
            self.assertEqual(residue._atoms, set(atoms[index * 2: index * 2 + 2]))
            self.assertEqual(residue._id, str(index))
            self.assertEqual(residue._name, str(index) + "nm")
            self.assertEqual(residue.temp_chain_id, "A" if index < 3 else "B")
        for atom in atoms:
            with self.assertRaises(AttributeError):
                atom.temp_residue_id
            with self.assertRaises(AttributeError):
                atom.temp_residue_name
            with self.assertRaises(AttributeError):
                atom.temp_chain_id



class ResiduesToChainsTests(TestCase):

    @patch("atomium.converters.pdbdatafile2model.Chain")
    def test_can_chains_from_residues(self, mock_chain):
        chain1, chain2 = Mock(), Mock()
        mock_chain.side_effect = [chain1, chain2]
        residues = []
        for n in range(6):
            residues.append(Mock(Residue))
            residues[-1].temp_chain_id = "A" if n < 3 else "B"
            residues[-1].next = MagicMock()
        chains = residues_to_chains(residues)
        self.assertEqual(chains, [chain1, chain2])
        residues[0].next.assert_called_with(residues[1])
        residues[1].next.assert_called_with(residues[2])
        residues[3].next.assert_called_with(residues[4])
        residues[4].next.assert_called_with(residues[5])
        mock_chain.assert_any_call(*residues[:3], chain_id="A")
        mock_chain.assert_any_call(*residues[3:], chain_id="B")
        for residue in residues:
            with self.assertRaises(AttributeError):
                residue.temp_chain_id
