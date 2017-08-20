from unittest import TestCase
from unittest.mock import patch, Mock, MagicMock
from atomium.converters.pdbdatafile2model import *
from atomium.files.pdbdatafile import PdbDataFile
from atomium.structures.models import Model
from atomium.structures.molecules import Residue, Molecule
from atomium.structures.atoms import Atom
from atomium.structures.reference import bonds

class PdbDataFileToModelTests(TestCase):

    @patch("atomium.converters.pdbdatafile2model.load_chains")
    @patch("atomium.converters.pdbdatafile2model.load_molecules")
    @patch("atomium.converters.pdbdatafile2model.bond_atoms")
    def test_can_get_model_from_data_file(self, mock_bond, mock_mol, mock_chains):
        data_file = Mock(PdbDataFile)
        data_file.atoms = "aaa"
        data_file.heteroatoms = "hhh"
        model = pdb_data_file_to_model(data_file)
        self.assertIsInstance(model, Model)
        mock_chains.assert_called_with("aaa", model)
        mock_mol.assert_called_with("hhh", model)
        mock_bond.assert_called_with(model)



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



class MoleculeLoadingTests(TestCase):

    @patch("atomium.converters.pdbdatafile2model.atoms_to_residues")
    @patch("atomium.converters.pdbdatafile2model.atom_dict_to_atom")
    def test_can_load_molecules(self, mock_atom, mock_residues):
        atom1, atom2, atom3 = Mock(), Mock(), Mock()
        mock_atom.side_effect = [atom1, atom2, atom3]
        mol1, mol2 = Mock(), Mock()
        mock_residues.return_value = [mol1, mol2]
        atomdict1, atomdict2, atomdict3 = Mock(), Mock(), Mock()
        model = Mock()
        model.add_molecule = MagicMock()
        load_molecules([atomdict1, atomdict2, atomdict3], model)
        mock_atom.assert_any_call(atomdict1)
        mock_atom.assert_any_call(atomdict2)
        mock_atom.assert_any_call(atomdict3)
        mock_residues.assert_called_with([atom1, atom2, atom3], molecule=True)
        model.add_molecule.assert_any_call(mol1)
        model.add_molecule.assert_any_call(mol2)



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


    def test_can_get_molecules_from_atoms(self):
        atoms = []
        for n in range(12):
            atoms.append(Mock(Atom))
            atoms[-1].temp_residue_id = str(n // 2)
            atoms[-1].temp_residue_name = str(n // 2) + "nm"
            atoms[-1].temp_chain_id = "A" if n < 6 else "B"
        molecules = atoms_to_residues(atoms, molecule=True)
        self.assertEqual(len(molecules), 6)
        for index, molecule in enumerate(molecules):
            self.assertIsInstance(molecule, Molecule)
            self.assertNotIsInstance(molecule, Residue)
            self.assertEqual(molecule._atoms, set(atoms[index * 2: index * 2 + 2]))
            self.assertEqual(molecule._id, str(index))
            self.assertEqual(molecule._name, str(index) + "nm")
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



class AtomBondingTests(TestCase):

    @patch("atomium.converters.pdbdatafile2model.make_intra_residue_bonds")
    def test_can_bond_atoms(self, mock_intra):
        residues = [Mock(), Mock()]
        model = Mock()
        model.residues.return_value = set(residues)
        bond_atoms(model)
        mock_intra.assert_called_with(set(residues), bonds)



class IntraResidueConnectionTests(TestCase):

    def test_can_connect_residue(self):
        residues = [Mock(), Mock(), Mock(), Mock()]
        atoms = [Mock(), Mock(), Mock(), Mock(), Mock(), Mock()]
        atoms += [Mock(), Mock(), Mock(), Mock(), Mock(), Mock()]
        for atom in atoms:
            atom.bond = MagicMock()
        for i, residue in enumerate(residues):
            residue.atoms.return_value = set(atoms[i * 3: i * 3 + 3])
            residue.name.return_value = ["CYS", "TYR", "MET", "VAL"][i]
            atoms[i * 3].name.return_value = "A"
            atoms[i * 3 + 1].name.return_value = "B"
            atoms[i * 3 + 2].name.return_value = ["C", "P", "U", "D"][i]
        d = {
         "CYS": {"A": ["B"], "B": ["A"]},
         "TYR": {"A": ["B", "P"], "B": ["A"], "P": ["A"]},
         "MET": {"A": ["B", "X"], "X": ["A"], "B": ["A"]}
        }
        make_intra_residue_bonds(residues, d)
        atoms[0].bond.assert_called_with(atoms[1])
        atoms[1].bond.assert_called_with(atoms[0])
        atoms[3].bond.assert_any_call(atoms[4])
        atoms[3].bond.assert_any_call(atoms[5])
        atoms[4].bond.assert_called_with(atoms[3])
        atoms[5].bond.assert_called_with(atoms[3])
        atoms[6].bond.assert_called_with(atoms[7])
        atoms[7].bond.assert_called_with(atoms[6])
        self.assertFalse(atoms[2].bond.called)
        self.assertFalse(atoms[8].bond.called)
        self.assertFalse(atoms[9].bond.called)
        self.assertFalse(atoms[10].bond.called)
        self.assertFalse(atoms[11].bond.called)
