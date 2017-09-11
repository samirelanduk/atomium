from unittest import TestCase
from unittest.mock import patch, Mock, MagicMock
from atomium.converters.pdbdatafile2models import *
from atomium.files.pdbdatafile import PdbDataFile
from atomium.structures.models import Model
from atomium.structures.molecules import Residue, Molecule
from atomium.structures.atoms import Atom
from atomium.structures.reference import bonds

class PdbDataFileToModelsTests(TestCase):

    @patch("atomium.converters.pdbdatafile2models.load_chains")
    @patch("atomium.converters.pdbdatafile2models.load_molecules")
    @patch("atomium.converters.pdbdatafile2models.bond_atoms")
    def test_can_get_models_from_data_file(self, mock_bond, mock_mol, mock_chains):
        data_file = Mock(PdbDataFile)
        data_file.atoms = [{"model": 1}, {"model": 2}, {"model": 3}]
        data_file.heteroatoms = "hhh"
        data_file.connections = "ccc"
        models = pdb_data_file_to_models(data_file)
        self.assertEqual(len(models), 3)
        self.assertIsInstance(models[0], Model)
        self.assertIsInstance(models[1], Model)
        self.assertIsInstance(models[2], Model)
        mock_chains.assert_any_call(
         [{"model": 1}, {"model": 2}, {"model": 3}], models[0], 1
        )
        mock_chains.assert_any_call(
         [{"model": 1}, {"model": 2}, {"model": 3}], models[1], 2
        )
        mock_chains.assert_any_call(
         [{"model": 1}, {"model": 2}, {"model": 3}], models[2], 3
        )
        mock_mol.assert_any_call("hhh", models[0], 1)
        mock_mol.assert_any_call("hhh", models[1], 2)
        mock_mol.assert_any_call("hhh", models[2], 3)
        mock_bond.assert_any_call("ccc", models[0])
        mock_bond.assert_any_call("ccc", models[1])
        mock_bond.assert_any_call("ccc", models[2])



class ChainLoadingTests(TestCase):

    @patch("atomium.converters.pdbdatafile2models.residues_to_chains")
    @patch("atomium.converters.pdbdatafile2models.atoms_to_residues")
    @patch("atomium.converters.pdbdatafile2models.atom_dict_to_atom")
    def test_can_load_chains(self, mock_atom, mock_residues, mock_chains):
        atom1, atom2, atom3 = Mock(), Mock(), Mock()
        mock_atom.side_effect = [atom1, atom2, atom3]
        residue1, residue2 = Mock(), Mock()
        mock_residues.return_value = [residue1, residue2]
        chain1, chain2 = Mock(), Mock()
        mock_chains.return_value = [chain1, chain2]
        atomdict1, atomdict2, atomdict3 = [{"model": 2}, {"model": 2}, {"model": 3}]
        model = Mock()
        model.add_chain = MagicMock()
        load_chains([atomdict1, atomdict2, atomdict3], model, 2)
        mock_atom.assert_any_call(atomdict1)
        mock_atom.assert_any_call(atomdict2)
        mock_residues.assert_called_with([atom1, atom2])
        mock_chains.assert_called_with([residue1, residue2])
        model.add_chain.assert_any_call(chain1)
        model.add_chain.assert_any_call(chain2)



class MoleculeLoadingTests(TestCase):

    @patch("atomium.converters.pdbdatafile2models.atoms_to_residues")
    @patch("atomium.converters.pdbdatafile2models.atom_dict_to_atom")
    def test_can_load_molecules(self, mock_atom, mock_residues):
        atom1, atom2, atom3 = Mock(), Mock(), Mock()
        mock_atom.side_effect = [atom1, atom2, atom3]
        mol1, mol2 = Mock(), Mock()
        mock_residues.return_value = [mol1, mol2]
        atomdict1, atomdict2, atomdict3 = [{"model": 2}, {"model": 2}, {"model": 3}]
        model = Mock()
        model.add_molecule = MagicMock()
        load_molecules([atomdict1, atomdict2, atomdict3], model, 2)
        mock_atom.assert_any_call(atomdict1)
        mock_atom.assert_any_call(atomdict2)
        mock_residues.assert_called_with([atom1, atom2], molecule=True)
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
        self.assertEqual(atom._bfactor, 15.56)
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
        self.assertEqual(atom._bfactor, 0)
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
            self.assertEqual(
             set(residue._id_atoms.values()), set(atoms[index * 2: index * 2 + 2])
            )
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
            self.assertEqual(
             set(molecule._id_atoms.values()), set(atoms[index * 2: index * 2 + 2])
            )
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

    @patch("atomium.converters.pdbdatafile2models.Chain")
    def test_can_get_chains_from_residues(self, mock_chain):
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

    @patch("atomium.converters.pdbdatafile2models.make_intra_residue_bonds")
    @patch("atomium.converters.pdbdatafile2models.make_inter_residue_bonds")
    @patch("atomium.converters.pdbdatafile2models.make_connections_bonds")
    def test_can_bond_atoms(self, mock_con, mock_inter, mock_intra):
        residues = [Mock(), Mock()]
        model = Mock()
        model.residues.return_value = set(residues)
        connections = "ccc"
        bond_atoms(connections, model)
        mock_intra.assert_called_with(set(residues), bonds)
        mock_inter.assert_called_with(set(residues))
        mock_con.assert_called_with(model, "ccc")



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



class InterResidueConnectionTests(TestCase):

    def test_can_connect_residues(self):
        residues = [Mock(), Mock(), Mock(), Mock()]
        atoms = [Mock(), Mock(), Mock(), Mock(), Mock(), Mock()]
        atoms += [Mock(), Mock(), Mock(), Mock(), Mock(), Mock()]
        for atom in atoms:
            atom.bond = MagicMock()
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
            residue.next.return_value = residues[i + 1] if i != 3 else None
        make_inter_residue_bonds(residues)
        atoms[2].bond.assert_called_with(atoms[3])
        atoms[5].bond.assert_called_with(atoms[6])
        atoms[8].bond.assert_called_with(atoms[9])
        self.assertFalse(atoms[0].bond.called)
        self.assertFalse(atoms[1].bond.called)
        self.assertFalse(atoms[3].bond.called)
        self.assertFalse(atoms[4].bond.called)
        self.assertFalse(atoms[6].bond.called)
        self.assertFalse(atoms[7].bond.called)
        self.assertFalse(atoms[9].bond.called)
        self.assertFalse(atoms[10].bond.called)
        self.assertFalse(atoms[11].bond.called)


    def test_can_skip_bond_where_atom_not_present(self):
        residue1, residue2 = Mock(), Mock()
        residue1.next.return_value = residue2
        residue2.next.return_value = None
        atom = Mock()
        atom.bond = MagicMock()
        atom.distance_to.return_value = 4.9
        residue1.atom.return_value = atom
        residue2.atom.return_value = None
        make_inter_residue_bonds([residue1, residue2])
        self.assertFalse(atom.bond.called)
        residue2.atom.return_value = atom
        residue1.atom.return_value = None
        make_inter_residue_bonds([residue1, residue2])
        self.assertFalse(atom.bond.called)


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
            atom.atom_id.return_value = index
            atom.bond = MagicMock()
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
        atoms[1].bond.assert_any_call(atoms[2])
        atoms[1].bond.assert_any_call(atoms[3])
        atoms[2].bond.assert_called_with(atoms[1])
        atoms[3].bond.assert_called_with(atoms[1])
        self.assertFalse(atoms[0].called)
