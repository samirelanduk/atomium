from unittest import TestCase
from unittest.mock import Mock, patch
from atomium.structures.molecules import Molecule, AtomicStructure
from atomium.structures.atoms import Atom

class MoleculeTest(TestCase):

    def setUp(self):
        self.atom1, self.atom2, self.atom3 = Mock(Atom), Mock(Atom), Mock(Atom)
        self.atoms = [self.atom1, self.atom2, self.atom3]



class MoleculeCreationTests(MoleculeTest):

    @patch("atomium.structures.molecules.AtomicStructure.__init__")
    def test_can_create_molecule(self, mock_init):
        mock_init.return_value = None
        mol = Molecule(self.atom1, self.atom2, self.atom3)
        self.assertIsInstance(mol, AtomicStructure)
        mock_init.assert_called_with(mol, self.atom1, self.atom2, self.atom3)
        self.assertEqual(mol._id, None)
        self.assertEqual(mol._name, None)


    def test_can_create_molecule_with_id(self):
        mol = Molecule(self.atom1, self.atom2, self.atom3, molecule_id="A100")
        self.assertEqual(mol._id, "A100")


    def test_molecule_id_must_be_str(self):
        with self.assertRaises(TypeError):
            Molecule(self.atom1, self.atom2, self.atom3, molecule_id=1000)


    def test_molecule_id_must_be_unique(self):
        mol = Molecule(self.atom1, self.atom2, self.atom3, molecule_id="A110")
        with self.assertRaises(ValueError):
            Molecule(self.atom1, self.atom2, self.atom3, molecule_id="A110")


    def test_can_create_molecule_with_name(self):
        mol = Molecule(self.atom1, self.atom2, self.atom3, name="HIS")
        self.assertEqual(mol._name, "HIS")


    def test_molecule_name_must_be_str(self):
        with self.assertRaises(TypeError):
            Molecule(self.atom1, self.atom2, self.atom3, name=1000)



class MoleculeIdChangesTests(MoleculeTest):

    def test_changing_id_updates_known_ids(self):
        mol = Molecule(self.atom1, self.atom2, self.atom3, molecule_id="A120")
        self.assertIn("A120", Molecule.known_ids)
        mol._id = "A121"
        self.assertIn("A121", Molecule.known_ids)
        self.assertNotIn("A120", Molecule.known_ids)


    def test_deleting_a_molecule_frees_up_its_id(self):
        mol = Molecule(self.atom1, self.atom2, self.atom3, molecule_id="A130")
        self.assertIn("A130", Molecule.known_ids)
        mol = Molecule(self.atom1, self.atom2, self.atom3, molecule_id="A131")
        self.assertNotIn("A130", Molecule.known_ids)
        self.assertIn("A131", Molecule.known_ids)



class MoleculeReprTests(MoleculeTest):

    def test_molecule_repr_no_id_or_name(self):
        mol = Molecule(self.atom1, self.atom2, self.atom3)
        self.assertEqual(str(mol), "<Molecule (3 atoms)>")


    def test_molecule_repr_id_no_name(self):
        mol = Molecule(self.atom1, self.atom2, self.atom3, molecule_id="B10A")
        self.assertEqual(str(mol), "<Molecule B10A (3 atoms)>")


    def test_molecule_repr_name_no_id(self):
        mol = Molecule(self.atom1, self.atom2, self.atom3, name="GLY")
        self.assertEqual(str(mol), "<Molecule (GLY, 3 atoms)>")


    def test_molecule_repr_id_and_name(self):
        mol = Molecule(self.atom1, self.atom2, molecule_id="B10A", name="GLY")
        self.assertEqual(str(mol), "<Molecule B10A (GLY, 2 atoms)>")



class MoleculeIdTests(MoleculeTest):

    def test_molecule_id_property(self):
        mol = Molecule(self.atom1, self.atom2, self.atom3, molecule_id="B10A")
        self.assertIs(mol._id, mol.molecule_id())


    def test_can_update_molecule_id(self):
        mol = Molecule(self.atom1, self.atom2, self.atom3, molecule_id="B10A")
        mol.molecule_id("A20")
        self.assertEqual(mol._id, "A20")


    def test_molecule_id_must_be_str(self):
        mol = Molecule(self.atom1, self.atom2, self.atom3, molecule_id="B10A")
        with self.assertRaises(TypeError):
            mol.molecule_id(10)



class MoleculeNameTests(MoleculeTest):

    def test_molecule_name_property(self):
        mol = Molecule(self.atom1, self.atom2, self.atom3, name="VAL")
        self.assertIs(mol._name, mol.name())


    def test_can_update_molecule_name(self):
        mol = Molecule(self.atom1, self.atom2, self.atom3, name="VAL")
        mol.name("HIS")
        self.assertEqual(mol._name, "HIS")


    def test_molecule_name_must_be_str(self):
        mol = Molecule(self.atom1, self.atom2, self.atom3, name="VAL")
        with self.assertRaises(TypeError):
            mol.name(10)
