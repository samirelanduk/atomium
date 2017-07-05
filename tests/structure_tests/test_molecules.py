from unittest import TestCase
from unittest.mock import Mock, patch
from atomium.structures.molecules import Molecule
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
        self.assertIsInstance(mol, Molecule)
        mock_init.assert_called_with(mol, self.atom1, self.atom2, self.atom3)
        self.assertEqual(mol._id, None)


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
