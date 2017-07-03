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
        mock_init.assert_called_with(self.atom1, self.atom2, self.atom3)
