from unittest import TestCase
from unittest.mock import patch, Mock, MagicMock
from atomium.structures.models import MoleculeStructure
from atomium.structures.molecules import Molecule, Residue
from atomium.structures.atoms import Atom

class MoleculeStructureTest(TestCase):

    def setUp(self):
        self.structure = MoleculeStructure()
        self.atom1, self.atom2 = Mock(Atom), Mock(Atom)
        self.atom3, self.atom4 = Mock(Atom), Mock(Atom)
        self.molecule1, self.molecule2 = Mock(Molecule), Mock(Residue)
        self.molecule1.molecule_id.return_value = "C"
        self.molecule2.molecule_id.return_value = "D"
        self.molecule1.name.return_value = "III"
        self.molecule2.name.return_value = "JJJ"
        self.atom1.molecule.return_value = self.molecule1
        self.atom2.molecule.return_value = self.molecule1
        self.atom3.molecule.return_value = self.molecule2
        self.atom4.molecule.return_value = self.molecule2
        self.molecule1.atoms.return_value = set([self.atom1, self.atom2])
        self.molecule2.atoms.return_value = set([self.atom3, self.atom4])
        self.structure.atoms = lambda: set([
         self.atom1, self.atom2, self.atom3, self.atom4
        ])
        self.structure.add_atom = MagicMock()
        self.structure.remove_atom = MagicMock()


class MoleculeStructureMoleculesTests(MoleculeStructureTest):

    def test_can_get_molecules(self):
        self.assertEqual(
         self.structure.molecules(),
         set([self.molecule1, self.molecule2])
        )


    def test_can_filter_none_from_molecules(self):
        self.atom4.molecule.return_value = None
        self.assertEqual(
         self.structure.molecules(),
         set([self.molecule1, self.molecule2])
        )


    def test_can_get_molecules_by_id(self):
        self.assertEqual(
         self.structure.molecules(molecule_id="C"), set([self.molecule1])
        )
        self.assertEqual(
         self.structure.molecules(molecule_id="D"), set([self.molecule2])
        )
        self.assertEqual(self.structure.molecules(molecule_id="E"), set())


    def test_can_get_molecules_by_name(self):
        self.assertEqual(
         self.structure.molecules(name="III"), set([self.molecule1])
        )
        self.assertEqual(
         self.structure.molecules(name="JJJ"), set([self.molecule2])
        )
        self.assertEqual(self.structure.molecules(name="GLY"), set())


    def test_can_filter_out_specific_molecule_types(self):
        self.assertEqual(
         self.structure.molecules(generic=True), set([self.molecule1])
        )


    def test_can_filter_out_water(self):
        self.molecule2.name.return_value = "WAT"
        self.assertEqual(
         self.structure.molecules(water=False), set([self.molecule1])
        )
        self.molecule2.name.return_value = "HOH"
        self.assertEqual(
         self.structure.molecules(water=False), set([self.molecule1])
        )



class MoleculeStructureMoleculeTests(MoleculeStructureTest):

    @patch("atomium.structures.models.MoleculeStructure.molecules")
    def test_molecule_calls_molecules(self, mock_molecules):
        mock_molecules.return_value = set([self.molecule2])
        molecule = self.structure.molecule(name="A")
        mock_molecules.assert_called_with(name="A")
        self.assertIs(molecule, self.molecule2)


    @patch("atomium.structures.models.MoleculeStructure.molecules")
    def test_molecule_can_return_none(self, mock_molecules):
        mock_molecules.return_value = set()
        self.assertIs(self.structure.molecule(name="E"), None)


    @patch("atomium.structures.models.MoleculeStructure.molecules")
    def test_molecule_can_get_molecule_by_id_and_name(self, mock_molecules):
        mock_molecules.return_value = set([self.molecule1])
        molecule = self.structure.molecule(molecule_id="A1", name="A")
        mock_molecules.assert_called_with(molecule_id="A1", name="A")
        self.assertIs(molecule, self.molecule1)



class MoleculeStructureMoleculeAdditionTests(MoleculeStructureTest):

    def test_can_add_molecule(self):
        self.structure._atoms = set([
         self.atom1, self.atom2
        ])
        self.structure.add_molecule(self.molecule2)
        self.structure.add_atom.assert_any_call(self.atom3)
        self.structure.add_atom.assert_any_call(self.atom4)


    def test_can_only_add_molecules(self):
        with self.assertRaises(TypeError):
            self.structure.add_molecule("self.molecule4")



class MoleculeStructureMoleculeRemovalTests(MoleculeStructureTest):

    def test_can_remove_molecule(self):
        self.structure.remove_molecule(self.molecule2)
        self.structure.remove_atom.assert_any_call(self.atom3)
        self.structure.remove_atom.assert_any_call(self.atom4)


    def test_can_only_remove_molecules(self):
        with self.assertRaises(TypeError):
            self.structure.remove_molecule("self.molecule4")
