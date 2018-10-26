from unittest import TestCase
from unittest.mock import Mock, patch, PropertyMock, MagicMock
from atomium.structures import Ligand, Atom, Molecule, Het, AtomStructure

class LigandTest(TestCase):

    def setUp(self):
        self.atoms = [Mock(Atom), Mock(Atom)]



class LigandCreationTests(LigandTest):

    @patch("atomium.structures.Het.__init__")
    def test_can_create_ligand(self, mock_init):
        ligand = Ligand(*self.atoms)
        self.assertIsNone(ligand._id)
        self.assertIsNone(ligand._name)
        self.assertIsNone(ligand._internal_id)
        self.assertIsNone(ligand._chain)
        self.assertFalse(ligand._water)
        mock_init.assert_called_with(ligand, *self.atoms)


    @patch("atomium.structures.Het.__init__")
    def test_can_create_ligand_with_attributes(self, mock_init):
        ligand = Ligand(*self.atoms, id="A1", name="XMP", chain="CH", internal_id="A", water=True)
        self.assertEqual(ligand._id, "A1")
        self.assertEqual(ligand._name, "XMP")
        self.assertEqual(ligand._chain, "CH")
        self.assertEqual(ligand._internal_id, "A")
        self.assertTrue(ligand._water)
        mock_init.assert_called_with(ligand, *self.atoms)



class LigandReprTests(LigandTest):

    def test_ligand_repr(self):
        ligand = Ligand(*self.atoms, id="A1", name="XMP", chain="CH", water=False)
        self.assertEqual(repr(ligand), "<Ligand XMP (A1)>")


    def test_water_ligand_repr(self):
        ligand = Ligand(*self.atoms, id="A1", name="HOH", chain="CH", water=True)
        self.assertEqual(repr(ligand), "<Water HOH (A1)>")



class LigandAtomsTests(LigandTest):

    def test_can_get_atoms(self):
        ligand = Ligand(*self.atoms, id="A1", name="XMP", chain="CH", water=False)
        self.assertEqual(ligand.atoms(), set(self.atoms))



class LigandWaterTests(LigandTest):

    def test_can_get_water_marker(self):
        ligand = Ligand(*self.atoms, id="A1", name="XMP", chain="CH", water=True)
        self.assertTrue(ligand.water)
        ligand._water = False
        self.assertFalse(ligand.water)



class LigandCopyingTests(LigandTest):

    @patch("atomium.structures.Ligand.atoms")
    def test_can_copy_ligand(self, mock_atoms):
        mock_atoms.return_value = self.atoms
        ligand = Ligand(*self.atoms, id="A1", name="XMP", internal_id="B", chain="C", water=True)
        p = patch("atomium.structures.Ligand")
        mock_ligand = p.start()
        try:
            copy = ligand.copy()
            mock_ligand.assert_called_with(
             self.atoms[0].copy.return_value, self.atoms[1].copy.return_value,
             id="A1", name="XMP", internal_id="B", water=True
            )
        finally: p.stop()
