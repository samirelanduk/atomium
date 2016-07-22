
from unittest import TestCase
import unittest.mock
from molecupy.structures import SmallMolecule, AtomicStructure, PdbAtom

class SmallMoleculeTest(TestCase):

    def setUp(self):
        self.atoms = [unittest.mock.Mock(spec=PdbAtom) for _ in range(10)]



class SmallMoleculeCreationTests(SmallMoleculeTest):

    def test_can_create_small_molecule(self):
        small_molecule = SmallMolecule("A500", "MOL", *self.atoms)
        self.assertIsInstance(small_molecule, AtomicStructure)
        self.assertEqual(small_molecule._molecule_id, "A500")
        self.assertEqual(small_molecule._molecule_name, "MOL")
        self.assertEqual(small_molecule._atoms, set(self.atoms))


    def test_molecule_id_must_be_str(self):
        with self.assertRaises(TypeError):
            SmallMolecule(1.1, "HET", *self.atoms)


    def test_molecule_name_must_be_str(self):
        with self.assertRaises(TypeError):
            SmallMolecule("A1", 1, *self.atoms)


    def test_small_molecule_repr(self):
        small_molecule = SmallMolecule("A500", "MOL", *self.atoms)
        self.assertEqual(str(small_molecule), "<SmallMolecule A500 (MOL)>")



class SmallMoleculePropertyTests(SmallMoleculeTest):

    def test_small_molecule_properties(self):
        small_molecule = SmallMolecule("A500", "MOL", *self.atoms)
        self.assertEqual(small_molecule.molecule_id(), "A500")
        self.assertEqual(small_molecule.molecule_name(), "MOL")


    def test_can_change_molecule_name(self):
        small_molecule = SmallMolecule("A500", "MOL", *self.atoms)
        self.assertEqual(small_molecule.molecule_name(), "MOL")
        small_molecule.molecule_name("LOM")
        self.assertEqual(small_molecule.molecule_name(), "LOM")
