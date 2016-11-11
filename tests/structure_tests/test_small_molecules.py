from unittest import TestCase
import unittest.mock
from molecupy.structures import SmallMolecule, AtomicStructure, Atom, BindSite

class SmallMoleculeTest(TestCase):

    def setUp(self):
        self.atoms = [unittest.mock.Mock(Atom) for _ in range(10)]
        for atom in self.atoms:
            atom._molecule = None



class SmallMoleculeCreationTests(SmallMoleculeTest):

    def test_can_create_small_molecule(self):
        small_molecule = SmallMolecule("A500", "MOL", *self.atoms)
        self.assertIsInstance(small_molecule, AtomicStructure)
        self.assertEqual(small_molecule._molecule_id, "A500")
        self.assertEqual(small_molecule._molecule_name, "MOL")
        self.assertEqual(small_molecule._bind_site, None)
        self.assertEqual(small_molecule._model, None)
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


    def test_small_molecules_update_atoms(self):
        for atom in self.atoms:
            self.assertIs(atom._molecule, None)
        small_molecule = SmallMolecule("A500", "MOL", *self.atoms)
        for atom in self.atoms:
            self.assertIs(atom._molecule, small_molecule)



class SmallMoleculePropertyTests(SmallMoleculeTest):

    def test_small_molecule_properties(self):
        small_molecule = SmallMolecule("A500", "MOL", *self.atoms)
        self.assertEqual(small_molecule.molecule_id(), "A500")
        self.assertEqual(small_molecule.molecule_name(), "MOL")
        self.assertEqual(small_molecule.bind_site(), None)
        self.assertEqual(small_molecule.model(), None)


    def test_can_change_molecule_name(self):
        small_molecule = SmallMolecule("A500", "MOL", *self.atoms)
        self.assertEqual(small_molecule.molecule_name(), "MOL")
        small_molecule.molecule_name("LOM")
        self.assertEqual(small_molecule.molecule_name(), "LOM")


    def test_molecule_name_can_only_be_changed_to_str(self):
        small_molecule = SmallMolecule("A500", "MOL", *self.atoms)
        with self.assertRaises(TypeError):
            small_molecule.molecule_name(100)


    def test_can_add_bind_site(self):
        bind_site = unittest.mock.Mock(spec=BindSite)
        bind_site._ligand = None
        small_molecule = SmallMolecule("A500", "MOL", *self.atoms)
        self.assertEqual(small_molecule.bind_site(), None)
        small_molecule.bind_site(bind_site)
        self.assertEqual(small_molecule.bind_site(), bind_site)
        self.assertEqual(bind_site._ligand, small_molecule)


    def test_small_molecule_can_only_add_bind_site(self):
        small_molecule = SmallMolecule("A500", "MOL", *self.atoms)
        with self.assertRaises(TypeError):
            small_molecule.bind_site("site")


    def test_small_molecule_integrates_added_atoms(self):
        small_molecule = SmallMolecule("A500", "MOL", *self.atoms[:-1])
        self.assertIs(self.atoms[-1]._molecule, None)
        small_molecule.add_atom(self.atoms[-1])
        self.assertIs(self.atoms[-1]._molecule, small_molecule)
        self.assertEqual(len(small_molecule.atoms()), 10)


    def test_small_molecule_deintegrates_removed_atoms(self):
        small_molecule = SmallMolecule("A500", "MOL", *self.atoms)
        self.assertIs(self.atoms[-1]._molecule, small_molecule)
        small_molecule.remove_atom(self.atoms[-1])
        self.assertIs(self.atoms[-1]._molecule, None)
        self.assertEqual(len(small_molecule.atoms()), 9)
