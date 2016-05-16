from unittest import TestCase
from molecupy import exceptions
from molecupy.structures import PdbAtom, AtomicStructure, PdbModel, PdbSmallMolecule

class ModelTest(TestCase):

    def check_valid_model(self, model):
        self.assertIsInstance(model, PdbModel)
        self.assertIsInstance(model, AtomicStructure)
        self.assertIsInstance(model.atoms, set)
        self.assertIsInstance(model.small_molecules, set)
        self.assertRegex(str(model), r"<Model \((\d+) atoms\)>")



class ModelCreationTests(ModelTest):

    def test_can_make_model(self):
        model = PdbModel()
        self.check_valid_model(model)



class SmallMoleculeAdditionTests(ModelTest):

    def setUp(self):
        self.atom1 = PdbAtom(1.0, 1.0, 1.0, "H", 1, "H1")
        self.atom2 = PdbAtom(1.0, 1.0, 2.0, "C", 2, "CA")
        self.molecule1 = PdbSmallMolecule(1, "MOL1", self.atom1, self.atom2)
        self.atom3 = PdbAtom(1.0, 1.0, 3.0, "H", 3, "H1")
        self.atom4 = PdbAtom(1.0, 1.0, 4.0, "C", 4, "CA")
        self.molecule2 = PdbSmallMolecule(2, "MOL2", self.atom1, self.atom2)
        self.model = PdbModel()


    def test_can_add_small_molecules(self):
        self.model.add_small_molecule(self.molecule1)
        self.assertIn(self.molecule1, self.model.small_molecules)


    def test_cannot_add_two_small_molecules_with_same_id(self):
        self.model.add_small_molecule(self.molecule1)
        self.molecule2.molecule_id = 1
        with self.assertRaises(exceptions.DuplicateSmallMoleculeError):
            self.model.add_small_molecule(self.molecule2)
        self.molecule2.molecule_id = 2
        self.model.add_small_molecule(self.molecule2)
        self.assertIn(self.molecule1, self.model.small_molecules)
        self.assertIn(self.molecule2, self.model.small_molecules)


    def test_add_small_molecule_atoms_to_model(self):
        self.model.add_small_molecule(self.molecule1)
        self.assertEqual(self.model.atoms, self.molecule1.atoms)
        self.model.add_small_molecule(self.molecule2)
        self.assertEqual(
         self.model.atoms,
         set(list(self.molecule1.atoms) + list(self.molecule2.atoms))
        )


    def test_small_molecules_must_be_small_molecules(self):
        with self.assertRaises(TypeError):
            self.model.add_small_molecule("mol")


    def test_can_get_small_molecules_by_id(self):
        self.model.add_small_molecule(self.molecule1)
        self.model.add_small_molecule(self.molecule2)
        self.assertEqual(
         self.model.get_small_molecule_by_id(1),
         self.molecule1
        )
        self.assertEqual(
         self.model.get_small_molecule_by_id(2),
         self.molecule2
        )
        self.assertEqual(
         self.model.get_small_molecule_by_id(3),
         None
        )


    def test_can_only_search_by_numeric_id(self):
        with self.assertRaises(TypeError):
            self.model.get_small_molecule_by_id(None)


    def test_can_get_small_molecules_by_name(self):
        self.model.add_small_molecule(self.molecule1)
        self.model.add_small_molecule(self.molecule2)
        self.assertEqual(
         self.model.get_small_molecule_by_name("MOL1"),
         self.molecule1
        )
        self.assertEqual(
         self.model.get_small_molecule_by_name("MOLX"),
         None
        )


    def test_can_get_multiple_small_molecules_by_name(self):
        self.model.add_small_molecule(self.molecule1)
        self.model.add_small_molecule(self.molecule2)
        self.molecule2.molecule_name = "MOL1"
        self.assertEqual(
         self.model.get_small_molecules_by_name("MOL1"),
         set([self.molecule1, self.molecule2])
        )
        self.assertEqual(self.model.get_small_molecules_by_name("XXX"), set())


    def test_can_only_search_by_string_name(self):
        with self.assertRaises(TypeError):
            self.model.get_small_molecule_by_name(None)
        with self.assertRaises(TypeError):
            self.model.get_small_molecules_by_name(None)
