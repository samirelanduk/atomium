from unittest import TestCase
from molecupy import exceptions
from molecupy.molecules import Atom, Molecule, AtomicStructure, Model

class ModelTest(TestCase):

    def setUp(self):
        self.atom1 = Atom(1.0, 1.0, 1.0, "H", atom_id=1, atom_name="H1")
        self.atom2 = Atom(1.0, 1.0, 2.0, "C", atom_id=2, atom_name="CA")
        self.atom3 = Atom(1.0, 1.0, 3.0, "O", atom_id=3, atom_name="OX1")
        self.atom2.covalent_bond_to(self.atom1)
        self.atom2.covalent_bond_to(self.atom3)
        self.molecule1 = Molecule(self.atom1, self.atom2, self.atom3)
        self.atom4 = Atom(1.0, 1.0, 4.0, "H", atom_id=4, atom_name="H1")
        self.atom5 = Atom(1.0, 1.0, 5.0, "C", atom_id=5, atom_name="CA")
        self.atom6 = Atom(1.0, 1.0, 6.0, "O", atom_id=6, atom_name="OX1")
        self.atom5.covalent_bond_to(self.atom4)
        self.atom5.covalent_bond_to(self.atom6)
        self.molecule2 = Molecule(self.atom4, self.atom5, self.atom6)
        self.atom7 = Atom(1.0, 1.0, 7.0, "H", atom_id=7, atom_name="H1")
        self.atom8 = Atom(1.0, 1.0, 8.0, "C", atom_id=8, atom_name="CA")
        self.atom9 = Atom(1.0, 1.0, 9.0, "O", atom_id=9, atom_name="OX1")
        self.atom8.covalent_bond_to(self.atom7)
        self.atom8.covalent_bond_to(self.atom9)
        self.molecule3 = Molecule(self.atom7, self.atom8, self.atom9)


    def check_valid_model(self, model):
        self.assertIsInstance(model, Model)
        self.assertIsInstance(model, AtomicStructure)
        self.assertIsInstance(model._molecules, set)
        self.assertIsInstance(model.get_molecules(), set)
        self.assertIsInstance(model.atoms, set)
        self.assertRegex(str(model), r"<Model \((\d+) atoms\)>")



class ModelCreationTests(ModelTest):

    def test_can_make_model(self):
        model = Model()
        self.check_valid_model(model)



class ModelPopulationTests(ModelTest):

    def test_can_add_molecules_to_model(self):
        model = Model()
        model.add_molecule(self.molecule1)
        self.assertEqual(len(model.get_molecules()), 1)
        self.assertEqual(len(model.atoms), 3)
        self.assertEqual(self.molecule1.model, model)
        model.add_molecule(self.molecule2)
        self.assertEqual(len(model.get_molecules()), 2)
        self.assertEqual(len(model.atoms), 6)
        self.assertEqual(self.molecule2.model, model)
        model.add_molecule(self.molecule3)
        self.assertEqual(len(model.get_molecules()), 3)
        self.assertEqual(len(model.atoms), 9)
        self.assertEqual(self.molecule3.model, model)


    def test_can_only_add_molecules_to_model(self):
        model = Model()
        with self.assertRaises(TypeError):
            model.add_molecule("molecule")
        atomic_structure = AtomicStructure(self.atom1, self.atom2, self.atom3)
        with self.assertRaises(TypeError):
            model.add_molecule(atomic_structure)


    def test_molecules_must_have_unique_IDs(self):
        model = Model()
        model.add_molecule(self.molecule1)
        self.molecule1.molecule_id = 1
        self.molecule2.molecule_id = 1
        with self.assertRaises(exceptions.DuplicateMoleculeError):
            model.add_molecule(self.molecule2)




class MoleculeRetrievalTests(ModelTest):

    def test_can_get_molecule_by_id(self):
        self.molecule1.molecule_id = 1
        self.molecule2.molecule_id = 2
        self.molecule3.molecule_id = 3
        model = Model()
        model.add_molecule(self.molecule1)
        model.add_molecule(self.molecule2)
        model.add_molecule(self.molecule3)
        self.assertIs(model.get_molecule_by_id(1), self.molecule1)
        self.assertIs(model.get_molecule_by_id(2), self.molecule2)
        self.assertIs(model.get_molecule_by_id(3), self.molecule3)
        self.assertIs(model.get_molecule_by_id(4), None)


    def test_can_only_search_by_numeric_id(self):
        model = Model()
        with self.assertRaises(TypeError):
            model.get_molecule_by_id(None)
