from unittest import TestCase
from molecupy import exceptions
from molecupy.molecules import Atom, Molecule, AtomicStructure, Model

class ModelTest(TestCase):

    def check_valid_model(self, model):
        self.assertIsInstance(model, Model)
        self.assertIsInstance(model, AtomicStructure)
        self.assertIsInstance(model.molecules, set)
        self.assertIsInstance(model.atoms, set)
        self.assertRegex(str(model), r"<Model \((\d+) atoms\)>")



class ModelCreationTests(ModelTest):

    def test_can_make_model(self):
        model = Model()
        self.check_valid_model(model)



class ModelPopulationTests(ModelTest):

    def test_can_add_molecules_to_model(self):
        atom1 = Atom(1.0, 1.0, 1.0, "H", atom_id=1, atom_name="H1")
        atom2 = Atom(1.0, 1.0, 2.0, "C", atom_id=2, atom_name="CA")
        atom3 = Atom(1.0, 1.0, 3.0, "O", atom_id=3, atom_name="OX1")
        atom2.covalent_bond_to(atom1)
        atom2.covalent_bond_to(atom3)
        molecule1 = Molecule(atom1, atom2, atom3)
        atom4 = Atom(1.0, 1.0, 4.0, "H", atom_id=4, atom_name="H1")
        atom5 = Atom(1.0, 1.0, 5.0, "C", atom_id=5, atom_name="CA")
        atom6 = Atom(1.0, 1.0, 6.0, "O", atom_id=6, atom_name="OX1")
        atom5.covalent_bond_to(atom4)
        atom5.covalent_bond_to(atom6)
        molecule2 = Molecule(atom4, atom5, atom6)
        atom7 = Atom(1.0, 1.0, 7.0, "H", atom_id=7, atom_name="H1")
        atom8 = Atom(1.0, 1.0, 8.0, "C", atom_id=8, atom_name="CA")
        atom9 = Atom(1.0, 1.0, 9.0, "O", atom_id=9, atom_name="OX1")
        atom8.covalent_bond_to(atom7)
        atom8.covalent_bond_to(atom9)
        molecule3 = Molecule(atom7, atom8, atom9)
        model = Model()

        model.add_molecule(molecule1)
        self.assertEqual(len(model.molecules), 1)
        self.assertEqual(len(model.atoms), 3)
        self.assertEqual(molecule1.model, model)
        model.add_molecule(molecule2)
        self.assertEqual(len(model.molecules), 2)
        self.assertEqual(len(model.atoms), 6)
        self.assertEqual(molecule2.model, model)
        model.add_molecule(molecule3)
        self.assertEqual(len(model.molecules), 3)
        self.assertEqual(len(model.atoms), 9)
        self.assertEqual(molecule3.model, model)


    def test_can_only_add_molecules_to_model(self):
        model = Model()
        with self.assertRaises(TypeError):
            model.add_molecule("molecule")
        atom1 = Atom(1.0, 1.0, 1.0, "H", atom_id=1, atom_name="H1")
        atom2 = Atom(1.0, 1.0, 2.0, "C", atom_id=2, atom_name="CA")
        atom3 = Atom(1.0, 1.0, 3.0, "O", atom_id=3, atom_name="OX1")
        atomic_structure = AtomicStructure(atom1, atom2, atom3)
        with self.assertRaises(TypeError):
            model.add_molecule(atomic_structure)
