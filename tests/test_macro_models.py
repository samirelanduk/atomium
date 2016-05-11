from unittest import TestCase
from molecupy import exceptions
from molecupy.molecules import Atom, Molecule, Model
from molecupy.macromolecules import Residue, ResiduicStructure, Chain, MacroModel, Complex

class MacroModelTest(TestCase):

    def setUp(self):
        self.atom1 = Atom(1.0, 1.0, 1.0, "H", atom_id=1, atom_name="H1")
        self.atom2 = Atom(1.0, 1.0, 2.0, "C", atom_id=2, atom_name="CA")
        self.atom3 = Atom(1.0, 1.0, 3.0, "O", atom_id=3, atom_name="OX1")
        self.atom2.covalent_bond_to(self.atom1)
        self.atom2.covalent_bond_to(self.atom3)
        self.residue1 = Residue(1, "MON1", self.atom1, self.atom2, self.atom3)
        self.atom4 = Atom(1.0, 1.0, 4.0, "H", atom_id=4, atom_name="H1")
        self.atom5 = Atom(1.0, 1.0, 5.0, "C", atom_id=5, atom_name="CA")
        self.atom6 = Atom(1.0, 1.0, 6.0, "O", atom_id=6, atom_name="OX1")
        self.atom5.covalent_bond_to(self.atom4)
        self.atom5.covalent_bond_to(self.atom6)
        self.residue2 = Residue(2, "MON2", self.atom4, self.atom5, self.atom6)
        self.atom7 = Atom(1.0, 1.0, 7.0, "H", atom_id=7, atom_name="H1")
        self.atom8 = Atom(1.0, 1.0, 8.0, "C", atom_id=8, atom_name="CA")
        self.atom9 = Atom(1.0, 1.0, 9.0, "O", atom_id=9, atom_name="OX1")
        self.atom8.covalent_bond_to(self.atom7)
        self.atom8.covalent_bond_to(self.atom9)
        self.residue3 = Residue(3, "MON3", self.atom7, self.atom8, self.atom9)
        self.residue1.connect_to(self.residue2, self.atom3, self.atom4)
        self.residue2.connect_to(self.residue3, self.atom6, self.atom7)


    def check_valid_macro_model(self, macro_model):
        self.assertIsInstance(macro_model, MacroModel)
        self.assertIsInstance(macro_model, Model)
        self.assertIsInstance(macro_model._molecules, set)
        self.assertIsInstance(macro_model.get_molecules(), set)
        self.assertIsInstance(macro_model.atoms, set)
        self.assertIsInstance(macro_model._chains, set)
        self.assertIsInstance(macro_model.get_chains(), set)
        self.assertIsInstance(macro_model._small_molecules, set)
        self.assertIsInstance(macro_model.get_small_molecules(), set)
        self.assertIsInstance(macro_model._complexes, set)
        self.assertIsInstance(macro_model.get_complexes(), set)
        self.assertRegex(str(macro_model), r"<MacroModel \((\d+) atoms\)>")



class MacroModelCreationTests(MacroModelTest):

    def test_can_make_macro_model(self):
        macro_model = MacroModel()
        self.check_valid_macro_model(macro_model)



class ChainAdditionTests(MacroModelTest):

    def test_can_add_chains(self):
        chain = Chain(self.residue1, self.residue2, self.residue3)
        macro_model = MacroModel()
        macro_model.add_chain(chain)
        self.assertEqual(len(macro_model.get_chains()), 1)
        self.assertEqual(len(macro_model.get_molecules()), 1)
        self.assertEqual(chain.model, macro_model)
        self.assertEqual(len(macro_model.atoms), 9)


    def test_chains_must_be_chains(self):
        macro_model = MacroModel()
        with self.assertRaises(TypeError):
            macro_model.add_chain("chain")



class SmallMoleculeAdditionTests(MacroModelTest):

    def test_can_add_small_molecules(self):
        small_molecule = Molecule(self.atom4, self.atom5, self.atom6)
        macro_model = MacroModel()
        macro_model.add_small_molecule(small_molecule)
        self.assertEqual(len(macro_model.get_small_molecules()), 1)
        self.assertEqual(len(macro_model.get_molecules()), 1)
        self.assertEqual(small_molecule.model, macro_model)
        self.assertEqual(len(macro_model.atoms), 3)


    def test_chains_must_be_chains(self):
        macro_model = MacroModel()
        with self.assertRaises(TypeError):
            macro_model.add_chain("mol")



class ComplexAdditionTests(MacroModelTest):

    def test_can_add_complexes(self):
        chain1 = Chain(self.residue1, self.residue2)
        chain2 = Chain(self.residue3)
        complex_ = Complex(chain1, chain2)
        macro_model = MacroModel()
        macro_model.add_complex(complex_)
        self.assertEqual(len(macro_model.get_complexes()), 1)
        self.assertEqual(len(macro_model.get_chains()), 2)
        self.assertEqual(len(macro_model.get_molecules()), 2)
        self.assertEqual(complex_.model, macro_model)
        self.assertEqual(chain1.model, macro_model)
        self.assertEqual(chain2.model, macro_model)
        self.assertEqual(len(macro_model.atoms), 9)


    def test_complexes_must_be_complexes(self):
        macro_model = MacroModel()
        with self.assertRaises(TypeError):
            macro_model.add_complex("complex")
