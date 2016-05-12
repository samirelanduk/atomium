import copy
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
        self.atom10 = Atom(1.0, 1.0, 10.0, "H", atom_id=10, atom_name="H1")
        self.atom11 = Atom(1.0, 1.0, 11.0, "C", atom_id=11, atom_name="CA")
        self.atom12 = Atom(1.0, 1.0, 12.0, "O", atom_id=12, atom_name="OX1")
        self.atom11.covalent_bond_to(self.atom10)
        self.atom11.covalent_bond_to(self.atom12)
        self.residue4 = Residue(4, "MON4", self.atom10, self.atom11, self.atom12)
        self.residue1.connect_to(self.residue2, self.atom3, self.atom4)
        self.residue3.connect_to(self.residue4, self.atom9, self.atom10)


    def check_valid_macro_model(self, macro_model):
        self.assertIsInstance(macro_model, MacroModel)
        self.assertIsInstance(macro_model, Model)
        self.assertIsInstance(macro_model.atoms, set)
        self.assertIsInstance(macro_model._molecules, set)
        self.assertIsInstance(macro_model._small_molecules, set)
        self.assertIsInstance(macro_model._chains, set)
        self.assertIsInstance(macro_model._complexes, set)
        self.assertRegex(str(macro_model), r"<MacroModel \((\d+) atoms\)>")



class MacroModelCreationTests(MacroModelTest):

    def test_can_make_macro_model(self):
        macro_model = MacroModel()
        self.check_valid_macro_model(macro_model)



class SmallMoleculeAdditionTests(MacroModelTest):

    def setUp(self):
        MacroModelTest.setUp(self)
        self.molecule1 = Molecule(self.atom1, self.atom2, self.atom3)
        self.molecule2 = Molecule(self.atom4, self.atom5, self.atom6)


    def test_can_add_small_molecules(self):
        macro_model = MacroModel()
        macro_model.add_small_molecule(self.molecule1)
        self.assertEqual(len(macro_model.get_small_molecules()), 1)
        self.assertEqual(len(macro_model.get_molecules()), 1)
        self.assertEqual(self.molecule1.model, macro_model)
        self.assertEqual(len(macro_model.atoms), 3)
        macro_model.add_small_molecule(self.molecule2)
        self.assertEqual(len(macro_model.get_small_molecules()), 2)
        self.assertEqual(len(macro_model.get_molecules()), 2)
        self.assertEqual(self.molecule2.model, macro_model)
        self.assertEqual(len(macro_model.atoms), 6)


    def test_molecules_must_be_molecules(self):
        macro_model = MacroModel()
        with self.assertRaises(TypeError):
            macro_model.add_chain("mol")


    def test_molecules_must_be_unique(self):
        macro_model = MacroModel()
        self.molecule1.molecule_id = 23
        self.molecule2.molecule_id = 23
        macro_model.add_small_molecule(self.molecule1)
        with self.assertRaises(exceptions.DuplicateMoleculeError):
            macro_model.add_small_molecule(self.molecule2)



class ChainAdditionTests(MacroModelTest):

    def setUp(self):
        MacroModelTest.setUp(self)
        self.chain1 = Chain(self.residue1, self.residue2)
        self.chain2 = Chain(self.residue3, self.residue4)


    def test_can_add_chains(self):
        macro_model = MacroModel()
        macro_model.add_chain(self.chain1)
        self.assertEqual(len(macro_model.get_chains()), 1)
        self.assertEqual(len(macro_model.get_molecules()), 1)
        self.assertEqual(self.chain1.model, macro_model)
        self.assertEqual(len(macro_model.atoms), 6)
        macro_model.add_chain(self.chain2)
        self.assertEqual(len(macro_model.get_chains()), 2)
        self.assertEqual(len(macro_model.get_molecules()), 2)
        self.assertEqual(self.chain2.model, macro_model)
        self.assertEqual(len(macro_model.atoms), 12)


    def test_chains_must_be_chains(self):
        macro_model = MacroModel()
        with self.assertRaises(TypeError):
            macro_model.add_chain("chain")


    def test_chains_must_be_unique(self):
        macro_model = MacroModel()
        self.chain1.chain_id = "A"
        self.chain2.chain_id = "A"
        macro_model.add_chain(self.chain1)
        with self.assertRaises(exceptions.DuplicateChainError):
            macro_model.add_chain(self.chain2)



class ComplexAdditionTests(MacroModelTest):

    def setUp(self):
        MacroModelTest.setUp(self)
        self.chain1 = Chain(self.residue1, self.residue2)
        self.chain2 = Chain(self.residue3, self.residue4)
        self.complex = Complex(self.chain1, self.chain2)


    def test_can_add_complexes(self):
        macro_model = MacroModel()
        macro_model.add_complex(self.complex)
        self.assertEqual(len(macro_model.get_complexes()), 1)
        self.assertEqual(len(macro_model.get_chains()), 2)
        self.assertEqual(len(macro_model.get_molecules()), 2)
        self.assertEqual(len(macro_model.get_small_molecules()), 0)
        self.assertEqual(self.complex.model, macro_model)
        self.assertEqual(self.chain1.model, macro_model)
        self.assertEqual(self.chain2.model, macro_model)
        self.assertEqual(len(macro_model.atoms), 12)


    def test_complexes_must_be_complexes(self):
        macro_model = MacroModel()
        with self.assertRaises(TypeError):
            macro_model.add_complex("complex")


    def test_complexes_must_be_unique(self):
        macro_model = MacroModel()
        self.chain1.chain_id = "A"
        self.chain2.chain_id = "A"
        macro_model.add_chain(self.chain1)
        with self.assertRaises(exceptions.DuplicateChainError):
            macro_model.add_chain(self.chain2)


    def test_complexes_must_be_unique(self):
        self.complex.complex_id = 1
        complex2 = copy.deepcopy(self.complex)
        macro_model = MacroModel()
        macro_model.add_complex(self.complex)
        with self.assertRaises(exceptions.DuplicateComplexError):
            macro_model.add_complex(complex2)
