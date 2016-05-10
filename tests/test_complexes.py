from unittest import TestCase
from molecupy import exceptions
from molecupy.molecules import Atom, AtomicStructure
from molecupy.macromolecules import Residue, Chain, MacroModel, Complex

class ComplexTest(TestCase):

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
        self.chain1 = Chain(self.residue1, self.residue2)
        self.chain2 = Chain(self.residue3, self.residue4)


    def check_valid_complex(self, complex_, check_complex_id=False, check_complex_name=False):
        self.assertIsInstance(complex_, Complex)
        self.assertIsInstance(complex_, AtomicStructure)
        self.assertIsInstance(complex_.chains, set)
        if complex_.model is not None:
            self.assertIsInstance(complex_.model, MacroModel)
        if check_complex_id:
            self.assertIsInstance(complex_.complex_id, int)
        if check_complex_name:
            self.assertIsInstance(complex_.complex_name, str)
        self.assertRegex(
         str(complex_),
         r"<Complex \((\d+) chains\)>"
        )



class ComplexCreationTests(ComplexTest):

    def test_can_create_complex(self):
        complex_ = Complex(self.chain1, self.chain2)
        self.check_valid_complex(complex_)


    def test_can_create_complex_with_id(self):
        complex_ = Complex(self.chain1, self.chain2, complex_id=10)
        self.check_valid_complex(complex_, check_complex_id=True)


    def test_complex_id_must_be_int(self):
        with self.assertRaises(TypeError):
            complex_ = Complex(self.chain1, self.chain2, complex_id=1.1)
        with self.assertRaises(TypeError):
            complex_ = Complex(self.chain1, self.chain2, complex_id="10")


    def test_can_create_complex_with_name(self):
        complex_ = Complex(self.chain1, self.chain2, complex_name="MOL")
        self.check_valid_complex(complex_, check_complex_name=True)


    def test_complex_name_must_be_str(self):
        with self.assertRaises(TypeError):
            complex_ = Complex(self.chain1, self.chain2, complex_name=1)
